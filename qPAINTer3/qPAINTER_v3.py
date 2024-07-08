"""
    qPAINTer v3 (simplified)
    ~~~~~~~~~~~~~~~~
    This version simply uses the number of binding events rather than mean dark times to estimate the number of binding sites present. 

    :author: Philipp Steen, 2024
"""

import numpy as np
import pandas as pd
from itertools import groupby
from operator import itemgetter
import multiprocessing
from pandarallel import pandarallel
import yaml
import matplotlib as mpl
import matplotlib.pyplot as plt


def Import(in_path):
    """Imports clustered .hdf5 from Picasso. Requires group info."""
    table = pd.read_hdf(in_path, key = 'locs')
    table = table.sort_values(by=['group', 'frame'])
    return(table)

def BindingEvents(frames):
    """Links localizations and returns all bright and dark times from a trace."""
    ranges =[]
    for k,g in groupby(enumerate(frames),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append((group[0],group[-1]))
    events = np.asarray(ranges)
    n_events = len(events)
    beginnings = (events[:,0])
    dark_times = beginnings[1:] - beginnings[:-1]
    bright_times = events[:,1] - events[:,0] + 1
    return(n_events, events, bright_times, dark_times)

def KineticsCalcs(group_df):
    """Finds consecutive binding events."""
    n_events, events, bright_times, dark_times = BindingEvents(np.asarray(group_df["frame"]))
    x = group_df['x'].mean()
    y = group_df['y'].mean()
    lpx = group_df['lpx'].mean()
    lpy = group_df['lpy'].mean()
    return(pd.Series(data=(n_events, events, bright_times,dark_times,x,y,lpx,lpy),index=["n_events","events","bright_times","dark_times","x","y","lpx","lpy"]))

def AssembleTable(table):
    pandarallel.initialize(nb_workers=min(30, multiprocessing.cpu_count()))
    table_k = table.groupby('group').parallel_apply(KineticsCalcs)
    return(table_k)

def MeanFrame(values):
    first_values = [sublist[0] for sublist in values]
    first_values_std = np.std(first_values)
    mean_frame = np.mean(first_values)
    return([first_values_std, mean_frame])

def ApplyFilter(df, filter_cutoff, ranges):
    df = df.drop(df[df['mf_std'] <= filter_cutoff].index)
    df = df.drop(df[df['mf'] <= ranges[0]].index)
    df = df.drop(df[df['mf'] >= ranges[1]].index)
    return(df)

def PlotFiltering(fig, axs, calib_table, count_table, xmax_calib, xmax_count, color, alpha):
    #Allows for visual inspection of the filtering process
    axs[0, 0].set_ylabel("Calibration")
    axs[0, 0].hist(calib_table["mf_std"], color = color, bins = np.linspace(0, xmax_calib[0], 20), alpha=alpha)
    axs[0, 0].set_xlabel("Mean Frame STD")
    axs[0, 1].hist(calib_table["mf"], color = color, bins = np.linspace(0, xmax_calib[1], 20), alpha=alpha)
    axs[0, 1].set_xlabel("Mean Frame")
    axs[1, 0].set_ylabel("Measurement")
    axs[1, 0].hist(count_table["mf_std"], color = color, bins = np.linspace(0, xmax_count[0], 20), alpha=alpha)
    axs[1, 0].set_xlabel("Mean Frame STD")
    axs[1, 1].hist(count_table["mf"], color = color, bins = np.linspace(0, xmax_count[1], 20), alpha=alpha)
    axs[1, 1].set_xlabel("Mean Frame")
    return(fig, axs)

class qPAINT:
    def __init__(self,
                 calib_path,
                 count_path):
        self.calib_path = calib_path
        self.count_path = count_path
        self.calib_table = []
        self.count_table = []

    def Load(self):
        self.calib_table = AssembleTable(Import(self.calib_path))
        self.count_table = AssembleTable(Import(self.count_path))

    def Filter(self, filter_cutoff, factor, showplot):
        #Filters sticking events from data
        #Determine the "mean frame" and standard deviation of the "mean frame" for each group
        self.calib_table['filter_results'] = self.calib_table['events'].parallel_apply(MeanFrame)
        self.calib_table[['mf_std', 'mf']] = pd.DataFrame(self.calib_table['filter_results'].tolist(), index=self.calib_table.index)
        self.calib_table = self.calib_table.drop(columns=['filter_results'])
        self.count_table['filter_results'] = self.count_table['events'].parallel_apply(MeanFrame)
        self.count_table[['mf_std', 'mf']] = pd.DataFrame(self.count_table['filter_results'].tolist(), index=self.count_table.index)
        self.count_table = self.count_table.drop(columns=['filter_results'])
        xmax_calib = self.calib_table['mf_std'].max(), self.calib_table['mf'].max()
        xmax_count = self.count_table['mf_std'].max(), self.count_table['mf'].max()
        if showplot:
            fig, axs = plt.subplots(2, 2, figsize=(10,8))
            fig, axs = PlotFiltering(fig, axs, self.calib_table, self.count_table, xmax_calib, xmax_count, "blue", 1)
        #Calculate necessary values
        calib_mean_frame_std = self.calib_table["mf"].std()
        calib_mean_frame = self.calib_table["mf"].mean()
        count_mean_frame_std = self.count_table["mf"].std()
        count_mean_frame = self.count_table["mf"].mean()
        #Apply filtering
        self.calib_table = ApplyFilter(self.calib_table, filter_cutoff, [calib_mean_frame-factor*calib_mean_frame_std, calib_mean_frame+factor*calib_mean_frame_std])
        self.count_table = ApplyFilter(self.count_table, filter_cutoff, [count_mean_frame-factor*count_mean_frame_std, count_mean_frame+factor*count_mean_frame_std])
        if showplot:
            fig, axs = PlotFiltering(fig, axs, self.calib_table, self.count_table, xmax_calib, xmax_count, "orange", 1)
        plt.show()

    def Calculate(self):
        #Calculate cluster sizes
        calib_n = np.mean(self.calib_table["n_events"])
        #calib_std = np.std(self.calib_table["n_events"]) 
        self.calib_table['cluster_size'] = self.calib_table["n_events"]/calib_n
        self.count_table['cluster_size'] = self.count_table["n_events"]/calib_n
    
    def SaveLargeClusters(self, threshold, destination):
        #Save pick file with large clusters for investigation in Picasso Render
        large = self.count_table[self.count_table['cluster_size']>=threshold]
        x = np.asarray(large['x'])
        y = np.asarray(large['y'])
        pick_diameter = 1.0
        def save_picks(name):
            yaml.default_flow_style = None
            picks = {"Centers": [[float(x[i]),float(y[i])] for i in range(0,len(y))],"Diameter": pick_diameter,"Shape": "Circle"}
            with open(str(name) + '_jup.yaml', "w") as f:
                yaml.dump(picks, f,default_flow_style=None)
        SaveEverything = save_picks(destination)

    def SaveCSV(self, destination):
        #Save a CSV file with the results
        reduced_calib_table = self.calib_table[['n_events', 'x', 'y', 'lpx', 'lpy', 'mf_std', 'mf', 'cluster_size']].copy()
        reduced_calib_table.to_csv(destination+'_qPAINT_calibration_results.csv', index = True)
        reduced_count_table = self.count_table[['n_events', 'x', 'y', 'lpx', 'lpy', 'mf_std', 'mf', 'cluster_size']].copy()
        reduced_count_table.to_csv(destination+'_qPAINT_counting_results.csv', index = True)


