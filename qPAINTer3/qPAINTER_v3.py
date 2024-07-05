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

class qPAINT:
    def __init__(self,
                 calib_path,
                 count_path):
        self.calib_path = calib_path
        self.count_path = count_path
        self.calib_table = []
        self.count_table = []

    def Start(self):
        self.calib_table = AssembleTable(Import(self.calib_path))
        self.count_table = AssembleTable(Import(self.count_path))
        calib_n = np.mean(self.calib_table["n_events"])
        #calib_std = np.std(self.calib_table["n_events"]) 
        self.calib_table['cluster_size'] = self.calib_table["n_events"]/calib_n
        self.count_table['cluster_size'] = self.count_table["n_events"]/calib_n
    
    def SaveLargeClusters(self, threshold, destination):
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
        reduced_calib_table = self.calib_table[['n_events', 'x', 'y', 'lpx', 'lpy', 'cluster_size']].copy()
        reduced_calib_table.to_csv(destination+'_qPAINT_calibration_results.csv', index = True)
        reduced_count_table = self.count_table[['n_events', 'x', 'y', 'lpx', 'lpy', 'cluster_size']].copy()
        reduced_count_table.to_csv(destination+'_qPAINT_counting_results.csv', index = True)


