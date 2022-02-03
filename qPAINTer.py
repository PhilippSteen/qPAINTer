"""
    qPAINTer
    ~~~~~~~~~~~~~~~~
    
    :author: Philipp Steen, 2021
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
import math
import sys
import os
import yaml
import statistics
from scipy.stats import norm

from matplotlib import cm
from matplotlib.colors import Normalize

import warnings
warnings.filterwarnings("ignore")

#Link Localizations
def DetermineBursts(group, ignore_dark): #Identifies the bright segments, ignoring dark times below the specified threshold
    bursts = []
    ignore_dark = ignore_dark+1
    burst = []
    previous_frame = -10000
    for index, row in group.iterrows():
        frame = row["frame"]
        if frame - previous_frame <= ignore_dark:
            burst.append(frame)
        else:
            bursts.append(burst)
            burst = []
            burst.append(frame)
        previous_frame = frame
    bursts.append(burst)
    return(bursts[1:]) #The first entry of bursts is always [], so it can be discarded

def RemoveShortBrights(bursts, ignore_bright): #Removes bright segments shorter than a specified duration
    new_bursts = []
    for entry in bursts:
        if entry[-1]-entry[0]>=ignore_bright:
            new_bursts.append(entry)
    return(new_bursts)

def CalculateBright(bursts): #Calculates the lengths of the bright times
    brights = []
    for entry in bursts:
        brights.append(entry[-1]-entry[0]+1)
    return(brights)

def CalculateDark(bursts): #Calculates the lengths of the dark times
    darks = []
    last = 0
    for entry in bursts:
        first = entry[0]
        difference = first-last
        last = entry[-1]
        darks.append(difference-1)
    return(darks[1:]) #The first entry is discarded, as it is just the time from the start of the measurement until the first bright segment
        
def LinkLocalizations(group, ignore_dark, ignore_bright):
    bursts = DetermineBursts(group, ignore_dark)
    bursts = RemoveShortBrights(bursts, ignore_bright)
    bright = CalculateBright(bursts)
    dark = CalculateDark(bursts)
    return(bright, dark, bursts) #Returns bright times and dark times as arrays

def MeanFrame(bursts, lower=0, upper=1e100, std_cutoff=0):
    burst_locations = []
    for entry in bursts:
        burst_locations.append(np.mean(entry))
    all_frames = np.sum(burst_locations)
    mean_frame = all_frames/len(burst_locations)
    mean_frame_std = np.std(burst_locations)
    if (mean_frame >= lower) and (mean_frame <= upper) and (mean_frame_std >= std_cutoff):
        return(mean_frame, mean_frame_std)
    else:
        return(-1000, -1000)

def ScreenLocs(file_name, bursts_table, show_screening, std_factor = 2, std_cutoff = 0):
    mean_frames = []
    mean_frame_stds = []
    good_picks = []
    
    for i in range(len(bursts_table)):
        mean_frame, mean_frame_std = MeanFrame(bursts_table[i])
        mean_frames.append(mean_frame)
        mean_frame_stds.append(mean_frame_std)
    
    mean,std=norm.fit(mean_frames)
    lower = mean-std_factor*std
    upper = mean+std_factor*std
    
    if show_screening == True:
        mpl.style.use('seaborn-poster')
        fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', figsize = (8,6))
        fig.tight_layout()

        #The following makes the histogram and gauss fit have the same y-scaling (the gauss is normalized otherwise)
        n, bins, patches = axs[0, 0].hist(mean_frames, bins=20, density = True)
        n1, bins1, patches1 = axs[0, 0].hist(mean_frames, bins=20, density = False, color = "grey")
        scalefactors = n1/n
        good_scalefactors = []
        for i in scalefactors:
            if math.isnan(i) == False:
                good_scalefactors.append(i)
        scalefactor = np.mean(good_scalefactors)
        
        xmin = np.min(mean_frames)
        xmax = np.max(mean_frames)
        x = np.linspace(xmin, xmax, 100)
        y = norm.pdf(x, mean, std)
        
        axs[0, 0].set_title('Mean frame analysis')
        axs[0, 0].plot(x, y*scalefactor, color = "red", label = "Gaussian fit")
        axs[0, 0].axvline(x=mean, linewidth = 0.8, color = 'black', label = "Mean")
        axs[0, 0].axvline(x=lower, linewidth = 0.8, color = 'teal', label = "Cut-off")
        axs[0, 0].axvline(x=upper, linewidth = 0.8, color = 'teal')
        axs[0, 0].set_ylabel('Counts')
        axs[0, 0].legend(loc = "lower left")
        
        axs[0, 1].set_title('Mean frame standard\ndeviation analysis')
        axs[0, 1].hist(mean_frame_stds, bins = 20, density = False, color = "grey")
        axs[0, 1].axvline(x=std_cutoff, linewidth = 0.8, color = 'teal', label = "Cut-off")
        axs[0, 1].legend(loc = "lower right")
    
    mean_frames_clean = []
    mean_frame_stds_clean = []

    for i in range(len(bursts_table)):
        mean_frame, mean_frame_std = MeanFrame(bursts_table[i], lower, upper, std_cutoff)
        if (mean_frame > 0) and (mean_frame_std > 0):
            mean_frames_clean.append(mean_frame)
            mean_frame_stds_clean.append(mean_frame_std)
            good_picks.append(i)
    
    if show_screening == True:
        axs[1, 0].hist(mean_frames_clean, bins=20, density = False, color = "grey")
        axs[1, 0].set_ylabel('Counts')
        axs[1, 0].set_xlabel('Frames')
        axs[1, 1].hist(mean_frame_stds_clean, bins = 20, density = False, color = "grey")
        axs[1, 1].set_xlabel('Frames')
        plt.subplots_adjust(wspace=0.02, hspace=0.02)
        #plt.tight_layout()
        plt.show()
        
        if file_name != False:
            fig.savefig(file_name+' - filtering.png', transparent=False, bbox_inches='tight')
        
    return(good_picks)

def expfunc(t,tau):
    return((1 - np.exp(-(t/tau)))) #This function is used for fitting


class qPAINT:
    def __init__(self, 
                 exposure_time = 0.1, 
                 calibration_path = "",
                 calibration_sites = 1,
                 data_path = "",
                 dbclusters_path = "",
                 title_0 = "Title",
                 title_1 = "Calibration",
                 title_2 = "Data",
                 ignore_dark = 1.0,
                 ignore_bright = 0,
                 prompt_cutoff = False,
                 save_everything = False,
                 output_location = "./outputs/"): #Specify default path
        
        self.exposure_time = exposure_time
        self.calibration_path = calibration_path
        self.calibration_sites = calibration_sites
        self.data_path = data_path
        self.dbclusters_path = dbclusters_path
        self.prompt_cutoff = prompt_cutoff
        self.title_0 = title_0
        self.title_1 = title_1
        self.title_2 = title_2
        self.ignore_dark = ignore_dark #Dark times of x frames or less are ignored
        self.ignore_bright = ignore_bright #Bright times of x frames and less are ignored
        self.save_everything = save_everything
        self.output_location = output_location
        
        self.data_fulltable = 0
        self.actual_numbers_df = 0
        
        self.mean_dark = 0
        self.mean_dark_std = 0
        
        self.actual_numbers = []
        self.counted_indices = []

    # Calculates the mean dark time for each group by conducting exponential fits for the CDF of the dark times of each group
    def CalcDark(self, input_path, title, screen_picks, std_factor, std_cutoff, show_screening, lowercutoff = 0, uppercutoff = 1e100):
        fulltable = pd.read_hdf(input_path, key = 'locs')
        fulltable.sort_values(by=['group', 'frame'])
        bright_time_table = []
        dark_time_table = []
        bursts_table = []
        index_table = []
        
        for i in range(max(fulltable["group"])+1): 
            bright_time, dark_time, bursts = LinkLocalizations((fulltable[fulltable["group"] == (i) ]), self.ignore_dark, self.ignore_bright)
            bright_time_table.append(bright_time)
            dark_time_table.append(dark_time)
            bursts_table.append(bursts)
        
        if self.save_everything == True:
            file_name = self.output_location+title
        else:
            file_name = False
        
        if screen_picks == True:
            good_picks = ScreenLocs(file_name, bursts_table, show_screening, std_factor, std_cutoff)
        else:
            good_picks = fulltable["group"]
        

        indices = []
        dark_times_by_spot = []
        
        mpl.style.use('seaborn-poster')
        fig, ax = plt.subplots(1, figsize = (8,6))
        fig.tight_layout()
        ax.set_xscale('log')
        ax.set_xlabel("Dark time (s)")
        ax.set_ylabel("Probability")   
        ax.set_title(title+"\nExponential dark time fits for all calculated spots")
        
        for i in range(len(dark_time_table)): 
            if i in good_picks:
                bright_time = bright_time_table[i]
                dark_time = dark_time_table[i]
                if len(dark_time) >= 2:
                    dark_times_s = [element * self.exposure_time for element in dark_time]
                    dark_sorted = np.sort(dark_times_s)
                    p = 1. * np.arange(len(dark_sorted)) / (len(dark_sorted) - 1)
                    popt, pcov = curve_fit(expfunc, dark_sorted, p)
                    if (popt[0] > self.ignore_dark) and (popt[0] > lowercutoff) and (popt[0] < uppercutoff):
                        dark_times_by_spot.append(popt[0])
                        indices.append(i)
                        plottingaxis = np.linspace(0, dark_sorted[-1], 5000)
                        ax.plot(plottingaxis, expfunc(plottingaxis, *popt), linewidth = "0.8", color = 'darkorange')
                        ax.axvline(x=popt[0], linewidth = 0.6, color = 'purple') 
                        ax.plot(dark_sorted, p, linewidth = "0.5", color = "black")
        
        final_table = pd.DataFrame({'Pick':indices, 'CDF-fitted dark time (s)':dark_times_by_spot})
        
        mean_mean_dark = np.mean(dark_times_by_spot)
        
        ax.axvline(x=mean_mean_dark, linewidth = 0.6, color = 'green', label = r'$\overline{\tau_{d*}}=$' + str(np.round(mean_mean_dark, decimals=2)) + " s")
        ax.legend(loc = "upper left")
        
        plt.tight_layout()
        plt.show()
        
        if self.save_everything == True:
            final_table.to_csv(file_name+' - dark_times.csv', index = False)
            fig.savefig(file_name+' - dark_times_fits.png', transparent=False, bbox_inches='tight')
        
        return(final_table)

    def Calibrate(self, screen_picks = False, std_factor = 2, std_cutoff = 1000, show_screening = False):
        path = self.calibration_path
        title = self.title_0 + ", " + self.title_1
        calibration_table = self.CalcDark(path, title, screen_picks, std_factor, std_cutoff, show_screening)
        
        if self.prompt_cutoff == True:
            print("Perhaps you would like to specify lower and upper cut-offs?")
            yes_or_no = input("Yes or No: ")
            if yes_or_no == "Yes" or yes_or_no == "yes" or yes_or_no == "Y" or yes_or_no == "y":
                lowercutoff = int(input("Lower cut-off: "))
                uppercutoff = int(input("Upper cut-off: "))
                calibration_table = self.CalcDark(path, title, screen_picks, std_factor, std_cutoff, show_screening, lowercutoff, uppercutoff)
            
        self.mean_dark = calibration_table["CDF-fitted dark time (s)"].mean()
        self.mean_dark_std = calibration_table["CDF-fitted dark time (s)"].std()
        
        mpl.style.use('seaborn-poster')
        fig, ax = plt.subplots(1, figsize = (8,6))
        fig.tight_layout()
        ax.hist(calibration_table["CDF-fitted dark time (s)"], bins = 60)
        ax.set_xlabel("Dark time (s)")
        ax.set_ylabel("Counts")   
        ax.set_title(title+"\nAll fitted dark times")
        ax.axvline(x=self.mean_dark, linewidth = 0.6, color = 'green', label = r'$\overline{\tau_{d*}}=$' + str(np.round(self.mean_dark, decimals=2)) + " s")
        ax.legend(loc = "upper right")
        plt.show()
        if self.save_everything == True:
            file_name = self.output_location+self.title_0+", "+self.title_1
            fig.savefig(file_name+' dark_times_hist.png', transparent=False, bbox_inches='tight')
        
    def Count(self, screen_picks = False, std_factor = 2, std_cutoff = 1000, show_screening = False):
        path = self.data_path
        title = self.title_0 + ", " + self.title_2
        
        #To find the x-y positions for each counted pick
        fulltable_xy = pd.read_hdf(path, key = 'locs')
        fulltable_xy.sort_values(by=['group', 'frame'])
        positions_table = fulltable_xy.groupby('group').mean()
        
        self.data_fulltable = self.CalcDark(path, title, screen_picks, std_factor, std_cutoff, show_screening)
        
        if self.prompt_cutoff == True:
            print("Perhaps you would like to specify lower and upper cut-offs?")
            yes_or_no = input("Yes or No: ")
            if yes_or_no == "Yes" or yes_or_no == "yes" or yes_or_no == "Y" or yes_or_no == "y":
                lowercutoff = int(input("Lower cut-off: "))
                uppercutoff = int(input("Upper cut-off: "))
                self.data_fulltable = self.CalcDark(path, title, screen_picks, std_factor, std_cutoff, show_screening, lowercutoff, uppercutoff)
            
        numbers_of_binding_sites = []
        x = []
        y = []
        checkgroup = []
        
        for index, row in self.data_fulltable.iterrows():
            numbers_of_binding_sites.append(np.round((self.mean_dark*self.calibration_sites)/row["CDF-fitted dark time (s)"]))
            
            desired_info = positions_table.loc[[row["Pick"]]]
            
            x.append(float(desired_info["x"]))
            y.append(float(desired_info["y"]))
            
        self.data_fulltable['number of binding sites'] = numbers_of_binding_sites
        self.data_fulltable['x'] = x
        self.data_fulltable['y'] = y

    def Discuss(self):
        mpl.style.use('seaborn-poster')
        fig, ax = plt.subplots(1, figsize = (8,6))
        plt.tight_layout()
        
        bins = np.arange(np.min(self.data_fulltable['number of binding sites']), np.max(self.data_fulltable['number of binding sites']) + 1.5) - 0.5
        
        ax.hist(self.data_fulltable['number of binding sites'], bins, density = False, rwidth=0.8)
        ax.set_xticks(bins + 0.5)
        
        ax.set_title(self.title_0 + ",\n " + self.title_2+"\nqPAINT-calculated number of binding sites")
        ax.set_xlabel("Calculated number of binding sites per cluster")
        ax.set_ylabel("Number of clusters")
        
        plt.tight_layout()
        plt.show()
        
        print(self.data_fulltable)
        
        if self.save_everything == True:
            file_name = self.output_location+self.title_0+", "+self.title_2
            self.data_fulltable.to_csv(file_name+' - qPAINT-results.csv', index = False)
            fig.savefig(file_name+' - qPAINT-results.png', transparent=False, bbox_inches='tight')
    
    def DiscussWithSizes(self):
        mpl.style.use('seaborn-poster')
        fig, ax = plt.subplots(1, figsize = (8,6))
        plt.tight_layout()
        
        bins = np.arange(np.min(self.data_fulltable['number of binding sites']), np.max(self.data_fulltable['number of binding sites']) + 1.5) - 0.5
        
        ax.hist(self.data_fulltable['number of binding sites'], bins, density = False, rwidth=0.8)
        ax.set_xticks(bins + 0.5)
        
        ax.set_title(self.title_0 + ",\n " + self.title_2+"\nqPAINT-calculated number of binding sites")
        ax.set_xlabel("Calculated number of binding sites per cluster")
        ax.set_ylabel("Number of clusters")
        
        plt.tight_layout()
        plt.show()
        
        #Cluster size section
        clusters_table = pd.read_hdf(self.dbclusters_path, key = 'clusters')
        convex_hull = []
        for index, row in self.data_fulltable.iterrows():
            desired_info = clusters_table.loc[row["Pick"]]
            convex_hull.append(float(desired_info["convex_hull"]))      
        
        self.data_fulltable['convex_hull'] = convex_hull
        print(self.data_fulltable)
        
        if self.save_everything == True:
            file_name = self.output_location+self.title_0+", "+self.title_2
            self.data_fulltable.to_csv(file_name+' - qPAINT-results.csv', index = False)
            fig.savefig(file_name+' - qPAINT-results.png', transparent=False, bbox_inches='tight')
    
    def SpatialPlot(self, display=False, cluster_size_cutoff=0):
        sums = []
        xvals = []
        yvals = []
        
        for index, row in self.data_fulltable.iterrows():
            sums.append(row["number of binding sites"])
            xvals.append(row["x"])
            yvals.append(row["y"])
            
        if display == True:
            dx = np.full(len(xvals), 0.3)
            dy = np.full(len(xvals), 0.3)
            starting_z = np.zeros(len(xvals))
            fig=plt.figure(figsize = (10,8))
            cmap = cm.get_cmap('plasma')
            norm = Normalize(vmin=min(sums), vmax=max(sums))
            colors = cmap(norm(sums))
            ax = fig.add_subplot(111, projection='3d')
            ax.bar3d(xvals, yvals, starting_z, dx, dy, sums, color=colors)
            ax.set_xlabel('X-position')
            ax.set_ylabel('Y-position')
            ax.set_zlabel('Number of binding sites per cluster')
            plt.tight_layout()
            plt.show()
            
            if self.save_everything == True:
                file_name = self.output_location+self.title_0 + ", " + self.title_2 + " - Centers_of_Clusters"
                fig.savefig(file_name+'.png', transparent=False, bbox_inches='tight')
        
        if self.save_everything == True:
            very_large_x = []
            very_large_y = []
            very_large_s = []
            
            #Find only the very large clusters
            for index, value in enumerate(sums):
                if value >= cluster_size_cutoff:
                    very_large_s.append(value)
                    very_large_x.append(xvals[index])
                    very_large_y.append(yvals[index])
            
            file_name = self.output_location+self.title_0 + ", " + self.title_2 + " - Centers_of_Clusters"
            
            pick_diameter = 1.0
            def save_picks_1(name):
                yaml.default_flow_style = None
                picks = {"Centers": [[float(very_large_x[i]),float(very_large_y[i])] for i in range(0,len(very_large_y))],"Diameter": pick_diameter,"Shape": "Circle"}
                with open(str(name) + '_jup.yaml', "w") as f:
                    yaml.dump(picks, f,default_flow_style=None)
            SaveEverything = save_picks_1(file_name)


##############################################
##########     USER INPUT BELOW     ##########
##############################################

Cell = qPAINT(exposure_time = 0.1,
                    calibration_path = "./testfiles/testfile_1.hdf5",
                    data_path = "./testfiles/testfile_2.hdf5",
                    dbclusters_path = "./testfiles/testfile_3.hdf5",
                    ignore_dark = 5, #In units of frames
                    prompt_cutoff = False,
                    title_0 = "Sample Cell",
                    title_1 = "Calibration",
                    title_2 = "Counting",
                    output_location = "./outputs/",
                    save_everything = True) #Set to True to save the output files

Cell.Calibrate(screen_picks = True, std_factor = 2, std_cutoff = 4000, show_screening = True)
Cell.Count(screen_picks = True, std_factor = 2, std_cutoff = 4000, show_screening = True)
#Cell.Discuss()
Cell.DiscussWithSizes()
Cell.SpatialPlot(display = True, cluster_size_cutoff=10)
