from qPAINTER_v3 import *

Measurement = qPAINT('/.../...cell1_dbscan_SBS.hdf5', # file containing single binding sites.hdf5
                   '/.../...cell1_dbscan.hdf5') # file containing clusters to be counted.hdf5
Measurement.Load()
Measurement.Filter(10000, 2, showplot=True)
Measurement.Calculate()
Measurement.SaveLargeClusters(10, '/.../name') # path to where results should be saved/name
Measurement.SaveCSV('/.../name') # path to where results should be saved/name

#Optional: The results can be visualized here
'''
import matplotlib as mpl
import matplotlib.pyplot as plt

fix, ax = plt.subplots(1, figsize=(10,10))
ax.scatter(Measurement.count_table['x'], Measurement.count_table['y'],s=Measurement.count_table['cluster_size'])
ax.set_aspect('equal')
plt.show()

fix, ax = plt.subplots(1, figsize=(10,10))
ax.hist(Measurement.calib_table["cluster_size"],bins=5)
plt.show()

fix, ax = plt.subplots(1, figsize=(10,10))
ax.hist(Measurement.count_table["cluster_size"],bins=100)
plt.show()

'''
#print(Measurement.count_table.head(50))