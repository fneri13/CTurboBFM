import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
# from Utils.styles import *
import os 

axialLocationIdx = -5

intputFiles = [
    '../axisymmetricChannel_2DSolver/Volume_CSV/results.csv',
    '../axisymmetricChannel_3DSolver/Volume_CSV/results.csv',
]

labels = [
    '2D',
    '3D',
]

linestyles = [
    'solid',
    'dashed',
]

datas = []
for inputFile in intputFiles:
    # Read the first three lines to extract grid sizes
    with open(inputFile, 'r') as f:
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])

    # Read the rest of the CSV data into a DataFrame
    df = pd.read_csv(inputFile, skiprows=3)
    data = {}

    for field in df.columns:
        data[field] = np.reshape(df[field].values, (ni, nj, nk))
    datas.append(data)


# PLOTS
axialLocationIdx = ni//2
for field in ['Velocity X', 'Velocity Y', 'Velocity Z', 'Density', 'Pressure', 'Entropy']:
    plt.figure()
    for i, data in enumerate(datas):
        plt.plot(data[field][axialLocationIdx,:,0], data['y'][axialLocationIdx,:,0], label=labels[i], linestyle=linestyles[i])
    plt.ylabel('r [m]')
    plt.xlabel(field)
    plt.legend()
    plt.savefig(field + '.pdf', bbox_inches='tight')


plt.show()