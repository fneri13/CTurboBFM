#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

cwd = os.getcwd()
files = os.listdir(cwd+'/monitorPoints')
inputFiles = [cwd + '/monitorPoints/' + file for file in files if file.endswith('.csv')]
inputFiles = sorted(inputFiles)

dfs = []
for inputFile in inputFiles:
    df = pd.read_csv(inputFile)
    df["Velocity_Magnitude[m/s]"] = np.sqrt(df["Velocity_X[m/s]"]**2 + df["Velocity_Y[m/s]"]**2 + df["Velocity_Z[m/s]"]**2)
    dfs.append(df)


fieldList = [
            'Velocity_X[m/s]',
            'Velocity_Y[m/s]',
            'Velocity_Z[m/s]',
            'Pressure[Pa]',
            'Velocity_Magnitude[m/s]',
             ]

for field in fieldList:
    plt.figure()
    for i in range(1,len(dfs)):
        plt.plot(dfs[i]["Time[s]"], dfs[i][field], label=f"Point {i}")

    plt.xlabel("Time [s]")
    plt.ylabel(field)
    plt.legend()
    plt.grid(alpha=0.3)
plt.show()
