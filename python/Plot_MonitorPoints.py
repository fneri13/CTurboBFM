#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from styles import *

os.makedirs("pictures", exist_ok=True)

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

for ii,field in enumerate(fieldList):
    plt.figure()
    
    if len(dfs) > 1:
        for i in range(1,len(dfs)):
            plt.plot(dfs[i]["Time[s]"], dfs[i][field], label=f"Point {i}")
    else:
        plt.plot(dfs[0]["Time[s]"], dfs[0][field])

    plt.xlabel("Time [s]")
    plt.ylabel(field)
    
    # if len(dfs) > 1:
    #     plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(f"pictures/monitorPoint_{ii}.pdf", bbox_inches='tight')
plt.show()
