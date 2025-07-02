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

pressureSensors = []
for df in dfs:
    pressureSensors.append(df["Pressure[Pa]"])
time = dfs[0]["Time[s]"]

print()


def scale_oscillations(pressure, deltaTheta, scaling=1.1):
    y = (pressure-pressure.min())/(pressure.max()-pressure.min())*deltaTheta/scaling
    y = y-y[0]
    return y

nSensors = len(pressureSensors)
plt.figure()
for i in range(0, len(pressureSensors)):
    deltaTheta = 360 / (nSensors)
    theta = i * 360 / (nSensors)
    plt.plot(time, scale_oscillations(pressureSensors[i], deltaTheta) + theta, label=f"probe {i}")

plt.xlabel(r"$t$ [s]")
plt.ylabel(r"$\theta$ [deg]")
plt.legend()
plt.grid(alpha=.3)

plt.show()