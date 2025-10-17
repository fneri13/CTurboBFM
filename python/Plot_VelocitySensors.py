#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
# from styles import *
import argparse

os.makedirs("Pictures", exist_ok=True)

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Plot monitor point data with a given time per revolution. Use 1 to plot in seconds.")
parser.add_argument("time_per_revolution", type=float, help="Time per revolution in seconds")
args = parser.parse_args()
time_per_revolution = args.time_per_revolution

def extractOrder(filepath):
    """Extracts the integer index from filenames like 'monitorPoint_12.csv'."""
    match = re.search(r'monitorPoint_(\d+)\.csv$', filepath)
    return int(match.group(1)) if match else -1

cwd = os.getcwd()
files = os.listdir(cwd + '/Monitor_Points')

inputFiles = [cwd + '/Monitor_Points/' + file for file in files if file.endswith('.csv')]
inputFiles = sorted(inputFiles, key=extractOrder)

dfs = []
for inputFile in inputFiles:
    df = pd.read_csv(inputFile)
    df["Velocity_Magnitude[m/s]"] = np.sqrt(df["Velocity_X[m/s]"]**2 + df["Velocity_Y[m/s]"]**2 + df["Velocity_Z[m/s]"]**2)
    dfs.append(df)

pressureSensors = []
for df in dfs:
    pressureSensors.append(df["Velocity_Magnitude[m/s]"])
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
    plt.plot(time/time_per_revolution, scale_oscillations(pressureSensors[i], deltaTheta) + theta, label=f"probe {i}")

if time_per_revolution != 1:
    plt.xlabel("Revs [-]")
else:
    plt.xlabel(r"$t$ [s]")
plt.ylabel(r"$\theta$ [deg]")
# plt.legend()
plt.grid(alpha=.3)
try:
    plt.savefig("Pictures/VelocitySensors.pdf", bbox_inches='tight')
except:
    print("Could not save figure. Pictures folder not existing.")

plt.show()