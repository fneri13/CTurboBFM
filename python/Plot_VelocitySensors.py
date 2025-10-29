#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
# from styles import *
import argparse
from Plot_MonitorPoints import read_monitor_points_data
from Plot_PressureSensors import plot_normalized_sensors_oscillation


def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Plot monitor point data with a given time per revolution. Use 1 to plot in seconds.")
    parser.add_argument("time_per_revolution", type=float, help="Time per revolution in seconds")
    args = parser.parse_args()
    time_per_revolution = args.time_per_revolution

    dfs = read_monitor_points_data()

    velocitySensors = []
    for df in dfs:
        velocitySensors.append(df["Velocity_Magnitude[m/s]"])
    time = dfs[0]["Time[s]"]

    plot_normalized_sensors_oscillation(velocitySensors, time, time_per_revolution, "Velocity")


if __name__ == "__main__":
    main()

plt.show()