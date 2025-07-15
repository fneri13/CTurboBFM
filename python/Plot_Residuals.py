#!/usr/bin/env python3

#### Plot the residuals of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from styles import *

def read_residuals(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def shift_to_zero(res):
    return res - res[0]

def plot_residuals(data):
    plt.figure()
    for key in data:
        plt.plot(shift_to_zero(data[key]), label=key)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xlabel("Iteration [-]")
    plt.ylabel("Residuals Drop [-]")
    plt.savefig("pictures/residuals.pdf", bbox_inches='tight')
    plt.show()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the residuals of a CTurboBFMsimulation.")
    parser.add_argument("input_file", help="Path to the inresiduals CSV file")
    args = parser.parse_args()

    input_filename = args.input_file
    
    data = read_residuals(input_filename)
    os.makedirs("pictures", exist_ok=True)
    plot_residuals(data)
