#!/usr/bin/env python3

#### Plot the turbo performance of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from styles import *

def read_file(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def shift_to_zero(res):
    return res - res[0]

def plot_greitzer(data):
    
    # discard the first iteration, which was the dumb initialization
    time = data["Time[s]"][1:]
    plenumPressure = data["PlenumPressure[Pa]"][1:]
    plenumInletMassflow = data["PlenumInletMassflow[kg/s]"][1:]
    plenumOutletMassflow = data["PlenumOutletMassflow[kg/s]"][1:]
    
    
    fig, axes = plt.subplots(1, 2, figsize=(9,4))
    
    axes[0].plot(time*1E3, plenumPressure/1E5)
    
    axes[1].plot(time*1E3, plenumInletMassflow, label="Plenum Inlet")
    axes[1].plot(time*1E3, plenumOutletMassflow, label="Plenum Outlet")
    axes[1].legend(ncol=1)
    
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_xlabel(r'$t \ \rm{[ms]}$')
    axes[0].set_ylabel(r'$p$ [bar]')
    axes[1].set_ylabel(r'$\dot{{m}}$ [kg/s]')
    
    plt.tight_layout()
    plt.savefig("Pictures/greitzer_dynamics.pdf", bbox_inches='tight')
        
    plt.show()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the Greitzer dynamics of a CTurboBFMsimulation.")
    parser.add_argument("input_file", help="Path to the turbo performance CSV file")
    args = parser.parse_args()

    input_filename = args.input_file
    
    data = read_file(input_filename)
    os.makedirs("Pictures", exist_ok=True)
    plot_greitzer(data)
