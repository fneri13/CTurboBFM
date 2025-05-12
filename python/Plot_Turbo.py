#!/usr/bin/env python3

#### Plot the turbo performance of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

def read_file(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def shift_to_zero(res):
    return res - res[0]

def plot_turbo(data):
    
    # mass flow plot
    plt.figure()
    plt.plot(data["Massflow"])
    plt.grid(alpha=0.3)
    plt.xlabel("Iteration [-]")
    plt.ylabel(r"$\dot{m}$ [kg/s]")
    plt.savefig("pictures/massflow.pdf", bbox_inches='tight')
    
    # beta plot
    plt.figure()
    plt.plot(data["PRtt"], label=r"$\beta_{tt}$")
    plt.plot(data["TRtt"], label=r"$\tau_{tt}$")
    plt.plot(data["ETAtt"], label=r"$\eta_{tt}$")
    plt.grid(alpha=0.3)
    plt.xlabel("Iteration [-]")
    plt.legend()
    plt.savefig("pictures/performance.pdf", bbox_inches='tight')
    
    # beta evolution plot
    plt.figure()
    plt.plot(data["Massflow"], data["PRtt"], '-k', lw=0.5)
    plt.scatter(data["Massflow"], data["PRtt"], c=np.linspace(1, len(data["Massflow"]), len(data["Massflow"])))
    plt.grid(alpha=0.3)
    plt.xlabel(r"$\dot{m}$ [kg/s]")
    plt.ylabel(r"$\beta_{tt}$ [-]")
    plt.colorbar()
    plt.savefig("pictures/betaEvolution.pdf", bbox_inches='tight')
    
    # eta evolution plot
    plt.figure()
    plt.plot(data["Massflow"], data["ETAtt"], '-k', lw=0.5)
    plt.scatter(data["Massflow"], data["ETAtt"], c=np.linspace(1, len(data["Massflow"]), len(data["Massflow"])))
    plt.grid(alpha=0.3)
    plt.xlabel(r"$\dot{m}$ [kg/s]")
    plt.ylabel(r"$\eta_{tt}$ [-]")
    plt.colorbar()
    plt.savefig("pictures/etaEvolution.pdf", bbox_inches='tight')
        
    plt.show()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the performance of a CTurboBFMsimulation.")
    parser.add_argument("input_file", help="Path to the turbo performance CSV file")
    args = parser.parse_args()

    input_filename = args.input_file
    
    data = read_file(input_filename)
    os.makedirs("pictures", exist_ok=True)
    plot_turbo(data)
