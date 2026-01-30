#!/usr/bin/env python3

#### Extract characteristics from the folders denominated *kPa*

import numpy as np
import pandas as pd
import argparse
import os
import subprocess
import matplotlib.pyplot as plt

def read_file(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def run_cases(inputFolder):
    """
    For each folder named like:
        KT_XXX
        PREQ_XXX
        P_XXX
        MF_XXX
    run the command:
        CTurboBFM input.ini
    """
    possibleNames = ['KT', 'PREQ', 'P', 'MF']
    
    folders = os.listdir(inputFolder)
    runFolders = [f for f in folders if f.split('_')[0] in possibleNames]

    # Sort using the number after the underscore
    runFolders.sort(key=lambda f: float(f.split('_')[1]))
    print("Run folders found:", runFolders)

    all_data = {
        "Massflow": [],
        "PRtt": [],
        "TRtt": [],
        "ETAtt": []
    }

    processes = []

    # Launch all runs in parallel (if < 10)
    if len(runFolders) < 10:
        print(f"Launching all {len(runFolders)} cases in parallel...")

        for run in runFolders:
            run_path = os.path.join(inputFolder, run)
            ini_file = os.path.join(run_path, "input.ini")

            if not os.path.exists(ini_file):
                print(f"WARNING: input.ini not found in {run_path}, skipping.")
                continue

            print(f"Starting CTurboBFM in {run_path} ...")

            p = subprocess.Popen(
                ["turbobfm", "input.ini"],
                cwd=run_path
            )
            processes.append((run, p))

        # Wait for all to finish
        for run, p in processes:
            ret = p.wait()
            if ret == 0:
                print(f"Run {run} completed successfully.")
            else:
                print(f"Run {run} FAILED with return code {ret}.")

    # Otherwise: run sequentially
    else:
        print("More than 10 runs — executing sequentially.")
        for run in runFolders:
            run_path = os.path.join(inputFolder, run)
            ini_file = os.path.join(run_path, "input.ini")

            if not os.path.exists(ini_file):
                print(f"WARNING: input.ini not found in {run_path}, skipping.")
                continue

            print(f"Running CTurboBFM in {run_path} ...")
            try:
                subprocess.run(
                    ["CTurboBFM", "input.ini"],
                    cwd=run_path,
                    check=True
                )
            except Exception as e:
                print(f"ERROR running CTurboBFM in {run}: {e}")

    return all_data



def write_characteristics(data):
    with open('characteristics.csv', 'w') as f:
        f.write("MassFlow [kg/s],PRtt [-],TRtt [-],ETAtt [-]\n")
        for i in range(len(data['Massflow'])):
            f.write(f"{data['Massflow'][i]:.6f},{data['PRtt'][i]:.6f},{data['TRtt'][i]:.6f},{data['ETAtt'][i]:.6f}\n")
    
    print("Characteristics saved to characteristics.csv")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract the characteristics from the runs located in input folder."
    )
    parser.add_argument("input_folder", help="Path to the folder containing kPa-named runs")
    args = parser.parse_args()

    input_folder = args.input_folder
    
    data = run_cases(input_folder)
    write_characteristics(data)
