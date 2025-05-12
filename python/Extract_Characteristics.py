#!/usr/bin/env python3

#### Extract characteristics from the folders denominated *kPa*

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

def read_file(input_filename):
    df = pd.read_csv(input_filename)
    data_dict = {col: df[col].to_numpy() for col in df.columns}
    return data_dict


def read_solutions(inputFolder):
    folders = os.listdir(inputFolder)
    runFolders = [folder for folder in folders if 'kPa' in folder]
    runFolders.sort()
    print("Run folders found: ", runFolders)
    
    massflow = []
    PRtt = []
    TRtt = []
    ETAtt = []
    
    for run in runFolders:
        input_filename = inputFolder + '/' + run + '/turbo.csv'
        data = read_file(input_filename)
        massflow.append(data['Massflow'][-1])
        PRtt.append(data['PRtt'][-1])
        TRtt.append(data['TRtt'][-1])
        ETAtt.append(data['ETAtt'][-1])
    
    return {'Massflow':np.array(massflow), 'PRtt':np.array(PRtt), 'TRtt':np.array(TRtt), 'ETAtt':np.array(ETAtt)}



def write_characteristics(data):
    with open('characteristics.csv', 'w') as f:
        f.write("MassFlow [kg/s],PRtt [-],TRtt [-],ETAtt [-]\n")
        for i in range(len(data['Massflow'])):
            f.write(f"{data['Massflow'][i]:.6f},{data['PRtt'][i]:.6f},{data['TRtt'][i]:.6f},{data['ETAtt'][i]:.6f}\n")
    
    print("Characteristics saved to characteristics.csv")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the characteristics from the runs located in input folder.")
    parser.add_argument("input_folder", help="Path to the folder containing kPa nominated Runs")
    args = parser.parse_args()

    input_folder = args.input_folder
    
    data = read_solutions(input_folder)
    write_characteristics(data)
    
    
