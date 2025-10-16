#!/usr/bin/env python3


import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

def delete_vts(inputFolder):
    folders = os.listdir(inputFolder)
    folders = [folder for folder in folders if os.path.isdir(inputFolder + '/' + folder)]
    for folder in folders:
        subfolders = os.listdir(inputFolder + '/' + folder)
        for subfolder in subfolders:
            if subfolder=='vts':
                os.system('rm -rf ' + inputFolder + '/' + folder + '/vts')
                print('Deleted ' + inputFolder + '/' + folder + '/vts')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the characteristics from the runs located in input folder.")
    parser.add_argument("input_folder", help="Path to the folder containing the simulations whose vts folder you want to remove")
    args = parser.parse_args()

    input_folder = args.input_folder
    
    data = delete_vts(input_folder)
    
    
