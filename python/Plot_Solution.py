#!/usr/bin/env python3

#### Plot the turbo performance of the simulation in csv format

import numpy as np
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from styles import *



def GetDataDict(input_filename, field):
    # Read the first three lines to extract grid sizes
    with open(input_filename, 'r') as f:
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])
    
    
    df = pd.read_csv(input_filename, skiprows=3)
    data_dict = {col: df[col].to_numpy().reshape((ni, nj, nk)) for col in df.columns}
    
    if field not in data_dict.keys():
        # check if the field name with spaces exists in the dataset
        try:
            fieldName = field.replace('_', ' ')
            if fieldName not in data_dict.keys():
                raise ValueError
        except:
            raise ValueError(f"Field '{field}' not found in the data columns. Available fields: {list(data_dict.keys())}. For composed names, instead of spaces use underscores '_' in the field name.")
    else:
        fieldName = field
    new_dict = {}
    new_dict['x'] = data_dict['x']
    new_dict['y'] = data_dict['y']
    new_dict['z'] = data_dict['z']
    new_dict[field] = data_dict[fieldName]
    return new_dict



def PlotData(data, iSlice, jSlice, kSlice, fieldName):
    
    # get the dimensionality
    nDim = 0
    
    if iSlice==':':
        nDim += 1
    if jSlice==':':
        nDim += 1
    if kSlice==':':
        nDim += 1
    
    if (nDim !=1 and nDim !=2):
        raise ValueError("Only 1D or 2D plots are supported. Please provide exactly one ':' slice.")
    
    if nDim==1:
        Plot1D(data, iSlice, jSlice, kSlice, fieldName)
    elif nDim==2:
        Plot2D(data, iSlice, jSlice, kSlice, fieldName)




def Plot1D(data, iSlice, jSlice, kSlice, fieldName):
    if iSlice ==':':
        j = int(jSlice)
        k = int(kSlice)
        x = data['x'][:,j,k]
        values = data[fieldName][:,j,k]
        xlabel = 'X [m]'
        title = f'alongI_j:{j}_k:{k}'
    elif jSlice ==':':
        i = int(iSlice)
        k = int(kSlice)
        x = data['y'][i,:,k]
        values = data[fieldName][i,:,k]
        xlabel = 'Y [m]'
        title = f'alongJ_i:{i}_k:{k}'
    elif kSlice ==':':
        i = int(iSlice)
        j = int(jSlice)
        x = data['z'][i,j,:]*180/np.pi
        values = data[fieldName][i,j,:]
        xlabel = 'Z [m]'
        title = f'alongK_i:{i}_j:{j}'
    
    plt.figure()
    plt.plot(x, values, '-o', mfc='none')
    plt.xlabel(xlabel)
    plt.ylabel(fieldName)
    plt.grid(alpha=0.2)
    # plt.title(title)
    plt.savefig(f"pictures/{fieldName}_1D_plot_{title}.pdf", bbox_inches='tight')




def Plot2D(data, iSlice, jSlice, kSlice, fieldName):
    if iSlice ==':':
        if jSlice ==':':
            k = int(kSlice)
            x = data['x'][:,:,k]
            xlabel = 'X [m]'
            y = data['y'][:,:,k]
            ylabel = 'Y [m]'
            values = data[fieldName][:,:,k]
            title = 'kPlane:'+str(k)
        elif kSlice ==':':
            j = int(jSlice)
            x = data['x'][:,j,:]
            xlabel = 'X [m]'
            y = data['z'][:,j,:]
            ylabel = 'Z [m]'
            values = data[fieldName][:,j,:]
            title = 'jPlane:'+str(j)
    elif jSlice ==':':
        i = int(iSlice)
        x = data['y'][i,:,:]
        xlabel = 'Y [m]'
        y = data['z'][i,:,:]
        ylabel = 'Z [m]'
        values = data[fieldName][i,:,:]
        title = 'iPlane:'+str(i)
    else:
        raise ValueError("Something wrong in the selection process of the 2D contour")
    
    plt.figure()
    plt.contourf(x, y, values, cmap = 'turbo', levels=30)
    plt.colorbar()
    plt.title(fieldName)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f"pictures/{fieldName}_2D_plot_{title}.pdf", bbox_inches='tight')
        
    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the solution field of a CTurboBFM result.")
    
    parser.add_argument("inputFile", type=str, help="Path to the CSV file containing the simulation results.")
    parser.add_argument("iSlice", type=str, help="Slice along i-index")
    parser.add_argument("jSlice", type=str, help="Slice along j-index")
    parser.add_argument("kSlice", type=str, help="Slice along k-index")
    parser.add_argument("fieldName", type=str, help="Field name to plot")
    
    args = parser.parse_args()

    inputFile = args.inputFile
    iSlice = args.iSlice
    jSlice = args.jSlice
    kSlice = args.kSlice
    fieldName = args.fieldName
    
    os.makedirs("pictures", exist_ok=True)
    data = GetDataDict(inputFile, fieldName)
    PlotData(data, iSlice, jSlice, kSlice, fieldName)
    plt.show()
    