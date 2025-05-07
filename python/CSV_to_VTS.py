#!/usr/bin/env python3

#### Convert structured CSV data to VTS format. Given the path to the csv results file (obtained with CTurboBFM), outputs a VTS file with the same name, in the same location.

import numpy as np
import pandas as pd
from pyevtk.hl import gridToVTK
import argparse
import os

def read_structured_csv(filename):
    with open(filename, 'r') as f:
        # Read grid dimensions
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])

    # Skip only the first 3 lines (dimensions), keep header line
    df = pd.read_csv(filename, skiprows=3)

    # Check consistency
    n_points = ni * nj * nk
    assert len(df) == n_points, f"Expected {n_points} points but found {len(df)}"

    shape = (ni, nj, nk)

    data = {
        col: df[col].to_numpy().reshape(shape, order='C')  # i-fastest
        for col in df.columns
    }

    return data


def writeVTK(data, filename):
    
    pointsData = {
        "Velocity": (
            np.ascontiguousarray(data['Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Velocity Z'], dtype=np.float64)
        )
    }

    try:
        pointsData["Grid Velocity"] = (
            np.ascontiguousarray(data['Grid Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Grid Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Grid Velocity Z'], dtype=np.float64)
        )

        pointsData["Relative Velocity"] = (
            np.ascontiguousarray(data['Relative Velocity X'], dtype=np.float64),
            np.ascontiguousarray(data['Relative Velocity Y'], dtype=np.float64),
            np.ascontiguousarray(data['Relative Velocity Z'], dtype=np.float64)
        )
    except KeyError:
        pass
        
    
    for key in data.keys():
        if key not in ['Velocity X', 'Velocity Y', 'Velocity Z', 'x', 'y', 'z',
                       'Grid Velocity X', 'Grid Velocity Y', 'Grid Velocity Z',
                       'Relative Velocity X', 'Relative Velocity Y', 'Relative Velocity Z']:
            pointsData[key] = np.ascontiguousarray(data[key], dtype=np.float64)

    gridToVTK(
        filename,
        np.ascontiguousarray(data['x']),
        np.ascontiguousarray(data['y']),
        np.ascontiguousarray(data['z']),
        pointData=pointsData
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert structured CSV data to VTK format.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    args = parser.parse_args()

    input_filename = args.input_file
    output_filename = os.path.splitext(input_filename)[0]  # Strip extension

    data = read_structured_csv(input_filename)
    writeVTK(data, output_filename)
    print("Written VTS file to", output_filename)
