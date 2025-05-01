import numpy as np
import pandas as pd
from pyevtk.hl import gridToVTK

def read_structured_csv(filename):
    with open(filename, 'r') as f:
        # Read grid dimensions
        ni = int(f.readline().strip().split('=')[1])
        nj = int(f.readline().strip().split('=')[1])
        nk = int(f.readline().strip().split('=')[1])

    # Now skip only the first 3 lines (dimensions), keep header line
    df = pd.read_csv(filename, skiprows=3)

    # Check consistency
    n_points = ni * nj * nk
    assert len(df) == n_points, f"Expected {n_points} points but found {len(df)}"

    shape = (ni, nj, nk)

    data = {
        col: df[col].to_numpy().reshape(shape, order='C')  # Fortran order: i-fastest
        for col in df.columns
    }

    return data


def writeVTK(data, filename):
    # Compute velocity components
    u = data['rhoU'] / data['rho']
    v = data['rhoV'] / data['rho']
    w = data['rhoW'] / data['rho']

    pointsData = {
        "Density": np.ascontiguousarray(data['rho']),
        "Velocity": (
            np.ascontiguousarray(u),
            np.ascontiguousarray(v),
            np.ascontiguousarray(w)
        ),
        "Total Energy": np.ascontiguousarray(data['rhoE'] / data['rho'])
    }

    gridToVTK(
        filename,
        np.ascontiguousarray(data['x']),
        np.ascontiguousarray(data['y']),
        np.ascontiguousarray(data['z']),
        pointData=pointsData
    )



fileName = "results_124_48_10.csv"
outputFilename = fileName.split(".")[0]
data = read_structured_csv(fileName)
writeVTK(data, outputFilename)


