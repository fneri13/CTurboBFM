import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


# File path
filename = "volumeData/results.csv"  # <-- Update this

# Read the header (first 3 lines)
with open(filename, 'r') as f:
    ni_line = f.readline().strip()
    nj_line = f.readline().strip()
    nk_line = f.readline().strip()
    header_line = f.readline().strip()

# Parse grid dimensions
NI = int(ni_line.split('=')[1])
NJ = int(nj_line.split('=')[1])
NK = int(nk_line.split('=')[1])

# Read the CSV data skipping the first 3 lines
data = pd.read_csv(filename, skiprows=3)

# Convert each column to NumPy arrays
arrays = {col: data[col].to_numpy() for col in data.columns}

# Optionally reshape to 3D arrays
shape = (NI, NJ)
arrays_3d = {key: val.reshape(shape, order='F') for key, val in arrays.items()}  # Fortran order (column-major) common in CFD data


plt.figure()
plt.contourf(arrays_3d['x'], arrays_3d['y'], arrays_3d['Pressure'], levels=32, cmap='turbo')
plt.colorbar(label='Pressure')
levels = np.linspace(0.16, 1.71, 32)
plt.contour(arrays_3d['x'], arrays_3d['y'], arrays_3d['Density'], levels=levels, colors='black', lw=0.1)
plt.title('Schulz Benchmark 3')
plt.xticks([])
plt.yticks([])
plt.gca().set_aspect('equal', adjustable='box')
# plt.savefig('results.pdf', bbox_inches='tight')
plt.show()