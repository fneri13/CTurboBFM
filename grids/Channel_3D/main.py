import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pickle
from TurboBFM.Preprocess.grid_generation import transfinite_grid_generation




"""
Generate a 2D rectangular geometry, that will be used as verification
"""
OUTPUT_FOLDER = 'Grid'
LX = 3
R_IN = 1
R_OUT = 2
NX = 96
NY = 32
NK = 100
THETA_MAX = 360*np.pi/180

x = np.linspace(0, LX, NX)
y = np.linspace(R_IN, R_OUT, NY)
theta = np.linspace(0, THETA_MAX, NK)



X, Y = np.meshgrid(x, y, indexing='ij')
NI, NJ = X.shape

# revolve in 3D
x3D = np.zeros((NI, NJ, NK))
y3D = np.zeros((NI, NJ, NK))
z3D = np.zeros((NI, NJ, NK))
for i in range(NI):
    for j in range(NJ):
        for k in range(NK):
            x3D[i,j,k] = X[i,j]
            y3D[i,j,k] = Y[i,j]*np.cos(theta[k])
            z3D[i,j,k] = Y[i,j]*np.sin(theta[k])


# Create a 3D scatter plots
mesh = pv.StructuredGrid(x3D, y3D, z3D)
plotter = pv.Plotter()
plotter.add_mesh(mesh, cmap='viridis', show_edges=True)
plotter.show_axes()


# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_%02i.csv' %(NI, NJ, NK), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK={NK}\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                file.write(f"{x3D[i,j,k]:.6e},{y3D[i,j,k]:.6e},{z3D[i,j,k]:.6e}\n")


plotter.show()
plt.show()







