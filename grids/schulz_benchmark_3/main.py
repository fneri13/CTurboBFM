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
LX = 1
LY = 1
NX = 400
NY = 400

x = np.linspace(0, LX, NX)
y = np.linspace(0, LY, NY)

xInterface = 0.5
yInterface = 0.5

X, Y = np.meshgrid(x, y, indexing='ij')
NI, NJ = X.shape


# Create a 3D scatter plots
mesh = pv.StructuredGrid(X, Y, np.zeros_like(X))
plotter = pv.Plotter()
plotter.add_mesh(mesh, cmap='viridis', show_edges=True)
plotter.show_axes()



density = np.zeros_like(X)
velocity_x = np.zeros_like(X)
velocity_y = np.zeros_like(X)
velocity_z = np.zeros_like(X)
temperature = np.zeros_like(X)

rho = [1.5, 0.5323, 0.138, 0.5323]
p = [1.5, 0.3, 0.029, 0.3]
ux = [0,1.206,1.206,0]
uy = [0,0,1.206,1.206]

Rgas = 287

for i in range(NI):
    for j in range(NJ):
        if X[i,j] > xInterface and Y[i,j] > yInterface:
            density[i,j] = rho[0]
            velocity_x[i,j] = ux[0]
            velocity_y[i,j] = uy[0]
            temperature[i,j] = p[0] / (Rgas * rho[0])
        elif X[i,j] < xInterface and Y[i,j] > yInterface:
            density[i,j] = rho[1]
            velocity_x[i,j] = ux[1]
            velocity_y[i,j] = uy[1]
            temperature[i,j] = p[1] / (Rgas * rho[1])
        elif X[i,j] < xInterface and Y[i,j] < yInterface:
            density[i,j] = rho[2]
            velocity_x[i,j] = ux[2]
            velocity_y[i,j] = uy[2]
            temperature[i,j] = p[2] / (Rgas * rho[2])
        elif X[i,j] > xInterface and Y[i,j] < yInterface:
            density[i,j] = rho[3]
            velocity_x[i,j] = ux[3]
            velocity_y[i,j] = uy[3]
            temperature[i,j] = p[3] / (Rgas * rho[3])



grid = {'X': X, 'Y': Y}
# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_fullData.csv' %(NX, NY), 'w') as file:
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK=1\n")
    file.write("x,y,z,Density,Velocity X,Velocity Y,Velocity Z,Temperature\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(1):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f},{density[i,j]:.6f},{velocity_x[i,j]:.6f},{velocity_y[i,j]:.6f},{velocity_z[i,j]:.6f},{temperature[i,j]:.6f}\n")
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_onlyCoords.csv' %(NX, NY), 'w') as file:
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK=1\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(1):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f}\n")


plotter.show()
plt.show()







