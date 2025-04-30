import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pickle
from scipy.optimize import fsolve
from TurboBFM.Preprocess.grid_generation import transfinite_grid_generation


"""
Test case taken from section 11.5.2 of Numerical Computation of Internal and External Flows: 
The Fundamentals of Computational Fluid Dynamics (Second Edition) by Charles Hirsch.
"""

OUTPUT_FOLDER = 'Grid'
L = 1
NX = 128
NY = 48
NK = 10
SPAN = L/8
STREAMWISE_COEFF = 1.1
SPANWISE_COEFF = 1.1





def func(alpha):
    y = np.cos(alpha)/2/np.sin(alpha) + 0.1-1/2/np.sin(alpha)
    return y
alpha = fsolve(func, 1)[0]
r_bump = 1/2/np.sin(alpha)

NX_bump = NX//3

x1 = np.linspace(0, L, NX_bump)
y1 = np.zeros_like(x1)

theta = np.linspace(0, 2*alpha, NX_bump)
x2 = 1.5*L + r_bump*np.cos(np.pi/2+alpha-theta)
y2 = -(r_bump-0.1*L)+r_bump*np.sin(np.pi/2+alpha-theta)

x3 = np.linspace(2*L, 3*L, NX_bump)
y3 = np.zeros_like(x3)

# three blocks
x_wall = [x1, x2, x3]
y_wall = [y1, y2, y3]

x_inlet = [np.zeros(NY),
           np.zeros(NY)+x1[-1],
           np.zeros(NY)+x2[-1]]
y_inlet = [np.linspace(0, L, NY),
           np.linspace(0, L, NY),
           np.linspace(0, L, NY)]

x_outlet = [np.zeros(NY)+x1[-1],
           np.zeros(NY)+x2[-1],
           np.zeros(NY)+x3[-1]]
y_outlet = [np.linspace(0, L, NY),
           np.linspace(0, L, NY),
           np.linspace(0, L, NY)]

x_up = [np.linspace(0, L, NX_bump),
        np.linspace(L, 2*L, NX_bump),
        np.linspace(2*L, 3*L, NX_bump)]
y_up = [np.zeros(NX_bump)+L,
        np.zeros(NX_bump)+L,
        np.zeros(NX_bump)+L]

Xmulti, Ymulti = [], []
stretch_stream = ['right', 'both', 'left']
stretch_span = ['bottom', 'bottom', 'bottom']
for i in range(3):
    xgrid, ygrid = transfinite_grid_generation(np.vstack((x_inlet[i], y_inlet[i])), 
                                               np.vstack((x_wall[i], y_wall[i])), 
                                               np.vstack((x_outlet[i], y_outlet[i])), 
                                               np.vstack((x_up[i], y_up[i])),
                                               stretch_type_stream=stretch_stream[i], stretch_type_span=stretch_span[i],
                                               streamwise_coeff=STREAMWISE_COEFF, spanwise_coeff=SPANWISE_COEFF)
    Xmulti.append(xgrid)
    Ymulti.append(ygrid)

# aseemble a single block
X = np.concatenate((Xmulti[0], Xmulti[1][1:,:], Xmulti[2][1:,:]), axis=0)
Y = np.concatenate((Ymulti[0], Ymulti[1][1:,:], Ymulti[2][1:,:]), axis=0)

NI, NJ = X.shape

X3d = np.zeros((NI, NJ, NK))
Y3d = np.zeros((NI, NJ, NK))
Z3d = np.zeros((NI, NJ, NK))
for k in range(NK):
    X3d[:,:,k] = X
    Y3d[:,:,k] = Y
    Z3d[:,:,k] = SPAN*k/(NK-1)

# Create a 3D scatter plots
mesh = pv.StructuredGrid(X3d, Y3d, Z3d)
plotter = pv.Plotter()
plotter.add_mesh(mesh, cmap='viridis', show_edges=True)
plotter.show_axes()


# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i_%02i.csv' %(NI, NJ, NK), 'w') as file:
    file.write(f"NDIMENSIONS=3\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK={NK}\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                file.write(f"{X3d[i,j,k]:.6f},{Y3d[i,j,k]:.6f},{Z3d[i,j,k]:.6f}\n")
    print(f"Grid saved to {OUTPUT_FOLDER}/grid_%02i_%02i_%02i.csv" %(NI, NJ, NK))


plotter.show()
plt.show()







