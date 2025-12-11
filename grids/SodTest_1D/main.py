import numpy as np
import matplotlib.pyplot as plt
import os

# input parameters
NX = 1500
L = 1
X_INTERFACE = 0.5
RHO_L, RHO_R = 1.0, 0.125
UX_L , UX_R = 0.0, 0.0
UY_L , UY_R = 0.0, 0.0
UZ_L , UZ_R = 0.0, 0.0
T_L, T_R = 1.0/(RHO_L* 287), 0.1/(RHO_R* 287)


# generate grid and initial conditions
x = np.linspace(0, L, NX)
rho = np.where(x < X_INTERFACE, RHO_L, RHO_R)
ux = np.where(x < X_INTERFACE, UX_L, UX_R)
uy = np.where(x < X_INTERFACE, UY_L, UY_R)
uz = np.where(x < X_INTERFACE, UZ_L, UZ_R)
T = np.where(x < X_INTERFACE, T_L, T_R)
y = np.zeros_like(x)
z = np.zeros_like(x)




os.makedirs('Grid', exist_ok=True)

# save the restart data
with open('Grid' + '/grid_%02i_RestartData.csv' %(NX), 'w') as file:
    file.write(f"NI={NX}\n")
    file.write(f"NJ=1\n")
    file.write(f"NK=1\n")
    file.write("x,y,z,Density,Velocity X,Velocity Y,Velocity Z,Temperature\n")
    for i in range(NX):
        for j in range(1):
            for k in range(1):
                file.write(f"{x[i]:.6f},{y[i]:.6f},{z[i]:.6f},{rho[i]:.6f},{ux[i]:.6f},{uy[i]:.6f},{uz[i]:.6f},{T[i]:.6f}\n")


# save only the grid
with open('Grid' + '/grid_%02i_OnlyCoords.csv' %(NX), 'w') as file:
    file.write(f"NDIMENSIONS=1\n")
    file.write(f"NI={NX}\n")
    file.write(f"NJ=1\n")
    file.write(f"NK=1\n")
    file.write("x,y,z\n")
    for i in range(NX):
        for j in range(1):
            for k in range(1):
                file.write(f"{x[i]:.6f},{y[i]:.6f},{z[i]:.6f}\n")
                
print("Grid and initial conditions generated successfully.")