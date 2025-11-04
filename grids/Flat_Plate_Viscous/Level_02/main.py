import numpy as np
import matplotlib.pyplot as plt
import os
import pickle



"""
Generate a 2D rectangular geometry, that will be used as verification
"""
OUTPUT_FOLDER = 'Grid'
LX = 0.2
H = 0.02
NX = 65
NY = 113
STRETCH = 1.083317311 # from hirsch

y = np.zeros(NY)
y[1] = 0.2E-5
for j in range(2, NY):
    dy = STRETCH * (y[j-1] - y[j-2])
    y[j] = y[j-1] + dy

x = np.zeros(NX)
x[1] = 0.1E-3
for i in range(2, NX):
    dx = STRETCH * (x[i-1] - x[i-2])
    x[i] = x[i-1] + dx


X, Y = np.meshgrid(x, y, indexing='ij')
NI, NJ = X.shape


plt.figure()
for i in range(NX):
    plt.plot(X[i, :], Y[i, :], 'k', lw=0.5)
for j in range(NY):
    plt.plot(X[:, j], Y[:, j], 'k', lw=0.5)

plt.xlabel(r'$x \ \rm{[m]}$')
plt.ylabel(r'$y \ \rm{[m]}$')



grid = {'X': X, 'Y': Y}
# Create output directory
if os.path.exists(OUTPUT_FOLDER):
    print('Output Folder already present')
else:
    os.mkdir(OUTPUT_FOLDER)
with open(OUTPUT_FOLDER + '/grid_%02i_%02i.csv' %(NX, NY), 'w') as file:
    file.write(f"NDIMENSIONS=2\n")
    file.write(f"NI={NI}\n")
    file.write(f"NJ={NJ}\n")
    file.write(f"NK=1\n")
    file.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(1):
                file.write(f"{X[i,j]:.6f},{Y[i,j]:.6f},{0:.6f}\n")

plt.show()







