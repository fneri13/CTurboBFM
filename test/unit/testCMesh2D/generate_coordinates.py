import numpy as np

NI, NJ, NK = 6, 6, 1
x = np.linspace(0,5,NI)
y = np.linspace(0,5,NJ)
z = np.linspace(0,1,NK)

with open("coordinates.csv", "w") as f:
    f.write(f"NDIMENSIONS=2\n")
    f.write(f"NI={NI}\n")
    f.write(f"NJ={NJ}\n")
    f.write(f"NK={NK}\n")
    f.write("x,y,z\n")
    for i in range(NI):
        for j in range(NJ):
            for k in range(NK):
                f.write(f"{x[i]},{y[j]},{z[k]}\n")