import numpy as np

NI, NJ, NK = 100, 200, 1
x = np.linspace(1,2,NI)
y = np.linspace(1,2,NJ)
z = np.linspace(1,2,NK)

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