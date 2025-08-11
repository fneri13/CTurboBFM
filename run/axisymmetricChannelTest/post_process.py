import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from Utils.styles import *
import os 


intputFile = 'volumeData/results.csv'
os.makedirs('processedData', exist_ok=True)
outputFile = 'processedData/axialProfile.pkl'


# Read the first three lines to extract grid sizes
with open(intputFile, 'r') as f:
    ni = int(f.readline().strip().split('=')[1])
    nj = int(f.readline().strip().split('=')[1])
    nk = int(f.readline().strip().split('=')[1])

# Read the rest of the CSV data into a DataFrame
df = pd.read_csv(intputFile, skiprows=3)
data = {}

for field in df.columns:
    data[field] = np.reshape(df[field].values, (ni, nj, nk))


radius = data['y'][ni//2, :, 0]
utheta = data['Velocity Z'][ni//2, :, 0]
density = data['Density'][ni//2, :, 0]
pressure = data['Pressure'][ni//2, :, 0]

plt.figure()
plt.plot(radius, utheta, label='Azimuthal Velocity', color='blue')
plt.xlabel('Radius')
plt.ylabel('Azimuthal Velocity')



    
plt.show()
        
        


