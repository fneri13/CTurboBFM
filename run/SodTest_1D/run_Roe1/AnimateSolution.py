import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

folder = 'Volume_CSV'
files = os.listdir(folder)
files = sorted([folder + '/' + file for file in files if file.endswith('.csv')])


for file in files:
    df = pd.read_csv(file, skiprows=3)
    plt.clf()
    plt.plot(df['x'], df['Density'])
    plt.xlabel('x')
    plt.ylabel('P')
    plt.xlim([0.1, 0.9])
    plt.pause(0.1)
    

plt.show()