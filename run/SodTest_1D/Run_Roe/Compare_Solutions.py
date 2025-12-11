import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from Utils.styles import *


files = [
    'Volume_CSV/results_250.csv',
    'Volume_CSV/results_1500.csv',
    'Volume_CSV/results_250_2nd.csv',
]

labels = [
    'N:250',
    'N:1500',
    'N:250 + MUSCL',
]

fields = [
    'Density',
    'Velocity X',
    'Pressure',
]

for field in fields:
    plt.figure()
    for file in files:
        df = pd.read_csv(file, skiprows=3)
        plt.plot(df['x'], df[field], label=labels[files.index(file)])
    plt.xlabel(r'$x$ [-]')
    plt.ylabel(r'%s [-]' % field)
    plt.grid(alpha=0.2)
    plt.legend()
    plt.savefig('pictures/%s.pdf' % field, bbox_inches='tight')
plt.show()