#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
# from styles import *
import argparse

os.makedirs("pictures", exist_ok=True)

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Plot mesh quality statistics.")
parser.add_argument("data_folder", type=str, help="Path of the folder containing the mesh quality data")
args = parser.parse_args()



folder = "./Mesh_Quality/"
skewness = np.loadtxt(folder + "skewness.csv", delimiter=",")
orthogonality = np.loadtxt(folder + "orthogonality.csv", delimiter=",")
aratio = np.loadtxt(folder + "aspect_ratio.csv", delimiter=",")


def plot_normalized_hist(data, color, title):
    N = len(data)
    weights = np.ones(N) / N  # Each entry contributes 1/N
    plt.hist(data, bins=30, weights=weights, color=color, edgecolor='black', alpha=0.7)
    plt.title(f"{title} (N={N})")
    plt.xlabel("Value")
    # plt.ylabel("Probability")
    plt.grid(alpha=0.3)

plt.figure(figsize=(15, 4))

plt.subplot(1, 3, 1)
plot_normalized_hist(skewness, 'steelblue', 'Skewness')

plt.subplot(1, 3, 2)
plot_normalized_hist(orthogonality, 'darkorange', 'Orthogonality')

plt.subplot(1, 3, 3)
plot_normalized_hist(aratio, 'seagreen', 'Aspect Ratio')

plt.savefig("pictures/mesh_quality.pdf", bbox_inches='tight')
plt.show()