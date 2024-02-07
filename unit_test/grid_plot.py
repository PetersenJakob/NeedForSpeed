import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("grid.csv", sep=",", header=None)

plt.plot(data[0], 1 + 0 * data[0], ".b", label="Uniform")
plt.plot(data[1], 2 + 0 * data[1], ".r", label="Exponential")
plt.plot(data[2], 3 + 0 * data[2], ".k", label="Hyperbolic")

plt.xlabel("Grid")
plt.legend()

plt.show()
