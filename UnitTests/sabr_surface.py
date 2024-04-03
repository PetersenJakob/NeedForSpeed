import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.ticker import LinearLocator

s = pd.read_csv("sabr_surface_s.csv", sep=",", header=None)
v = pd.read_csv("sabr_surface_v.csv", sep=",", header=None)
data = pd.read_csv("sabr_surface.csv", sep=",", header=None)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

s, v = np.meshgrid(s, v)

# Plot the surface.
surf = ax.plot_surface(s, v, data.transpose(), cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_xlabel('S')
ax.set_ylabel('v')
ax.set_zlabel('Price')

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
