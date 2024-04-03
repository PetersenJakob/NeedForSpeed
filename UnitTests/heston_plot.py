import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("heston.csv", sep=",", header=None)

plt.plot(data[2], data[3], "-b")
plt.plot(data[2], data[4], "-r")
plt.plot(data[2], data[5], "-g")
plt.xlabel("Strike")
plt.ylabel("Option price diff")
plt.legend()
plt.show()


plt.plot(data[2], data[6], "-b")
plt.plot(data[2], data[7], "-r")
plt.xlabel("Strike")
plt.ylabel("Implied vol")
plt.legend()
plt.show()
