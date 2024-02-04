import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("slr.txt", sep=",", header=None)

print(data)

plt.plot(data[1], data[3], "-g")
plt.plot(data[1], data[5], "-r")

coef_max = np.polyfit(data[1], data[3],1)
print(coef_max)
poly1d_fn_max = np.poly1d(coef_max)
plt.plot(data[1], poly1d_fn_max(data[1]), '-.k')

coef_l2 = np.polyfit(data[1], data[5],1)
print(coef_l2)
poly1d_fn_l2 = np.poly1d(coef_l2)
plt.plot(data[1], poly1d_fn_l2(data[1]), '-.k')

plt.xlabel("log(dx)")
plt.ylabel("log(error)")

plt.show()
