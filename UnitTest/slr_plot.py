import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("slr.csv", sep=",", header=None)

# Maximum norm.
plt.plot(data[1], data[3], ".b")
# l1 vector norm.
plt.plot(data[1], data[5], ".r")
# l2 vector norm.
plt.plot(data[1], data[7], ".g")
# L1 function norm.
plt.plot(data[1], data[9], ".k")
# L2 function norm.
plt.plot(data[1], data[11], ".y")

# Linear regression of maximum norm data.
coef_max = np.polyfit(data[1], data[3],1)
print("Maximum norm: ", coef_max)
poly1d_fn_max = np.poly1d(coef_max)
plt.plot(data[1], poly1d_fn_max(data[1]), '-b', label="Max")

# Linear regression of l1 vector norm.
coef_l1_vec = np.polyfit(data[1], data[5],1)
print("l1 vector norm: ", coef_l1_vec)
poly1d_fn_l1_vec = np.poly1d(coef_l1_vec)
plt.plot(data[1], poly1d_fn_l1_vec(data[1]), '-r', label="l1 vec")

# Linear regression of l2 vector norm.
coef_l2_vec = np.polyfit(data[1], data[7],1)
print("l2 vector norm: ", coef_l2_vec)
poly1d_fn_l2_vec = np.poly1d(coef_l2_vec)
plt.plot(data[1], poly1d_fn_l2_vec(data[1]), '-g', label="l2 vec")

# Linear regression of L1 function norm.
coef_l1_func = np.polyfit(data[1], data[9],1)
print("L1 function norm: ", coef_l1_func)
poly1d_fn_l1_func = np.poly1d(coef_l1_func)
plt.plot(data[1], poly1d_fn_l1_func(data[1]), '-k', label="L1 func")

# Linear regression of L2 function norm.
coef_l2_func = np.polyfit(data[1], data[11],1)
print("L2 function norm: ", coef_l2_func)
poly1d_fn_l2_func = np.poly1d(coef_l2_func)
plt.plot(data[1], poly1d_fn_l2_func(data[1]), '-y', label="L2 func")

plt.xlabel("log(step size)")
plt.ylabel("log(error)")
plt.legend()

plt.show()
