import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv("black_scholes_call.csv", sep=",", header=None)

plt.plot(data[0], data[1], "-b", label="Numerical")
plt.plot(data[0], data[2], "-r", label="Analytical")
plt.plot(data[0], data[3], "--k", label="Payoff")

plt.xlabel("Stock price")
plt.ylabel("Call price")
plt.legend()

plt.show()
