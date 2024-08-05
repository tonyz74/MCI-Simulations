import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from scipy.integrate import quad


fname = argv[1]
file = open(fname, "r")

lines = file.readlines()
steps = int(lines[0])
lb, rb = [float(i) for i in lines[2].split(" ")]
data = np.array([float(i) for i in lines[3:]])


x = np.linspace(lb, rb, steps)
fig, ax = plt.subplots()
ax.plot(x, data, "o", color="r")

ax.set_xlabel("x")
ax.set_ylabel("y")

func = eval(lines[1])
if func is not None:
    func = np.vectorize(func)
    true_integ = quad(func, lb, rb)[0]
    ax.plot(x, func(x), color="b")
    ax.plot(x, func(x) / true_integ, color="g")
    ax.plot(x, data * true_integ, "+", color="aqua")

plt.show()
