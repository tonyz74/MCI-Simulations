import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
file = open(fname, "r")
lines = file.readlines()

steps = int(lines[0])
dims = int(lines[1])

bounds = []
i = 2
for k in range(dims):
    a, b = [float(j) for j in lines[i].split()]
    bounds.append((a, b))
    i += 1

vals = []
for k in range(dims):
    now = []
    for _ in range(steps):
        now.append(float(lines[i]))
        i += 1
    vals.append(np.array(now))

fig, axs = plt.subplots(dims, figsize=(10.0, 8.0))
fig.suptitle("PDFs in First Dimension")

assert(dims == 1)
if dims == 1:
    axs = [axs]

for k in range(dims):
    x = np.linspace(bounds[k][0], bounds[k][1], steps)

    axs[k].set_xlabel(f"x{k}")
    axs[k].set_ylabel("p")
    axs[k].set_xlim(bounds[k][0], bounds[k][1])
    axs[k].plot(x, vals[k], "o")

plt.show()
