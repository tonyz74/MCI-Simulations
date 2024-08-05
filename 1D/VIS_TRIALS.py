import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad
import arviz as az
import sys

ln = list(open(sys.argv[1]).readlines())
cnt = int(ln[0])
a, b = [float(i) for i in ln[2].split()]
vals = [float(i) for i in ln[3:]]


fig, ax = plt.subplots(figsize=(10, 6.2), layout="constrained")
az.plot_dist(
    vals, kind="hist", ax=ax,
    hist_kwargs={"bins": 50, "density": True}
)

ax.set_xlabel("Estimated Integral")
ax.set_ylabel("Relative Frequency")

func = eval(ln[1])
if func is not None:
    true_integ = quad(func, a, b)[0]
    ax.axvline(x=true_integ, color="r", linewidth=2)

mu, stddev = norm.fit(vals)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, stddev)
plt.plot(x, p, 'k', linewidth=2)

calc_std = np.std(vals)
title = "mu = %.10f, var = %.10f" % (mu, stddev ** 2)
plt.title(title)

plt.locator_params(axis='x', nbins=8)
plt.show()
