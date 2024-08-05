import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad
import arviz as az
import sys

ln = list(open(sys.argv[1]).readlines())
cnt = int(ln[0])
vals = [float(i) for i in ln[3:]]


fig, ax = plt.subplots(figsize=(10, 6.2), layout="constrained")
fig.suptitle("Distribution of Integration Results")

az.plot_dist(
    vals, kind="hist", ax=ax,
    hist_kwargs={"bins": 50, "density": True}
)

ax.set_xlabel("Estimated Integral")
ax.set_ylabel("Relative Frequency")

true_val = 4.87304976542
# true_val = 0.927295218002
ax.axvline(x=true_val, linewidth=4, color="r")

mu, stddev = norm.fit(vals)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, stddev)
plt.plot(x, p, 'k', linewidth=2)

ferr = abs(true_val - mu) / true_val

calc_std = np.std(vals)
fit = (
    "Normal Fit Results:\n" +
    "Mean: %.10f\n" +
    "Variance: %.10f\n" +
    "Standard Deviation: %.10f\n" +
    "Percentage Error: %.10f%%"
) % (mu, stddev ** 2, stddev, ferr * 100.0)

ax.text(0.05, 0.8, fit, transform=ax.transAxes)

plt.locator_params(axis='x', nbins=8)
plt.show()
