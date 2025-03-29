import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw


# Define the function y(x) using the -1 branch
def y_of_x(x):
    ln_x = np.log(x)
    return -1 / ln_x * lambertw(-1 / ln_x, k=-1).real  # Use W_{-1}


# Define the domain: (e^(-e), e^(1/e))
x_min = np.exp(-np.e)
x_max = np.exp(1 / np.e)
x = np.linspace(x_min, x_max, 1000)

# Evaluate y(x)
y = y_of_x(x)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(
    x, y, label=r"$y(x) = -\frac{1}{\ln x} \cdot W_{-1}\left(-\frac{1}{\ln x} \right)$"
)
plt.title("Correct Real-Valued Solution Using $W_{-1}$ Branch")
plt.xlabel("$x$")
plt.ylabel("$y(x)$")
plt.grid(True)
plt.legend()
plt.xlim(x_min, x_max)
plt.ylim(bottom=0)
plt.show()
