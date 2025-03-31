import numpy as np
from scipy import integrate
import time

# Parameters
NUM_SAMPLES = 100
X_MIN = 1e-3
X_MAX = np.exp(1 / np.e) - X_MIN
A, B = 0, np.pi
CHEB_N = 128

# Storage
times = {"romberg": [], "quad": [], "clenshaw": []}
errors = {"romberg": [], "quad": [], "clenshaw": []}


# Clenshaw-Curtis Quadrature Implementation
def clenshaw_curtis_integrate(f, a, b, N=128):
    k = np.arange(0, N + 1)
    x = np.cos(np.pi * k / N)  # Chebyshev nodes on [-1, 1]
    x_mapped = 0.5 * (b - a) * x + 0.5 * (b + a)

    fx = f(x_mapped)

    # Fast cosine transform weights (barycentric weights)
    c = np.zeros(N + 1)
    c[0] = c[-1] = 2
    c[1:-1] = 1
    w = np.real(np.fft.irfft(np.concatenate([c, c[-2:0:-1]])))[: N + 1]
    w *= 2 / N

    return 0.5 * (b - a) * np.dot(w, fx)


# Define integrand
def make_f(x_val):
    def f(t):
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.asarray(t)
            denom = t
            numer = np.log(x_val) * np.sin(t) * np.exp(t / np.tan(t))
            bad_mask = (denom == 0) | ~np.isfinite(numer)
            argument = 1 - (numer / denom)
            # Ensure log argument > 0
            argument = np.where((argument <= 0) | bad_mask, np.nan, argument)
            return np.log(argument)

    return f


# Sampling loop
x_vals = np.random.uniform(X_MIN, X_MAX, size=NUM_SAMPLES)

for x in x_vals:
    f = make_f(x)

    # Reference via QUAD with tight tolerance
    try:
        ref_val, _ = integrate.quad(f, A, B, epsabs=1e-12, epsrel=1e-12)
    except Exception:
        continue

    # Romberg
    try:
        t0 = time.time()
        romb_val = integrate.romberg(f, A, B, divmax=10, show=False)
        t1 = time.time()
        times["romberg"].append(t1 - t0)
        errors["romberg"].append(abs(romb_val - ref_val))
    except Exception:
        continue

    # Quad
    try:
        t0 = time.time()
        quad_val, _ = integrate.quad(f, A, B)
        t1 = time.time()
        times["quad"].append(t1 - t0)
        errors["quad"].append(abs(quad_val - ref_val))
    except Exception:
        continue

    # Clenshaw-Curtis
    try:
        t0 = time.time()
        clenshaw_val = clenshaw_curtis_integrate(f, A, B, N=CHEB_N)
        t1 = time.time()
        times["clenshaw"].append(t1 - t0)
        errors["clenshaw"].append(abs(clenshaw_val - ref_val))
    except Exception:
        continue

# Print summary
print("\n=== RESULTS ===")
for method in times:
    if times[method]:
        print(f"\n{method.upper()}:")
        print(
            f"  Avg Time:  {np.mean(times[method]):.6f} s ± {np.std(times[method]):.6f} s"
        )
        print(
            f"  Avg Error: {np.mean(errors[method]):.2e} ± {np.std(errors[method]):.2e}"
        )
    else:
        print(f"\n{method.upper()}: No successful runs.")
