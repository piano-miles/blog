import math

x = 0.04
L = math.log(x)

z1 = (
    1 + (math.e - 1) * math.exp(math.e * 0.5) * math.sqrt(math.exp(-math.e) - x)
) / math.e
z2 = math.exp(0.15 * math.e - 1) * x**0.15
z3 = (1 - math.sqrt(1 - math.exp(2 * math.e) * x * x)) / math.e

for i in range(10):
    print(z1, z2, z3)
    Lxz1 = L * x**z1
    Lxz2 = L * x**z2
    Lxz3 = L * x**z3

    z1 *= 1 + (math.log(z1) - Lxz1) / (z1 * L * Lxz1 - 1)
    z2 *= 1 + (math.log(z2) - Lxz2) / (z2 * L * Lxz2 - 1)
    z3 *= 1 + (math.log(z3) - Lxz3) / (z3 * L * Lxz3 - 1)
