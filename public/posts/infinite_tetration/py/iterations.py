import math

BRANCH_POINT: float = math.exp(-math.e)
MAX_VALUE: float = math.exp(1 / math.e)


def compute_infinite_tetration(
    x: float,
    MAX_ITER: int = 128,
    tolerance: float = 1e-9,
    verbose: bool = False,
) -> tuple[float, ...]:
    """Computes the infinite tetration of non-negative real numbers using the Newton-Raphson method.

    Args:
        x (float): non-negative float to infinitely tetrate.
        MAX_ITER (int, optional): Number of maximum iterations allowed. Defaults to 32.
        tolerance (float, optional): Convergence tolerance. Defaults to 1e-9.

    Returns:
        tuple[float, ...]: The converged value, or the branches of the solution.
    """

    if x < 0.0:
        if math.isclose(x, -1.0):
            return (-1.0,)
        return (float("nan"),)
    elif math.isclose(x, 0.0):
        return (0.0, 1.0)
    elif math.isclose(x, 1.0):
        return (1.0,)

    elif x < BRANCH_POINT:
        z1 = (
            1 + (math.e - 1) * math.exp(math.e * 0.5) * math.sqrt(math.exp(-math.e) - x)
        ) / math.e
        z2 = math.exp(0.15 * math.e - 1) * x**0.15
        z3 = (1 - math.sqrt(1 - math.exp(2 * math.e) * x * x)) / math.e
        L = math.log(x)

        # Flags to check for convergence
        c1 = c2 = c3 = False

        # These will be L * x^z for each branch.
        Lxz1 = Lxz2 = Lxz3 = 1

        # Iterate for x^x^y = y until convergence
        for i in range(MAX_ITER):
            # Keeping track of previous z's to check for convergence
            zp1, zp2, zp3 = z1, z2, z3
            if not c1:
                Lxz1 = L * x**z1
            if not c2:
                Lxz2 = L * x**z2
            if not c3:
                Lxz3 = L * x**z3

            # These cases should NEVER occur. It is here just in case.
            if z1 < 0:
                print(f"WARNING: iteration {i + 1} encountered negative z1: {z1}")
                z1 = x**z1
            if z2 < 0:
                print(f"WARNING: iteration {i + 1} encountered negative z1: {z2}")
                z2 = x**z2
            if z3 < 0:
                print(f"WARNING: iteration {i + 1} encountered negative z1: {z3}")
                z3 = x**z3

            try:
                # Newton-Raphson iterations
                if not c1:
                    z1 *= 1 + (math.log(z1) - Lxz1) / (z1 * L * Lxz1 - 1)
                if not c2:
                    z2 *= 1 + (math.log(z2) - Lxz2) / (z2 * L * Lxz2 - 1)
                if not c3:
                    z3 *= 1 + (math.log(z3) - Lxz3) / (z3 * L * Lxz3 - 1)

            except ValueError as Err:
                print(f"Iteration {i + 1} failed with error: {Err}")
                print(f"z1 = {z1}, Lxz1 = {Lxz1}")
                print(f"z1 = {z2}, Lxz1 = {Lxz2}")
                print(f"z1 = {z3}, Lxz1 = {Lxz3}")
                return (z1, z2, z3)
            except ZeroDivisionError as Err:
                print(f"Iteration {i + 1} failed with error: {Err}")
                print(f"z1 = {z1}, Lxz1 = {Lxz1}")
                print(f"z1 = {z2}, Lxz1 = {Lxz2}")
                print(f"z1 = {z3}, Lxz1 = {Lxz3}")
                return (z1, z2, z3)

            # Check for convergence
            if not c1:
                c1 = abs(z1 - zp1) < tolerance * max(abs(z1), 1.0)
            if not c2:
                c2 = abs(z2 - zp2) < tolerance * max(abs(z2), 1.0)
            if not c3:
                c3 = abs(z3 - zp3) < tolerance * max(abs(z3), 1.0)
            if c1 and c2 and c3:
                if verbose:
                    print(f"Converged in {i} iterations.")
                return (z1, z2, z3)

        print(
            f"WARNING: iteration {i + 1} failed to converge after {MAX_ITER} iterations."
        )
        print(f"x: {x}")
        return (z1, z2, z3)

    elif x < MAX_VALUE:
        z = ((1 + 2.051 * x) * x - 6.037) / ((16.183 - 4.006 * x) * x - 15.134)
        L = math.log(x)

        # Iterate for x^y = y until convergence
        for i in range(MAX_ITER):
            zp = z  # previous z

            if z < 0:
                # This should never occur. It is here just in case.
                print(f"WARNING: iteration {i + 1} encountered negative z: {z}")
                z = x**z

            try:
                # Newton-Raphson iteration
                z *= (math.log(z) - 1) / (L * z - 1)

            except ValueError as Err:
                print(f"Iteration {i + 1} failed with ValueError: {Err}")
                print(f"z = {z}")
                return (z,)
            except ZeroDivisionError as Err:
                print(f"Iteration {i + 1} failed with ZeroDivisionError: {Err}")
                print(f"z = {z}")
                return (z,)

            if abs(z - zp) < tolerance * max(abs(z), 1.0):
                if verbose:
                    print(f"Converged in {i} iterations.")
                return (z,)

        print(
            f"WARNING: iteration {i + 1} failed to converge after {MAX_ITER} iterations."
        )
        print(f"x: {x}")
        return (z,)

    elif math.isclose(x, MAX_VALUE):
        return (math.e,)
    else:
        return (float("inf"),)
