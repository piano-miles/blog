import math

z = 3.0
print(f"z_0: {z:.11f}")

for i in range(4):
    z += math.sin(z)
    print(f"z_{i+1}: {z:.11f}")

print(f"π:   {math.pi:.11f}")
