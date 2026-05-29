import numpy as np

from Affine_ArithmeticClassV3 import AffineArray


X = AffineArray.from_intervals([
    (-2, 2),
    (-1, 1),
    (-1, 1)
])

x = X[0]
r = X[1]
s = X[2]

print(x)
print(r)
print(s)

print(type(x))

add1 = x + 3
add2 = x + r
add3 = -5 + s
add4 = s - r

addx1 = x + add1
addx2 = x - add2 + s

print(add1)
print(add2)
print(add3)
print(add4)
print(addx2)


prod1 = s*r
print("prod1 = ", prod1)

