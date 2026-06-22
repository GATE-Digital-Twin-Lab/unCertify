import numpy as np

from Affine_ArithmeticClassV3 import AffineArray


X = AffineArray.from_intervals([
    (-2, 2),
    (-1, 1),
    (-0.5, 2.5)
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

Y = AffineArray.from_intervals([
    (1, 2),
    (1, 1.5),
    (0.5, 1)
])

y1 = Y[0]
y2 = Y[1]
y3 = Y[2]

div1 = x/y1
print("div1 =", div1)
print("div1 =", div1.interval)

expVal1 = y3.exp(cheb=True)
expVal2 = (-y3).exp(cheb=True)

print("expVal1 =", expVal1)
print("expVal1 =", expVal1.interval)
print("expVal2 =", expVal2)
print("expVal2 =", expVal2.interval)

lnVal1 = y3.log()
print("logVal1 =", lnVal1)
print("logVal1 =", lnVal1.interval)
lnVal1 = y3.log(cheb=True)
print("logVal1 =", lnVal1)
print("logVal1 =", lnVal1.interval)

absVal = s.abs()
print("absVal =", absVal)
print("absVal =", absVal.interval)

absVal = (s).abs(cheb=True)
print("absVal =", absVal)
print("absVal =", absVal.interval)

sqrtVal = y1.sqrt()
print("sqrtVal =", sqrtVal)
print("sqrtVal =", sqrtVal.interval)

sqrtVal = y1.sqrt(cheb=True)
print("sqrtVal =", sqrtVal)
print("sqrtVal =", sqrtVal.interval)

Theta = AffineArray.from_intervals([
    (-np.pi/2, np.pi),
    (0.25*np.pi, 0.8*np.pi),
    (0.4*np.pi, 2.4*np.pi)
])

th1 = Theta[0]
th2 = Theta[1]
th3 = Theta[2]

sineVal = th2.sin()
print("sineVal = ", sineVal)
print("sineVal = ", sineVal.interval)

sineVal = th2.sin(cheb=True)
print("sineVal = ", sineVal)
print("sineVal = ", sineVal.interval)

cosVal = th1.cos()
print("cosVal = ", cosVal)
print("cosVal = ", cosVal.interval)

cosVal = th1.cos(cheb=True)
print("cosVal = ", cosVal)
print("cosVal = ", cosVal.interval)

# tanVal = th3.tan()
# print("tanVal = ", tanVal)
# print("tanVal = ", tanVal.interval)

# tanVal = th3.tan(cheb=True)
# print("tanVal = ", tanVal)
# print("tanVal = ", tanVal.interval)

funVal = (th3.log()).exp()
print(th3)
print("funVal = ", funVal)
print("funVal = ", funVal.interval)

Theta = AffineArray.from_intervals([
    (-np.pi/2, np.pi/2),
    (0.25*np.pi, 0.8*np.pi),
    (0.4*np.pi, 2.4*np.pi)
])

th1 = Theta[0]
th2 = Theta[1]
th3 = Theta[2]

cotVal = th2.cotan()
print("cotVal = ", cotVal)
print("cotVal = ", cotVal.interval)

cotVal = th2.cotan(cheb=True)
print("cotVal = ", cotVal)
print("cotVal = ", cotVal.interval)

sinhVal = th1.sinh()
print("sinhVal = ", sinhVal)
print("sinhVal = ", sinhVal.interval)

sinhVal = th1.sinh(cheb=True)
print("sinhVal = ", sinhVal)
print("sinhVal = ", sinhVal.interval)

coshVal = th2.cosh()
print("coshVal = ", coshVal)
print("coshVal = ", coshVal.interval)

coshVal = th2.cosh(cheb=True)
print("coshVal = ", coshVal)
print("coshVal = ", coshVal.interval)

tanhVal = th1.tanh()
print("tanhVal = ", tanhVal)
print("tanhVal = ", tanhVal.interval)

tanhVal = th1.tanh(cheb=True)
print("tanhVal = ", tanhVal)
print("tanhVal = ", tanhVal.interval)

Theta = AffineArray.from_intervals([
    (-1.0, 1.0),
    (0.0, 0.8),
    (-0.4, 0.25)
])

th1 = Theta[0]
th2 = Theta[1]
th3 = Theta[2]

arcsinVal = th2.arcsin()
print("arcsinVal = ", arcsinVal)
print("arcsinVal = ", arcsinVal.interval)

arcsinVal = th2.arcsin(cheb=True)
print("arcsinVal = ", arcsinVal)
print("arcsinVal = ", arcsinVal.interval)

arccosVal = th3.arccos()
print("arccosVal = ", arccosVal)
print("arccosVal = ", arccosVal.interval)

arccosVal = th3.arccos(cheb=True)
print("arccosVal = ", arccosVal)
print("arccosVal = ", arccosVal.interval)


arctanVal = th1.arctan()
print("arctanVal = ", arctanVal)
print("arctanVal = ", arctanVal.interval)

arctanVal = th1.arctan(cheb=True)
print("arctanVal = ", arctanVal)
print("arctanVal = ", arctanVal.interval)

Theta = AffineArray.from_intervals([
    (-1.0, 1.0),
    (-0.2, -0.8),
    (-0.4, 0.2)
])

th1 = Theta[0]
th2 = Theta[1]
th3 = Theta[2]

powVal = th3.pow(-3)
print("powVal = ", powVal)
print("powVal = ", powVal.interval)

powVal = th3.pow(-3, cheb=True)
print("powVal = ", powVal)
print("powVal = ", powVal.interval)





