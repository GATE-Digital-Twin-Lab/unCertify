import numpy as np

def scale(x, lo, hi):
    return lo + (hi - lo)*x

#---------------------------------------------------
#  Numpy functions
#---------------------------------------------------


def Branin_numpy(x):

    a = 1.0
    b = 5.1 / (4.0 * np.pi**2)
    c = 5.0 / np.pi
    r = 6.0
    s = 10.0
    t = 1.0 / (8.0 * np.pi)

    x1 = x[0]
    x2 = x[1]

    # Map normalized variables:
    #
    # x1 in [0,1] -> [-5,10]
    # x2 in [0,1] -> [0,15]
    x1 = 15.0 * x1 - 5.0
    x2 = 15.0 * x2

    funVal = (
        a * (x2 - b * x1**2 + c * x1 - r)**2
        + s * (1.0 - t) * np.cos(x1)
        + s
        + 5.0 * x1
    )

    return funVal


def Ackley_numpy(x, d=2):

    a = 20.0
    b = 0.2
    c = 2 * np.pi

    funVal = a + np.exp(1)

    sum1 = 0.0
    sum2 = 0.0

    for i in range(d):

        xi = scale(x[i], -32.768, 32.768)

        sum1 += xi**2
        sum2 += np.cos(c * xi)

    term1 = np.exp(
        -b * np.sqrt(sum1 / d)
    )

    term2 = np.exp(
        sum2 / d
    )

    funVal += -a * term1
    funVal += -term2

    return funVal


def Eggholder_numpy(x):

    x1 = scale(x[0], -512, 512)
    x2 = scale(x[1], -512, 512)

    funVal = (
        -(x2 + 47)
        * np.sin(
            np.sqrt(
                np.abs(
                    x2 + 47 + x1 / 2
                )
            )
        )

        - x1
        * np.sin(
            np.sqrt(
                np.abs(
                    x1 - (x2 + 47)
                )
            )
        )
    )

    return funVal


#---------------------------------------------------
#  AA functions
#---------------------------------------------------


def test_BraninAA(x, cheb=False):

    a = 1
    b = 5.1 / (4 * np.pi**2)
    c = 5  / np.pi
    r = 6
    s = 10
    t = 1 / (8 * np.pi)

    x1 = x[0]
    x2 = x[1]

    x1 = 15 * x1 - 5
    x2 = 15 * x2

    funVal = (
        a * (x2 - b * x1.pow(2, cheb=cheb) + c * x1 - r).pow(2, cheb=cheb) 
        + s * (1 - t) * (x1).cos(cheb=cheb) + s + 5*x1
    )

    return funVal



def test_AckleyAA(x, d=2, cheb=False):

    a = 20.0
    b = 0.2
    c = 2*np.pi

    funVal = a + np.exp(1)

    sum1 = 0.0
    sum2 = 0.0

    for i in range(d):

        # xi = 2*32.768*x[i] - 32.768
        xi = scale(x[i], -32.768, 32.768)

        sum1 += xi.pow(2, cheb=False)
        sum2 += (c * xi).cos(cheb=cheb)

    term1 = (-b * (sum1/d).sqrt(cheb=cheb)).exp(cheb=cheb)
    term2 = ((sum2/d)).exp(cheb=cheb)

    funVal += -a * term1
    funVal += -term2

    return funVal

def test_EggholderAA(x, cheb=False):

    # x1 = 512 * x[0] - 512
    # x2 = 512 * x[1] - 512

    x1 = scale(x[0], -512, 512)
    x2 = scale(x[1], -512, 512)

    funVal = (
        -(x2 + 47)
        * (
            (x2 + 47 + x1/2)
            .abs(cheb=cheb)
            .sqrt(cheb=cheb)
        ).sin(cheb=cheb)

        - x1
        * (
            (x1 - (x2 + 47))
            .abs(cheb=cheb)
            .sqrt(cheb=cheb)
        ).sin(cheb=cheb)
    )

    return funVal