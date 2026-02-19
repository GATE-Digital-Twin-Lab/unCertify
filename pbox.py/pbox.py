# -*- coding: utf-8 -*-
# @author: Scott Ferson

"""
pbox: Interval and probability-box (p-box) library
-----------------------------------------------------------------

pbox provides two core uncertainty objects:

  • Interval  – an uncertain number known only to lie between bounds
  • Pbox      – a distribution or family of distributions represented
                 by lower and upper cumulative distribution functions

The library supports:
  - basic interval arithmetic
  - elementary p-box arithmetic assuming independence
  - confidence boxes (c-boxes)
  - logical and set operations
  - plotting intervals and p-boxes, and automatic plotting for p-boxes

Quick Start
-----------

# Intervals

x = I(2, 5)        # [2, 5]
x.left()           # 2
x.right()          # 5

y = I(4, 3)        # automatically corrected to [3, 4]
z = I("7.5")       # interpreted via significant digits as [7.45, 7.55]

v = x + y          # [5.0, 9.0]
w = x * y * z      # [44.7, 151.0]

# P-boxes and c-boxes

u = U(0, 1)        # uniform(0,1)
n = N(0, 1)        # normal(0,1)
b = B(2, 3)        # beta(2,3)

m = N(I(5,6), 1)   # p-box of normals with mean in [5,6]
p = KN(2, 10)      # c-box for 2 successes in 10 trials
o = MMM(0, 10, 1)  # p-box of all distributions with mean 1 over range [0,10]

# General uncertainty arithmetic and logic

s = u + n + Interval(0,1) # Pbox(range=[-3.09, 5.09], mean=[0.5, 1.5])
e = env(u, n)             # Pbox(range=[-3.09, 3.09], mean=[0.0, 0.5])
a = KN(2,10) * KN(2,100)  # Pbox(range=[2.18e-06, 0.077], mean=[0.0036, 0.0081])

# Plotting

plot(x)            # interval
plot(x,form='e')   # ellipse, 't' for triangle
plot(a)            # p-box
plot(e,fmt='b:')   # blue with dotted lines

# MIT License (c) 2026 Scott Ferson
"""

"""
RECONCILATION OF THE PYTHON AND R VERSIONS OF THE PBOX LIBRARY

The libraries pbox.r and pba.r define gamma() and gamma2(), but their meanings 
are exactly swapped in this pbox.py library. I am currently thinking that the R 
code should be updated to be consistent with the Python versions.

The massreassignexamples() function in the R version is a bit more elaborate.
It has an extra line of functions, including the Schmitt trigger and barbell.

"""


"""
FUNCTIONS AVAILABLE IN THIS LIBRARY

* Natural inputs

The library is designed to be useful for analysts who aren’t sure what a normal 
distribution is, and have never heard of a p-box. Although it employs state-of-
the-art probabilistic and non-Laplacian uncertainty, it can make use of very 
coarse inputs, including those expressed verbally such as 

    “between 50 and 100”
    “0 out of 10”
    “less than 25%”
    “about 9.3” 
    “268 count”
    “1 in 10,000” (expert opinion)
    
The library can make use of sample data too, even very small data sets that 
other analyses cannot handle.  And it can also use qualitative information, 
summaries from reports, and design constraints or specifications


The library fashions its inputs from whatever information a user has and assumptions 
they are comfortable making. 

orderof
O
plusminus
plusminuspercent

KN, km

format_sigdigs
format_decimals
sgnumber
sigorder
lastsigfig
about
    
around
almost
elicited
exactly
<<above>>  downto
<<below>>   upto
nearly

It can use the inputs in a wide variety of arithmetic and logical calculations 
and models to reveal what can be known without requiring extra assumptions.

The library can also make use of sample data, including data with significant 
imprecision, censoring, missingness and other kinds of measurement uncertainty.

* Constructors to specify intervals, distributions and p-boxes

You can specify intervals explicitly or implicitly.
There are over seventy named distributions.  You can give their parameters as 
scalars or as intervals (in which case you get a p-box).
There are also distribution p-box constructors that don't require you to decide
on a 

* Estimators that make use of sample data
(MM, ML, ME, RB, CB)

* Arithmetic operations among uncertain numbers
This library supports these binary convolutions, binary nonconvolutional 
aggregations, unary transformations, statistical accessors,

+
-
*
/
**
minimum               # # # # # # # #
maximum               # # # # # # # #

* New operations with uncertain numbers
env
imp (intersection)
least
greatest
mixture               # # # # # # # #

* Standard and new statistics and accessors
mean
sd
var
prob
ci
support
straddles
mignitude
breadth
width
iqr
IQR

* Standard transformations and some new ones
negate
reciprocate
complement
exp
log
sqrt
square
round
floor
ceil
trunc
abs
sign
sigilium
mignitude
sin
cos
tan
atan
asin
acos
asec
acsc
lambertw
fatten
widen

* New mass reassignment functions
above
below
between
lowest
highest
rescale
truncate
constrainedto
on01
censor

* Logical functions and operations
The library offers several logical functions used in fault tree analysis and
logical modeling.
                                              exclusive-
     negation        conjunction  disjunction  disjunction  logical[s]    
     complement      andc         orc         xorc          Correlated   
     cond            andf         orf         xorf          Dependence   
     equivalence     andi         ori         xori          like_INDEPENDENCE   
     imply           andn         orn         xorn          permissible_lucas   
     modusponens[2]  andp         orp         xorp          check_lucas   
     modustollens[2] AND          OR          XOR           lucas_from_frank   

Functions with capitalized names assume independence among the logical inputs.
These include AND, BUFFER, EQUIVALENCE, INHIBIT, MOON (M out of N), NAND, NOR, 
NOT, OR, PAND, SCHMITT, XNOR, and XOR. Other than NOT, BUFFER and SCHMITT, they 
all accept multiple probabilistic inputs, which may be variously characterized 
as Booleans, scalars, intervals, distributions, p-boxes or c-boxes so long as 
they are constrained to the range [0,1]. The value 0 is understood to represent 
False, the value 1 stands for True, and the dunno interval [0,1] represents 
ignorance about the truth value.

<<Some logical functions are not exposed: FMP1R, FMP2R, LMP1R, Sandc, 
normalize_dependence, m_out_of_n_prob, splay>>   

"""

print("Immediate needs:")
print("simple inputs, IP version of count")
print("powerfunction, muth, lomax, cantor, burr, etc.")
print("MEquantiles, MEdiscretemean")
print('mixture')
print("sawinconrad bug"); 
if 0:
    """
    U(1,4)
    U(I(1,2), I(4,5))
    U(I(1,10), I(4,5))       # should be BIGGER than previous, but it's samller
    a=U(I(1,2), I(4,5))
    b=U(I(1,4), I(4,5))
    c=U(I(1,10), I(4,5))
    plot(a,lw=3); plot(b,fmt='b'); plot(c,fmt='c'); plt.show()
    
    a=sawinconrad(I(0,3),I(3,4),10); a    
    b=sawinconrad(I(0,4),I(3,4),10); b     # should be bigger, but it's smaller
    plot(a); plot(b,fmt='b'); #plot(c,fmt='c')
    
    # triangular, trapezoidal, maybe minmax... constructors may be wrong too
    """
print("Consider numpy's ufuncs")
# ...then NumPy’s ufuncs do the right thing with your interval type if you’ve 
# implemented __array_ufunc__.
# Exponentiation sd**0.5 may or may not behave the same depending on your overloads.
# So using np.sqrt everywhere keeps you aligned with NumPy’s ufunc machinery, 
# which is the right long‑term design for interval arithmetic.
print("Loose ends:  FLEX Frechet, MOON for p-boxes, if x isn't sorted, Pbox(x) should maybe return EDF(x)")


import math
import numpy as np
from numbers import Number       
import matplotlib.pyplot as plt
import scipy.stats as sps

# ============================================================
# Global options
# ============================================================

class IvO:
    ordered = True                      # set to FALSE to allow Kaucher objects
    autocorrect = True                # set to FALSE to support Kaucher objects
    digits = None             # set default decimal digits for repr() and str()
    tol = 1e-8       # width below which an interval is interpreted as a scalar 
    asscalar = 0           # if sigdigs don't work, 0=midpoint, 1=left, 2=right
    form = "b"                     # display form: b=box, e=ellipse, t=triangle
    ylab = ""                            # y-axis label when plotting intervals
    suppress_np_warnings = True       # quiets warnings from divby0 and invalid
    # or use  with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
    quieterrors = True       # domain errors, sqrt(negs), log(negs), asin, etc.
    # prolly need an error-behavior option for powers and exponentiations
    why_negative      = 'undefined for negative values'
    why_nonpositive   = 'undefined for nonpositive values' 
    why_outofunitdisk = 'undefined for values outside [-1,+1]' 
    why_outofrange    = 'cannot compute'  # e.g., exp(1000) OverflowError

class PbO:
    steps = 200
    bOt = 0.001 # p-value of lowest quantile for unbounded p-box (unless 1/steps is smaller); set to 0 for infinite tail
    tOp = 0.999 # p-value of highest quantile for unbounded p-box (unless 1-1/steps is larger); set to 1 for infinite tail
    Bzero = 1e-6
    Bone = 1 - 1e-6
    allowNA = True # permit missing values in p
    cumulative = True
    print_digits = 6

    @staticmethod
    def ii(): return np.arange(0, PbO.steps) / PbO.steps

    @staticmethod
    def jj(): return np.arange(1, PbO.steps+1) / PbO.steps

    @staticmethod 
    def iii(): return np.concatenate((np.array([min(PbO.bOt, 1/PbO.steps)]), np.arange(1, PbO.steps) / PbO.steps))

    @staticmethod
    def jjj(): return np.concatenate((np.arange(1, PbO.steps) / PbO.steps, np.array([max(PbO.tOp, 1-1/PbO.steps)])))
      
def setoption(**kwargs):   # use comme ça:  setoption(digits=3)
    for key, value in kwargs.items():
        matched = False
        if hasattr(IvO, key): 
            setattr(IvO, key, value)
            matched = True
        if hasattr(PbO, key):
            setattr(PbO, key, value)
            matched = True
        if not matched: raise KeyError(f"Unknown option: {key}")        
        
# ============================================================
# NA type
# ============================================================

class NAType:
    def __repr__(self): return "NA"           # essentially a directed infinity

NA = NAType()

def is_na(x): return isinstance(x, NAType)

def is_missing(x): return x is None or is_na(x)

"""
Although R has NaN (IEEE corruption), and NULL (absence), and NA (missing), 
Python has only float('nan') (IEEE NaN) and None (absence of value), and no 
native NA or coherent way to indicate missingness.  So we introduce one:
    
   Interval(NA, 5)                     # [-inf, 5]
   Interval(4, NA)                     # [4, inf]
   Interval(NA, NA)                    # [-inf, inf]     # fully vacuous
   Pbox( [1,2,3,NA,5], [4,NA,6,7,8] )  # Pbox( [1,2,3,3,5], [4,6,6,7,8] )  
   
Note that the outward-directed rounding heals the NA to propagates uncertainty.
"""

# ============================================================
# Interval class
# ============================================================

def format_decimals(a, d):           # format_decimals(100/3,3) yields '33.333'
    if d is None: d = IvO.digits
    if d is None: return f"{a}"
    fmt = f"{{:.{d}f}}"
    return f"{fmt.format(a)}"

def format_sigdigs(a, D):               # format_sigdigs(100/3,3) yields '33.3'
    fmt = f"{{:.{D}g}}"
    return f"{fmt.format(a)}"

def sgnumber(user_input: str):     # number ± its significant-digit imprecision
    user_input = user_input.strip().lower()
    tens = '0'
    if 'e' in user_input: mantissa, tens = user_input.split('e', 1)
    else: mantissa = user_input
    if '.' in mantissa: j = len(mantissa.split('.')[1])
    #else: j = len(mantissa.split('0', 1)) - len(mantissa) + 1           
    else: j = len(mantissa.rstrip('0')) - len(mantissa) 
    pm = 10**(-j) * 10**int(tens) / 2
    #print('input:',user_input,', mantissa:',mantissa, ', j:',j, ', tens:',tens, ', pm:',pm)
    return([float(user_input)-pm, float(user_input)+pm])

class Interval:
    def __init__(self, lo, hi=None, auto=True):
        """
        Interval(12,22)                                   # [12.0, 22.0]
        Interval(22,12)                                   # [12.0, 22.0]
        Interval(12)                                      # 12.0
        Interval(12,)                                     # 12.0
        Interval(12,None)                                 # 12.0
        Interval(12,NA)                                   # [12.0, inf]
        Interval(12,np.inf)                               # [12.0, inf]
        Interval(12,float('inf'))                         # [12.0, inf]
        Interval(NA,NA)                                   # [-inf, inf]
        Interval('12')                                    # [11.5, 12.5]
        Interval(Interval(12,13), Interval(21,22))        # [12, 22]
        Interval(Interval(12,50), Interval(9,22))         # [12, 22]
        """
        def _scalarize(x):          # numPy scalar (np.float64, np.int64, etc.)
            if isinstance(x, np.generic): return float(x)
            if isinstance(x, np.ndarray) and x.size == 1: return float(x)
            return x
        lo = _scalarize(lo)
        hi = _scalarize(hi)
        if is_na(lo): lo = float("-inf")
        if is_na(hi): hi = float("inf")
        if isinstance(lo, str) and hi is None:     # significant-digit interval
            lo_val, hi_val = sgnumber(lo)
            self._lo = float(lo_val)
            self._hi = float(hi_val)
        elif hi is None and isinstance(lo, Interval):       # existing interval
            self._lo = float(lo.lo)
            self._hi = float(lo.hi)
        elif hi is None and np.isscalar(lo):              # degenerate interval
            self._lo = float(lo)
            self._hi = float(lo)
        elif hi is None and isinstance(lo, np.ndarray) and lo.size == 2:
            self._lo = float(lo[0])
            self._hi = float(lo[1])
        elif hi is None and hasattr(lo, "__len__") and len(lo) == 2:
            self._lo = float(lo[0])
            self._hi = float(lo[1])
        elif hi is not None:                                  # explicit bounds
            self._lo = float(left(lo))
            self._hi = float(right(hi))
        else: raise ValueError("Bad interval input")
        if self._hi < self._lo and IvO.autocorrect and auto:
            self._lo, self._hi = self._hi, self._lo

    @property
    def lo(self):
        return self._lo

    @property
    def hi(self):
        return self._hi
    
    def setlo(self, new_lo): # 
        """We cannot support direct user changes to an existing interval such as with 
        A = Interval(1,2)
        A.hi = 3
    to make A into the interval [1,3].  In Python, if we allowed such post hoc 
    changes, we could not prevent spooky change-at-a-distance behavior like 
        A = Interval(1,2)
        B = A
        B.lo = 0
    which makes B into [0,2] but also changes A so it becomes [0,2] too.  To 
    prevent this behavior, we must make Intervals immutable, so we have to use
    special functional setter-like methods like setlo() which copy the object 
    and change its contents on the fly.  The instruction
        B = B.setlo(0)
    is the value-semantic analogue of B.lo=0 that alters the value of B, but 
    not anything that it used to point to, or anything that may have previously
    been used to previously create it.
        A = Interval(1,2)
        B = A
        B = B.setlo(0)
        B   # [0.0, 2.0]
        A   # [1.0, 2.0]  no spooky action at a distance"""
        return Interval(new_lo, self._hi)
    
    def sethi(self, new_hi):
        return Interval(self._lo, new_hi)

    def copy(self):
        return Interval(self.lo, self.hi)

    def left(self): return self.lo

    def right(self): return self.hi

    def range(self): return (self.lo, self.hi)
    
    support = range
    
    def width(self): return self.hi - self.lo
       
    def rad(self): return (self.hi - self.lo)/2
       
    def mid(self): return (self.hi + self.lo)/2

    def __iter__(self):
        if self.lo == self.hi: yield self.lo              # degenerate interval
        else:
            yield self.lo
            yield self.hi
        
    def __len__(self): return 1 if self.lo == self.hi else 2 
       
    # def __repr__(self,d=None):
    #     if abs(self.lo - self.hi) <= IvO.tol:
    #         return f"{self.lo}"
    #     return f"[{self.lo}, {self.hi}]"

    def __repr__(self,d=None):
        if abs(self.lo - self.hi) <= IvO.tol:
            return format_decimals(self.lo,d)
        return '['+format_decimals(self.lo,d)+', '+format_decimals(self.hi,d)+']'

    def __neg__(self): return Interval(-self.hi, -self.lo)
    
    def __abs__(self):                     # does not connect to Interval.abs()
        if 0 <= self.lo: return self
        if self.hi <= 0: return -self
        return Interval(0, max(abs(self.lo),abs(self.hi)))
            
    def mignitude(self):               # least distance between values and zero
        if self.contains(0): return 0   
        return min(abs(self.lo), abs(self.hi))

    def square(self): 
        if 0 <= self.lo: return Interval(left(self)**2, right(self)**2)
        if self.hi <= 0: return Interval(right(self)**2, left(self)**2)
        return Interval(0, max(self.lo**2, self.hi**2))
        
    def __contains__(self, other): # associated with 'in'; see also contains()
        if isinstance(other, Interval): return (other.lo >= self.lo) and (self.hi >= other.hi)
        return self.lo <= other and self.hi >= other    
    
    def subsetof(self, other): # subset or equal to
        other = Interval(other)
        return((other.lo <= self.lo) and (self.hi <= other.hi))
    
    def contains(self, other): # reverse of subsetof, or is it converse, inverse, obverse??
        other = Interval(other)
        return(other.subsetof(self))

    def __add__(self, other):
        if is_pbox(other): return NotImplemented
        if isinstance(other, Interval): return Interval(self.lo + other.lo, self.hi + other.hi)
        return Interval(self.lo + other, self.hi + other)
    
    __radd__ = __add__

    def __sub__(self, other):
        if is_pbox(other): return NotImplemented
        if isinstance(other, Interval): return Interval(self.lo - other.hi, self.hi - other.lo)
        return Interval(self.lo - other, self.hi - other)

    def __rsub__(self, other):
        if is_pbox(other): return NotImplemented
        return Interval(other - self.hi, other - self.lo)

    def __mul__(self, other):
        if is_pbox(other): return NotImplemented
        if isinstance(other, Interval):
            vals = [
                self.lo * other.lo, self.lo * other.hi,
                self.hi * other.lo, self.hi * other.hi
            ]
            return Interval(min(vals), max(vals))
        if other >= 0:
            return Interval(self.lo * other, self.hi * other)
        return Interval(self.hi * other, self.lo * other)

    __rmul__ = __mul__

    def __truediv__(self, other): # self / other
        if straddles(other): raise ZeroDivisionError("Division by interval containing zero")
        if is_pbox(other): return NotImplemented
        if isinstance(other, Interval): return self * Interval(1 / other.hi, 1 / other.lo)
        if isscalar(other):
            if 0 < other: return Interval(self.lo / other, self.hi / other)
            else: return Interval(self.hi / other, self.lo / other)
    
    def __rtruediv__(self, other): # other / self
        if straddles(self): raise ZeroDivisionError("Division by interval containing zero")
        if is_pbox(other): return NotImplemented
        #return Interval(other/self.hi,  other/self.lo)
        return env(other / self.hi,  other / self.lo)

    def env(self, other):                                         # convex hull 
        return Interval(min(self.lo,left(other)), max(self.hi,right(other)))

    def imp(self, other):                                        # intersection
        return Interval(max(self.lo,left(other)), min(self.hi,right(other)))

    def smin(self, other):                     # endpoint or 'sidewise' minimum 
        return Interval(min(self.lo,left(other)), min(self.hi,right(other)))

    def smax(self, other):                     # endpoint or 'sidewise' maximum 
        return Interval(max(self.lo,left(other)), max(self.hi,right(other)))

    def dot(self, other) :                                 # hyperbolic product
        return Interval(self.lo * left(other), self.hi * right(other)) 
     
    def identical(self,other): 
        if not isinstance(other,Interval): return False
        return self.lo==other.lo and self.hi==other.hi
   
    def __eq__(self, other):
        # Copilot insists the next line is better than the following one, but I think it is delusional
        if left(self)==right(self) and left(other)==right(other) and left(self)==left(other): return True
        #if left(self)==right(self)==left(other)==right(other): return True
        if right(self) < left(other) or right(other) < left(self): return False
        return Interval(0,1) # dunno
    
    def __ne__(self, other):
        eq = self.__eq__(other)
        if eq is True:  return False
        if eq is False: return True
        return Interval(0,1)

    def __lt__(self, other):
        if right(self) < left(other): return True
        if left(self) >= right(other): return False
        return Interval(0,1)

    def __le__(self, other):
        if right(self) <= left(other): return True
        if left(self) > right(other): return False
        return Interval(0,1)

    def __gt__(self, other):
        if left(self) > right(other): return True
        if right(self) <= left(other): return False
        return Interval(0,1)
    
    def __ge__(self, other):
        if left(self) >= right(other): return True
        if right(self) < left(other): return False
        return Interval(0,1)
    
    def __rlt__(self, other): return Interval(other).__lt__(self)
    
    def __rle__(self, other): return Interval(other).__le__(self)
    
    def __rgt__(self, other): return Interval(other).__gt__(self)
    
    def __rge__(self, other): return Interval(other).__ge__(self)

def plusminus(c, r): return Interval(c - r, c + r)

def plusminuspercent(x, p): return plusminus(x, abs(x)*p/100)

def midwid(m, w) : return plusminus(m, w/2)

i = I = Interval
PM = plusminus
PMP = plusminuspercent

# ============================================================
# Pbox class
# ============================================================

def long(u):
    try : return len(u)
    except TypeError : return 1
    
def is_monotone(a, incr=True):
    if len(a) == 1: return True
    for i in range(1, len(a)):
        if a[i] >= a[i-1]: continue
        else: return False
    return True 

class Pbox:
    
    def __init__(self, u, d=None, ml=None, mh=None):
        if is_interval(u) and d is None: u,d = left(u), right(u)
        if isscalar(u) or isscalar(d): 
            many = max(long(u), long(d))
            #many = PbO.steps    # force p-boxes to all have the same n
            if isscalar(u): u = [u] * many
            if isscalar(d): d = [d] * many
        uu = list(u)
        dd = list(d if d is not None else u)
        if len(uu) != len(dd): raise ValueError("Left and right sides must have the same length")
        for i in range(len(uu)):
            if is_missing(uu[i]): # NA or None
                if i==0 : uu[i] = float("-inf")
                else:     uu[i] = uu[i-1]           # outward rounding downward
            #print(type(uu[i]), uu[i])
        for i in reversed(range(len(dd))):
            if is_missing(dd[i]): # NA or None
                if i==len(dd)-1: dd[i] = float("inf")             
                else:            dd[i] = dd[i+1]      # outward rounding upward
        u = np.asarray(uu, dtype=float)
        d = np.asarray(dd, dtype=float)      
        #if not is_monotone(u): print('NOT MONOTONE:\n',u)
        if not is_monotone(u): raise ValueError("Left side nonmonotonic")
        if not is_monotone(d): raise ValueError("Right side nonmonotonic")
        self._u = u
        self._d = d
        self._n = len(u)
        self._ml = float(np.mean(u) if ml is None else ml)
        self._mh = float(np.mean(d) if mh is None else mh)

    @property                                              # read-only property
    def u(self): return self._u.copy()     # .copy() to avoid external mutation

    @property
    def d(self): return self._d.copy()

    @property
    def n(self): return self._n

    @property
    def ml(self): return self._ml

    @property
    def mh(self): return self._mh

    def left(self):
        return float(self._u[0])

    def right(self):
        return float(self._d[-1])

    def leftside(self):
        return self._u.copy()

    def rightside(self):
        return self._d.copy()

    def steps(self):
        return self._n

    def mean(self):
        return Interval(self._ml, self._mh)

    def copy(self):
        return Pbox(self._u.copy(), self._d.copy(), self._ml, self._mh)
       
    def setleftside(self, value, index=None):            # see Interval.setlo()    
        if index is None:
            if isinstance(value, (list, tuple, np.ndarray)): new_u = np.asarray(value, dtype=float)
            else: new_u = np.full(self._n, float(value))
        else:
            new_u = self._u.copy()
            new_u[index] = float(value)
        return Pbox(new_u, self._d.copy(), self._ml, self._mh)

    def setrightside(self, value, index=None):
        new_d = self._d.copy()
        if index is not None: new_d[index] = float(value)
        else:    
            if isinstance(value, (list, tuple, np.ndarray)): new_d = np.asarray(value, dtype=float)
            else: new_d = np.full(self._n, float(value))
        return Pbox(self._u.copy(), new_d, self._ml, self._mh)

    def setmean(self, value):
        """Makes a copy of the p-box and reset its mean, but does not check 
        that the new mean is consistent with the p-box edges.  
        
        To change p-box B's mean, you must use B.setmean(newmean)."""
        return Pbox(self._u.copy(), self._d.copy(), left(value), right(value))

    def computemean(self):
        return Interval(mean(self.leftside()),mean(self.rightside()))
            
    def mean(self):    # use ends(mean()) if you need to iterate both endpoints
        return Interval(self.ml, self.mh)
    
    def std(self): return (np.sd(self.u), np.sd(self.d)) # sidewise standard deviations; N.B. this is NOT the standard deviation of a p-box   
    
    def cut(self, p, tight=True):
        if p < 0 or p > 1: raise ValueError("Second argument for cut must be a probability between zero and one")
        n = self.n  
        if tight:   
            p_long = p * n
            fractional = (p_long % 1) == 0
            idx_u = min(n, (1 if fractional else 0) + math.ceil(p_long))
            idx_d = max(1, math.ceil(p_long))
            return Interval(self.u[idx_u-1], self.d[idx_d-1])
        if p == 1: lower = n
        else:
            if (p % (1/n)) == 0: lower = round(p * n)
            else: lower = math.ceil(p * n)
        if p == 0: upper = 1
        else:
            if (p % (1/n)) == 0: upper = round(p * n) + 1
            else: upper = math.floor(p * n) + 1
        return Interval(self.u[max(lower, 1)-1], self.d[min(upper, n)-1])
        
    def median(self):  # conservative w.r.t. discretization  
        """Returns median conservative with respect to discretization."""
        return Interval(self.u[self.n // 2 - (1-self.n % 2)], self.d[self.n // 2])
    # # the conservative median=cut(0.5,False) rather than the optimistic one
    # setoption(steps=10)  # even number of discretization steps
    # a = N(5,1)    
    # m = median(a) # the optimistic median would be 5, rather than [4.75,5.25]
    # print(m)
    # plt.plot([left(m),left(m),right(m),right(m)],[0,1,1,0])
    # plt.plot([1,7],[0.5,0.5]);  plot(a)    
    # plt.show()   
    # setoption(steps=9)  # odd number of discretization steps   
    # a = N(5,1)    
    # m = median(a) # same as the optimistic median
    # print(m)
    # plt.plot([left(m),left(m),right(m),right(m)],[0,1,1,0])
    # plt.plot([1,7],[0.5,0.5]);  plot(a)    
    
    def iqr(self,tight=False):
        return Interval(left(self.cut(0.25,tight)), right(self.cut(0.75,tight)))
    
    def __repr__(self):
        return (f"Pbox(range=[{self.left()}, {self.right()}], "
                f"mean="+Interval(self.ml, self.mh,auto=False).__repr__()+")")
     
    def identical(self, other):
        if not isinstance(other, Pbox): other = Pbox(other)
        return all(same(x,y) for x,y in zip(self.u, other.u)) \
               and all(same(x,y) for x,y in zip(self.d, other.d)) \
               and same(self.ml, other.ml) and same(self.mh, other.mh)

    def __neg__(self):
        return Pbox(-self.d[::-1], -self.u[::-1], ml=-self.mh, mh=-self.ml)

    #def __abs__(self): return abs(self) # this doesn't create Pbox.abs()

    def __add__(self, other):
        if isinstance(other, Pbox):
            return conv_pbox(self, other, op="+")
        if isinstance(other, Interval):
            return self + as_pbox(other)
        return Pbox(self.u + other, self.d + other,
                    ml=self.ml + other, mh=self.mh + other)

    __radd__ = __add__
  
    def __sub__(self, other):
        other = as_pbox(other)
        return conv_pbox(self, -other, op="+")

    def __rsub__(self, other): return as_pbox(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, Pbox):
            return conv_pbox(self, other, op="*")
        if isinstance(other, Interval):
            return self * as_pbox(other)
        if other >= 0:
            return Pbox(self.u * other, self.d * other,
                        ml=self.ml * other, mh=self.mh * other)
        return - (self * (-other))

    __rmul__ = __mul__
    
    def __truediv__(self, other): # self/other
        if is_pbox(other): return self * reciprocate(other)
        if is_interval(other) or is_scalar(other):
            return self * reciprocate(Pbox(other))
        return NotImplemented  
    
    def __rtruediv__(self, other): return reciprocate(self) * other # other/self

    def __eq__(self, other): # Equality is not generally meaningful for distributions or p-boxes
        raise ValueError("Equality comparisons are not meaningful for distributions or p-boxes")

    def __ne__(self, other): # Equality is not generally meaningful for distributions or p-boxes
        raise ValueError("Equality comparisons are not meaningful for distributions or p-boxes")
     
    def prob(self, s=0):
        return Interval( len(self.d[self.d<s])/self.n, len(self.u[self.u<=s])/self.n)

    def xprob(self, s=0):              # required by the inequality comparisons
        return Interval( len(self.d[self.d<=s])/self.n, len(self.u[self.u<s])/self.n)

    # A < B = Pr(A − B < 0)
    # A ≤ B = Pr(A − B ≤ 0)
    # A > B = Pr(B − A < 0)
    # A ≥ B = Pr(B − A ≤ 0)        
     
    def __lt__(self, other): return (self-other).prob()
    def __le__(self, other): return (self-other).xprob()   
    def __gt__(self, other): return (other-self).xprob()
    def __ge__(self, other): return (other-self).prob()   
    def __rlt__(self, other): return Pbox(other).__lt__(self)
    def __rle__(self, other): return Pbox(other).__le__(self)
    def __rgt__(self, other): return Pbox(other).__gt__(self)
    def __rge__(self, other): return Pbox(other).__ge__(self)
    
    def summary(self):
        return {     
            'name':       '',
            'units':      '',
            'shape':      '',
            'mean':       mean(self),
            'check':      self.computemean(),
            #'sd':         sd(self),
            #'var':        var(self),
            'breadth':    breadth(self),
            'iqwidth':    self.iqr().width(),
            'iqr':        self.iqr(),
            'support':    support(self),
            'left':       left(self),
            'pc01':       cut(self, 0.01),
            'pc05':       cut(self, 0.05),
            'pc25':       cut(self, 0.25),
            'median':     cut(self, 0.50),
            'pc75':       cut(self, 0.75),
            'pc95':       cut(self, 0.95),
            'pc99':       cut(self, 0.99),
            'right':      right(self),
            'steps':      steps(self) }

def summarize(x: Pbox): 
    labs = {     
        'line0':      'Summary',
        'name':       '  Name:                ',
        'units':      '  Units:               ',
        'shape':      '  Shape:               ',
        'mean':       '  Average:             ',
        'check':      '  Checked:             ',
        'sd':         '  Std dev:             ',
        'var':        '  Variance:            ',
        'breadth' :   '  Breadth:             ',
        'iqwidth':    '  Interquartile width: ',
        'iqr':        '  Interquartile range: ',
        'support':    '  Support range:       ',
        'line1':      '  Order statistics',
        'left':       '    Left (min) value:  ',
        'pc01':       '    1st percentile:    ',
        'pc05':       '    5th percentile:    ',
        'pc25':       '    25th percentile:   ',
        'median':     '    Median (50th%ile): ',
        'pc75':       '    75th percentile:   ',
        'pc95':       '    95th percentile:   ',
        'pc99':       '    99th percentile:   ',
        'right':      '    Right (max) value: ',
        'steps':      '  Discretizations:     ',         }
    xs = x.summary()
    for k,v in labs.items(): print(v, xs.get(k, '')) 
    
# ============================================================
# Detection
# ============================================================

def isscalar(x):                       # isscalar(I(2,2))                 False    
    return np.isscalar(x)              # isscalar(float('NaN'))            True
                                                    
def is_scalar(x):                      # is_scalar(I(2,2))                 True    
    return left(x) == right(x)         # is_scalar(float('NaN'))          False  
                                       
def is_logical(a): return 0 <= left(a) <= right(a) <= 1

def is_interval(x): return isinstance(x, Interval)

def is_pbox(x): return isinstance(x, Pbox)

def is_zero(x): return is_scalar(x) and same(left(x),0)

"""----------------------------------------------------------------------------
             value   repr   isscalar  is_scalar  .lo   .hi   list   len  ends
a = I(2,3)   [2,3]   [2,3]  False     False      2     3     [2,3]  2    (2,3)
b = I(4,4)   [4,4]   4      False     True       4     4     [4]    1    (4,4)
c = 5        5       5      True      True       *     *     **     ***  (5,5)
-------------------------------------------------------------------------------
* AttributeError: 'int' object has no attribute...  
** TypeError: 'int' object is not iterable
*** TypeError: object of type 'int' has no len()
-------------------------------------------------------------------------------
for A in a: print(A)     # yields 2 and 3
for A in b: print(A)     # just yields 4
for A in c: print(A)     # precipitates TypeError: 'int' object is not iterable
-------------------------------------------------------------------------------
for A in ends(a): print(A)  # yields 2 and 3
for A in ends(b): print(A)  # yields 4 and 4 redundantly
for A in ends(c): print(A)  # yields 5 and 5 redundantly
-------------------------------------------------------------------------------
To convert an interval a into a (short) np.array, use np.array(list(a)). If you 
forget to use list(), you'll get a length-one array consisting of an interval.
-------------------------------------------------------------------------------
"""

# ============================================================
# Coercion helpers
# ============================================================

def as_interval(x):
    if isinstance(x, Interval):            return x
    if isinstance(x, Pbox):                return Interval(x.left(), x.right())
    if isscalar(x):                        return Interval(x, x)
    if hasattr(x,"__len__") and len(x)==2: return Interval(x[0], x[1])
    raise TypeError("Cannot coerce to an interval")

def as_pbox(x):
    if is_pbox(x): return x
    #if is_interval(x): return Pbox(left(x), right(x))
    if is_interval(x): return Pbox(x)
    if isscalar(x): return Pbox(x)
    return Pbox(x) # punt

def ends(x): return (left(x), right(x))

def as_scalar(x): # also see _scalarize()
    if isscalar(x):               return float(x)
    if is_scalar(x):              return left(x)
    if isinstance(x, Number):     return float(x)
    if isinstance(x, Interval):
        y = x.mid()
        good = str([y, x.left(), x.right()][IvO.asscalar])   
        for d in range(1,14):
            st = format_sigdigs(y,d)
            looks = I(st)
            print(d,looks,looks.contains(x) )
            if looks.contains(x) : good = st
            else: return float(good)
        return float(good)
    if isinstance(x, Pbox):       return as_scalar(x.median())
    #if hasattr(x,"__len__") and len(x)==2: return Interval(x[0], x[1])
    raise TypeError("Cannot coerce to an scalar") 
    """   
              IvO.asscalar = 0 =          = 1 =          = 2 =
    as_scalar(I(2,2))               #    2              2              2
    as_scalar(I(1,3))               #    2              1              3
    as_scalar(MMM(0,9,1))           #    1.0101         0              2.0202
    as_scalar(I(1/3, 1/3))          #    0.3333333333   0.3333333333   0.3333333333
    as_scalar(I('1.23'))            #    1.23           1.23           1.23   
    as_scalar(I(14.36184,14.36245)) #    14.362         14.362         14.362     
    as_scalar(True)                 #    1              1              1
    as_scalar(dunno)                #    0.5            0              1
    as_scalar(empty)                #    nan            inf            -inf  
    as_scalar(Interval(NA, 5))      #    -inf           -inf           5   
    as_scalar(Interval(4, NA))      #    inf            4              inf
    as_scalar(Interval(NA, NA))     #    inf            -inf           inf   
    as_scalar(float('NaN'))         #    nan            nan            nan  
    as_scalar(float('inf'))         #    inf            inf            inf  
    #as_scalar(None)                #    TypeError: 'NoneType' object is not iterable
    """      

def as_vectors(data): return np.array([left(datum) for datum in data]), np.array([right(datum) for datum in data])

def as_intervals(lo,hi): return [Interval(L,R) for L,R in zip(lo,hi)]
    
# ============================================================
# Equality, identicality, sameness and similarity
# ============================================================

"""
There are four kinds of 'equality' comparisons among uncertain numbers A and B:
    
    Operation           Assesses
    A == B              epistemic possible equality 
    A.identical(B)      exact structural equality
    same(A,B)           tolerance‑based equality of structure
    A === B             degree of similarity

The equality operator == works as you would expect for scalars, as do the other
various inequality comparison operators.  For intervals, however, == assesses 
whether the underlying numeric quantities COULD BE equal.  It asks about their
possible equality.  Intervals that do not overlap always yield False when
compared with ==.  But, unless both intervals happen to be degenerate scalars, 
if they overlap at all, the == operator will always return the dunno interval 
[0,1] which means 'either True or False' rather than simply True or False.
Python is sensitive to even very small differences in floating point numbers,
but it has limits to its ability to represent differences in floating point 
numbers.  These examples illustrate the behavior of the == equality operator:
    2 == 2                                  # True 
    2 == 2.000000000000001                  # False   (tiny differences matter)
    2 == 2.0000000000000001                 # True  (Python's numerical limits)
    I(1,2) == I(1,2)                        # [0,1]
    I(1,2) == I(1,  2.0001)                 # [0,1]  (precision doesn't matter)
    I(1,2) == I(1,  2.000000000000000000001)# [0,1]
    I(1,2) == I(3,4)                        # False   (intervals don't overlap)
The == operator is NOT defined for probability distributions or p-boxes because
the probability of any two random values being equal is always zero.

The identical() method compares not the underlying numerical values represented 
by the uncertain numbers, but the forms of the uncertain numbers themselves.
The identical() method checks for structural equality.  It returns True only 
if the structure and the numerical values that define that structure exactly 
agree to within machine precision. Again, even very tiny differences can break
equality, but there are limits to floating point representation in Python:
    Interval.identical(I(1,2),I(1,2.000000000000001))  # False
    Interval.identical(I(1,2),I(1,2.0000000000000001)) # True (Python's limits)

The function same(a,b) compares values using a numerical tolerance that is tied 
to the interval display convention.  If values differ by a magnitude no larger 
than this value, they are considered the same.  For example,
    same(2,          2.0000001)       # False
    same(2,          2.00000001)      # True    (if IvO.tol = 1e-8)
    same(I(1,2), I(1,2))              # True
    same(I(1,2), I(1,2.0001))         # False
    same(I(1,2), I(1,2.00000001))     # True    (if IvO.tol = 1e-8)
    a = N(5,1);  same(a,a)            # True
    same(a,N(5,1))                    # True    
The same tolerance also determines how intervals are displayed.  For example, 
    I(   2,          2.00000001)      # 2.0     (if IvO.tol = 1e-8)
    I(   2,          2.0000001)       # [2.0, 2.0000001]
You can set the tolerance value with an instruction like "setoption(val=1e-6)".

The similarity operation === has not yet been implemented, but it will behave
like these examples:
    [1,2] === [1,2]                #    1 
    [1,2] === [1,2.001]            #    0.9995002 
    [1,2] === minmax(1,2)          #    1 
    [1,2] === minmaxmean(1,2,1.3)  #    0.7622963 
This operator must be implemented with a method or function as Python does not
allow defining or overloading a suitable infix operator.          
"""
   
def identical(a,b):
    if type(a) != type(b): return False
    if isinstance(a, Interval): return a.identical(b)
    if isinstance(a, Pbox): return a.identical(b)
    return a == b  # including scalars, strings and everything else
 
def same(a,b):
    if type(a) != type(b): return False
    if isscalar(a): return abs(a-b) < IvO.tol
    if isinstance(a, Interval): return same(a.lo, b.lo) and same(a.hi, b.hi)
    if isinstance(a, Pbox):
        return all(same(x,y) for x,y in zip(a.u, b.u)) \
           and all(same(x,y) for x,y in zip(a.d, b.d)) \
           and same(a.ml, b.ml) and same(a.mh, b.mh)
    return a == b
 
# ============================================================
# Envelope and imposition
# ============================================================

# def Env(*objs):
#     if all(is_scalar(x) or is_interval(x) for x in objs):
#         return Interval(min(left(x) for x in objs), max(right(x) for x in objs))
#     pboxes = [as_pbox(x) for x in objs]
#     if len(pboxes) == 1: return pboxes[0]
#     # flex-envelope across pboxes with arbitrary step counts
#     N = PbO.steps
#     sizes = [pb.n for pb in pboxes]
#     M = max(N, min(200, np.prod(sizes)))
#     Ls = [M // pb.n for pb in pboxes]
#     expanded_us = [np.repeat(pb.u, L) for pb, L in zip(pboxes, Ls)]
#     expanded_ds = [np.repeat(pb.d, L) for pb, L in zip(pboxes, Ls)]
#     # envelope on expanded grid
#     U = np.minimum.reduce(expanded_us)
#     D = np.maximum.reduce(expanded_ds)
#     # downsample to canonical step count, with outward-directed rounding
#     step = M // N
#     U2 = np.empty(N)
#     D2 = np.empty(N)
#     for k in range(N):
#         block_slice = slice(k*step, (k+1)*step)
#         U2[k] = U[block_slice].min()
#         D2[k] = D[block_slice].max()
#     return Pbox(U2, D2, min(pb.ml for pb in pboxes), max(pb.mh for pb in pboxes))

# def Imp(*objs):
#     if all(is_scalar(x) or is_interval(x) for x in objs):
#         lo = max(left(x) for x in objs)
#         hi = min(right(x) for x in objs)
#         if hi < lo: raise ValueError("Imposition empty")
#         return Interval(lo, hi)
#     if all(is_pbox(x) for x in objs):
#         u = np.maximum.reduce([x.leftside() for x in objs])
#         d = np.minimum.reduce([x.rightside() for x in objs])
#         if np.any(d < u): raise ValueError("Imposition does not exist")
#         ml = max(x.ml for x in objs)
#         mh = min(x.mh for x in objs)
#         return Pbox(u, d, ml, mh)
#     raise TypeError("imp expects all Interval or all Pbox")
    
def combine(objs, L, R, means=False):
    """
    Generic side-wise combination operator for scalars, intervals, and p-boxes.
    L  : function that reduces a list/array of left-sides (e.g., np.minimum.reduce)
    R : function that reduces a list/array of right-sides (e.g., np.maximum.reduce)
    """
    if all(is_scalar(x) or is_interval(x) for x in objs): # all scalars,intervals
        return Interval(L([left(x) for x in objs]), R([right(x) for x in objs]))
    pboxes = [as_pbox(x) for x in objs]         # promote everything to p-boxes 
    if len(pboxes) == 1: return pboxes[0]
    N = PbO.steps
    sizes = [pb.n for pb in pboxes]
    M = max(N, min(200, np.prod(sizes)))                # size of the fine grid
    Ls = [M // pb.n for pb in pboxes]                       # expansion factors
    expanded_us = [np.repeat(pb.u, L) for pb, L in zip(pboxes, Ls)]
    expanded_ds = [np.repeat(pb.d, L) for pb, L in zip(pboxes, Ls)]
    U,D = L(expanded_us), R(expanded_ds)   # apply aggregators on expanded grid
    step = M // N
    U2,D2 = np.empty(N), np.empty(N) # conservative downsample to default steps
    for k in range(N):                              # outward-directed rounding
        sl = slice(k*step, (k+1)*step)
        U2[k] = U[sl].min()                   # left side: conservative minimum
        D2[k] = D[sl].max()                  # right side: conservative maximum
    if means: ml,mh = L([b.ml for b in pboxes]), R([b.mh for b in pboxes])
    else: ml,mh = float(np.mean(U2)), float(np.mean(D2))
    return Pbox(U2, D2, ml, mh)

def env(*objs, means=False): 
    '''Compute the envelope of the inputs, which may be a combination of uncertain
    numbers, that is, p-boxes or c-boxes, intervals, distributions, scalars (real
    numbers), or Booleans.

    When the argument 'means' is False, which is the default, the operation merely
    bounds the uncertain numbers. The mean of the output structure is computed anew 
    from these bounds derived from the inputs. To use env() as an aggregation, set 
    means=True. Aggregation implies that at least one of the inputs is surely 
    bounding the quantity's true distribution, so their envelope must also surely 
    bound this quantity's distribution even though we don't know which input does. 
    Because this argument also applies to the means, the mean of the distribution 
    must be inside the envelope of the means of the inputs. The result is a proper 
    characterization of the quantity's distribution and mean given that we don't 
    know which input describes it.'''  
    return combine(objs, np.minimum.reduce, np.maximum.reduce, means)

#- add optional check_empty=True logic for imp

def imp(*objs, means=False): 
    '''Compute the intersection of the uncertain inputs, which may be a combination 
    of intervals, p-boxes, and c-boxes.
    
    When the argument 'means' is False, which is the default, the operation merely
    intersects the uncertain numbers. The mean of the output structure is computed 
    anew from these bounds derived from the inputs. To use imp() as an aggregation, 
    set means=True. Aggregation implies that each of the inputs surely bounds the 
    true distribution of a quantity of interest. If we know that each input surely 
    bounds the quantity's distribution, their intersection must also surely bound 
    this quantity's distribution. Likewise, the mean of the quantity must be inside 
    the intersection of the means of the inputs. If the intersection of inputs does 
    not exist, or its mean is empty (or its endpoints are inverted), the assumption 
    that all the inputs are sure characterizations of the distribution is false.'''
    return combine(objs, np.maximum.reduce, np.minimum.reduce, means)

def constrainedto(a,LO,HI):      # truncates (or shoves) to the support [LO,HI]
    if is_pbox(a) : 
        u = np.minimum(np.maximum(LO,leftside(a)),HI)
        d = np.maximum(np.minimum(rightside(a),HI),LO)
        # what if old mean conflicts with new mean, e.g., in01(N(1,2))?
        return Pbox(u, d, ml=max(a.ml, np.mean(u)), mh=min(np.mean(d), a.mh))
    #return Interval(max(left(a), 0), min(1, right(a)), auto=False)  # on01(2) yields [2, 1]
    return Interval(min(max(left(a), LO),HI), max(min(HI, right(a)),LO), auto=False)

def on01(a): return constrainedto(a,0,1)     # truncates (or shoves) onto [0,1]

# ============================================================
# Sidewise min and max (nonconvolutional n-ary operations)
# ============================================================

def smin(*objs, means=False): return combine(objs, np.minimum.reduce, np.minimum.reduce,means)

def smax(*objs, means=False): return combine(objs, np.maximum.reduce, np.maximum.reduce,means)
   
least = smin
greatest = smax

# see also the convolutional binary functions minimum() and maximum()

# ============================================================
# Convolution (assuming independence or not)
# ============================================================

def negate(x):
    return -x

def reciprocate(x): 
    if left(x) <= 0 and right(x) >= 0: raise ZeroDivisionError("Division by zero")
    if isscalar(x): return 1/x
    if is_interval(x): return 1/x
    u = 1/(rightside(x)[::-1])
    d = 1/(leftside(x)[::-1])    
    return Pbox(u, d, np.mean(u), np.mean(d)) 		 

def zbuffer(b: Pbox):  # buffer the (nonnegative) numeric values away from zero
    if not isinstance(b, Pbox): raise TypeError("zbuffer() requires a Pbox")
    b.u = np.maximum(b.u, PbO.Bzero)
    b.d = np.minimum(b.d, PbO.Bone)
    return b

def conv_pbox(x, y, op="+"):
    Xu, Xd = x.u, x.d
    Yu, Yd = y.u, y.d
    m, p = x.n, y.n
    n = min(200, m * p)
    L = (m * p) // n
    Xu_rep = np.repeat(Xu, p)
    Xd_rep = np.repeat(Xd, p)
    Yu_rep = np.tile(Yu, m)
    Yd_rep = np.tile(Yd, m)
    if op == "+":
        cu = Xu_rep + Yu_rep
        cd = Xd_rep + Yd_rep
    elif op == "*":
        cu = Xu_rep * Yu_rep
        cd = Xd_rep * Yd_rep
        if (x.left() <= 0 <= x.right()) or (y.left() <= 0 <= y.right()):
            c2 = Xu_rep * Yd_rep
            c3 = Xd_rep * Yu_rep
            cu = np.minimum.reduce([cu, c2, c3])
            cd = np.maximum.reduce([cd, c2, c3])
            
            
    # # case RandomNbr::plus:
        # for (i=0; i<m; i++) for (j=0; j<p; j++)	big[i*p+j] = x.d[i] + y.d[j];
        # break;            
        #
    # # case RandomNbr::times:
        # for (i=0; i<m; i++) for (j=0; j<p; j++)	big[i*p+j] = x.d[i] * y.d[j];
        # // Dave Myers' suggestion for making sigma convolution work with straddling factors
        # for (i=0; i<m; i++) for (j=0; j<p; j++)	big[i*p+j] = max(big[i*p+j], x.d[i] * y.u[j]);
        # for (i=0; i<m; i++) for (j=0; j<p; j++)	big[i*p+j] = max(big[i*p+j], x.u[i] * y.d[j]);
        # for (i=0; i<m; i++) for (j=0; j<p; j++)	big[i*p+j] = max(big[i*p+j], x.u[i] * y.u[j]);
        #
    # # case RandomNbr::maximum:
        # for (i=0; i<m; i++) for (j=0; j<p; j++)
        # {
        #     if (x.d[i] > y.d[j]) big[i * p + j] = x.d[i];
        #     else big[i * p + j] = y.d[j];
        # }
        # break;
    # # case RandomNbr::minimum:
        # for (i=0; i<m; i++) for (j=0; j<p; j++)
        # {
        #     if (x.d[i] < y.d[j]) big[i * p + j] = x.d[i];
        #     else big[i * p + j] = y.d[j];
        # }
        # break;
            
            
    else: raise ValueError("Unsupported op")
    cu = np.sort(cu)
    cd = np.sort(cd)
    Zu = cu[::L][:n]
    Zd = cd[L-1::L][:n]
    if   op == "+": ml, mh = ends(x.mean() + y.mean())
    elif op == "*": ml, mh = ends(x.mean() * y.mean())
    else: raise ValueError("Unsupported operation "+op)
    return Pbox(Zu, Zd, ml, mh)

'''
def frechetconv(x,y,op='+'):
    if op=='-': return frechetconv(x,(-y),'+')
    if op=='/': return frechetconv(x,reciprocate.pbox(y),'*')
    if op=='*':
        if straddles(x) or straddles(y): raise NotImplemented('Frechet straddling zero') #return(imp(balchprod(x,y),naivefrechetconv(x,y,'*')))
        if is_zero(x) or is_zero(y): return 0 # prevents an infinite loop
        if (right(x) <= 0) and (right(y)<=0): return frechetconv(-x,-y,'*')
        if right(x) <= 0: return -frechetconv(-x,y,'*')
        if right(y) <= 0: return -frechetconv(x,-y,'*')
    n = PbO.steps    
    zu = np.zeros(n, dtype=float);  zd = np.zeros(n, dtype=float)
    for i in range(n):
        j = np.arange(i, n)
        k = np.arange(n-1, i-1, -1)
        if   op=='+': zd[i] = np.min(x.d[j] + y.d[k])          # lower envelope
        elif op=='*': zd[i] = np.min(x.d[j] * y.d[k])          # lower envelope
        elif op==min: zd[i] = np.min(np.minimum(x.d[j],y.d[k]))# lower envelope
        elif op==max: zd[i] = np.min(np.maximum(x.d[j],y.d[k]))# lower envelope
        j = np.arange(0, i+1)
        k = np.arange(i, -1, -1)       
        if   op=='+': zu[i] = np.max(x.u[j] + y.u[k])          # upper envelope
        elif op=='*': zu[i] = np.max(x.u[j] * y.u[k])          # upper envelope
        elif op==min: zu[i] = np.max(np.minimum(x.u[j],y.u[k]))# upper envelope
        elif op==max: zu[i] = np.max(np.maximum(x.u[j],y.u[k]))# upper envelope
    ml = -np.inf
    mh = np.inf
    if   op == "+": ml, mh = ends(x.mean() + y.mean())
    elif op == "*": ml, mh = ends(x.mean() * y.mean())
    elif op == min:
        if right(y) < left(x): ml, mh = mean(y)
        elif right(x) < left(y): ml, mh = mean(x)
        else: ml, mh = imp(env(smin(mean(x),mean(y)), smin(right(x),right(y))), mean(x)+mean(y)-env(smax(left(x),left(y)), smax(mean(x),mean(y))))
        ml, mh = imp(Interval(np.mean(zu),np.mean(zd)),Interval(ml,mh))
        # ml, mh = imp(z.mymean, VKmeanmaximum(x, y, RandomNbr::dw))   # causes the max(3,N(5,1)) bug
    elif op == max:
        if right(y) < left(x): ml, mh = mean(x)
        elif right(x) < left(y): ml, mh = mean(y)
        else: ml, mh = imp(env(smax(mean(x),mean(y)), smax(right(x),right(y))), mean(x)+mean(y)-env(smin(left(x),left(y)), smin(mean(x),mean(y))))
        ml, mh = imp(Interval(np.mean(zu),np.mean(zd)),Interval(ml,mh))
        # ml, mh = imp(z.mymean, VKmeanmaximum(x, y, RandomNbr::dw))   # causes the max(3,N(5,1)) bug
    else: raise ValueError("Unsupported operation "+op)
    return Pbox(u=zu.copy(), d=zd.copy(), ml=ml, mh=mh )    
''' 
    
  
    
 
# Assumed available in your environment:
#   class Pbox(u, d, ml=None, mh=None)
#   class Interval(lo, hi)
#   PbO.steps
#   left(x), right(x), mean(x), ends(x), smin(a,b), smax(a,b)
#   imp(a, b), env(a, b)
#   straddles(x), is_zero(x)
#   reciprocate.pbox(y)


def _frechet_raw(x, y, op):
    """
    Mixed-resolution Fréchet envelopes for +, *, min, max.
    This is a faithful generalization of your for-loop algorithms.
    x, y are p-boxes with possibly different lengths.
    Returns zd_raw, zu_raw (no resampling).
    """
    xd, xu = np.asarray(x.d, float), np.asarray(x.u, float)
    yd, yu = np.asarray(y.d, float), np.asarray(y.u, float)

    nx, ny = len(xd), len(yd)
    n = min(nx, ny)  # keep the same "i" semantics as your original code

    zd = np.empty(n, dtype=float)
    zu = np.empty(n, dtype=float)

    for i in range(n):
        # LOWER ENVELOPE (d)
        j = np.arange(i, nx)
        k = np.arange(ny - 1, ny - 1 - len(j), -1)  # mirror as in original
        if op == '+':
            vals = xd[j] + yd[k]
            zd[i] = np.min(vals)
        elif op == '*':
            vals = xd[j] * yd[k]
            zd[i] = np.min(vals)
        elif op is min:
            vals = np.minimum(xd[j], yd[k])
            zd[i] = np.min(vals)
        elif op is max:
            vals = np.maximum(xd[j], yd[k])
            zd[i] = np.min(vals)
        else:
            raise ValueError(f"Unsupported op {op}")

        # UPPER ENVELOPE (u)
        j = np.arange(0, i + 1)
        k = np.arange(i, -1, -1)
        j = j[j < nx]
        k = k[k < ny]
        m = min(len(j), len(k))
        j = j[:m]
        k = k[:m]

        if op == '+':
            vals = xu[j] + yu[k]
            zu[i] = np.max(vals)
        elif op == '*':
            vals = xu[j] * yu[k]
            zu[i] = np.max(vals)
        elif op is min:
            vals = np.minimum(xu[j], yu[k])
            zu[i] = np.max(vals)
        elif op is max:
            vals = np.maximum(xu[j], yu[k])
            zu[i] = np.max(vals)

    return zd, zu


def _resample_monotone(zd_raw, zu_raw, n_target):
    """
    Resample zd, zu to length n_target, preserving:
      - zu nondecreasing
      - zd nonincreasing
      - zu <= zd pointwise
    """
    m = len(zd_raw)
    if m == n_target:
        zd = zd_raw.copy()
        zu = zu_raw.copy()
    else:
        old_idx = np.linspace(0.0, 1.0, m)
        new_idx = np.linspace(0.0, 1.0, n_target)
        zd = np.interp(new_idx, old_idx, zd_raw)
        zu = np.interp(new_idx, old_idx, zu_raw)

    # enforce monotonicity
    zu = np.maximum.accumulate(zu)          # nondecreasing
    zd = np.minimum.accumulate(zd[::-1])[::-1]  # nonincreasing

    # enforce zu <= zd
    zu = np.minimum(zu, zd)

    return zd, zu


def frechetconv(x, y, op='+'):
    """
    Fréchet convolution for +, -, *, /, min, max.
    Mixed-resolution safe, resampled to PbO.steps.
    """
    # Reduce to +, *, min, max
    if op == '-':
        return frechetconv(x, -y, '+')
    if op == '/':
        return frechetconv(x, reciprocate.pbox(y), '*')

    # Sign logic for multiplication
    if op == '*':
        if straddles(x) or straddles(y):
            raise NotImplementedError("Fréchet product straddling zero not implemented")
        if is_zero(x) or is_zero(y):
            return 0  # or a zero Pbox in your style
        if (right(x) <= 0) and (right(y) <= 0):
            return frechetconv(-x, -y, '*')
        if right(x) <= 0:
            return -frechetconv(-x, y, '*')
        if right(y) <= 0:
            return -frechetconv(x, -y, '*')

    # Core envelopes (faithful to your for-loop semantics)
    zd_raw, zu_raw = _frechet_raw(x, y, op)

    # Resample to PbO.steps with monotonicity enforced
    n_target = PbO.steps
    zd, zu = _resample_monotone(zd_raw, zu_raw, n_target)

    # Mean / variance logic
    ml, mh = -np.inf, np.inf

    if op == '+':
        ml, mh = ends(mean(x) + mean(y))

    elif op == '*':
        ml, mh = ends(mean(x) * mean(y))

    elif op is min:
        if right(y) < left(x):
            ml, mh = mean(y)
        elif right(x) < left(y):
            ml, mh = mean(x)
        else:
            ml, mh = imp(
                env(smin(mean(x), mean(y)), smin(right(x), right(y))),
                mean(x) + mean(y) - env(smax(left(x), left(y)), smax(mean(x), mean(y)))
            )
        ml, mh = imp(Interval(np.mean(zu), np.mean(zd)), Interval(ml, mh))

    elif op is max:
        if right(y) < left(x):
            ml, mh = mean(x)
        elif right(x) < left(y):
            ml, mh = mean(y)
        else:
            ml, mh = imp(
                env(smax(mean(x), mean(y)), smax(right(x), right(y))),
                mean(x) + mean(y) - env(smin(left(x), left(y)), smin(mean(x), mean(y)))
            )
        ml, mh = imp(Interval(np.mean(zu), np.mean(zd)), Interval(ml, mh))

    else:
        raise ValueError(f"Unsupported operation {op}")

    return Pbox(u=zu, d=zd, ml=ml, mh=mh)  
    
    
  
 



'''
x = U(0,5)
y = U(3,4)
fm = frechetconv(x,y,min);  plot(fm)  
fa = frechetconv(x,y);  plot(fa)

PbO.steps = 10
x = U(0,5)
PbO.steps = 20
y = U(3,4)
PbO.steps = 30
fm2 = frechetconv(x,y,min);  plot(fm2)  
fa2 = frechetconv(x,y);  plot(fa2)
  
plot(fa); plot(fa2,fmt='b')     # Hooray! Copilot made a FLEX Frechet algorithm
plot(fm); plot(fm2,fmt='b')     # Hooray! Copilot made a FLEX Frechet algorithm

x = U(0,5)
y = U(3,4)
z = frechet_min(x,y)
z

x = U(0,5)
y = U(3,4)
z = frechet_max(x,y)
z

#gray(frechetconv(x,y,min),lw=4);  plot(perfectconv(x,y,min),fmt='m:')
#gray(frechetconv(x,y,max),lw=4);  plot(perfectconv(x,y,max),fmt='m:')

'''

   
  
    
  
    
  
    
def frechetconvSLOW(x, y, op='+'):
    n = len(x.d)
    assert len(y.d) == n == len(x.u) == len(y.u)
    zd = np.empty(n, dtype=float)
    zu = np.empty(n, dtype=float)
    for i in range(n):
        outlier = inf
        for j in range(i, n):
            k = i - j + n - 1
            here = x.d[j] + y.d[k]
            if here < outlier: outlier = here
        zd[i] = outlier

        outlier = -inf
        for j in range(0, i + 1):
            k = i - j
            here = x.u[j] + y.u[k]
            if here > outlier: outlier = here
        zu[i] = outlier
    return Pbox(zu, zd)    
  

import operator

ops = {
    '+': operator.add,
    '*': operator.mul,
    '-': operator.sub,
    '/': operator.truediv,
    '**':operator.pow,
#    'v': pmax,        # join
#    '^': pmin,        # meet
    } # add more as needed

def docall(op, args): return ops[op](*args)

def perfectconv(a, b, op='+'):
    def protectedlog(a): return(math.log(zbuffer(a)))
    #if op == ('^','**'): return(exp(perfectconv.pbox(protectedlog(a),b,'*')))  # prolly doesn't work for ^ in interesting cases
    if op in ('-','/'): cu,cd = docall(op, [a.u, b.d]), docall(op, [a.d, b.u])
    elif op == max:     cu,cd = np.maximum(a.u, b.u), np.maximum(a.d, b.d)
    elif op == min:     cu,cd = np.minimum(a.u, b.u), np.minimum(a.d, b.d)
    else:               cu,cd = docall(op, [a.u, b.u]), docall(op, [a.d, b.d])
    return Pbox(np.sort(cu), np.sort(cd))
    #scu = np.sort(cu + 0); scd = np.sort(cd + 0); if (all(cu == scu) and all(cd == scd)) return pbox(u=scu, d=scd,  dids=paste(a@dids,b@dids), bob=a@bob) else return pbox(u=scu, d=scd,  dids=paste(a@dids,b@dids)) 
  
def oppositeconv(a, b, op='+'): # prolly doesn't work for ^ in interesting cases
    if (op in ('-', '/')):
        cu = docall(op, [a.u, b.d[::-1]])
        cd = docall(op, [a.d, b.u[::-1]])
    else:
        cu = docall(op, [a.u, b.u[::-1]])
        cd = docall(op, [a.d, b.d[::-1]])
    return Pbox(np.sort(cu),np.sort(cd))
    #pbox(u=sort(cu), d=ort(cd),  dids=paste(a@dids,b@dids)) # neither bob nor not-bob
  
"""
positiveconv = positiveconv.pbox = function(a,b,op='+') {                # positive (PQD) dependence
  if (op %in% c('-')) return(negativeconv.pbox(a, -b, '+'))
  if (op %in% c('/')) return(negativeconv.pbox(a, reciprocate(b), '*'))
  if (op %in% c('^')) if ((straddles(a-1)) || (left(a) < 0)) stop('say what?')
  n = Pbox$steps
  cu = cd = rep(0,n)
  if (op %in% c('+','*','^')) {
        for (i in 1:n) {
                infimum = inf
                for (j in i:n) {                                                           # convert the for loops to increase the speed of this function
                        if (op == '+') here = a@d[[j]] + b@d[[n*i/j]] else 
                        if (op == '*') here =  a@d[[j]] * b@d[[n*i/j]]
                        if (op == '^') here =  a@d[[j]] ^ b@d[[n*i/j]]
                        if (here<infimum) infimum = here
                        }
                cd[i] = infimum
                supremum = -inf
                for (j in 1:i) {
                        kk = floor(1 + n * ((i-1)/n-(j-1)/n)/(1-(j-1)/n))
                        if (op == '+') here = a@u[[j]] + b@u[[kk]] else
                        if (op == '*') here = a@u[[j]] * b@u[[kk]]
                        if (op == '^') here = a@u[[j]] ^ b@u[[kk]]
                        if (here>supremum) supremum = here
                        }
                cu[i] = supremum
                }
        if (op=='+') v = env(var(a)+var(b), var(a)+var(b)+2*sqrt(var(a)*var(b))) else v=interval(0,Inf)
		if (op=='^') if (right(a)<=1) {safe=cu; cu=sort(cd); cd=sort(safe)} # not sure this is correct
        return(pbox(cu, cd, ml=a@ml+b@ml, mh=a@mh+b@mh, vl=left(v), vh=right(v), dids=paste(a@dids,b@dids)))
        }               
  }

# Debugging the ^ operator in positiveconv
#par(mfrow=c(2,2))
#checkem = function(a,b) {
#  plot(frechetconv.pbox(a,b,'^'),col='blue')  # frechet (blue) should enclose all
#  lines(positiveconv(a,b,'^'),lw=3)             # should enclose independent (red) and perfect (green)
#  lines(conv.pbox(a,b,'^'),col='red')
#  lines(perfectconv.pbox(a,b,'^'),col='green')
#  }
#a = U(2,4)             # bases all above one
#b = N(5,0.5)           # exponents all positive
#checkem(a,b)
#
#a = U(2,4)             # bases all above one
#b = N(2,1)             # positive and negative exponents
#checkem(a,b)
#
#a = U(0.2,0.4)         # bases all below one but positive
#b = N(5,0.5)           # exponents all positive
#checkem(a,b)
#
#a = U(0.2,0.4)         # bases all below one but positive
#b = N(2,1)             # positive and negative exponents
#checkem(a,b)

negativeconv = negativeconv.pbox = function(a,b,op='+') {                # negative (PQD) dependence
  if (op %in% c('-')) return(positiveconv.pbox(a, -b, '+'))
  if (op %in% c('/')) return(positiveconv.pbox(a, reciprocate(b), '*'))
  if (op %in% c('^')) stop('say what?')
  n = Pbox$steps
  cu = cd = rep(0,n)
#  stop('say what?')
  if (op %in% c('-','/')) {
	for (i in 1:n) {
		infimum = inf
		for (j in i:n) {                                                           # convert the for loops to increase the speed of this function
			if (op == '-') here = a@d[[j]] - b@d[[n*i/j]] else 
			if (op == '/') here =  a@d[[j]] / b@d[[n*i/j]]
			if (here<infimum) infimum = here
			}
		cu[i] = infimum
		supremum = -inf
		for (j in 1:i) {
			kk = floor(1 + n * ((i-1)/n-(j-1)/n)/(1-(j-1)/n))
			if (op == '-') here = a@u[[j]] - b@u[[kk]] else
			if (op == '/') here = a@u[[j]] / b@u[[kk]]
			if (here>supremum) supremum = here
			}
		cd[i] = supremum
		}
	if (op=='-') {m = interval(a@ml-b@mh, a@mh-b@ml); v = env(var(a)+var(b), var(a)+var(b)+2*sqrt(var(a)*var(b))) } else {m=interval(-Inf,Inf); v=interval(0,Inf) }
	return(pbox(cu, cd, ml=left(m), mh=right(m), vl=left(v), vh=right(v), dids=paste(a@dids,b@dids)))
	}
  }
"""

def pmin(a,b):
    return frechetconv(a,b,min) 
   
def pmax(a,b):
    return frechetconv(a,b,max) 
   
# def pmin(a,b,d=Dependence.UNKNOWN):
#     d = normalize_dependence(d)
#     if d == Dependence.INDEPENDENCE:return conv_pbox(a,b,min)
#     if d == Dependence.PERFECT:     return perfectconv(a,b,min)
#     if d == Dependence.OPPOSITE:    return oppositeconv(a,b,min)
#     #if d == Dependence.POSITIVE:    return <<>>
#     #if d == Dependence.NEGATIVE:    return <<>>
#     return frechetconv(a,b,min) 
     
# ============================================================
# Quantile functions
# ============================================================

# These quantile functions assume that p is an array of probability 
# values, and any distribution parameters are scalars.  The function 
# qpbox() is used to assemble p-boxes from interval distribution 
# parameters via enveloping.  This strategy works with most of these
# quantile functions, but it depends on the reasonableness of the
# distribution parameterizations.  It does not work well for some
# distributions such as uniform, triangular, trapezoidal, and the
# sawinconrad distribution, which need a special function to make
# p-boxes.

# In scipy (pronounced skippy), the expression qnbinom(0,123,) yields 
# np.float64(-1.0), which is a very busy way to say -1.  Of course, in R, 
# the expression returns 0, which is correct, as the negative binomial 
# distribution is defined on {0,1,2,3,...}, not on -1.  This 'quirk' of
# returning -1 for quantiles of nonnegative discrete distributions is 
# actually just a bug, but it pervades all the discrete distributions 
# that scipy (pronounced skippy!) supports. We fix them here by max-ing 
# the outcome with 0.
 
def max0_ppf(dist):
    def q(p, *params): return np.maximum(0, dist.ppf(p, *params))
    return q

qbernoulli        = max0_ppf(sps.bernoulli)
qbetabinomial     = max0_ppf(sps.betabinom)
qbinomial         = max0_ppf(sps.binom)
qgeometric        = max0_ppf(sps.geom)
qhypergeometric   = max0_ppf(sps.betabinom)
qnegativebinomial = max0_ppf(sps.nbinom)
qpascal           = max0_ppf(sps.geom)
qpoisson          = max0_ppf(sps.poisson)

def qbeta(p, a, b):
    #if (straddles(a) and straddles(b)): return(Pbox(0,1))
    if a==0: return np.zeros_like(p)
    if b==0: return np.ones_like(p)    
    return sps.beta(a=a, b=b).ppf(p)

def qbeta1(p, m, s):
    tmp = m*(1-m)/(s*s) - 1   # tmp has repeated variables  
    return sps.beta(a=m*tmp, b=(1-m)*tmp).ppf(p)

#def qgamma(p, a, t): return sps.gamma.ppf(p, a, scale=t)

def qgamma(p, shape, scale=1, rate=None): 
    if rate is not None: scale = 1/rate    
    return sps.gamma.ppf(p, a=shape, scale=scale)

def qlogistic(p, mu, s): return sps.logistic.ppf(p, loc=mu, scale=s)

def qlognormal(p, m, s):
    # Convert mean/sd -> mu/sigma
    sigma = np.sqrt(np.log(1 + (s*s)/(m*m)))
    mu = np.log(m) - 0.5*sigma*sigma
    return sps.lognorm.ppf(p, s=sigma, scale=np.exp(mu))

def qMMML(p, lo, hi, m):
    mid = (hi - m) / (hi - lo)
    out = np.empty_like(p)
    mask = (p <= mid)
    out[mask]  = lo
    out[~mask] = (m - hi) / p[~mask] + hi
    return out                          # u=ifelse(p<=mid,min,(mean-max)/p+max)

def qMMMR(p, lo, hi, m):
    mid = (hi - m) / (hi - lo)
    out = np.empty_like(p)
    mask = (mid <= p)
    out[mask]  = hi
    out[~mask] = (m - lo * p[~mask]) / (1 - p[~mask])
    return out                        # d=ifelse(mid<=p,max,(mean-min*p)/(1-p))

def qnormal(p, mu=0, sigma=1): return sps.norm(loc=mu, scale=sigma).ppf(p)

def qsawinconrad(p, min, mu, max):  
    from scipy.optimize import root_scalar
    def f(alpha):
        num = max * np.exp(alpha * max) - min * np.exp(alpha * min)
        den = np.exp(alpha * max) - np.exp(alpha * min)
        return num / den - 1/alpha - mu
    sol = root_scalar(f, bracket=[-15, 15])
    alpha = sol.root
    return np.log(1 + p * (np.exp(alpha * (max - min)) - 1)) / alpha + min

def qsmirnov(p, n): return sps.ksone.ppf(p,n) # OMG supposed to be the same as R's smirnov
    
def qtrapezoidal(p, a, b, c, d):
    if same(d,a): return np.full_like(p,a)
    if same(c,b): return qtriangular(p,a,b,d)
    h = 2 / (c+d-b-a)
    p1 = h * (b-a)/2
    p2 = p1 + h * (c-b)
    r = np.where(p <= p2, (p - p1) / h + b, d - np.sqrt(2 * (1-p) * (d-c) / h))
    r[p<=p1] = a + np.sqrt(2 * p[p<=p1] * (b-a)/h)   
    return r

def qtriangular(p, a, c, b):
    Fc = (c - a) / (b - a)
    return np.where(
        p < Fc,
        a + np.sqrt(p * (b - a) * (c - a)),
        b - np.sqrt((1 - p) * (b - a) * (b - c)))

def quniform(p, a, b): return sps.uniform(loc=a, scale=b - a).ppf(p)

def qweibull(p, k, lam): return sps.weibull_min.ppf(p, c=k, scale=lam)

# ============================================================
# Distribution constructors
# ============================================================

# def qpbox(qfun, i, j, m, *params):                 # recursive implementation
#     if isinstance(qfun, list):
#         qfunL = qfun[0]
#         qfunR = qfun[1]
#     else: qfunL = qfunR = qfun
#     args = [as_interval(p) for p in params]
#     if all(is_scalar(iv) for iv in args):
#         scalars = [iv.left() for iv in args]
#         u = qfunL(i, *scalars)
#         d = qfunR(j, *scalars)
#         ml = left(m)
#         mh = right(m)
#         return Pbox(u, d, ml=ml, mh=mh)
#     corners = [(iv.left(), iv.right()) for iv in args]
#     import itertools
#     combos = list(set(list(itertools.product(*corners)))) # unique combinations
#     results = []
#     for combo in combos:
#         try: # recursively call qpbox with scalar parameters
#             pb = qpbox(qfun, i, j, m, *combo)
#             if validate_pbox(pb): results.append(pb) # Validate p-box
#         except Exception: pass  # Skip invalid combinations
#     if not results: raise ValueError("No parameters yield a valid distribution")
#     env_pb = env(*results)
#     return Pbox(env_pb.leftside(), env_pb.rightside(), ml=left(m), mh=right(m))
    
import itertools

def qpbox(qfun, i, j, m, *params):    # nonrecursive version using new iterator
    if isinstance(qfun, (list, tuple)): qfunL, qfunR = qfun
    else: qfunL = qfunR = qfun
    args = [as_interval(p) for p in params]
    corners = [(iv.left(), iv.right()) for iv in args]
    combos = itertools.product(*corners)
    left_vecs, right_vecs = [], []
    for combo in combos:
        try:
            uvals = np.asarray(qfunL(i, *combo))
            dvals = np.asarray(qfunR(j, *combo))
            if np.isnan(uvals).any() or np.isnan(dvals).any(): continue # skip this corner combination
            left_vecs.append(uvals)
            right_vecs.append(dvals)
        except Exception: pass
    if not left_vecs or not right_vecs: raise ValueError("No parameters yield a valid distribution")  
    left_stack  = np.vstack(left_vecs) # stack into 2D arrays (n_combos,n_points)
    right_stack = np.vstack(right_vecs)
    u_env = left_stack.min(axis=0)   # elementwise envelope across combinations
    d_env = right_stack.max(axis=0)
    return Pbox(u_env, d_env, ml=left(m), mh=right(m))

def validate_pbox(pb):
    if not (hasattr(pb,'u') and hasattr(pb,'d')): return False, 'missing u or d'
    u = np.asarray(pb.u);    d = np.asarray(pb.d)
    if u.ndim != 1 or d.ndim != 1 or len(u) != len(d): return False, 'edges must be matching vectors'
    if np.isnan(u).any() or np.isnan(d).any(): return False, 'NaN not allowed'    
    #if not np.all(np.isfinite(u)) or not np.all(np.isfinite(d)): return False, 'infinities not allowed'  
    if not(is_monotone(u) and is_monotone(d)): return False, 'edges must be nondecreasing'
    if np.any(u > d): return False, 'leftside cannot exceed rightside'
    return True, 'ok'

# def beta(v, w):
#     m = 1/(1 + w/v)
#     return qpbox(sps.beta.ppf, PbO.ii(), PbO.jj(), m, v, w)

def beta(v, w):
    if (straddles(v) and straddles(w)): return(Pbox(0,1))
    if v==0: return Pbox(0,0)
    if w==0: return Pbox(1,1)
    else: return qpbox(qbeta, PbO.ii(), PbO.jj(), Interval(0,1/(1+w/right(v))) if straddles(v) else 1/(1+w/v), v, w)
    #else: return qpbox(qbeta, PbO.ii(), PbO.jj(), 1/(1+w/v), v, w)
    
def beta1(m, s):
    return qpbox(qbeta, PbO.ii(), PbO.jj(), m, m, s)

def chisquared(df):
    return qpbox(sps.chi2.ppf, PbO.ii(), PbO.jjj(), df, df)

def exponential(m=1, rate=None):
    if rate is not None: m = 1/rate
    return qpbox(sps.expon.ppf, PbO.ii(), PbO.jjj(), m, 0, m)

def F(df1, df2):
    m = df2/(df2 - 2) if df2 > 2 else np.inf  # mean exists for df2 > 2: df2/(df2-2)
    return qpbox(sps.f.ppf, PbO.ii(), PbO.jjj(), m, df1, df2)

# Whenever you make a gamma distribution, we recommend NAMING the parameters.
# Confusion in gamma distribution parameterization is the worst in probability:
#    Evans et al. (1993) says gamma(scale, shape) and
#    Forbes et al. (2011) still says gamma(scale, shape), so
#    Risk Calc says gamma(scale, shape), but
#    Scipy says gamma(shape, loc, scale) and
#    R says gamma(shape, rate = 1, scale = 1/rate),
#    Wikipedia says gamma(shape, scale) or gamma(shape, rate),
# and the conventions are often not even consistent within a package for the 
# related distributions such as the inverse-gamma & normal-gamma distributions.

def gamma(shape, scale=1, rate=None):
    m = shape * scale
    #return qpbox(sps.gamma.ppf, PbO.ii(), PbO.jjj(), m, shape, 0, scale) # the zero is loc
    return qpbox(qgamma, PbO.ii(), PbO.jjj(), m, shape, scale)

def gamma1(mean, sd):
    shape = (mean/sd)**2
    scale = sd*sd/mean
    return gamma(shape, scale)

def gamma2(rate, shape):                          # DEPRECATED parameterization
    return gamma(shape=shape, scale=1/rate)

def gamma3(scale, shape): # Risk Calc             # DEPRECATED parameterization
    return gamma(shape=shape, scale=scale)

"""
A = gamma(2,13)
B = gamma1(2,13)
C = gamma2(2,13)
D = gamma3(2,13)              
cyan(A); blue(B); green(C); red(D)


# this is what's happening in R 

a = gamma(2,13)
b = gamma1(2,13)
c = gamma2(2,13)
c; blue(b); cyan(c); red(a)


"""

def gumbel(mu, beta_):
    # R's gumbel(mean, sd) mapping was: mean = mu + gamma*beta, sd = beta*pi/sqrt(6)
    # Here we assume parameters are location mu and scale beta_ directly.
    m = mu + 0.5772156649015329*beta_
    return qpbox(sps.gumbel_r.ppf, PbO.iii(), PbO.jjj(), m, mu, beta_)

def smirnov(n):   
    mean = sps.ksone.stats(n,moments='m')     # # # # what if n is an interval?
    return qpbox(qsmirnov, PbO.ii(),PbO.jjj(), mean, n) 

def laplace(mu, b):
    return qpbox(sps.laplace.ppf, PbO.iii(), PbO.jjj(), mu, mu, b)

def logistic(mu, s):
    return qpbox(qlogistic, PbO.iii(), PbO.jjj(), mu, mu, s)

def lognormal(m, s): # m, s are the mean and standard deviation of the resulting distribution
    sigma = sqrt(log(1 + square(s)/square(m)))
    mu = log(m) - 0.5*square(sigma)
    return qpbox(sps.lognorm.ppf, PbO.ii(), PbO.jjj(), m, sigma, exp(mu))

def lognormal2(mu, sigma): # mu, sigma are the mean and std of the underlying normal distribution (the distribution of the logs)
    mean = np.exp(mu + 0.5 * sigma * sigma)
    return qpbox(sps.lognorm.ppf, PbO.ii(), PbO.jjj(), mean, sigma, 0, np.exp(mu))

def lognormal3(gm, gsd): # geometric mean, geometric standard deviation
    mu = np.log(gm)
    sigma = np.log(gsd)
    mean = np.exp(mu + 0.5 * sigma * sigma)
    return qpbox(sps.lognorm.ppf, PbO.ii(), PbO.jjj(), mean, sigma, 0, np.exp(mu))

def loguniform1(mean, sd):
    # R's loguniform1(mean, sd) is a method-of-moments constructor;
    # here we assume you want a log-uniform on [a,b] with given mean/sd.
    # For a log-uniform on [a,b], X = exp(U), U~Uniform(log a, log b).
    # mean and sd equations are messy; if you already have a,b, better to
    # call a direct constructor. For now, we leave this as a placeholder
    # or simple uniform in log-space with given mean/sd of X.
    # You may want to replace this with your exact R logic if needed.
    raise NotImplementedError("loguniform1(mean, sd) needs your specific moment-matching formula.")

def minmaxmean(lo, hi, mu):
    return qpbox([qMMML,qMMMR], PbO.ii(), PbO.jj(), mu, lo, hi, mu)

def negativebinomial(r, p, n=200):
    """Negative binomial given r = shape (number of failures) and p = success probability."""
    return qpbox(qnegativebinomial, PbO.iii(), PbO.jjj(), r*(1/p-1), r, p)

def normal(mu, sigma):
    return qpbox(qnormal, PbO.iii(), PbO.jjj(), mu, mu, sigma)

normal1 = normal

def pareto(xm, alpha):
    m = xm*alpha/(alpha - 1) if alpha > 1 else np.inf
    return qpbox(sps.pareto.ppf, PbO.ii(), PbO.jjj(), m, alpha, xm)

def powerfunction(a, b):
    # R's powerfunction(a,b): support [0,1] with shape b, then scaled/shifted.
    # Here we assume X = a + (Y)*(something) was already handled in MMpowerfunction.
    # If you have a canonical powerfunction distribution, plug its ppf here.
    raise NotImplementedError("powerfunction(a,b) needs your chosen parameterization.")

# def sawinconrad(min, mu, max):
#     #a = left(min);      b = right(max)
#     c = left(mu);       d = right(mu)
#     if c < a: c = a
#     if d > b: d = b
#     return qpbox(qsawinconrad, PbO.ii(), PbO.jjj(), mu, min, mu, max)

def sawinconrad(min, mu, max):
    mu = I(min,max).imp(mu)
    print(min,mu,max)
    return qpbox(qsawinconrad, PbO.ii(), PbO.jjj(), mu, min, mu, max)

# a = sawinconrad(0,3,10)
# b = sawinconrad(0,4,10)
# ab = env(a,b)
# c = sawinconrad(0,I(3,4),10)
# C = Sawinconrad(0,I(3,4),10)
# ab; plot(c,fmt='g'); plot(C,fmt='b')


# a = sawinconrad(0,3,12)
# b = sawinconrad(0,4,10)
# ab = env(a,b)
# c = sawinconrad(0,I(3,4),I(10,12))
# C = Sawinconrad(0,I(3,4),I(10,12))
# ab; plot(c,fmt='g'); plot(C,fmt='b');plot(c,fmt='g')
# 1


def student(df):
    return qpbox(sps.t.ppf, PbO.iii(), PbO.jjj(), 0, df) # mean = 0 for df > 1



# trapezoidal
# function(min, lmode, rmode, max, ...){
#   a <- left(min);   b <- right(min)
#   c <- left(lmode); d <- right(lmode)
#   e <- left(rmode); f <- right(rmode)
#   g <- left(max);   h <- right(max)
#   if (c<a) c <- a   # implicit constraints
#   if (e<c) e <- c
#   if (g<e) g <- e
#   if (h<f) f <- h
#   if (f<d) d <- f
#   if (d<b) b <- d
#   # moments
#   Strapezoidalmean <- function(a, b, c, d){
#     ab <- a + b
#     cd <- c + d
#     if (nothing(cd-ab)) h <- 1 else h <- 1.0 / (3.0 * (cd - ab))
#     return(h * (c * cd + d*d - (a * ab + b*b)))
#     }
#   Strapezoidalvar <- function(a, b, c, d){
#     ab <- a + b
#     cd <- c + d
#     if (nothing(cd-ab)) h <- 1 else h <- 1.0 / (3.0 * (cd - ab))
#     m <- h * (c * cd + d*d - (a * ab + b*b))
#     if (nothing(d - a)) {
#                   m <- a   
#                   v <- 0.0
#            } else v <- 0.5 * h * (cd * (c*c + d*d) -  ab * (a*a + b*b)) - m*m
#     return(v)
#     }
#   ml <- Strapezoidalmean(a,c,e,g)
#   mh <- Strapezoidalmean(b,d,f,h)
#   if (g<=b) vl <- 0 else vl <- Strapezoidalvar(b,closest(lmode,(b+g)/2),closest(rmode,(b+g)/2),g)
#   vh <- Strapezoidalvar(a,closest(interval(c,d),a),closest(interval(e,f),h),h)
#   pbox(u=qtrapezoidal(ii(), a, c, e, g), d=qtrapezoidal(jj(), b, d, f, h), shape='trapezoidal', ml=ml, mh=mh, vl=vl, vh=vh, ...)
#   }



def trapezoidal(min, lmode, rmode, max):
    return qpbox(qtrapezoidal, min, lmode, rmode, max)
  

def triangular(a, c, b):
    loc = a
    scale = b - a
    shape = (c - a)/scale
    m = (a + b + c)/3
    return qpbox(qtriangular, PbO.ii(), PbO.jj(), m, shape, loc, scale)

def uniform(a, b):
    return qpbox(quniform, PbO.ii(), PbO.jj(), (a + b) / 2, a, b)

def uniform1(mean, sd):
    r = sd * np.sqrt(3)
    return qpbox(quniform, PbO.ii(), PbO.jj(), mean, mean - r, mean + r)

def weibull(scale, shape):
    # scipy: weibull_min(c=shape, scale=scale)
    # mean = scale * Gamma(1 + 1/shape)
    
    from scipy.special import gamma as Gamma
    
    m = scale * Gamma(1 + 1/shape)
    return qpbox(sps.weibull_min.ppf, PbO.ii(), PbO.jjj(), m, shape, 0, scale)


##########################################################################
# Discrete distributions
##########################################################################

def bernoulli(p):
    return qpbox(qbernoulli, PbO.ii(), PbO.jj(), p, p)

def binomial(n, p):
    return qpbox(qbinomial, PbO.ii(), PbO.jj(), n*p, n, p)

def Betabinomial(n, alpha, beta_):
    m = n*alpha/(alpha + beta_)
    return qpbox(qbetabinomial, PbO.ii(), PbO.jj(), m, n, alpha, beta_)


def betabinomial(n, alpha, beta):
    # Handle interval or scalar inputs uniformly
    aL, aU = left(alpha), right(alpha)
    bL, bU = left(beta), right(beta)
    nL, nU = left(n), right(n)

    # If n includes 0, the support includes 0
    if nU == 0: return Pbox(0, 0)

    # Degenerate cases from alpha or beta hitting zero
    degenerate_min = None
    degenerate_max = None

    # alpha = 0 forces X = 0
    if aL == 0: degenerate_min = 0

    # beta = 0 forces X = n
    if bL == 0: degenerate_max = nU

    # If both alpha and beta include 0, envelope is [0, n]
    if degenerate_min == 0 and degenerate_max == nU: return Pbox(0, nU)

    # If only alpha includes 0, include degenerate-at-0 in envelope
    tiny = 1e-8

    if degenerate_min == 0 and degenerate_max is None:
        core = betabinomial(n, max(aL, tiny), beta)   # interior part
        return env(Pbox(0,0), core)

    # If only beta includes 0, include degenerate-at-n in envelope
    if degenerate_max == nU and degenerate_min is None:
        core = betabinomial(n, alpha, max(bL, tiny))
        return env(core, Pbox(nU, nU))

    # Otherwise: standard Beta-Binomial
    m = n * 1/(1 + beta/alpha)
    return qpbox(qbetabinomial, PbO.ii(), PbO.jj(), m, n, alpha, beta)

def poisson(lam):
    return qpbox(qpoisson, PbO.ii(), PbO.jjj(), lam, lam)

def geometric(p):  # number of failures before first success, mean = (1-p)/p
    return qpbox(qgeometric, PbO.ii(), PbO.jjj(), (1 - p)/p, p)


B = beta
MMM = minmaxmean
N = normal
T = triangular
U = uniform




# ============================================================
# Histogram functions
# ============================================================
  
def EDF(x):
    """Construct an empirical p-box from data x, which may contain intervals."""
    def SEDF(x): # empirical distribution from a vector of scalars
        x = np.asarray(x)
        n = len(x)
        sx = np.argsort(x)
        steps = PbO.steps
        u = x[sx[np.arange(steps) * n // steps]]
        d = x[sx[np.arange(steps) * n // steps]]
        if steps % n == 0: d = u.copy()
        return Pbox(u, d)
    x, y = as_vectors(x)
    if np.all(x==y): return SEDF(x)
    else: return env(SEDF(x),SEDF(y))
    
def fatten(a, dm=0, leftbound=None, rightbound=None):
    """Raise and lower edges of p-box a by dm and widen the support to [leftbound, rightbound]."""
    if leftbound is None: leftbound = left(a)
    if rightbound is None: rightbound = right(a)
    Au, Ad = leftside(a), rightside(a)
    ii_vals = PbO.ii() >= (1 - dm)
    jj_vals = PbO.jjj() <= dm
    Au = np.concatenate([np.full(np.sum(ii_vals), leftbound),Au[~ii_vals]])
    Ad = np.concatenate([Ad[~jj_vals],np.full(np.sum(jj_vals), rightbound)])
    # # update variance bounds
    # Avl, Avh = dwVariance(A)
    return Pbox(Au,Ad, float(np.mean(Au)),float(np.mean(Ad)))

def KS_critical(n, conf=0.95, two_sided=True):
    """The critical value d is such that P(D_n<=d)=conf (two-sided by default).
    You almost surely want two_sided set to True. Use False for a stochastic 
    dominance test, or you need lower bounds on quantiles, or to assert that 
    the true distribution is not stochastically larger than the EDF.  We use 
    SciPy's kstwo, which implements the Marsaglia–Tsang–Wang algorithm."""
    alpha = 1.0 - conf
    if two_sided: return sps.kstwo.isf(alpha, n) # two-sided Kolmogorov–Smirnov
    else: return np.sqrt(-0.5 * np.log(alpha) / n) # asympototic approximation; P(D_n^+ <= d) ≈ 1-exp(-2 n d^2) => d ≈ sqrt(-0.5*log(alpha)/n)

def histogram(x, y=None, mn=None, mx=None, conf=0.95):
    """Construct a p-box representing the Kolmogorov-Smirnov confidence band 
    around the empirical distribution of possibly imprecise values."""
    if y is None:  x, y = as_vectors(x)
    if np.all(x==y): H = EDF(x)
    else: H = env(EDF(x), EDF(y))
    if mn is None: mn = min(np.min(x), np.min(y))
    if mx is None: mx = max(np.max(x), np.max(y))
    return fatten(H, KS_critical(len(x), conf), leftbound=mn, rightbound=mx) #return Pbox(fatten(H, KSDmax(n, conf), leftbound=mn, rightbound=mx), shape='histogram')


# Conal Brown code for normal confidence bands

# single-point confidence bands



# ============================================================
# Scott's "simple inputs" 
# ============================================================

"""
Uncertain numbers are WAY too complex. We must be useful to analysts who aren’t 
sure what a normal distribution is, and have never heard of a p-box.  These are 
the analysts who most need a UC tool that maybe employs fancy probabilistic and
non-Laplacian uncertainty, but demands very coarse inputs, possibly expressed 
verbally such as “between 50 and 100”, “1 out of 10”, “less than 25%” or “about
9.3” and fashions inputs from whatever information a user has and assumptions 
they're comfortable making.

orderof
O
plusminus
plusminuspercent

KN, km

format_sigdigs
format_decimals
sgnumber
sigorder
lastsigfig
about
    
around
almost
elicited
exactly
<<above>>  downto
<<below>>   upto
nearly
"""

# Fermi seat-of-the-pants estimates ===========================================

"""
certain kinds of high-low data
temperature extreme guages
extreme waterlevel gauges 
but also published data without detailed records but recalling the range and sample size

J.D. Bowman from NIOSH suggests a non-rigorous but possibly useful estimation 
strategy for probability distributions in the style--or at least brashness--of 
Enrico Fermi who was famous for his seat-of-the-pants estimates.  It works like
this: assume the shape, e.g., lognormal, and guess its parameters crudely from
the range of n randomly sampled values from it (you don't need the actual data, 
just their range) using a transformation of the mean deviate range for such 
distributions. Alternatively, they can be estimated from upper and lower bounds 
on a specified percentage of random values drawn from the distribution. These
bounds could be established from expert elicitation.

This approach makes many assumptions, including that the distribution shape is 
known, that the range of sample values is known but the values themselves are 
not, that one can specify their coverage probability or the sample size, and 
that coverage is symmetric or sampling was random. The approach could work with 
many distribution shapes, although it's not clear it’ll work with 3-parameter 
shapes. In theory one could pool input from several experts although it's not 
clear how to aggregate or weight their respective opinions.


<< Copilot says:
    
    
- distribution‑agnostic
- exact quantile inversion
- distribution‑specific range correction

    
- If the expert gives you the range of a sample,
you should shrink it to the expected central range before fitting a distribution.
This is the only part of Bowman’s method that is not distribution‑specific — 
and therefore the only part that belongs in a unified engine.

✔ Why this matters
If you remove this correction, then:
- elicited ranges from data become too wide
- fitted distributions become too dispersed
- gamma and Weibull become badly biased
- normal/lognormal become inconsistent with sample‑based elicitation
So yes — the mean‑normal‑range idea is still essential.
It’s the one piece of Bowman’s method that generalizes cleanly.

⭐ 4. So what was “clever” about Bowman’s approach?
Not the two‑quantile idea — that’s standard.
The clever part was:
- noticing that experts often give sample extremes, not quantiles
- realizing that the expected range of a normal sample is computable
- using that to shrink the elicited range to a central interval
- and then mapping that to lognormal parameters via geometric mean/sd
He found a way to turn messy expert input into something mathematically tractable.
Your unified engine keeps the good part (range correction) and replaces the brittle part (distribution‑specific algebra) with something principled (quantile inversion).

⭐ 5. Where this leaves you
You now have:
- a clean, extensible, mathematically honest elicitation engine
- exact quantile matching for all distributions
- consistent behavior across normal/lognormal/Weibull/gamma/logistic
- the mean‑normal‑range correction preserved
- no more accidental 95% ranges
- no more gamma distributions shifted to the right
This is the architecture you were always heading toward.
If you want, I can help you:
- add 3‑quantile elicitation
- add bounded distributions (beta, triangular, uniform)
- add skew‑normal or generalized gamma
- add automatic initial‑guess heuristics
- add diagnostics (e.g., check for impossible quantile pairs)
Just tell me where you want to take it next.
>>

References
Bowman, J.D., S.A. Shulman, and S. Sivagenesan. 2010. Expanded Fermi estimates 
for variances and correlations from expert judgments. [manuscript and slide 
presentation]. NIOSH and University of Cincinnati.

Pelega, M., M.D. Normanda, J. Horowitzb and M.G. Corradinia. 2007. An expanded 
Fermi solution for microbial risk assessment. International Journal of Food 
Microbiology 113(1): 92-101. 
"""

from scipy.integrate import quad

def mean_normal_range(n):
    """Expected range of n independent N(0,1) samples."""
    phi = sps.norm.pdf
    Phi = sps.norm.cdf
    def log_f(x): return np.log(n) + np.log(phi(x)) + (n-1)*np.log(Phi(x))
    m = np.max(log_f(np.linspace(-8, 8, 2001)))
    def integrand(x): return x * np.exp(log_f(x) - m)
    I = quad(integrand, -8, 8, limit=200)[0]
    return 2 * I * np.exp(m)

def mean_range_lognormal(m, s, n, sims=100000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = qlognormal(u, m, s)
    return np.mean(x.max(axis=1) - x.min(axis=1))

def mean_range(ppf, n, sims=200000):  # distribution-specific mean sample range 
    """Monte Carlo estimate of E[max - min] for n IID samples."""
    if n <= 1: return 1
    x = ppf(np.random.rand(sims, n))
    return np.mean(x.max(axis=1) - x.min(axis=1))

_mean_range_cache = {}

def mean_range_ppf(ppf, theta, n, sims=200000):
    """Monte Carlo estimate of E[max - min] for n independent indentically 
    distributed samples from distribution ppf with parameters theta."""
    if n <= 1: return 0.0
    key = (id(ppf), tuple(theta), n)
    if key in _mean_range_cache: return _mean_range_cache[key]
    u = np.random.rand(sims, n)
    x = ppf(u, *theta)
    mr = np.mean(x.max(axis=1) - x.min(axis=1))
    _mean_range_cache[key] = mr
    return mr

def fermi_normal(lower, upper, n):
    return normal((lower+upper)/2, (upper-lower)/mean_normal_range(n))

def fermi_lognormal(lower, upper, n, rigor=True, tol_mean=0.05, tol_range=0.10, m_grid=40, s_grid=40, sims_range=50000):
    """P-box of all lognormal distributions for which the mean is the average of
    upper and lower, and the mean range of n random samples is their difference."""
    m_target = (lower + upper) / 2.0
    R_target = upper - lower
    if n > 1 and R_target > 0:                 # crude normal-based scale guess
        mr_norm = 3.077                    # approx for n≈10; refine if desired
        s0 = R_target / mr_norm
    else: s0 = max(R_target, 1.0)
    m_vals = np.linspace(0.5*m_target, 1.5*m_target, m_grid)   # parameter grid
    s_vals = np.linspace(0.3*s0, 3.0*s0, s_grid)               # parameter grid
    admissible = []
    exterior = []
    for m in m_vals:
        for s in s_vals:
            if m <= 0 or s <= 0: continue
            mean_ok = abs(m - m_target)/m_target <= tol_mean
            R_hat = mean_range_lognormal(m, s, n, sims=sims_range)
            range_ok = abs(R_hat - R_target)/R_target <= tol_range
            if mean_ok and range_ok: admissible.append((m, s))
            else: exterior.append((m, s))
    if not admissible: raise RuntimeError("no admissible parameters found")
    admissible = np.array(admissible)
    exterior   = np.array(exterior)
    boundary = set()                   # build a fast lookup for grid adjacency
    m_step = m_vals[1] - m_vals[0]
    s_step = s_vals[1] - s_vals[0]
    admissible_set = set((round(m,12), round(s,12)) for (m,s) in admissible)
    exterior_set   = set((round(m,12), round(s,12)) for (m,s) in exterior)
    for (m, s) in admissible:
        boundary.add((m, s))                                 # admissible point
        if rigor:
            for dm in [-m_step, 0, m_step]:   # check neighbors in 8 directions
                for ds in [-s_step, 0, s_step]:
                    if dm == 0 and ds == 0:
                        continue
                    m2 = m + dm
                    s2 = s + ds
                    key = (round(m2,12), round(s2,12))
                    if key in exterior_set:
                        boundary.add((m2, s2))
    boundary = np.array(list(boundary))   # (m,s) pairs (admissible + exterior)
    m_list = boundary[:,0]
    s_list = boundary[:,1]
    m_interval = (m_list.min(), m_list.max())
    s_interval = (s_list.min(), s_list.max())
    pL = PbO.ii()   # left grid
    pU = PbO.jjj()  # right grid
    QL = np.full_like(pL, np.inf, dtype=float)
    QU = np.full_like(pU, -np.inf, dtype=float)
    for (m, s) in boundary:
        QL = np.minimum(QL, qlognormal(pL, m, s))     # lower quantile envelope 
        QU = np.maximum(QU, qlognormal(pU, m, s))     # upper quantile envelope
    return Pbox(QL, QU, *m_interval)  #, s_interval, boundary


"""
a = fermi_lognormal(0.2, 3.0, 10)
red(a)

def fermiF(f,a,b,n=None): #,pr=0.9):
    fab = fermi_normal(a,b,n) #,pr=pr)
    plot(fab)
    plt.title(str(n))
    plt.plot([a,a,b,b],[1,0,0,1],'xkcd:grey',ls=':')
    A,B = ends(support(fab))
    #pr = 0.9
    #alpha = (1-pr)/2
    #plt.plot([A,B,B,A],[1-alpha,1-alpha,alpha,alpha],'r',ls=':',lw=1)
    plt.show()
  

#init_splot(4,4,sharex=True,sharey=True)
fermiF(normal, a=0.2, b=3, n=1) 
fermiF(normal, a=0.2, b=3, n=2) 
fermiF(normal, a=0.2, b=3, n=3) 
fermiF(normal, a=0.2, b=3, n=4) 
fermiF(normal, a=0.2, b=3, n=5) 
fermiF(normal, a=0.2, b=3, n=6) 
fermiF(normal, a=0.2, b=3, n=7) 
fermiF(normal, a=0.2, b=3, n=8) 
fermiF(normal, a=0.2, b=3, n=9) 
fermiF(normal, a=0.2, b=3, n=10) 
fermiF(normal, a=0.2, b=3, n=25)
fermiF(normal, a=0.2, b=3, n=50)
fermiF(normal, a=0.2, b=3, n=100) 
fermiF(normal, a=0.2, b=3, n=250) 
fermiF(normal, a=0.2, b=3, n=500) 
fermiF(normal, a=0.2, b=3, n=1000) 


for n in [1,5,10,100,1000]: fermiF(normal, a=0.2, b=3, n=n) 
"""




    
"""

# confidence bands from range & n

fermi.norm.confband = function(x1, x2, n, pr=0.9, conf=0.95, bOt=0.005, tOp=0.995, long=500) {
  t = fermi.norm(x1,x2,n,pr)
  m = t[[1]]
  s = t[[2]]
  x = seq(qnorm(bOt,m,s),qnorm(tOp,m,s),length.out=long)
  if (conf != 0.95) stop('Cannot handle confidence level other than 95%')
  Dmax = approx.ksD95(n)
  fermip = pnorm(x,m,s)
  list(x=x,left=fermip+Dmax,right=fermip-Dmax,best=fermip)
  }

fermi.lnorm.confband = function(x1, x2, n, pr=0.9, conf=0.95, bOt=0.005, tOp=0.995, long=500) {
  t = fermi.lnorm(x1,x2,n,pr)
  m = t[[1]]
  s = t[[2]]
  x = seq(qlnorm(bOt,m,s),qlnorm(tOp,m,s),length.out=long)
  if (conf != 0.95) stop('Cannot handle confidence level other than 95%')
  Dmax = approx.ksD95(n)
  fermip = plnorm(x,m,s)
  list(x=x,left=fermip+Dmax,right=fermip-Dmax,best=fermip)
  }


par(mfrow=c(1,1))
a = fermi.lnorm.confband(.2, 2.3, 50)  # n = 50, lognormal
plot(a$x,a$best,type='l',xlim=c(0,10))
lines(a$x[a$left<=1],a$left[a$left<=1],type='l')
lines(a$x[0<=a$right],a$right[0<=a$right],type='l')

"""

"""
Naked expert elicitations of probabilities of rare events

There are basically two approaches to estimating probability without actual 
data: expert elicitation (i.e., guessing), and disaggregation into constituent 
components whose probabilities are easier to estimate (i.e., breaking into 
subproblems). When the latter approach is no longer workable, analysts must 
resort to the former and rely on expert opinion and estimation. But how should 
we characterize probabilities of events that are so rare that they have never 
been observed? By what principles can such characterizations be projected in 
probabilistic analyses? Sometimes elaborate elicitation strategies are employed 
to estimate rare-event probabilities, but the results are often expressed as 
probabilities with no indication about the uncertainty associated with the 
estimate. How might analysts model expert opinions about event probabilities of 
the form “1 in 10 million”, “about 1 in 1000”, or “it’s never seen in over 100 
years of observation”, so they can be used in calculations that account for 
rather than ignore epistemic uncertainty? Several strategies provide partial 
solutions, addressing significant digits, hedged expressions, estimating order 
of magnitude, precision overstatement bias, and uncertainty about the Bayesian 
prior. For instance, presumably the assertion that an event has a probability 
of 1 in 1,000 would include probability values as low as 0.5 in 1000, and as 
large as 1.5 in 1000. Linguistic analysis reveals a simple scheme to decode 
approximator words such as ‘about’, ‘around’, or ‘at least’ in natural-language 
expressions. Robust Bayes analysis can account for uncertainty about the prior. 
These strategies can be combined in a coherent probabilistic analysis that 
minimally captures the express epistemic uncertainty implied by common 
utterances from experts. The analysis is broadly acceptable under both Bayesian 
and frequentist interpretations of probability, and it distinguishes epistemic 
and aleatory uncertainties. The result represents a lower bound on the final 
uncertainty.
"""

"""
# Automatic intervals under significant-digit interpretation           

Rules for counting whether a digit is significant:               Example #    ± 
-All non-zero digits in a number are significant.                 64.2   3  0.1
-Zeros between two non-zero digits are significant.               1043   4    1
-Leading zeros are not significant; they're just place holders.   0.0021 2 1e-4
-Trailing zeros to the right of the decimal are significant.      52.00  4 0.01
-Trailing zeros in a number with a decimal point are significant. 390.   3    1
-Trailing zeros in a number with no decimal are not significant.  340    2   10
-All digits in a mantissa in scientific notation are significant. 1.02e4 3  100 
-Exact numbers such as definitions have infinitely many significant figures. 

These rules predate IEEE‑754 and remain standard in laboratory science. They
represent an important convention for recording measurements of various kinds, 
although they are not universally observed outside of science and engineering. 
Note, however, that they are not a comprehensive scheme for recording handling 
measurement imprecision. For example, it is possible for a measurement's 
imprecision to imply zero or even a negative number of significant digits, 
which can only be expressed by an interval with two explicit numbers rather 
than a single number with encoding or formatting.


"""

def format_sigdigs(x, D):
    """Format x with exactly D significant digits."""
    sign = "-" if x < 0 else ""
    x = abs(x)
    if x == 0: return sign + "0." + "0"*(D-1)
    fmt = f"{{:.{D}g}}"
    s = fmt.format(x)                        # get rounded value using g-format
    if "e" in s or "E" in s: # if scientific notation get mantissa and exponent 
        mantissa, exp = s.lower().split("e")
        exp = int(exp)
    else: mantissa, exp = s, None
    if exp is not None:                                  # try decimal notation
        # convert scientific notation to a decimal string
        digits = mantissa.replace(".", "")
        # Pad digits to D significant digits
        digits = digits + "0"*(D - len(digits))
        digits = digits[:D]
        decimal_pos = 1 + exp        # position of decimal point after shifting
        if 0 < decimal_pos < 20:                 # reasonable decimal expansion
            if decimal_pos >= D:                # all digits to left of decimal
                dec = digits + "0"*(decimal_pos - D)
                return sign + dec + "."
            else:                                        # insert decimal point
                dec = digits[:decimal_pos] + "." + digits[decimal_pos:]
                return sign + dec.rstrip("0").rstrip(".") + ("." if dec.endswith("0") else "")
        # Otherwise: scientific notation is appropriate    
    if "." in mantissa:    # if already has a decimal point, use trailing zeros
        whole, frac = mantissa.split(".")
        digits = whole + frac
        digits = digits + "0"*(D - len(digits))
        digits = digits[:D]
        if len(whole) >= D: return sign + whole[:D] + "."
        else: return sign + whole + "." + digits[len(whole):]  
    digits = mantissa + "0"*(D - len(mantissa))           # integer-like string
    digits = digits[:D]
    if len(digits) <= len(mantissa): return sign + digits + "."
    return sign + digits[:len(mantissa)] + "." + digits[len(mantissa):]

def format_decimals(a, d):           # format_decimals(100/3,3) yields '33.333'
    if d is None: d = IvO.digits
    if d is None: return f"{a}"
    fmt = f"{{:.{d}f}}"
    return f"{fmt.format(a)}"

# def sgnumber(user_input: str):     # number ± its significant-digit imprecision
#     user_input = user_input.strip().lower()
#     tens = '0'
#     if 'e' in user_input: mantissa, tens = user_input.split('e', 1)
#     else: mantissa = user_input
#     if '.' in mantissa: j = len(mantissa.split('.')[1])
#     #else: j = len(mantissa.split('0', 1)) - len(mantissa) + 1           
#     else: j = len(mantissa.rstrip('0')) - len(mantissa) 
#     pm = 10**(-j) * 10**int(tens) / 2
#     #print('input:',user_input,', mantissa:',mantissa, ', j:',j, ', tens:',tens, ', pm:',pm)
#     return([float(user_input)-pm, float(user_input)+pm])





def sigorder(s: str) -> int:
    """Return log10 order of the last significant digit of a significand string.
    The input s should represent the mantissa, i.e., lacking the E, e, and the 
    exponent parts of scientific notation. This function is not related to the 
    magnitude of a value, but rather to its implicit imprecision:
      Sigorder  Examples
        -3         0.001     0.009      1.234    87.164
        -2         0.01      0.09      34.50      0.02
        -1         0.1       0.9       34.5    4328.6  
         1         1         9     201000.        3
         2        10        90      14230        20. 
         3       100       900        700     62800
         4      1000      9000     201000   5959000"""
    s = str(s)
    n = len(s)
    i = n - 1
    if s.endswith('.'): return 1                   # Case: string ends with '.'
    while i >= 0 and s[i] != '.': i -= 1  # Scan backward to find decimal point
    if i >= 0: return -(n - 1 - i) # decimal point found
    # No decimal point: find last nonzero digit
    ch = ' '
    wh = 0
    i = n
    while i > 0 and ch in ('0', ' '):
        i -= 1
        wh += 1
        ch = s[i]
    if i == 0 and ch in ('0', ' '): return -1
    else: return wh

def lastsigfig(s: str) -> str:
    """Return the last significant digit (as a character) from a significand string."""
    s = str(s)
    n = len(s)
    i = n - 1
    if s.endswith('.'): return s[-2]   # If ends with '.', return previous char
    ch = s[i]
    # Scan backward until decimal point or start
    while i > 0:
        if s[i] == '.': return ch
        i -= 1
    # No decimal point: find last nonzero digit
    ch = ' '
    i = n
    while i > 0 and ch in ('0', ' '):
        i -= 1
        ch = s[i]
    if i == 0 and ch in ('0', ' '): return ' '
    else: return ch

def about(s=None, pbox=False, v=None, r=None, f=None):
    """Linguistic hedge characterizing uncertainty about an imprecise number.
    Specify the base value either as a string representation s (needed to infer 
    significant digits) or numerically by supplying v, r, and f:
       v: numeric magnitude,
       r: significance order, and
       f: whether the last significant digit is '5'.
    The limits of the hedge are returned a tuple unless pbox is True, in which 
    case the left bound, right bound, and an enclosing p-box of the hedge's 
    uncertainty are returned as a triple."""
    if v is None: v = float(s)
    if r is None: r = sigorder(s)
    if f is None: f = (lastsigfig(s) == '5')
    A = -0.2085;    B = 0.4285;   C = 0.2807;  D = 0.0940;      E = 0.0147;     
    Fp = -0.0640;   G = -0.0102;  H = 0.0404;  sigma = 0.5837;  r2 = 0.7412
    z = np.log10(v)
    L = (A + B*z + C*r + D*f + E*z*r + Fp*z*f + G*r*f + H*z*r*f)   # regression 
    w = 10**L
    a_low  = 10**z - w/2
    a_high = 10**z + w/2
    if not pbox: return (a_low, a_high)
    q = lognormal(10**(sigma**2/2), np.sqrt(10**(2*sigma**2) - 10**(sigma**2)))
    p = env(min(a_low, a_high) - q, q + max(a_low, a_high))
    return (a_low, a_high, p)


"""

def show_all_about(s):
    a = about(s, True)
    plot(a[2])   # P-box:  ~ ( range=[158,242], mean=[172,228], var=[2.6,3]) 
    plt.plot([a[0],a[0]], [0,1], 'xkcd:grey', ls=':')
    plt.plot([a[1],a[1]], [0,1], 'xkcd:grey', ls=':')
    plt.plot(about('200'), (0.5,0.5), 'b')

sgnumber('200')           # [150.0, 250.0]
about('200')              # 173.6951 226.3049
about(v=200, r=3, f=0)    # 173.6951 226.3049
show_all_about('200')

sgnumber('200.')          # [199.5, 200.5]
about('200.')    	      # 193.82 206.18
about(v=200, r=1, f=0)    # 193.8200 206.1800
about(v=200, r=0, f=0)    # 197.0046 202.9954
about(v=200, r=-1, f=0)   # 198.5481 201.4519
show_all_about('200.')

"""

def O(x): return env(x/10, x*10)

def orderof(x): return env(x/2, x*5)

# ============================================================
# Maximum entropy (and miscellaneous) distributions
# ============================================================

def antweiler_triangular(x):
    return triangular(Min=np.min(x), Mode=3*np.mean(x)-np.max(x)-np.min(x), Max=np.max(x))

def betapert(min, max, mode):
    mu = (min + max + 4*mode)/6
    if abs(mode - mu) < 1e-8: alpha1 = alpha2 = 3
    else:
        alpha1 = (mu - min)*(2*mode - min - max)/((mode - mu)*(max - min))
        alpha2 = alpha1*(max - mu)/(mu - min)
    return min + (max - min) * beta(alpha1, alpha2)

def MEminmax(min, max): return uniform(min, max)

def MEminmaxmean(min, max, mean): return sawinconrad(min, mean, max) #http://mathoverflow.net/questions/116667/whats-the-maximum-entropy-probability-distribution-given-bounds-a-b-and-mean, http://www.math.uconn.edu/~kconrad/blurbs/analysis/entropypost.pdf for discussion of this solution.

def MEmeansd(mean, sd): return normal(mean, sd)

def MEminmean(min, mean): return min + exponential(mean - min)

# MEmeangeomean omitted

# def MEdiscretemean(x, mu, steps=10, iterations=50):
#     # x is a list or numpy array
#     import numpy as np

#     def fixc(x, r):
#         return 1/np.sum(r**x)

#     r = br = 1
#     c = bc = fixc(x, r)
#     d = bd = (mu - np.sum((c * r**x) * x))**2

#     for j in range(1, steps+1):
#         step = 1/j
#         for i in range(iterations):
#             r = abs(br + (np.random.rand() - 0.5) * step)
#             c = fixc(x, r)
#             d = (mu - np.sum((c * r**x) * x))**2
#             if d < bd:
#                 br = r
#                 bc = c
#                 bd = d

#     w = bc * br**x
#     w = w / np.sum(w)

#     z = []
#     k = len(x)
#     for i in range(k):
#         z.extend([x[i]] * int(w[i] * MC.many))

#     if len(z) > MC.many:
#         z = z[:MC.many]
#     elif len(z) < MC.many:
#         extra = np.random.permutation(z)[:(MC.many - len(z))]
#         z = z + list(extra)

#     return mc(np.random.permutation(z))

# MEquantiles = quantiles

def MEdiscreteminmax(min, max): return np.minimum(trunc(uniform(min, max+1)), max)

def MEmeanvar(mean, var): return MEmeansd(mean, np.sqrt(var))

def MEminmaxmeansd(min, max, mean, sd): return beta1((mean - min)/(max - min), sd/(max - min)) * (max - min) + min

MEmmms = MEminmaxmeansd

def MEminmaxmeanvar(min, max, mean, var): return MEminmaxmeansd(min, max, mean, np.sqrt(var))

# ============================================================
# Method of matching moment distributions 
# ============================================================

def MMbernoulli(x): return bernoulli(np.mean(x))

def MMbeta(x): return beta1(np.mean(x), np.std(x,ddof=1))

def MMbetabinomial(n, x):
    x = np.asarray(x)
    m1 = np.mean(x)
    m2 = np.mean(x**2)
    d = n*(m2/m1 - m1 - 1) + m1
    return betabinomial(n, (n*m1 - m2)/d, (n - m1)*(n - m2/m1)/d)

def MMbinomial(x):
    a = np.mean(x)
    b = np.std(x,ddof=1)
    return binomial(round(a/(1 - b*b/a)), 1 - b*b/a)

def MMchisquared(x): return chisquared(round(np.mean(x)))

def MMexponential(x): return exponential(np.mean(x))

def MMF(x):
    w = 2/(1 - 1/np.mean(x))
    return F(round((2*w**3 - 4*w**2)/((w-2)**2 * (w-4) * np.std(x,ddof=1)**2 - 2*w**2)), round(w))

def MMgamma(x):
    a = np.mean(x)
    b = np.std(x,ddof=1)
    return gamma(b*b/a, (a/b)**2)

def MMgeometric(x): return geometric(1/(1 + np.mean(x)))

MMpascal = MMgeometric

def MMgumbel(x): return gumbel(np.mean(x) - 0.57721*np.std(x,ddof=1)*np.sqrt(6)/np.pi,np.std(x,ddof=1)*np.sqrt(6)/np.pi)

def MMlognormal(x): return lognormal(np.mean(x), np.std(x,ddof=1))

def MMlaplace(x): return laplace(np.mean(x), np.std(x,ddof=1)/np.sqrt(2))

def MMlogistic(x): return logistic(np.mean(x), np.std(x,ddof=1)*np.sqrt(3)/np.pi)

def MMloguniform(x): return loguniform1(np.mean(x), np.std(x,ddof=1))

def MMnormal(x): return normal(np.mean(x), np.std(x,ddof=1))

MMgaussian = MMnormal

def MMpareto(x):
    a = np.mean(x)
    b = np.std(x,ddof=1)
    return pareto(a/(1 + 1/np.sqrt(1 + a*a/b*b)), 1 + np.sqrt(1 + a*a/b*b))

def MMpoisson(x): return poisson(np.mean(x))

def MMpowerfunction(x):
    a = np.mean(x)
    b = np.std(x,ddof=1)
    return powerfunction(a/(1 - 1/np.sqrt(1 + (a/b)**2)), np.sqrt(1 + (a/b)**2) - 1)

def MMt(x):
    if 1 < np.std(x,ddof=1): return student(2/(1 - 1/np.std(x,ddof=1)**2))
    else: raise ValueError("Improper standard deviation for student distribution")

MMstudent = MMt

def MMuniform(x):
    a = np.mean(x)
    b = np.std(x,ddof=1)
    return uniform(a - np.sqrt(3)*b, a + np.sqrt(3)*b)

MMrectangular = MMuniform

def MMuniform1(w):
    mu1 = np.mean(w)
    mu2 = np.mean(w**2)
    m = np.sqrt(3*(mu2 - mu1*mu1))
    return uniform(mu1 - m, mu1 + m)

def MMtriangular(x, iters=100, dives=10):

    def skewness(x):
        m = np.mean(x)
        return np.sum((x - m)**3)/((len(x)-1)*np.std(x,ddof=1)**3)

    M = np.mean(x)
    V = np.var(x)
    S = skewness(x)

    a = aa = np.min(x)
    b = bb = np.max(x)
    c = cc = 3*M - a - b

    many = iters
    s1 = np.std(x,ddof=1)

    for k in range(dives):
        s1 = s2 = s3 = s1/2
        a = np.random.normal(aa, s1, many)
        b = np.random.normal(bb, s2, many)
        c = np.random.normal(cc, s3, many)

        m = (a + b + c)/3
        k2 = (a*a + b*b + c*c - a*b - a*c - b*c)
        v = k2/18
        s = (np.sqrt(2)*(a+b-2*c)*(2*a-b-c)*(a-2*b+c)) / (5 * k2**1.5)

        d = (M - m)**2 + (V - v)**2 + (S - s)**2
        i = np.argmin(d)

        aa = a[i]
        bb = b[i]
        cc = c[i]
    return triangular(a, c, b)

##########################################################################
# Maximum likelihood distributions
##########################################################################

def MLbernoulli(x): return bernoulli(np.mean(x))

MLnormal = MMnormal
MLgaussian = MMnormal
MLexponential = MMexponential
MLpoisson = MMpoisson
MLgeometric = MMgeometric
MLpascal = MMgeometric

def MLuniform(x): return uniform(np.min(x), np.max(x))

MLrectangular = MLuniform

def MLpareto(x):
    return pareto(np.min(x), len(x)/np.sum(np.log(x) - np.log(np.min(x))))

def MLlaplace(x):
    med = np.median(x)
    return laplace(med, np.sum(np.abs(x - med))/len(x))

def MLlognormal_(x):
    n = len(x)
    mu = np.sum(np.log(x))/n
    return lognormal(meanlog=mu, stdlog=np.sum((np.log(x)-mu)**2)/n)

def MLlognormal(x):
    return np.exp(MLnormal(np.log(x)))

def MLweibull(x, shapeinterval=(0.001, 500)):
    from scipy.optimize import root_scalar
    def f(k): return np.sum(x**k*np.log(x))/np.sum(x**k)-np.sum(np.log(x))/len(x)-1/k
    sol = root_scalar(f, bracket=shapeinterval)
    k = sol.root
    el = np.exp(np.log(np.sum(x**k)/len(x))/k)
    return weibull(scale=el, shape=k)

# ============================================================
# Bayesian distributions and robust Bayes p-boxes
# ============================================================

def km_bayes(k,m):
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given k successes and m failures 
  # randomly observed in k + m independent Bernoulli trials 
  if ((k < 0)  or (m < 0)): raise ValueError('Improper arguments to function km')
  return beta(k+0.5, m+0.5)

def KN_bayes(k,n):
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given only k successes out 
  # of n randomly observed independent Bernoulli trials 
  if ((k < 0)  or (n < k)): raise ValueError('Improper arguments to function KN')
  return beta(k+0.5 ,n-k+0.5)

def count_bayes(n, alpha=0, beta=0, m=None, s=None):
    """Bayesian characterization of uncertainty in count data.

    Counts are assumed to follow a poisson(λ) model, with a gamma prior
    on λ having shape `alpha` and rate `beta`. The posterior for λ is
    gamma(shape = alpha + sum(n), rate = beta + len(n)).

    This posterior induces a negative binomial posterior predictive
    distribution for future counts. You may use the objective prior
    (alpha=0, beta=0), a similar improper prior (alpha=1, beta=0), or
    the Jeffreys prior (alpha=0.5, beta=0). Alternatively, you may
    specify the prior for λ in terms of its mean `m` and standard
    deviation `s`.
    
    Usage:
    count_bayes(123)                 # Pbox(range=[79,176],mean=123)
    count_bayes([123,145])           # Pbox(range=[93,181],mean=134)
    count_bayes(I('100'))            # Pbox(range=[23,208],mean=[50,150])
    count_bayes([I('100'),I('100')]) # Pbox(range=[26,199],mean=[50,150])
    """    
    #count_bayes([123,145],alpha=I(0.25, 0.36), beta=I(0.0025, 0.003)) # Pbox(range=[93, 181], mean=[133.92, 134.01])
    #count_bayes([123,145],alpha=I(25, 36), beta=I(0.25, 0.3)) # Pbox(range=[88.0, 181], mean=[127.39, 135.11])    
    n = np.atleast_1d(np.asarray(n, dtype=object))
    k = len(n)
    if m is not None and s is not None: 
        alpha = (m/s)**2
        beta = m/s**2
    # nL = np.sum([left(x)  if isinstance(x, Interval) else x for x in n])
    # nU = np.sum([right(x) if isinstance(x, Interval) else x for x in n])
    nL = np.sum([left(x) for x in n])
    nU = np.sum([right(x) for x in n])
    p = 1 / ( 1 + 1/( beta + len(n)))  # p = (beta+len(n)) / (beta+len(n) + 1)
    return negativebinomial(Interval(alpha + nL, alpha + nU), p)

##########################################################################
# Additional compound and conjugate distribution constructors for Bayesian inference
##########################################################################

# Generally, the default values for the prior hyperparameters are 
# intended to yield the uninformative prior.  When there is no 
# data, the posterior returned is just the prior distribution.  
# As the sample size increases, the posterior tends to the data 
# manifested as the likelihood.  Typically, as the prior hyper-
# parameters (a,b) increase, the posterior grows away from the 
# likelihood and towards the prior.  Generally, the default values 
# for the prior hyperparameters yield the uninformative case. 

# compound distributions already defined
# BB, betabinomial, gammaexponential, NB, negativebinomial, poissonbinomial

def BCbernoulli(x, a=0.5, b=0.5, only=True):
    """BCbernoulli uses conjugacy to estimate the distribution of the next 
    binary outcome from 0s and 1s.
    
    # The data x is an array of zeros and ones (failures and successes).
    # Inputs a and b are the parameters of the prior beta distribution.
    # See also the km_bayes( ) and KN_bayes( ) functions.
    """
    s = sum(x)
    n = len(x)
    lk = beta(s + 1, n - s + 1)
    pp = betabinomial(1, s + a, n - s + b)
    if only: return pp
    else: return {
            "pr": beta(a, b),
            "po": beta(s + a, n - s + b),
            "pp": pp,
            "lk": lk   }


"""
km_bayes(5,15)      # Pbox(range=[0.0, 1.0], mean=0.2619047619047619)
b=BCbernoulli([0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],only=False)
for k,v in b.items(): print(k,v)
                    #pr Pbox(range=[0.0, 1.0], mean=0.5)
                    #po Pbox(range=[0.0, 1.0], mean=0.738095238095238)
                    #pp Pbox(range=[0.0, 1.0], mean=0.738095238095238)
                    #lk Pbox(range=[0.0, 1.0], mean=0.7272727272727273)
0.2619+0.73809   # 0.99999
"""

##############################################################################
# (All the plotting examples preserved exactly as comments)
# par(mfrow=c(3,4))
# options(digits=3)
# x = c(0,0,0,1,1)
#   BCbernoulli(x,only=FALSE) #Jeffreys is the default
#   b = BCbernoulli(x,0,0)  # Haldane prior
#   b; title(paste('Haldane', mean(b)))
#   abline(h=1-sum(x)/length(x),col='green')
#   b = BCbernoulli(x)  # Jeffreys prior
#   b; title(paste('Jeffreys prior', mean(b)))
#   b = BCbernoulli(x,1,1)  # Bayes-Laplace prior
#   b; title(paste('Bayes-Laplace', mean(b)))
#   b = BCbernoulli(x,2,2)  # Walley prior
#   b; title(paste('Walley', mean(b)))
# ...
##############################################################################


def BCbinomial(N, k, n, a=0.5, b=0.5, only=True):
    """BCbinomial uses conjugacy to estimate a probability from k successes out of N trials.
    The input k is the count of successes, and n is a corresponding number of trials.
    Both k and n may be arrays, but they must have the same length.
    Inputs a and b are the parameters of the prior beta.
    N is the number of trials to use for the posterior predictive distribution.
    See also the km_bayes( ) and KN_bayes( ) functions.
    """
    s = sum(k)
    sn = sum(n)
    lk = beta(s + 1, sn - s + 1)
    pp = betabinomial(N, s + a, sn - s + b)
    if only: return pp
    else: return {
            "pr": beta(a, b),
            "po": beta(s + a, sn - s + b),
            "pp": pp,
            "lk": lk   }


##############################################################################
# (All the long demonstration blocks preserved as comments)
# https://stats.stackexchange.com/questions/512148/beta-binomial-vs-updating-a-prior-beta-distribution
# ...
##############################################################################


def BCpoisson(x, a=0, b=0, r=None, only=True):
    """
    # the default hyperparameters seem to be the uninformative case
    """
    s = sum(x)
    n = len(x)

    lk = gamma(shape=s + 1, rate=n)
    pr = gamma(shape=a, rate=b)
    po = gamma(shape=a + s, rate=b + n)

    # posterior predictive (Negative Binomial form)
    pp = negativebinomial(a + s, 1 - 1 / (1 + b + n))

    if only:
        return pp
    else:
        return {
            "pr": pr,
            "po": po,
            "pp": pp,
            "lk": lk
        }


##############################################################################
# doit = function(x,a,b) {
#   pl(0,18)
#   points(x,rep(-0.017,length(x)),col='red')
#   bc = BCpoisson(x,a,b,only=FALSE)
#   blue(bc$pr)
#   gray(bc$po)
#   green(bc$pp)
#   edf(bc$ppr)
# }
# ...
##############################################################################

"""
def BCgeometric(x, a=0, b=0, r=runif(MC.many), only=True):
    # the default values for the hyperparameters a and b are the uninformative case,
    # although they will make the function crash if x=NULL
    s = sum(x)
    n = len(x)
    lk = beta(n + 1, s + 1)
    pr = beta(a, b)
    po = beta(a + n, b + s)
    por = rbeta(r_(r), a + n, b + s)
    pp = mc(rgeom(MC.many, por))   # do we know an analytical formula?
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp, "lk": lk}

##############################################################################
# doit = function(x,a,b) {
#   pl(0,18)
#   points(x,rep(-0.017,length(x)),col='red')
#   bc = BCgeometric(x,a,b,only=FALSE)
#   blue(bc$pr)
#   gray(bc$po)
#   green(bc$pp)
# }
# par(mfcol=c(4,2))
# a = 1/2; b = 1/2
# doit(rgeom(3000,0.5),a,b)   # should be the same as geometric(0.5)
# doit(rgeom(10,0.5),a,b)
# doit(c(3,14,1),a,b)
# doit(c(),a,b)                    # the posterior should equal the prior
# a = 1; b = 10
# doit(rgeom(3000,0.5),a,b)
# doit(rgeom(10,0.5),a,b)
# doit(c(3,14,1),a,b)
# doit(c(),a,b)                    # the posterior should equal the prior
##############################################################################


def BCuniform_knownmin(x, A, a=None, b=None, r=runif(MC.many), only=True):
    #  x_i ~ uniform(A,theta), that is, from a uniform distribution whose minimum is A
    #  and whose maximum needs to be established
    if a is None: a = A
    if b is None: b = A + 1
    lk = pareto(max(x), len(x) + 1)
    pr = pareto(a, b)
    po = pareto(max(a, max(x)), b + len(x))  # Masatoshi says max(x) is m
    por = qpareto(r_(r), max(a, max(x)), b + len(x))  # Masatoshi says max(x) is m
    pp = mc(runif(MC.many, A, por))
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp, "lk": lk}

# Alias to match R's BCuniform.knownmin = BCuniform
BCuniform = BCuniform_knownmin

# x = runif(25,5,13)
# bc = BCuniform(x,A=5,a=1,b=1,only=FALSE)
# bc


def BCnegativebinomial(x, R, a=0, b=0, r=runif(MC.many), only=True):
    # the default hyperparameters are the uninformative case
    s = sum(x)
    n = len(x)
    lk = beta(R * n + 1, s + 1)
    pr = beta(a, b)
    po = beta(a + R * n, b + s)
    por = rbeta(r_(r), a + R * n, b + s)
    pp = mc(rnbinom(MC.many, R, por))   # do we know an analytical formula?
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp, "lk": lk}


def BCnormal_knownsigma(x, sigma, m0=None, s0=None, r=runif(MC.many), only=True):
    # the default hyperparameters seem to be uninformative
    # (using s0=sd(x)/sqrt(length(x)) makes the posterior differ from the likelihood more strongly)
    # increasing s0 makes the prior more uninformative, unlike the typical behaviour of (a,b)
    # hyperparameters in other functions
    if m0 is None: m0 = mean(x)
    if s0 is None: s0 = sd(x)
    s = sum(x)
    n = len(x)
    pr = normal(m0, s0)
    lk = normal(s / n, sigma / sqrt(n))
    mprime = (m0 / s0**2 + s / sigma**2) / (1 / s0**2 + n / sigma**2)
    sprime = 1 / sqrt(1 / s0**2 + n / sigma**2)
    if abs(s0 < 1e-20):
        mprime = m0
        sprime = 0
    po = normal(mprime, sprime)
    por = rnorm(r_(r), mprime, sprime)
    pp = normal(mprime, sqrt(sprime**2 + sigma**2))   # pp = mc(rnorm(MC$many,por,sigma))
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp, "lk": lk}


def BCnormal_knownmu(x, mu, a=0, b=0, r=runif(MC.many), only=True):
    # the hyperparameter defaults are the uninformative case
    # (so the pr won't be defined as it would be improper)
    n = len(x)
    s = sum((xi - mu)**2 for xi in x)
    pr = sqrt(inversegamma(shape=a, rate=b))  # we parameterize N with sd, not var
    po = sqrt(inversegamma(shape=a + n/2, rate=b + s/2))
    por = sqrt(1 / rgamma(MC.many, shape=a + n/2, rate=b + s/2))
    pp = mc(rnorm(MC.many, mu, por))
    # ppa = mean(x)+sd(x)*sqrt(1+1/n)*student(n-1)
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp}

##############################################################################
# x = rnorm(1000,15,2)
# bc = BCnormal.knownmu(x,15,5,5)
# bc
# edf(x)
##############################################################################


def rnormgamma(n, mu, lambda_, alpha, beta):
    # normal-gamma deviates:
    # (1) Sample tau from a gamma distribution with parameters alpha and beta,
    # (2) Sample x from a normal distribution with mean mu and variance 1/(lambda * tau)
    # E(x) = mu;  E(tau) = alpha/beta
    if hasattr(n, "__len__"):
        n = len(n)
    tau = rgamma(n, alpha, beta)
    x = rnorm(n, mu, sqrt(1 / (lambda_ * tau)))  # tau and x are NOT independent
    return {"x": x, "tau": tau}


def seepr(m, l, a, b):
    pr = rnormgamma(MC.many, m, l, a, b)
    prx = mc(pr["x"])
    prt = mc(pr["tau"])
    pl(min(left(prx), left(prt)), max(right(prx), right(prt)))
    blue(prx)
    cyan(prt)
    title(f"{m} {l} {a} {b}")


def BCnormal(x, mu0=None, lambda0=1, alpha0=1, beta0=1, only=True):
    # x_i ~ N(mu, 1/tau)
    # The selection of the prior for a normal involves choosing values for 4 parameters.
    # mu0 is your guess about the mean of the normal data, and second is related to the
    # dispersion of this estimate about the mean.  You can use the seepr( ) function to
    # visualize the priors you select for BCnormal.  In practice, setting mu0 to mean(x)
    # and the other values to one seems to often give reasonable results, so they have
    # been specified as defaults, but I'm sure this violates some crucial Bayesian stricture.
    if mu0 is None: mu0 = mean(x)
    pr = rnormgamma(MC.many, mu0, lambda0, alpha0, beta0)
    n = len(x)
    xbar = mean(x)
    s = sd(x)
    po = rnormgamma(
        MC.many,
        (lambda0 * mu0 + n * xbar) / (lambda0 + n),
        lambda0 + n,
        alpha0 + n/2,
        beta0 + (n * s + (lambda0 * n * (xbar - mu0)**2) / (lambda0 + n)) / 2
    )
    pp = mc(rnorm(MC.many, po["x"], 1 / po["tau"]))
    pr = {"x": mc(pr["x"]), "tau": mc(pr["tau"])}
    po = {"x": mc(po["x"]), "tau": mc(po["tau"])}
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp}


def BCexponential(x, a=0, b=0, r=runif(MC.many), only=True):
    # the prior and posterior estimate the MEAN of the exponential
    # which is the reciprocal of its RATE parameter used by rexp()
    # the default hyperparameters are the uninformative case
    sm = "exponential(theta)"
    s = sum(x)
    n = len(x)
    lk = gamma(shape=n + 1, rate=s)     # reciprocated in the returned list
    pr = gamma(shape=a, rate=b)         # reciprocated in the returned list
    po = gamma(shape=a + n, rate=b + s) # reciprocated in the returned list
    por = rgamma(r_(r), shape=a + n, rate=b + s)
    pp = mc(rexp(MC.many, por))
    if only:
        return pp
    else:
        return {"pr": 1/pr, "po": 1/po, "pp": pp, "lk": 1/lk}


def BCpareto_knownmin(x, xm, a=0, b=0, r=runif(MC.many), only=True):
    # default hyperparameters seem to be the uninformative case
    s = sum(log(xi / xm) for xi in x)
    n = len(x)
    sm = "pareto(xm,theta)"
    lk = gamma(shape=n, rate=s)  # this is just a guess; prolly wrong
    pr = gamma(shape=a, rate=b)
    po = gamma(shape=a + n, rate=b + s)
    por = rgamma(r_(r), shape=a + n, rate=b + s)
    pp = mc(qpareto(r_(r), xm, por))
    if only:
        return pp
    else:
        return {"pr": pr, "po": po, "pp": pp, "lk": lk}

# Alias
BCpareto = BCpareto_knownmin

"""



# ============================================================
# Confidence boxes (c-boxes)
# ============================================================

def KN(k, n):
    if left(k) < 0 or right(n) < right(k): raise ValueError("Improper argument to KN")
    if straddles(k) and straddles(n): return Pbox(0,1)
    q = [lambda i, k, n: qbeta(i, k, n-k+1), lambda j, k, n: qbeta(j, k+1, n-k)]
    return qpbox(q, PbO.iii(), PbO.jjj(), Interval(k/(n+1), (k+1)/(n+1)), k, n)
    #return env(qpbox(qbeta, PbO.iii(), PbO.jjj(), k/(n+1), k, n-k+1), qpbox(qbeta, PbO.iii(), PbO.jjj(), (k+1)/(n+1), k+1, n-k))

def FKN(k, n): # binomial rate inference for trials designed with a fixed-K stopping rule;  this is negative-binomial sampling
    if left(k) < 0 or right(n) < right(k): raise ValueError("Improper argument to FKN")
    if straddles(k) and straddles(n): return Pbox(0,1)
    q = [lambda i, k, n: qbeta(i, k, n-k+1), lambda j, k, n: qbeta(j, k, n-k)]
    M = Interval(k/(n+1), Interval(0,1) if n==0 else k/n)
    return qpbox(q, PbO.iii(), PbO.jjj(), M, k, n)
    #return env(qpbox(qbeta, PbO.iii(), PbO.jjj(), k/(n+1), k, n-k+1), qpbox(qbeta, PbO.iii(), PbO.jjj(), Interval(0,1) if n==0 else k/n, k, n-k))

def km(k, m):
    if left(k) < 0 or left(m) < 0: raise ValueError("Improper argument to km")
    #if (is.pbox(k) || is.pbox(m)) return(uchenna(pbox(k), pbox(m)))
    if straddles(k) and straddles(m): return Pbox(0,1)
    q = [lambda i, k, m: qbeta(i, k, m+1), lambda j, k, m: qbeta(j, k+1, m)]
    M = Interval(0 if k == 0 else 1/(1 + (m+1)/k), 1/(1 + m/(k+1)))
    return qpbox(q, PbO.iii(), PbO.jjj(), M, k, m)
    #return env(qpbox(qbeta, PbO.iii(), PbO.jjj(), 0 if k==0 else 1 /(1+(m+1)/k), k, m+1), qpbox(qbeta, PbO.iii(), PbO.jjj(), 1 /(1+m/(k+1)), k+1, m))

# x[i] ~ Bernoulli(p), x[i] is either 0 or 1

def CBbernoulli(x):
    n = len(x)
    k = np.sum(x)
    return env(bernoulli(k/(n+1)), bernoulli((k+1)/(n+1)))

def CBbernoulli_p(x):
    n = len(x)
    k = np.sum(x)
    return env(beta(k, n-k+1), beta(k+1, n-k))

# x[i] ~ binomial(N, p), for known N, x[i] is a nonnegative integer ≤ N

def CBbinomial(N, x):
    n = len(x)
    k = np.sum(x)
    return env(betabinomial(N, k, n*N - k + 1),
               betabinomial(N, k+1, n*N - k))

def CBbinomial_p(N, x):
    n = len(x)
    k = np.sum(x)
    return env(beta(k, n*N - k + 1),
               beta(k+1, n*N - k))

# x[i] ~ binomial(N, p), for unknown N, x[i] is a nonnegative integer

def CBbinomialnp(x):
    raise NotImplementedError("see https://sites.google.com/site/cboxbinomialnp/")

def CBbinomialnp_n(x):
    raise NotImplementedError("see https://sites.google.com/site/cboxbinomialnp/")

def CBbinomialnp_p(x):
    raise NotImplementedError("see https://sites.google.com/site/cboxbinomialnp/")

# x[i] ~ Poisson(mean), x[i] is a nonnegative integer

def CBpoisson(x):
    n = len(x)
    k = np.sum(x)
    return env(
        negativebinomial(size=k,   prob=1 - 1/(n+1)),
        negativebinomial(size=k+1, prob=1 - 1/(n+1)))

def CBpoisson_mean(x):
    n = len(x)
    k = np.sum(x)
    return env( gamma(shape=k,   rate=n), gamma(shape=k+1, rate=n))

# x[i] ~ exponential(mean)

def CBexponential(x):
    return gammaexponential(shape=len(x), rate=np.sum(x))

def CBexponential_lambda(x):
    return gamma(shape=len(x), rate=np.sum(x))

def CBexponential_mean(x):
    return 1 / zbuff(CBexponential_lambda(x))

# x[i] ~ normal(mu, sigma)

def CBnormal(x):
    n = len(x)
    return np.mean(x) + np.std(x, ddof=1) * student(n-1) * np.sqrt(1 + 1/n)

def CBnormal_mu(x):
    n = len(x)
    return np.mean(x) + np.std(x, ddof=1) * student(n-1) / np.sqrt(n)

def CBnormal_sigma(x):
    n = len(x)
    return np.sqrt(np.var(x, ddof=1) * (n-1) * inversechisquared(n-1))

# x[i] ~ lognormal(mu, sigma), log(x[i]) ~ normal(mu, sigma)

def CBlognormal(x):
    lx = np.log(x)
    n = len(x)
    return np.exp(np.mean(lx) + np.std(lx, ddof=1) * student(n-1) * np.sqrt(1 + 1/n))

def CBlognormal_mu(x):
    lx = np.log(x)
    n = len(x)
    return np.mean(lx) + np.std(lx, ddof=1) * student(n-1) / np.sqrt(n)

def CBlognormal_sigma(x):
    lx = np.log(x)
    n = len(x)
    return np.sqrt(np.var(lx, ddof=1) * (n-1) * inversechisquared(n-1))

# x[i] ~ uniform(midpoint, width)  (Monte Carlo)

def CBuniform(x):
    n = len(x)
    r = np.max(x) - np.min(x)
    w = (r / np.random.beta(n-1, 2, size=10_000)) / 2
    m = (np.max(x) - w) + (2*w - r) * np.random.uniform(0, 1, size=10_000)
    samples = np.random.uniform(m - w, m + w)
    return histogram(samples, conf=0)

def CBuniform_midpoint(x):
    r = np.max(x) - np.min(x)
    w = r / np.random.beta(len(x)-1, 2, size=10_000)
    m = (np.max(x) - w/2) + (w - r) * np.random.uniform(0, 1, size=10_000)
    return histogram(m, conf=0)

def CBuniform_width(x):
    r = np.max(x) - np.min(x)
    w = r / np.random.beta(len(x)-1, 2, size=10_000)
    return histogram(w, conf=0)

# x[i] ~ uniform(minimum, maximum)

def CBuniform_minimum(x):
    n = len(x)
    r = np.max(x) - np.min(x)
    w = r / np.random.beta(n-1, 2, size=10_000)
    m = (np.max(x) - w/2) + (w - r) * np.random.uniform(0, 1, size=10_000)
    return histogram(m - w/2, conf=0)

def CBuniform_maximum(x):
    n = len(x)
    r = np.max(x) - np.min(x)
    w = r / np.random.beta(n-1, 2, size=10_000)
    m = (np.max(x) - w/2) + (w - r) * np.random.uniform(0, 1, size=10_000)
    return histogram(m + w/2, conf=0)

# x[i] ~ F, a continuous but unknown distribution

def CBnonparametric(x):
    return env(
        histogram(np.concatenate([x, [np.inf]]), conf=0),
        histogram(np.concatenate([x, [-np.inf]]), conf=0))

# x1[i] ~ normal(mu1, sigma1), x2[j] ~ normal(mu2, sigma2), independent

def CBnormal_meandifference(x1, x2):
    return CBnormal_mu(x2) - CBnormal_mu(x1)

# x[i] = Y + error[i], error[j] ~ F unknown, Y fixed

def CBnonparametric_deconvolution(y, error):
    z = []
    for e in error:
        z.extend(y - e)
    z = np.sort(np.array(z))
    Q = Get_Q(len(y), len(error))
    w = Q / np.sum(Q)
    return env(
        mixture(z, w),
        mixture(np.concatenate([z[1:], [np.inf]]), w))

def ci(b, c=0.95, alpha=None, beta=None):
    if alpha is None or beta is None:
        alpha = (1 - c) / 2
        beta = 1 - (1 - c) / 2
    return Interval(left(cut(b, alpha)), right(cut(b, beta)))

confidenceinterval = ci

# ============================================================
# Plotting
# ============================================================

"""
The plotting functions draw on the current axes if they exist, and create them 
only when needed. They never force figure creation or flushing.  Use plt.show() 
to flush a graph. If the caller does not specify the ax argument, it defaults 
to plt, i.e., matplotlib.pyplot's alias. You can recover axes after plotting 
with ax=plt.gca() if the figure is still open.
"""

def seq(start,stop,many=101) :
    return [start + i*(stop-start)/(many-1) for i in range(many)]

def plot_interval(x, ax=None, form=None, color="black", **kwargs): # kwargs might line lw, linewidth
    if (ax is None) or (ax == '') : ax = plt
    lo, hi = x.left(), x.right() 
    if form is None:  form = IvO.form
    if form=="e":  # ellipse
        tt = seq(0,3.14159,25)
        a = (hi - lo) / 2
        c = (lo + hi) / 2
        x = [c + a * np.cos(t) for t in tt]
        y = [np.sin(t) for t in tt]
        ax.plot(x,y, color=color, **kwargs)  
    elif form == "t":  # triangle
        xs = [lo, (lo + hi) / 2, hi]
        ys = [0, 1, 0]
        ax.plot(xs, ys, color=color, **kwargs)
    else: ax.plot([lo, lo, hi, hi], [0, 1, 1, 0], color=color, **kwargs) # box (default)

def plot_pbox(pb, ax=None, cumulative=None, fmt='', **kwargs):
    if (ax is None) or (ax == '') : ax = plt
    if cumulative is None: cumulative = PbO.cumulative
    if fmt == '' : fmt = ('r-','k-')         # don't just make this the default
    if isinstance(fmt, str) : fmt = (fmt,fmt)
    n = pb.steps()
    u = pb.leftside()
    d = pb.rightside()
    p_u = [i / n for i in range(0, n)] + [1, 1]
    x_u = list(u) + [u[-1], d[-1]]
    if not cumulative: p_u = [1 - v for v in p_u]
    ax.step(x_u, p_u, fmt[0], where="pre", **kwargs)
    p_d = [0, 0] + [i / n for i in range(1, n + 1)]
    x_d = [u[0], d[0]] + list(d)
    if not cumulative: p_d = [1 - v for v in p_d]
    ax.step(x_d, p_d, fmt[1], where="post", **kwargs)

def ecdf(d) :
    d = np.array(d)
    N = d.size
    pp= np.concatenate((np.arange(N),np.arange(1,N+1)))/N
    dd = np.concatenate((d,d))
    dd.sort()
    pp.sort()
    return dd,pp

def edf(d, ax=None, fmt=None, **kwargs) :
    if (ax is None) or (ax == '') : ax = plt
    z,p = ecdf(d)
    ax.plot(z,p, fmt, **kwargs)

def plot(x, ax=None, **kwargs):
    if is_pbox(x):     return plot_pbox(x, ax=ax, **kwargs)
    if is_interval(x): return plot_interval(x, ax=ax, **kwargs)
    return edf(x, ax=ax, **kwargs) 
    raise TypeError("don't know to plot "+str(type(x)))

def red(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'r', **kwargs)
def blue(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'b', **kwargs)
def gray(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'xkcd:grey', **kwargs)
def cyan(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'c', **kwargs)
def green(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'g', **kwargs)
def black(x, ax=None, fmt=None, **kwargs):  plot(x, ax=None, fmt=fmt if fmt is not None else 'k', **kwargs)

# ============================================================
# Subplotting     
# ============================================================

# Call init_splot(nr,nc) to make nr-by-nc subplots, then use splot() to advance 
# to a new plot or splot(i) to pick the ith graph, on which lines() will draw.

_splot_index = 0       # global state for the sequential graphing plot iterator
_splot_axes = None
_splot_current = None

def init_splot(nrows=1, ncols=1, **kwargs): # calls subplots and pass it kwargs, e.g. sharex=True, figsize=(4,5)
    global _splot_axes, _splot_index, _splot_current
    fig, ax = plt.subplots(nrows, ncols, **kwargs)
    if isinstance(ax, np.ndarray):  _splot_axes = ax.flatten()
    else:                           _splot_axes = np.array([ax])
    _splot_index = 0
    _splot_current = None
    return fig, ax

def splot_open(abend=False): 
    if _splot_current is None: 
        if abend: raise RuntimeError("no current plot")
        else: return False
    else: return True
    
def splot(x=None, y=None, i=None, **kwargs):  
    '''Updates the current subplot grid, or selects the ith one, and plots an 
    uncertain number x, or the (x,y) pair.'''
    global _splot_axes, _splot_index, _splot_current
    if _splot_axes is None: raise RuntimeError("call init_splot() before splot()")
    if i is not None:
        if not (0 <= i < len(_splot_axes)): raise IndexError(f"subplot index {i} out of range")
        _splot_current = _splot_axes[i]
    else:
        if _splot_index >= len(_splot_axes): raise RuntimeError("all subplots have been used")
        _splot_current = _splot_axes[_splot_index]
        _splot_index += 1
    if x is None: return _splot_current
    if is_interval(x) or is_pbox(x): return plot(x, ax=_splot_current, **kwargs)
    if y is not None: return _splot_current.plot(x, y, **kwargs)
    return _splot_current.plot(x, **kwargs)

def lines(*args, **kwargs):
    if _splot_current is None: raise RuntimeError("no current plot")
    if len(args)==1 and (is_interval(args[0]) or is_pbox(args[0])): return plot(args[0], ax=_splot_current, **kwargs)
    return _splot_current.plot(*args, **kwargs)      # ordinary matplotlib line

def title(s, **kwargs): _splot_current.set_title(s, **kwargs)

def xlabel(s, **kwargs): _splot_current.set_xlabel(s, **kwargs)

def ylabel(s, **kwargs): _splot_current.set_ylabel(s, **kwargs)

def textxy(x, y, s, **kwargs): _splot_current.text(x, y, s, **kwargs)

text = textxy  # alias

# ============================================================
# Generic structural accessor functions
# ============================================================

def left(x):
    if isinstance(x, Pbox):        return x.left()
    if isinstance(x, Interval):    return x.left()
    if isscalar(x):                return float(x)
    if isinstance(x, np.ndarray):
        if x.size == 1:            return float(x)
        if x.size == 2:            return float(x[0])
        return min(x)
    return min(x)   # lists, tuples, iterables, scalars

def right(x):
    if isinstance(x, Pbox):        return x.right()
    if isinstance(x, Interval):    return x.right()
    if isscalar(x):                return float(x)
    if isinstance(x, np.ndarray):
        if x.size == 1:            return float(x)
        if x.size == 2:            return float(x[1])
        return max(x)
    return max(x)

def leftside(x):
    if isinstance(x, Pbox):        return x.leftside()
    if isinstance(x, Interval):    return x.left()
    if isinstance(x, list):        return [left(elem) for elem in x]
    if isscalar(x):                return float(x)
    return min(x)

def rightside(x):
    if isinstance(x, Pbox):        return x.rightside()
    if isinstance(x, Interval):    return x.right()
    if isinstance(x, list):        return [right(elem) for elem in x]
    if isscalar(x):                return float(x)
    return max(x)

def steps(x):
    if isinstance(x, Pbox):        return x.n
    if isinstance(x, Interval):    return 1
    if isscalar(x):                return 1
    raise TypeError(f"steps() not defined for type {type(x)}")

def mean(x):
    if isinstance(x, Pbox):        return x.mean()
    if isinstance(x, Interval):    return x
    if isscalar(x):                return x
    if isinstance(x, np.ndarray):  return np.mean(x)
    raise TypeError(f"mean() not defined for type {type(x)}")

def sd(x, pop=True):
    if isinstance(x, Pbox):        raise NotImplementedError('sd for p-boxes') #return (None, None) # Placeholder for your extremal-variance algorithm
    if isinstance(x, Interval):    return Interval(0.0, abs(x.hi - x.lo) / 2)
    if isscalar(x):                return 0.0
    raise TypeError(f"sd() not defined for type {type(x)}")

def breadth(x):
    if isinstance(x, Pbox):        return np.sum(x.d - x.u) / x.n
    if isinstance(x, Interval):    return right(x) - left(x)
    if isscalar(x):                return 0
    raise TypeError(f"breadth() not defined for type {type(x)}")

def var(x, pop=True): return sd(x,pop)**2

def straddles(x, z = 0): return left(x) <= z <=right(x) 

def support(x, interval=True): return Interval(left(x), right(x)) if interval else [left(x), right(x)] 

def cut(x,s,tight=True):
    try: return x.cut(s,tight)
    except: return support(x,True)
    
def iqr(x,tight):
    try: return x.iqr(tight)
    except: return support(x,True)
    
def width(x):
    try: return x.width()
    except: return right(x) - left(x)

# ============================================================
# Mass reassignment functions
# ============================================================

# already have constrainedto, on01

truncate = constrainedto

def rev(a): return a[::-1]         # I can never remember this odd construction

def above(x, s):   # exclude all values from x but those above s
    if right(x) < s: raise ValueError('No values are above the threshold')
    if isinstance(x,Interval): return constrainedto(x,s,float('inf'))      
    if not is_pbox(x): raise ValueError('The argument to the above() truncation must be uncertain')
    p = 1-left(x.prob(s))
    n = PbO.steps
    ij = seq(0, n*p-1, n)
    zu = np.maximum(s,rev(rev(x.u)[np.ceil(ij).astype(int)]))
    zd = np.maximum(s,rev(rev(x.d)[np.floor(ij).astype(int)]))      
    return Pbox(u = zu, d = zd)

def below(x, s):   # exclude all values from x but those below s
    if s < left(x): raise ValueError('No values are below the threshold')
    if isinstance(x,Interval): return constrainedto(x,-float('inf'),s)      
    if not is_pbox(x): raise ValueError('The argument to the below() truncation must be uncertain')
    p = right(x.prob(s))
    n = PbO.steps
    ij = seq(0, n*p-1, n)
    zu = np.minimum(s,x.u[np.floor(ij).astype(int)])
    zd = np.minimum(s,x.d[np.ceil(ij).astype(int)])
    return Pbox(u = zu, d = zd)

def between(x, m, M):
    return below(above(x, m), M)

def lowest(x, p): # exclude all but the lowest p% of x, useful for Will Powley's rejection rescaling
  if not is_pbox(x): raise ValueError('The argument to the lowest() truncation must be uncertain')
  i = np.round(seq(0,PbO.steps*(p/100),PbO.steps)).astype(int)
  zu = x.u[i]
  zd = x.d[i]
  return Pbox(u = zu, d = zd)
  
def highest(x, p): # exclude all but the highest p% of x, useful for Will Powley's rejection rescaling
  if not is_pbox(x): raise ValueError('The argument to the highest() truncation must be uncertain')
  i = np.round(seq(0,PbO.steps*(p/100),PbO.steps)).astype(int)
  zu = rev(rev(x.u)[i])
  zd = rev(rev(x.d)[i])
  return Pbox(u = zu, d = zd)
  
def rescale(x, m, M): return m + (M-m) * (x - left(x))/width(x)  # linearly rescale x to the specified range [m,M], assumes scalar m

def censor(x, m, M):            # remove detail from x between m and M
  zu = x.u
  zu[(m < zu) & (zu < M)] = m
  zd = x.d
  zd[(m < zd) & (zd < M)] = M
  return Pbox(zu,zd)

def massreassignexamples(a=None, b=None, th=4, ith = Interval(4,5), pth = 25): 
    # tutorial display of various mass reassignment functions [the R version is more elaborate, with another line including the barbell transform]
    def c(*inputs):
        p = []
        for x in inputs : p.append(x)
        return np.array(p)    
    if a is None: a = N(5,1)
    if b is None: b = U(4,5)
    abelow = below(a,th)
    aabove = above(a,th)
    alowest = lowest(a,pth)
    ahighest = highest(a,pth)
    apmin = pmin(a,b)   
    apmax = pmax(a,b)  
#    apmini = a %|m|% b
#    apmaxi = a %|M|% b
    asmin = least(a,b)      # smin()
    asmax = greatest(a,b)   # smax()
    abetween = between(a,left(ith),right(ith))
    abetween2 = between(a, left(cut(a,pth/100)), right(cut(a,1-pth/100)))
    atruncate = truncate(a,left(ith),right(ith))
    arescale = rescale(a,left(ith),right(ith))
    aenv = env(a,b)
    atrunc = trunc(a)
    aceil = ceil(a)
    afloor = floor(a)
    around = round(a)
    acensor = censor(a,left(ith),right(ith))
    #old.par <- par(mfrow=c(4,5),mar=c(2,4,3,1))
    fig, ax = init_splot(4,5, sharex=True, sharey=True)
    plt.rcParams['axes.titlepad']=1 # alternatively, use title('t',size=8,pad=1)
    splot(a,c='xkcd:grey'); title(size=8,s='highest'); lines(c(0,left(ahighest),left(ahighest)),1-c(pth,pth,100)/100,c='b',ls=':'); lines(ahighest,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='above');   lines(c(th,th),c(0,1),ls=':',c='b'); lines(aabove,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='smax');    lines(b,c='b');   lines(asmax,lw=3)
    splot(); 
    #splot(a,c='xkcd:grey'); title(size=8,s='%|M|%');   lines(b,c='b');   lines(apmaxi,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='pmax');    lines(b,c='b');   lines(apmax,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='lowest');  lines(c(0,right(alowest),right(alowest)),c(pth,pth,100)/100,ls=':',c='b'); lines(alowest,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='below');   lines(c(th,th),c(0,1),ls=':',c='b'); lines(abelow,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='smin');    lines(b,c='b');   lines(asmin,lw=3)
    splot(); 
    #splot(a,c='xkcd:grey'); title(size=8,s='%|m|%');   lines(b,c='b');   lines(apmini,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='pmin');    lines(b,c='b');   lines(apmin,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='between(cut,cut)'); lines(abetween2,lw=3); lines(c(left(abetween2),left(abetween2),2.5,2.5,right(abetween2),right(abetween2)),c(0,pth/100,pth/100,1-pth/100,1-pth/100,1),c='b',ls=':')
    splot(a,c='xkcd:grey'); title(size=8,s='between'); lines(abetween,lw=3); lines(c(left(ith),right(ith)),c(0,0),c='b',ls=':')
    splot(a,c='xkcd:grey'); title(size=8,s='truncate');lines(atruncate,lw=3);lines(c(left(ith),right(ith)),c(0,0),c='b',ls=':')
    splot(a,c='xkcd:grey'); title(size=8,s='rescale'); lines(arescale,lw=3); lines(c(left(ith),right(ith)),c(0,0),c='b',ls=':')
    splot(a,c='xkcd:grey'); title(size=8,s='env');     lines(aenv,lw=3)
    splot(a,c='xkcd:grey'); title(size=8,s='ceil');    lines(aceil,lw=3); 
    splot(a,c='xkcd:grey'); title(size=8,s='floor');   lines(afloor,lw=3);     
    splot(a,c='xkcd:grey'); title(size=8,s='trunc');   lines(atrunc,lw=3); 
    splot(a,c='xkcd:grey'); title(size=8,s='round');   lines(around,lw=3); 
    splot(a,c='xkcd:grey'); title(size=8,s='censor');  lines(acensor,lw=3); 

# ============================================================
# Unary transformations
# ============================================================

# already have defined abs, negate, reciprocate

# recurrant functions (harder to implement):  sin, cos, tan, sec, csc, ctn

def complement(x):
    return 1-x
   
def tonp(x, f, *args, **kwargs):
    with np.errstate(divide='ignore', invalid='ignore') if IvO.suppress_np_warnings else np.errstate():
        if isinstance(x, Pbox): return Pbox(f(x.u, *args, **kwargs), f(x.d, *args, **kwargs))
        if isinstance(x, Interval): return Interval(f(left(x), *args, **kwargs), f(right(x), *args, **kwargs))
        return f(x, *args, **kwargs)

def exp(x, *args, **kwargs): return tonp(x,np.exp, *args, **kwargs)

def square(x, *args, **kwargs): return tonp(x,np.square, *args, **kwargs)

def atan(x, *args, **kwargs): return tonp(x,np.atan, *args, **kwargs)

def acot(x, *args, **kwargs): return tonp(x,lambda z: np.pi/2 - np.atan(z), *args, **kwargs)  # acot = lambda z: np.pi/2 - np.atan(z)

def round(x, *args, **kwargs): return tonp(x,np.round, *args, **kwargs)

def trunc(x, *args, **kwargs): return tonp(x,np.trunc, *args, **kwargs) # different from truncate()

def ceil(x, *args, **kwargs): return tonp(x,np.ceil, *args, **kwargs)

def floor(x, *args, **kwargs): return tonp(x,np.floor, *args, **kwargs)

def makefinite(x): 
    def finite_clip(arr):
        arr = np.asarray(arr, dtype=float)
        finite_vals = arr[np.isfinite(arr)]
        if finite_vals.size == 0: return arr  # nothing to clip against
        lo,hi = finite_vals.min(), finite_vals.max()
        arr = np.where(arr == np.inf, hi, arr)
        arr = np.where(arr == -np.inf, lo, arr)
        return arr
    u = finite_clip(x.u)
    d = finite_clip(x.d)
    return Pbox(u,d,np.mean(u),np.mean(d))

def retaindomain(x, endpoints):
    return between(x,*endpoints)
  
def domaintonp(x, domain, why, f, *args, **kwargs): # for domain-restricted unary functions like log, sqrt, asin, acos, lambertW
    if not (support(x) in domain):
        if IvO.quieterrors:
            try: x = retaindomain(x, domain)
            except ValueError:
                if isinstance(x, (Interval, Pbox)): raise ValueError(f.__name__+' '+why)
    result = tonp(x, f, *args, **kwargs)
    return float(result) if isscalar(result) else result

def log(x, base=math.e, *args, **kwargs):
    f = lambda z: np.log(z)/np.log(base)
    return domaintonp(x, positives, IvO.why_nonpositive, f, *args, **kwargs)

def sqrt(x, *args, **kwargs):
    return domaintonp(x, positifs, IvO.why_negative, np.sqrt, *args, **kwargs)

def asin(x, *args, **kwargs):
    return domaintonp(x, unitdisk, IvO.why_outofunitdisk, np.arcsin, *args, **kwargs)

def acos(x, *args, **kwargs):
    return domaintonp(x, unitdisk, IvO.why_outofunitdisk, np.arccos, *args, **kwargs)

def lambertw(x, k=0, tol=1e-8, *args, **kwargs):
    from scipy.special import lambertw
    f = lambda z: lambertw(z, k=0, tol=1e-8)
    return domaintonp(x, Interval(-1/np.e,inf), 'undefined for such negative values', f, *args, **kwargs)

def sign(x):
    if (right(x)<0): return -1 
    elif (0<left(x)): return +1 
    elif (is_scalar(x) and identical(x.u[[1]],0)): return 0 
    elif (right(x)<=0): return Interval(-1,0) 
    elif (0<=left(x)): return Interval(0,1) 
    else: return Interval(-1,1)

def sigilium(x):
    if right(x)<0: return -1 
    elif 0<left(x): return +1 
    else: return Pbox(u=np.sign(x.u), d=np.sign(x.d))

'''N.B. The function tonp_nn() is a template for 'lifting' calculation of some
function through a p-box by referencing a function that works for intervals.
This particular algorithm presumes the stacked intervals that comprise a p-box
are [u[k], d[k]], where u and d are the left and right quantiles of the p-box.
The cut() function, on the other hand, uses a much more complex for the index
k (and the index for the left side differs from that for the right), which is 
necessary to ensure conservative characterization of the interval quantiles.

PbO.steps = 1000  # reveals a strangeness, surely a bug somewhere
x = N(5,1)
gray(x)
n = steps(x)
for p in seq(0.001,0.999):

    # conservative
    if (p % (1/n)) == 0: lower = np.round(p * n)
    else: lower = np.ceil(p * n)
    if (p % (1/n)) == 0: upper = np.round(p * n) + 1
    else: upper = np.floor(p * n) + 1
    ku = int(max(lower, 1)) - 1
    kd = int(min(upper, n)) - 1
    plt.plot([x.u[ku], x.d[kd]], [p,p], 'b')     

    # tight
    pn = p * n
    fractional = (pn % 1) == 0
    ku = int(min(n, (1 if fractional else 0) + ceil(pn))) - 1
    kd = int(max(1, ceil(pn))) - 1
    plt.plot([x.u[ku], x.d[kd]], [p,p], 'c')
    
# what tonp_nm() is assuming
for k in range(n):
    po = (k+0.5)/n * p
    plt.plot([x.u[k], x.d[k]], [po,po], 'r')

plt.plot([2,8.6],[1,1],lw=0.25)



PbO.steps = 1001  
x = N(5,1)    # crashes



PbO.steps = 1001  # reveals strangeness, surely a bug in tonp_nm's idea
x = N(5,1)
gray(x)
n = steps(x)
for p in seq(0.1,0.9):

    # conservative
    plt.plot([*ends(cut(x,p))], [p,p], 'b')     

    # tight
    plt.plot([*ends(cut(x,p,tight=True))], [p,p], 'c')
    
# what tonp_nm() is assuming
for k in range(n):
    po = (k+0.5)/n * p
    plt.plot([x.u[k], x.d[k]], [po,po], 'r')

plt.plot([2,8.6],[1,1],lw=0.25)

'''
def tonp_nm(x, f_interval, *args, **kwargs): # for nonmonotone functions (doesn't work on list or nparrays)
    with np.errstate(divide='ignore', invalid='ignore') if IvO.suppress_np_warnings else np.errstate():
        if isinstance(x, Pbox):
            n = x.n
            u2 = np.empty(n)
            d2 = np.empty(n)
            for k in range(n):
                Ik = Interval(x.u[k], x.d[k])         # full interval at step k
                Jk = f_interval(Ik, *args, **kwargs)
                u2[k] = left(Jk)
                d2[k] = right(Jk)
            u2.sort()
            d2.sort()
            return Pbox(u2, d2) 
        if isinstance(x, Interval): return f_interval(x, *args, **kwargs)
        return f_interval(Interval(x, x), *args, **kwargs).lo

def abs(x, *args, **kwargs):
    return tonp_nm(x, Interval.__abs__, *args, **kwargs)
    
def square(x, *args, **kwargs):
    return tonp_nm(x, Interval.square, *args, **kwargs)
   
def interval_cos(I):
    a,b = I          
    #if 2*np.pi < b-a: return unitdisk
    A,B = a//np.pi, b//np.pi
    if 1 < B-A: return unitdisk
    fa,fb = np.cos(ends(I))
    if A==B: return Interval(fa,fb)
    if A % 2 == 0: return Interval(-1, max(fa,fb))
    else: return Interval(min(fa,fb),1)
    
def cos(x, *args, **kwargs): return tonp_nm(x, interval_cos, *args, **kwargs) # doesn't work on list or nparrays
   
def interval_sin(I): return interval_cos(I - np.pi/2)

def sin(x, *args, **kwargs): return tonp_nm(x, interval_sin, *args, **kwargs) # doesn't work on list or nparrays
    
def interval_tan(I):
    a,b = ends(I + np.pi/2)
    A,B = a//np.pi, b//np.pi
    if 0 < B-A: return anyreal
    fa,fb = np.tan(ends(I))
    return Interval(fa,fb)
       
def tan(x, *args, **kwargs): return tonp_nm(x, interval_tan, *args, **kwargs) # doesn't work on list or nparrays

def interval_asec(x): 
    #if isinstance(x,(Number,np.ndarray)): return np.acos(1/x)
    if unitdisk.contains(x): return float('nan')
    if x.contains(unitdisk): return Interval(0,np.pi)
    if straddles(x,-1): return Interval(np.acos(left(x)),np.pi)
    if straddles(x,1):  return Interval(0,np.acos(right(x)))
    return acos(I(1/right(x), 1/left(x)))




def interval_asec(I):
    a, b = I.lo, I.hi

    # left branch: (-inf, -1]
    JL = None
    if a <= -1:
        L_lo = a
        L_hi = min(b, -1.0)
        if L_lo <= L_hi:
            JL = Interval(asec_scalar(L_lo), asec_scalar(L_hi))

    # right branch: [1, inf)
    JR = None
    if b >= 1:
        R_lo = max(a, 1.0)
        R_hi = b
        if R_lo <= R_hi:
            JR = Interval(asec_scalar(R_lo), asec_scalar(R_hi))

    if JL is None and JR is None: raise ValueError("asec interval: no in-domain values in input")
        # no in-domain values at all
        # you can raise, return None, or a special "empty" interval
        

    if JL is None: return JR
    if JR is None: return JL

    # convex hull of both branches
    return Interval(min(JL.lo, JR.lo), max(JL.hi, JR.hi))




def asec(x, *args, **kwargs): return tonp_nm(censor(x,*unitdisk), interval_asec, *args, **kwargs) # doesn't work on list or nparrays


'''
for m in [5,4,3,2,1,0]:
    a = N(m,0.1)
    plot(unitdisk)
    red(censor(a,*unitdisk))
    plot(asec(a))
    plt.title(str(m))
    plt.show()

x = np.array(seq(-4,4))
plt.plot(x,np.acos(1/x))
'''





arcsin = asin
arccos = acos
arctan = atan
arccot = acot
arcsec = asec
#arccsc = acsc

# ----------------------------------------------------------------------
# Some constant structures
# ----------------------------------------------------------------------
    
inf = float('inf')
dunno = Interval(0,1)
anyreal = Interval(-inf,inf)
positives = Interval(np.nextafter(0, 1), inf)
negatives = Interval(-inf, -np.nextafter(0, 1))
unitdisk = Interval(-1,1)
positifs = Interval(0,inf)
negatifs = Interval(-inf,0)
empty = I(inf,-inf,auto=False)

# ----------------------------------------------------------------------
# Optional IPython auto-plotting support
# ----------------------------------------------------------------------

from IPython import get_ipython

def _pbox_display_hook(obj, p, cycle):
    plot(obj)
    p.text(repr(obj))
    return None                    # None tells IPython to also print repr(obj)

ip = get_ipython()
if ip is not None:
    ip.display_formatter.formatters['text/plain'].for_type(Pbox, _pbox_display_hook)
    
# -----------------------------------------------------------------------------
# End of Python uncertain number library pbox.py
# -----------------------------------------------------------------------------





""" CHAPTER

 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS 
 KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS       KNOWN BUGS      


PbO.steps = 2001
u = U(0,1)
n = N(5,1)
env(u,n)
PbO.steps = 200
env(u,n)        # ValueError: zero-size array to reduction operation minimum which has no identity


PbO.steps = 2001
u = U(0,1)
PbO.steps = 200
n = N(5,1)
env(u,n)        # ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape after 1 dimensions. The detected shape was (2,) + inhomogeneous part.


# -----------------------------------------------------------------------------
# End of KNOWN BUGS
# -----------------------------------------------------------------------------
"""



""" CHAPTER

  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY
  LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY   LOGIC LIBRARY

# -----------------------------------------------------------------------------
# Logic library should be physically incorporated into pbox.py here
# -----------------------------------------------------------------------------

#from logic import *   # Copilot freaked out about such an import


# -----------------------------------------------------------------------------
# End of LOGIC LIBRARY
# -----------------------------------------------------------------------------
"""


""" CHAPTER

  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
  IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    IOANNA    
"""
# ----------------------------------------------------------------------
# Abortive ioanna6.py Python library (maybe some stuff worth rescuing)
# ----------------------------------------------------------------------
'''

# There are about 3,000 lines of code here 


# -*- coding: utf-8 -*-
"""
A few Python algorithms for fitting precise distributions to data using

     maximum likelihood
     method of matching moments
     confidence boxes
     maximum entropy
     Bayesian inference
     maximum a posteriori
     PERT
     Fermi methods

I spent too much time on the ancillary functions that construct named probability
distributions.  I assume you already have most of those, and you'll want to swap 
out the functions I made on lines 160-328, or maybe everything on lines 22-328.

     
@author: Scott Ferson
Created starting 13:56:06 ET, Thursday, 12 November 2024
"""


###############################################################################
# Ancillary (infrastructural) functions 
###############################################################################
#
# MOST OF THE FUNCTIONS IN THIS SECTION SHOULD BE REPLACED BY YOUR OWN FUNCTIONS.
#
# These functions use a Monte Carlo assemblage of deviates, stored as a Numpy
# ndarray, to model a probability distribution, and a pair of them to model a 
# p-box.  Thus, if B is a p-box, the array B[0:many] is its left side, and the 
# array B[many:(2*many)] is its rightside.
#
# Perhaps the more serious mistake is using scipy.stats for the random deviate
# algorithms rather than numpy.random.  See https://stackoverflow.com/questions/4001577/difference-between-random-draws-from-scipy-stats-rvs-and-numpy-random
# for the difference: basically, scipy.stats is creating *distributions* from 
# which we draw random values with its rvs functions, whereas numpy.random 
# is generating random values directly, with a bit less overhead.  In essence, 
# scipy generates a random variable while numpy generates random numbers.  But
# it's not clear that we need this extra stuff.
#
# Are we using the divide-by-n formula or the formula that divides by n-1?
# The MLE and MoMM estimates both expect the population formulas for standard 
# deviation.  So, in R, we needed a correction psd = sqrt(((n-1)*sd(x)^2)/n) =
# sd(x)*sqrt(1-1/n), but this correction is not needed in Python's Numpy as
#   np.std(x)             computes the population standard deviation
#   np.std(x, ddof=1)     computes the sample standard deviation
#   np.var(x)             computes the population variance
#   np.var(x, ddof=1)     computes the sample variance
#
# Python, like R, can compose arguments to make compound distributions, e.g.,
# sps.norm.rvs(np.arange(100),2,size=100) makes a compound distribution
# from normals with increasing means.
#
# Default arguments in Python are a little bit clumsier than in R.  But it is
# possible to emulate them in Python.  For example, the Python function gg() 
# will behave like the R function g().
#
# g <- function(shape,rate=1,scale=1/rate) {rate = 1/scale; cat('rate:',rate,'  ','scale:',scale,'\n')}
# g(0)         # rate: 1    scale: 1 
# g(0,1)       # rate: 1    scale: 1 
# g(0,2)       # rate: 2    scale: 0.5 
# g(0,rate=1)  # rate: 1    scale: 1 
# g(0,rate=2)  # rate: 2    scale: 0.5 
# g(0,scale=1) # rate: 1    scale: 1 
# g(0,scale=2) # rate: 0.5  scale: 2 
#
#def gg(shape,rate=1,scale=None) :
#    if scale is None : scale = 1/rate
#    rate = 1/scale
#    print('rate:',rate,'  ','scale:',scale)
# gg(0)         # rate: 1.0    scale: 1.0
# gg(0,1)       # rate: 1.0    scale: 1.0
# gg(0,2)       # rate: 2.0    scale: 0.5
# gg(0,rate=1)  # rate: 1.0    scale: 1.0
# gg(0,rate=2)  # rate: 2.0    scale: 0.5
# gg(0,scale=1) # rate: 1.0    scale: 1
# gg(0,scale=2) # rate: 0.5    scale: 2
#

import sys
import traceback
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sps

many = 10000  # increase for more accuracy; decrease for speed

def stop(msg) :
    print(msg)
    print(traceback.format_exc()) # way too much: traceback.print_stack()
    sys.exit(1)

def left(x) : return(np.min(x))

def right(x) : return(np.max(x))

def env(x,y) : return(np.concatenate((x,y)))

def leftside(x) : return(x[0:many])

def rightside(x) : 
    if many < len(x) : return(x[(many):(2*many)])
    else :  return(x[0:many])

def ci(b, c=0.95, alpha=None, beta=None) :
    if alpha is None : alpha=(1-c)/2
    if beta is None : beta=1-(1-c)/2
    left = np.sort(b[0:many])[round(alpha*many)]
    if (many < len(b)) : right = np.sort(b[many:len(b)])[round(beta*many)] 
    else : right = np.sort(b[0:many])[round(beta*many)]
    return((left,right))

def pl(x,y=None) : 
    if not y is None : x = [x,y]
    plt.ylim(0, 1)
    plt.xlim(min(x), max(x))

def ecdf(d) :
    d = np.array(d)
    N = d.size
    pp= np.concatenate((np.arange(N),np.arange(1,N+1)))/N
    dd = np.concatenate((d,d))
    dd.sort()
    pp.sort()
    return dd,pp

def edf(d,c=None,lw=None,ls=None) :
    if d.size==(2*many) : # p-box
        z,p = ecdf(d[0:many])
        plt.plot(z,p,c=c,lw=lw,ls=ls)
        z,p = ecdf(d[(many):(2*many)])
        plt.plot(z,p,c=c,lw=lw,ls=ls)
    else : # distribution
        z,p = ecdf(d)
        plt.plot(z,p,c=c,lw=lw,ls=ls)
    #plt.ylabel('Cumulative probability') # just makes the graph smaller

def red(d,c=None,lw=None,ls=None) : edf(d,c='r',lw=lw,ls=ls)
def cyan(d,c=None,lw=None,ls=None) : edf(d,c='c',lw=lw,ls=ls)
def blue(d,c=None,lw=None,ls=None) : edf(d,c='b',lw=lw,ls=ls)
def green(d,c=None,lw=None,ls=None) : edf(d,c='g',lw=lw,ls=ls)
def black(d,c=None,lw=None,ls=None) : edf(d,c='k',lw=lw,ls=ls)
def yellow(d,c=None,lw=None,ls=None) : edf(d,c='y',lw=lw,ls=ls)
def orange(d,c=None,lw=None,ls=None) : edf(d,c='orange',lw=lw,ls=ls)
def purple(d,c=None,lw=None,ls=None) : edf(d,c='purple',lw=lw,ls=ls)
     
def digamma(x) : 
    # Python by John Burkardt based on Fortran by Jose Bernardo 
    # https://people.math.sc.edu/Burkardt/py_src/asa103/digamma.py
    if (x <= 0) : return -np.Inf
    if (x <= 0.000001) : return -0.57721566490153286060-1.0/x+1.6449340668482264365*x
    value = 0.0
    while (x < 8.5) :
        value = value - 1.0 / x
        x = x + 1.0
    r = 1.0 / x
    value = value + np.log (x) - 0.5 * r
    r = r * r
    return value - r*(1/12.0 - r*(1/120.0 - r*(1/252.0 - r*(1/240.0 - r*(1/132.0))))) 
    
def uniroot(f,a) : 
    # https://stackoverflow.com/questions/43271440/find-a-root-of-a-function-in-a-given-range
    # https://docs.scipy.org/doc/scipy-0.18.1/reference/optimize.html#root-finding
#    from scipy.optimize import brentq
#    return(brentq(f, min(a), max(a))) #,args=(t0)) # any function arguments beyond the varied parameter
    from scipy.optimize import fsolve
    return(fsolve(f, (min(a) + max(a))/2)) 

def zbuff(x) : return(x) # use 1/zbuff(x) if x touches zero; unnecessary with Monte Carlo distribution models



###############################################################################        
# Precise distribution constructors 
#
# Most of these functions should be replaced by better Python implementations. 
# These functions serve as placeholders so that the other constructors using 
# MLE, maxent, MoMM, and Bayes, etc. can be implemented and tested.  
#
# There are two problems that demand these functions be replaced.  (1) The 
# basic distribution constructors should yield p-boxes when any arguments are
# intervals (which these algorithms don't do). (2) These algorithms produce 
# distributions represented internally as collections of Monte Carlo deviates 
# rather than some semi-analytical or discrete representation used in Risk Calc.
# See the preamble to https://sites.google.com/site/confidenceboxes/software        

def bernoulli(p) : return(np.random.uniform(size=many) < p)

def beta(a,b) :
    #if (a==0) and (b==0) : return(env(np.repeat(0.0, many), np.repeat(1.0, many)))  # this is [0,1]
    if (a==0) and (b==0) : return(bernoulli(0.5))  # or should it be [0,1]?
    if (a==0) : return(np.repeat(0.0, many))
    if (b==0) : return(np.repeat(1.0, many))            
    return(sps.beta.rvs(a,b,size=many))

def beta1(m,s) : return(beta(m * (m * (1 - m) / (s**2) - 1), (m * (m * (1 - m) / (s**2) - 1)) * (1/m - 1)))

def betabinomial2(size,v,w) : return(sps.binom.rvs(size,beta(v,w),size=many))

def betabinomial(size,v,w) : return(sps.betabinom.rvs(size,v,w,size=many))

def binomial(size,p) : return(sps.binom.rvs(size,p,size=many))

def chisquared(v) : return(sps.chi2.rvs(v,size=many))

def delta(a) : return(np.repeat(a,many))

def exponential(rate=1,mean=None) :
    if mean is None : mean = 1/rate
    #rate = 1/mean
    return(sps.expon.rvs(scale=mean,size=many))

def exponential1(mean=1) :
    return(sps.expon.rvs(scale=mean,size=many))

def F(df1,df2) : return(sps.f.rvs(df1,df2,size=many))

def gamma(shape,rate=1,scale=None) :
    if scale is None : scale = 1/rate
    rate = 1/scale
    return(sps.gamma.rvs(a=shape,scale=1/rate,size=many))

def gammaexponential(shape,rate=1,scale=None) :
    if scale is None : scale = 1/rate
    rate = 1/scale
    #expon(scale=gamma(a=shape, scale=1/rate))
    return(sps.expon.rvs(scale=1/sps.gamma.rvs(a=shape,scale=scale,size=many),size=many))

def geometric(m) : return(sps.geom.rvs(m,size=many))

def gumbel(loc,scale) : return(sps.gumbel_r.rvs(loc,scale,size=many))

def inversechisquared(v) : return(1/chisquared(v))
    
def inversegamma(shape, scale=None, rate=None) : 
    if scale is None and not rate is None : scale = 1/rate
    return(sps.invgamma.rvs(a=shape,scale=scale,size=many))

def laplace(a,b) :  return(sps.laplace.rvs(a,b,size=many))

def logistic(loc,scale) : return(sps.logistic.rvs(loc,scale,size=many))

def lognormal(m,s) : 
    m2 = m**2; s2 = s**2
    mlog = np.log(m2/np.sqrt(m2+s2))
    slog = np.sqrt(np.log((m2+s2)/m2))
    return(sps.lognorm.rvs(s=slog,scale=np.exp(mlog),size=many))

def lognormal2(mlog,slog) : return(sps.lognorm.rvs(s=slog,scale=np.exp(mlog),size=many))

#lognormal = function(mean=NULL, std=NULL, meanlog=NULL, stdlog=NULL, median=NULL, cv=NULL, name='', ...){
#  if (is.null(meanlog) & !is.null(median)) meanlog = log(median)
#  if (is.null(stdlog) & !is.null(cv)) stdlog = sqrt(log(cv^2 + 1))
#  # lognormal(a, b) ~ lognormal2(log(a^2/sqrt(a^2+b^2)),sqrt(log((a^2+b^2)/a^2)))
#  if (is.null(meanlog) & (!is.null(mean)) & (!is.null(std))) meanlog = log(mean^2/sqrt(mean^2+std^2))
#  if (is.null(stdlog) & !is.null(mean) & !is.null(std)) stdlog = sqrt(log((mean^2+std^2)/mean^2))
#  if (!is.null(meanlog) & !is.null(stdlog)) Slognormal0(meanlog,stdlog,name) else stop('not enough information to specify the lognormal distribution')
#  }

def loguniform_solve(m,v) :
  def loguniform_f(a,m,v) : return(a*m*np.exp(2*(v/(m**2)+1)) + np.exp(2*a/m)*(a*m - 2*((m**2) + v)))
  def LUgrid(aa, w) : return(left(aa)+(right(aa)-left(aa))*w/100.0)
  aa = (m - np.sqrt(4*v), m)   # interval
  a = m
  ss = loguniform_f(a,m,v)
  for j in range(4) :
    for i in range(101) :  # 0:100 
      a = LUgrid( aa, i)
      s = abs(loguniform_f(a,m,v))
      if s < ss :
          ss = s
          si = i 
    a = LUgrid(aa, si)
    aa = (LUgrid(aa, si-1), LUgrid(aa, si+1))  # interval
  return(a)

def loguniform(min=None, max=None, minlog=None, maxlog=None, mean=None, std=None) :
    if (min is None) and (not (minlog is None)) : min = np.exp(minlog)
    if (max is None) and (not (maxlog is None)) : max = np.exp(maxlog)  
    if (max is None) and (not (mean is None)) and (not (std is None)) and (not (min is None)) : max = 2*(mean**2 +std**2)/mean - min
    if (min is None) and (max is None) and (not (mean is None)) and (not(std is None)) :
        min = loguniform_solve(mean,std**2)
        max = 2*(mean**2 +std**2)/mean - min
    return(sps.loguniform.rvs(min, max, size=many))

def loguniform1(m,s) : return(loguniform(mean=m, std=s))

def negativebinomial(size,prob) : return(sps.nbinom.rvs(size,prob,size=many))

def normal(m,s) : return(sps.norm.rvs(m,s, size=many))

def pareto(mode, c) : return(sps.pareto.rvs(c,scale=mode,size=many))

def poisson(m) : return(sps.poisson.rvs(m,size=many))

def powerfunction(b,c) : return(sps.powerlaw.rvs(c,scale=b,size=many))

# parameterisation of rayleigh differs from that in pba.r
def rayleigh(loc,scale) : return(sps.rayleigh.rvs(loc,scale,size=many))

def sawinconrad(min, mu, max) : # WHAT are the 'implicit constraints' doing?     
  def sawinconradalpha01(mu) :
      def f(alpha) : return(1/(1-1/np.exp(alpha)) - 1/alpha - mu)
      if np.abs(mu-0.5)<0.000001 : return(0)      
      return(uniroot(f,np.array((-500,500))))
  def qsawinconrad(p, min, mu, max) : 
        alpha = sawinconradalpha01((mu-min)/(max-min))
        if np.abs(alpha)<0.000001 : return(min+(max-min)*p) 
        else : min+(max-min)*((np.log(1+p*(np.exp(alpha)-1)))/alpha)
  a = left(min);   b = right(max)
  c = left(mu);    d = right(mu)
  if c<a : c = a   # implicit constraints
  if b<d : d = b
  #return(qsawinconrad(np.random.uniform(size=many), min, mu, max))
  return(qsawinconrad(np.random.uniform(size=many), min, mu, max))
  
def student(v) : return(sps.t.rvs(v,size=many))

def uniform(a,b) : return(sps.uniform.rvs(a,b-a,size=many)) # who parameterizes like this?!?!

def triangular(min,mode,max) : return(np.random.triangular(min,mode,max,size=many)) # cheating: uses random rather than sps

def histogram(x) : return(x[(np.trunc(sps.uniform.rvs(size=many)*len(x))).astype(int)])

def mixture(x,w=None) :
    if w is None : w = np.repeat(1,len(x))
    print(many)
    r = np.sort(sps.uniform.rvs(size=many))[::-1]
    x = np.concatenate(([x[0]],x))
    w = np.cumsum(np.concatenate(([0],w)))/np.sum(w)
    u = []
    j = len(x)-1
    for p in r : 
        while True :
            if w[j] <= p : break
            j = j - 1
        u = np.concatenate(([x[j+1]],u))
    return(u[np.argsort(sps.uniform.rvs(size=len(u)))])


"""
You mentioned that you'd prepared slides for Monday.  Do you wanna share those
by email, or just show them next Monday?

c-boxes

Characterising distribs

Thanks for sharing the pba.py and related files.
Shouldn't the test

     if args[i].__class__.__name__ != "Interval":

be isolated (wrapped) in a function?

    def isntinterval(a) :
        return a.__class__.__name__ != "Interval"
    
Also, the test seems fragile.  Could it be something more flexible such as

    def isntinterval(a) :
        return left(a) != right(a)
    
so entering a tuple or a list in the arguments would be interpreted robustly?


I am so impressed with Python that you can end an enumeration with a comma as in

    [1,2,3,]
    
That may be the second best thing about Python's syntax (after the asterisk).


Given what you were saying this afternoon, I was surprised that your pbox 
constructor norm() doesn't include the sps distribution object what it
returns:

 return Pbox(
        Left,
        Right,
        steps=steps,
        shape="norm",
        mean_left=mean.left,
        mean_right=mean.right,
        var_left=var.left,
        var_right=var.right,  # I am very impressed that you can have this comma here
        )








def Beta(*args) :
    o = sps.beta(*args)
    return o.rvs(size=many)



    new_args = itertools.product(*args)

    bounds = []

    mean_hi = -np.inf
    mean_lo = np.inf
    var_lo = np.inf
    var_hi = 0

    for a in new_args:
        bounds.append(dists[function_name].ppf(p, *a))




It seems like it would be helpful to be able to access the methods from within
a pbox.  Would this be bad or troublesome or create too much overhead?


Shouldn't the test

     if args[i].__class__.__name__ != "Interval":

be isolated (wrapped) in a function?

    def isntinterval(a) :
        return a.__class__.__name__ != "Interval"
    
Also, the test seems fragile.  Could it be something more flexible such as

    def isntinterval(a) :
        return left(a) != right(a)
    
so entering a tuple or a list in the arguments would be interpreted robustly?





def Normal(*args) : #, steps=Params.steps):
    args = list(args)  
    print(args)
    return sps.norm.rvs(*arg, size=many)

    Left, Right, mean, var = __get_bounds("norm", steps, *args)
    return Pbox(
        Left,
        Right,
        steps=steps,
        shape="norm",
        mean_left=mean.left,
        mean_right=mean.right,
        var_left=var.left,
        var_right=var.right,  # I am very impressed that you can have this comma here
        )

rvs(a, b, loc=0, scale=1, size=1, random_state=None)

cdf(x, a, b, loc=0, scale=1)



















import sps as s

a = s.lognorm.rvs(s=2,scale=np.exp(0.2), size=many)
la = s.lognorm.rvs(loc=2,s=2,scale=np.exp(0.2), size=many)

edf(a); red(la)
min(a)
min(la)



def beta(*args) : return s.beta(*args).rvs(size=many)
def normal(*args) : return s.norm(*args).rvs(size=many)

def lognormal(*args) : return s.norm(*args).rvs(size=many)
    
def uniform(*args) : return s.uniform(*args).rvs(size=many)  # U[loc, loc + scale]





    m2 = m**2; s2 = s**2
    mlog = np.log(m2/np.sqrt(m2+s2))
    slog = np.sqrt(np.log((m2+s2)/m2))
    return(sps.lognorm.rvs(s=slog,scale=np.exp(mlog),size=many))

def lognormal2(mlog,slog) : return(sps.lognorm.rvs(s=slog,scale=np.exp(mlog),size=many))

#lognormal = function(mean=NULL, std=NULL, meanlog=NULL, stdlog=NULL, median=NULL, cv=NULL, name='', ...){
#  if (is.null(meanlog) & !is.null(median)) meanlog = log(median)
#  if (is.null(stdlog) & !is.null(cv)) stdlog = sqrt(log(cv^2 + 1))
#  # lognormal(a, b) ~ lognormal2(log(a^2/sqrt(a^2+b^2)),sqrt(log((a^2+b^2)/a^2)))
#  if (is.null(meanlog) & (!is.null(mean)) & (!is.null(std))) meanlog = log(mean^2/sqrt(mean^2+std^2))
#  if (is.null(stdlog) & !is.null(mean) & !is.null(std)) stdlog = sqrt(log((mean^2+std^2)/mean^2))
#  if (!is.null(meanlog) & !is.null(stdlog)) Slognormal0(meanlog,stdlog,name) else stop('not enough information to specify the lognormal distribution')
#  }






    

def bernoulli(p=0.25);               sh(x,'bernoulli(p=0.25)') 
def beta(a=2,b=3) ;                  sh(x,'beta(a=2,b=3)') 
def betabinomial2(size=10,v=2,w=3);  sh(x,'betabinomial2(size=10,v=2,w=3)') 
def betabinomial(size=10,v=2,w=3);   sh(x,'betabinomial(size=10,v=2,w=3)')  
def binomial(12,0.4);                sh(x,'binomial(size=12,p=0.4)')
def chisquared(v=6);                 sh(x,'chisquared(v=6)')
def exponential(mean=2);             sh(x,'exponential(mean=2)') 
def F(6,11);                         sh(x,'F(df1=6,df2=11)')
def gamma(shape=4,rate=2);           sh(x,'gamma(shape=4,rate=2)')
def gammaexponential(shape=4,rate=2);sh(x,'gammaexponential(shape=4,rate=2)')
def geometric(m=0.3);                sh(x,'geometric(m=0.3)')
def gumbel(2,4);                     sh(x,'gumbel(loc=2,scale=4)')
def inversechisquared(14);           sh(x,'inversechisquared(df=14)') 
def inversegamma(shape=2,scale=4);   sh(x,'inversegamma(shape=2,scale=4)')
def laplace(a=4,b=5);                sh(x,'laplace(a=4,b=5)') 
def logistic(2,3);                   sh(x,'logistic(loc=2,scale=3)')
def lognormal(m=2,s=1);              sh(x,'lognormal(m=10,s=1)')
def lognormal2(mlog=-2,slog=1);      sh(x,'lognormal2(mlog=-2,slog=1)')
def loguniform(min=2, madef 6);        sh(x,'loguniform(min=2, madef 6)')
def negativebinomial(size=10,prob=0.25); sh(x,'negativebinomial(size=10,prob=0.25)') 
def normal(m=5,s=1) ;                sh(x,'normal(m=5,s=1)')
def pareto(mode=3, c=2);             sh(x,'pareto(mode=3, c=2)')
def poisson(m=4);                    sh(x,'poisson(m=4)')
#def rayleigh(4,3);                   sh(x,'rayleigh(4,3)')
def sawinconrad(2,4,9) ;             sh(x,'student(2,4,9)')
def student(v=5) ;                   sh(x,'student(v=5)')
def triangular(2,5,11);              sh(x,'triangular(2,5,11)')
def uniform(a=2,b=4) ;               sh(x,'uniform(a=2,b=4)') 



"""




        

##########################################################################
# Scott's maximum likelihood estimation constructors
##########################################################################


# Wikipedia [https://en.wikipedia.org/wiki/Maximum_likelihood_estimation] says
# "From the perspective of Bayesian inference, MLE is generally equivalent to 
# maximum a posteriori (MAP) estimation with a prior distribution that is uniform 
# in the region of interest. In frequentist inference, MLE is a special case of 
# an extremum estimator, with the objective function being the likelihood."

# Wikipedia [https://en.wikipedia.org/wiki/Maximum_likelihood_estimation] says
# "Maximum-likelihood estimators have no optimum properties for finite samples, 
# in the sense that (when evaluated on finite samples) other estimators may  
# have greater concentration around the true parameter-value. [Pfanzagl, Johann 
# (1994). Parametric Statistical Theory. Walter de Gruyter. pp. 207–208. 
# doi:10.1515/9783110889765. ISBN 978-3-11-013863-4. MR 1291393]"


# Some of these functions may not support intervals in the data x           #**

def sMLbernoulli(x) : return(bernoulli(x.mean())) # is this legit?

def sMLnormal(x) : return(normal(x.mean(),x.std()))                          #**

def sMLgaussian(x) : return(MLnormal(x))

def sMLexponential(x) : return(exponential(1/x.mean()))

def sMLpoisson(x) : return(poisson(x.mean())) # is this legit?

def sMLgeometric(x) : return(geometric(1/(1+x.mean())))

def sMLgumbel(x) :
    loc, scale = sps.gumbel_r.fit(x)
    return(gumbel(loc,scale))

def sMLpascal(x) : return(MLgeometric(x))

def sMLuniform(x) : return(uniform(min(x), max(x)))                         #**

def sMLrectangular(x) : return(MLuniform(x))

def sMLpareto(x) : return(pareto(min(x), len(x)/np.sum(np.log(x)-np.log(min(x)))))  #**

def sMLlaplace(x) : return(laplace(x.median(), np.sum(np.abs(x-x.median())/len(x))))  #**

def sMLdoubleexponential(x) : return(MLlaplace(x))
        
def sMLlognormal2(x) :                                                      #**
    n = len(x)
    mu = np.sum(np.log(x))/n
    return(lognormal2(mlog=mu, slog=np.sum((np.log(x)-mu)**2)/n))  # this function gives clearly poor results
    
def sMLlognormal(x) : return(np.exp(MLnormal(np.log(x)))) # just uses transformation, which seems unlikely to be true, but fitdistrplus package uses it too

def sMLloguniform(x) :
    a,b,_,_ = sps.loguniform.fit(x)
    return(loguniform(a,b))

def sMLweibull(x, shapeinterval=None) :                                     #**
    if shapeinterval is None : shapeinterval = np.array((0.001,500))
    def wf(k) : return(np.sum(x**k * np.log(x)) / np.sum(x**k) - np.sum(np.log(x)) / len(x) - 1/k)
    k = uniroot(wf, shapeinterval)
    el = np.exp(np.log(np.sum(x^k)/len(x))/k)
    return(sps.weibull_min.rvs(scale=el, shape=k))
  
def sMLgamma(data) :                                                        #**
    xbar = data.mean()
    shape=(xbar/data.std())**2  # initial estimate of shape from MoM
    logxbar = np.log(xbar)
    meanlog = np.log(data).mean()
    def f(x) : return(np.log(x) - digamma(x) - logxbar + meanlog)
    shape = uniroot(f,shape*np.array((0.5,1.5)))
    rate = shape/xbar
    return(gamma(shape=shape, rate=rate))



###############################################################################
# Alternative maximum likelihood estimation constructors using scipy.stats
###############################################################################

# Wikipedia [https://en.wikipedia.org/wiki/Maximum_likelihood_estimation] says
# From the perspective of Bayesian inference, MLE is generally equivalent to 
# maximum a posteriori (MAP) estimation with a prior distribution that is uniform 
# in the region of interest. In frequentist inference, MLE is a special case of 
# an extremum estimator, with the objective function being the likelihood.

# Some of these functions may support intervals in the data x.  See 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html#scipy.stats.rv_continuous.fit
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html#scipy.stats.CensoredData

def MLbernoulli(x) : return sps.bernoulli.rvs(x.mean(),size=many)
def MLbeta(x) : return sps.beta.rvs(*sps.beta.fit(x),size=many)
def MLbetabinomial(x) : return sps.betabinom.rvs(*sps.betabinom.fit(x),size=many)
def MLbinomial(x) : return sps.binom.rvs(*sps.binom.fit(x),size=many)
def MLchisquared(x) : return sps.chis.rvs(*sps.chi2.fit(x),size=many)
def MLexponential(x) : return sps.expon.rvs(*sps.expon.fit(x),size=many)
def MLF(x) : return sps.f.rvs(*sps.f.fit(x),size=many)
def MLgamma(x) : return sps.gamma.rvs(*sps.gamma.fit(x),size=many)
def MLgammaexponential(x) : return sps.gammaexponential(*sps.gammaexpon.fit(x),size=many)
def MLgeometric(x) : return sps.geom.rvs(*sps.geom.fit(x),size=many)
def MLgumbel(x) : return sps.gumbel_r.rvs(*sps.gumbel_r.fit(x),size=many)
def MLlaplace(x) :  return sps.laplace.rvs(*sps.laplace.fit(x),size=many)
def MLlogistic(x) : return sps.logistic.rvs(*sps.logistic.fit(x),size=many)
def MLlognormal(x) : return sps.lognorm.rvs(*sps.lognorm.fit(x),size=many)
def MLloguniform(x) : return sps.loguniform.rvs(*sps.loguniform.fit(x),size=many)
def MLnegativebinomial(x) : return sps.nbinom.rvs(*sps.nbinom.fit(x),size=many)
def MLnormal(x) : return sps.norm.rvs(*sps.norm.fit(x),size=many)
def MLpareto(x) : return sps.pareto.rvs(*sps.pareto.fit(x),size=many)
def MLpoisson(x) : return sps.poisson.rvs(*sps.poisson.fit(x),size=many)
def MLpowerfunction(x) : return sps.powerlaw.rvs(*sps.powerlaw.fit(x),size=many)
def MLrayleigh(x) : return sps.rayleigh.rvs(*sps.rayleigh.fit(x),size=many)
def MLstudent(x) : return sps.t.rvs(*sps.t.fit(x),size=many)
def MLtriangular(x) : return sps.triang.rvs(*sps.triang.fit(x),size=many)
def MLuniform(x) : return sps.uniform.rvs(*sps.uniform.fit(x),size=many)



###############################################################################
# Method-of-Moment distribution constructors (matching central moments of x)
###############################################################################

# Some of tThese functions may not support intervals in the data x          #**

def MMbernoulli(x) : return(bernoulli(x.mean())) # assumes x is zeros and ones

def MMbeta(x) : return(beta1(x.mean(), x.std()))                            #**

def MMbetabinomial(n,x) :                                                   #**
    # n must be provided; it's not estimated from data
    # https://en.wikipedia.org/wiki/Beta-binomial_distribution#Example:
    # MMbetabinomial(n=12,rep(0:12,c(3,24,104,286,670,1033,1343,1112,829,478,181,45,7))) 
    m1 = x.mean()
    m2 = (x**2).mean()
    d = n * (m2 / m1 - m1 - 1) + m1
    return(betabinomial(n, (n*m1 - m2) / d, (n-m1)*(n - m2/m1) / d))

def MMbinomial(x) :                                                         #**
    a = x.mean()
    b = x.std()
    return(binomial(int(np.abs(np.round(a/(1-b**2/a)))), np.abs(1-b**2/a)))

def MMchisquared(x) : return(chisquared(np.round(x.mean())))

def MMexponential(x) : return(exponential(mean=x.mean()))

def MMF(x) :                                                                #**
    w = 2/(1-1/x.mean())
    return(F(np.round((2*w**3 - 4*w**2) / ((w-2)**2 * (w-4) * x.std()**2 - 2*w**2)), np.round(w)))

def MMgamma(x) :                                                            #**
    a = x.mean()
    b = x.std()
    return(gamma(b**2/a, 1/(a/b)**2))  #gamma1(a, b) ~ gamma(b²/a, (a/b)²)

def MMgeometric(x) : return(geometric(1/(1+x.mean())))

def MMpascal(x) : return(geometric(1/(1+x.mean())))

#def MMgumbel0(x) : return(gumbel(x.mean() - 0.57721 * x.std() * np.sqrt(6)/ np.pi, x.std() * np.sqrt(6)/ np.pi))       #**  # https://stackoverflow.com/questions/51427764/using-method-of-moments-with-gumbel-r-in-python-scipy-stats-gumbel-r
def MMgumbel(x) :                                                           #**
    # https://stackoverflow.com/questions/51427764/using-method-of-moments-with-gumbel-r-in-python-scipy-stats-gumbel-r  
    scale = np.sqrt(6)/np.pi * np.std(x)
    loc = np.mean(x) - np.euler_gamma*scale
    return(gumbel(loc, scale))

def MMextremevalue(x) : return(gumbel(x.mean() - 0.57721 * x.std() * np.sqrt(6)/ np.pi, x.std() * np.sqrt(6)/ np.pi))       #**

def MMlognormal(x) : return(lognormal(x.mean(), x.std()))                   #**

def MMlaplace(x) : return(laplace(x.mean(), x.std()/np.sqrt(2)))            #**

def MMdoubleexponential(x) : return(laplace(x.mean(), x.std()/np.sqrt(2)))  #**

def MMlogistic(x) : return(logistic(x.mean(), x.std() * np.sqrt(3)/np.pi))  #**

def MMloguniform(x) : return(loguniform(mean=x.mean(), std=x.std()))        #**

def MMnormal(x) : return(normal(x.mean(), x.std()))                         #**

def MMgaussian(x) : return(normal(x.mean(), x.std()))                       #**

def MMpareto(x) :                                                           #**
    a = x.mean()
    b = x.std()
    return(pareto(a/(1+1/np.sqrt(1+a**2/b**2)), 1+np.sqrt(1+a**2/b**2)))

def MMpoisson(x) : return(poisson(x.mean()))

def MMpowerfunction(x) :                                                    #**
    a = x.mean()
    b = x.std()
    return(powerfunction(a/(1-1/np.sqrt(1+(a/b)**2)), np.sqrt(1+(a/b)**2)-1))

def MMt(x) : return(student(2/(1-1/x.std()**2)))

def MMstudent(x) : 
    if (1<x.std()) : return(student(2/(1-1/x.std()**2))) 
    else : stop('Improper standard deviation for student distribution')

def MMuniform(x) :                                                          #**
    a = x.mean()
    b = np.sqrt(3)*x.std()
    return(uniform(a-b, a+b))

def MMrectangular(x) : return(MMuniform(x))       

def MMtriangular(x,iters=100,dives=10) :                                    #**
  # iterative search for triangular distribution parameters using method of 
  # matching moments (you solve the thing analytically! too messy without help)
  # testing code indicated with #-#
  #-#some = 10
  #-#A = runif(1,0,10)
  #-#B = A + runif(1,0,10)
  #-#C = runif(1,A,B)
  #-#x = qtriangular(runif(some), A,C,B)
  def skewness(x) : 
      m = x.mean()
      return(np.sum((x-m)**3)/((len(x)-1)*x.std()**3))  # std uses the population formula, may need the sample formula
  M = np.mean(x)
  V = np.var(x)
  S = skewness(x)
  a = aa = min(x) # apparently double assignments work
  b = bb = max(x)
  c = cc = 3*M-a-b 
  many = iters
  s1 = np.std(x)
  for k in range(dives) :
    s1 = s2 = s3 = s1/2
    a = np.random.normal(aa,s1,many)
    b = np.random.normal(bb,s2,many)
    c = np.random.normal(cc,s3,many)
    m = (a+b+c)/3
    k = (a**2+b**2+c**2-a*b-a*c-b*c)
    v = k/18
    s = (np.sqrt(2)*(a+b-2*c)*(2*a-b-c)*(a-2*b+c)) / (5 * k ** (3/2))
    d = (M-m)**2 + (V-v)**2 + (S-s)**2
    i = np.argmin(d)  # which.min(d)
    aa = a[i]
    bb = b[i]
    cc = c[i]
  #-#gray(triangular(A,B,C), new = TRUE)
  #-#blue(x)
  #-#green(triangular(aa,bb,cc))
  #-#A;aa; B;bb; C;cc  # the order is min, max, mode
  print(aa,bb,cc)
  return(triangular(aa,cc,bb)) 
  


###############################################################################
#  Confidence boxes (c-boxes) for parameters and the next observable value
###############################################################################

def km(k,m) :
  # The formula env(beta(k,m+1),beta(k+1,m)) is intervalized and pared to [0,1].
  # If we're using Monte Carlo deviates to model distributions and p-boxes, 
  # we don't need Bone, Bzero, np.minimum or np.maximum 
  Bzero = 1e-6
  Bone = 1-Bzero
  if ((left(k) < 0)  or (left(m) < 0)) : stop('Improper arguments to function km')
  #if is.pbox(k) or is.pbox(m) : return(uchenna(pbox(k),pbox(m)))
  #else :
  return(np.minimum(np.maximum(env(beta(left(k),right(m)+1),beta(right(k)+1,left(m))),Bzero),Bone))

def KN(k,n) :
  # The formula env(beta(k,n-k+1),beta(k+1,max(0,n-k))) is intervalized and 
  # whittled down to [0,1].  If we're using Monte Carlo deviates to represent
  # distributions and p-boxes, we don't need Bone, Bzero, minimum or maximum 
  if ((left(k) < 0) or (right(n) < right(k))) : stop('Improper arguments to function KN')
  Bzero = 1e-6
  Bone = 1-Bzero
# return(np.minimum(np.maximum(env(beta(     k,       n -     k +1),beta(      k +1,np.maximum(0,     n -      k))) ,Bzero),Bone))
  return(np.minimum(np.maximum(env(beta(left(k),right(n)-left(k)+1),beta(right(k)+1,np.maximum(0,left(n)-right(k)))),Bzero),Bone))

def FKN(k,n) :  # binomial rate inference for trials designed with a fixed-K stopping rule
  if (left(k) < 0) or (right(n) < right(k)) : stop('Improper arguments to function KN')
  Bzero = 1e-6
  Bone = 1-Bzero
  return(np.minimum(np.maximum(env(beta(left(k),right(n)-left(k)+1),beta(right(k),np.maximum(0,left(n)-right(k)))),Bzero),Bone))

# the functionality of CBbernoulli and CBbinomial is condensed into km and KN

# x[i] ~ Bernoulli(p), x[i] is either 0 or 1
def CBbernoulli(x) : 
    n = len(x)
    k = sum(x)
    return(env(bernoulli(k/(n+1)), bernoulli((k+1)/(n+1))))
def CBbernoulli_p(x) :
    n = len(x)
    k = sum(x)
    return(env(beta(k, n-k+1), beta(k+1, n-k)))

# x[i] ~ binomial(N, p), for known N, x[i] is a nonnegative integer less than or equal to N
def CBbinomial(N,x) :
    n = len(x)
    k = sum(x)
    return(env(betabinomial(N,k,n*N-k+1),betabinomial(N,k+1, n*N-k)))
def CBbinomial_p(N,x) :
    n = len(x)
    k = sum(x)
    return(env(beta(k, n*N-k+1), beta(k+1, n*N-k)))

# x[i] ~ binomial(N, p), for unknown N, x[i] is a nonnegative integer
# see https://sites.google.com/site/cboxbinomialnp/
def CBbinomialnp(x) : stop('see https://sites.google.com/site/cboxbinomialnp/')
def CBbinomialnp_n(x) : stop('see https://sites.google.com/site/cboxbinomialnp/')
def CBbinomialnp_p(x) : stop('see https://sites.google.com/site/cboxbinomialnp/')

# x[i] ~ Poisson(mean), x[i] is a nonnegative integer
def CBpoisson(x) :
    n = len(x)
    k = sum(x)
    return(env(negativebinomial(size=k, prob=1-1/(n+1)),negativebinomial(size=k+1, prob=1-1/(n+1))))
def CBpoisson_mean(x) :
    n = len(x)
    k = sum(x)
    return(env(gamma(shape=k, rate=n),gamma(shape=k+1, rate=n)))

# x[i] ~ exponential(mean), x[i] is a nonnegative integer
def CBexponential(x) :
    n = len(x)
    k = sum(x)
    return(gammaexponential(shape=n, rate=k))
def CBexponential_mean(x) :
    n = len(x)
    k = sum(x)
    return(1/gamma(shape=n, rate=k))
def CBexponential_lambda(x) : return(1/CBexponential_mean(x))

# x[i] ~ normal(mu, sigma)
def CBnormal(x) : 
    n = len(x)
    return(np.mean(x) + np.std(x) * student(n - 1) * np.sqrt(1 + 1 / n))# pop or sample std?
def CBnormal_mu(x) :
    n = len(x)
    return(np.mean(x) + np.std(x) * student(n - 1) / np.sqrt(n))# pop or sample std?
def CBnormal_sigma(x) :
    n = len(x) 
    return(np.sqrt(np.var(x)*(n-1)*inversechisquared(n-1))) # pop or sample var?

# x[i] ~ lognormal(mu, sigma), x[i] is a positive value whose logarithm is distributed as normal(mu, sigma)
def CBlognormal(x) : 
    n = len(x)
    return(np.exp(np.mean(np.log(x)) + np.std(np.log(x)) * student(n - 1) * np.sqrt(1+1/n)))
def CBlognormal_mu(x) : 
    n = len(x)
    return(np.mean(np.log(x)) + np.std(np.log(x)) * student(n - 1) / np.sqrt(n))
def CBlognormal_sigma(x) : 
    n = len(x)
    return(np.sqrt(np.var(np.log(x))*(n-1)*inversechisquared(n-1)))

## x[i] ~ lognormal(mean, sd), where mean and sd are the mean and standard deviation of the x[i] values
#
## Would like a formula for the mean-stdev parameterization for lognormal, but this awkward strategy doesn't work:
#CBlognormal.mean <- function(x) {
#mu = CBlognormal.mu(x)
#sigma = CBlognormal.sigma(x)
#return(exp(mu %|+|% sigma^2)/2))
#}

# x[i] ~ uniform(midpoint, width)
# x[i] ~ uniform(minimum, maximum)
def CBuniform(x) :
    r=max(x)-min(x)
    w=(r/beta(len(x)-1,2))/2
    m=(max(x)-w)+(2*w-r)*uniform(0,1); 
    return(uniform(m-w, m+w))
def CBuniform_midpoint(x) : 
    r = max(x)-min(x) 
    w = r/beta(len(x)-1,2)
    m = (max(x)-w/2)+(w-(max(x)-min(x)))*uniform(0,1)
    return(m)
def CBuniform_width(x) : 
    r = max(x)-min(x) 
    return(r/beta(len(x)-1,2))
def CBuniform_minimum(x) : 
    r=max(x)-min(x); 
    w=r/beta(len(x)-1,2)
    m=(max(x)-w/2)+(w-r)*uniform(0,1)
    return(m-w/2)
def CBuniform_maximum(x) : 
    r=max(x)-min(x) 
    w=r/beta(len(x)-1,2)
    m=(max(x)-w/2)+(w-r)*uniform(0,1)
    return(m+w/2)
           
# x[i] ~ F, a continuous but unknown distribution
# N.B. the infinities don't plot, but they are there
def CBnonparametric(x) : return(env(histogram(np.concatenate((x, [-np.inf]))), histogram(np.concatenate((x, [np.inf])))))

# x1[i] ~ normal(mu1, sigma1), x2[j] ~ normal(mu2, sigma2), x1 and x2 are independent
def CBnormal_meandifference(x1, x2) : return(CBnormal_mu(x2) - CBnormal_mu(x1))

"""

# x[i] = Y + error[i],  error[j] ~ F,  F unknown,  Y fixed,  x[i] and error[j] are independent
def CBnonparametric_deconvolution(x, error) : # i.e., the c-box for Y

  def Get_Q( m_in , c_in , k = None) :
    if k is None : k = np.arange((m_in*c_in+1)) 
    def Q_size_GLBL( m ) : return(1 + m + m*(m+1)/2 + m*(m+1)*(m+2)*(3*m+1)/24)
    def Q_size_LoCL( m , c ) : return(1 + c + m*c*(c+1)/2 )
    def Grb_Q( m_in , c_in , Q_list ) : 
      m = max( m_in , c_in )
      c = min( m_in , c_in )
      i_min = Q_size_GLBL( m - 1 ) + Q_size_LoCL( m , c-1 ) + 1
      return(Q_list[i_min:(i_min + m*c)])
  
    def AddingQ( m , Q_list ) :
      Q_list[ Q_size_GLBL( m - 1 ) + 1 ] = 1       
      for c in range(m) :
          i_min = Q_size_GLBL( m - 1 ) + Q_size_LoCL( m , c ) + 1
          Q1 = np.concatenate(( Grb_Q( m-1 , c+1 , Q_list ) , np.repeat(0,(c+1))  ))
          Q2 = np.concatenate(( np.repeat(0,m), Grb_Q( m , c , Q_list )  ))
          Q_list[ i_min:(i_min + m*(c+1)) ] = Q1 + Q2
      return(Q_list[(Q_size_GLBL( m-1 ) + 1):Q_size_GLBL( m )])

    def Bld_Q( m_top ) :
      print('yo')
      Q_out = np.repeat(0,Q_size_GLBL( m_top ))
      Q_out[0] = 1
      for m in range(m_top) :
        Q_out[ (Q_size_GLBL( m ) + 1):(Q_size_GLBL( m+1 )) ] = AddingQ( m+1 , Q_out )
      return(Q_out)

    # body of Get_Q
    m = max( m_in , c_in )
    c = min( m_in , c_in )
    return(Grb_Q(m, c, Bld_Q(m))[k+1])

  
  # body of CBnonparametric_deconvolution
  z = []
  for err in error : z = np.append(z, [x - err])
  z.sort()
  Q = Get_Q(len(x), len(error))
  w = Q / sum( Q )
  return(env(mixture(z,w), mixture(np.append(z[1:],[np.inf]),w)))






def uchenna(kbox, mbox) { # computes the km(k,m) c-box when k and m are themselves c-boxes or p-boxes, as created by gilding for instance
  n = Pbox$steps  # this routine will take about  7 seconds if Pbox$steps is 100, but about a full minute if it is 200
  nk = steps(kbox) 
  nm = steps(mbox)
  Lk = Lm = Rk = Rm = NULL
  for (i in 1:nk) for (j in 1:nm) {  
    Lk = c(Lk, kbox@u[[i]])
    Rk = c(Rk, kbox@d[[i]])
    Lm = c(Lm, mbox@u[[j]])    
    Rm = c(Rm, mbox@d[[j]])   
  u = sort(qbeta(ii(), rep(Lk,      each=n), rep(Rm+1, each=n)))
  d = sort(qbeta(jj(), rep(Rk+1, each=n), rep(Lm,      each=n))) 
  u = u[ii() * n * nk * nm + 1]     
  d = d[jj() * n * nk * nm]         
  return(pbox(u,d))

  

################################################################
# HISTORICAL codes 
################################################################
# P-box for the confidence distribution on a proportion given k successes out of n trials;  error if n<k #
################################################################

balchbox.0 = function(trials, successes, name='') {               # original formulation
  if ((successes==0) && (trials==0)) return(pbox(0,1,shape='beta',name=name)) else
  if (successes==trials) return(pbox(u=qbeta(ii(), successes, trials-successes+1),d=1,shape='beta',name=name)) else
  if (successes==0) return(pbox(u=0,d=qbeta(jj(), successes+1,trials-successes),shape='beta',name=name)) else
  pbox(u=qbeta(ii(), successes, trials-successes+1),d=qbeta(jj(), successes+1,trials-successes),shape='beta',name=name)
  }

balchbox <- function(n, k, name='') {                             # better moments
  if ((k==0) && (n==0)) return(pbox(0,1,shape='beta',name=name))  # could be interval(0,1)
  if (k==n) return(pbox(env(beta(k, n-k+1),1),name=name))
  if (k==0) return(pbox(env(0,beta(k+1,n-k)),name=name))
  return(pbox(env(beta(k, n-k+1),beta(k+1,n-k)),name=name))
  }

balchbox = function(n, k, name='') {                              # slightly faster
  if (n<k) stop('The value of n (',n,') must be larger than or equal to k (',k,') in balchbox')
  uu = function() qbeta(ii(), k, n-k+1)
  dd = function() qbeta(jj(), k+1,n-k)
  if ((k==0) && (n==0)) {u=0;    d=1 }   else
  if (k==n)             {u=uu(); d=1}    else
  if (k==0)             {u=0;    d=dd()} else 
                        {u=uu(); d=dd()}
  pbox(u=u,d=d,shape='beta',name=name)
  }

balch.ci <- function(b,p1=0.025,p2=1-p1) interval(left(cut(b,p1)), right(cut(b,p2)))

#nk <- function(n, k, name='') {   # beta-binomial p-box implied by k successes out of n trials                         
#  if ((k==0) && (n==0)) return(pbox(0,1,shape='beta-binomial',name=name))  # could be interval(0,1)
#  if (k==n) return(pbox(env(BB(k, n-k+1, n), n),name=name))
#  if (k==0) return(pbox(env(0,BB(k+1, n-k, n)),name=name))
#  return(pbox(env(BB(k, n-k+1, n),BB(k+1, n-k, n)),name=name))
#  }
#
# KN <- function(k,n) return(env(beta(k, n-k+1), beta(k+1, n-k)))
# km <- function(k,m) return(env(beta(k, m+1),beta(k+1, m)))


########################################
# Distribution-free p-box constructors #
########################################

# must be able to accept interval arguments

minmax <- function(min, max, name=''){
  pbox(u=rep(min,Pbox$steps), d=rep(max,Pbox$steps), shape='{min, max}', name=name, ml=min, mh=max, vl=0, vh=(max-min)^2/4)
  }

# must be able to accept interval arguments

#minmaxmean <- function(a,b,c) return(mmms(a,b,c,(b-a)^2/4))

minmaxmean <- function(min, max, mean, name=''){
  mid = (max - mean) / (max - min)
  p <- ii()           
  u = ifelse(p <= mid, min, (mean - max) / p + max)
  p <- jj()
  d <- ifelse(mid <= p, max, (mean - min * p) / (1 - p))
  pbox(u=u, d=d, shape='{min, max, mean}', name=name, ml=mean, mh=mean, vl=0, vh=(max-min)*(max-mean)-(max-mean)*(max-mean))
  }

# must be able to accept interval arguments

minmean <- function(min, mean, name=''){
  p <- jjj()
  d <- ((mean - min) / (1-p)) + min
  pbox(u=rep(min,Pbox$steps), d=d, shape='{min, mean}', name=name, ml=mean, mh=mean, vl=0, vh=Inf)
  }

# must be able to accept interval arguments

meanstd <- function(mean, std, name=''){
  p <- iii()         
  u <- mean - std * sqrt(1/p - 1) 
  p <- jjj()
  d <- mean + std * sqrt(p / (1 - p))
  pbox(u=u, d=d, shape='{mean, std}', name=name, ml=mean, mh=mean, vl=std^2, vh=std^2)
  }

meanvar <- meanvariance <- function(mean, var, name='') meanstd(mean, sqrt(var), name)

# must be able to accept interval arguments

posmeanstd <- function(mean, std, name=''){
  p <- ii()          
  u <- pmax(0,mean - std * sqrt(1/p - 1))
  p <- jjj()
  d <- pmin(mean / (1 - p), mean + std * sqrt(p / (1 - p)))
  pbox(u=u, d=d, shape='{positive, mean, std}', name=name, ml=mean, mh=mean, vl=std^2, vh=std^2)
  }

# must be able to accept interval arguments

minmaxmode <- function(min, max, mode, name=''){
  p <- ii()
  u <- p * (mode - min) + min;
  p <- jjj()
  d <- p * (max - mode) + mode
  ml <- (min+mode)/2
  mh <- (mode+max)/2
  vl <- 0
  vh <- (max-min)*(max-min)/12
  pbox(u=u, d=d, shape='{min, max, mode}', name=name, ml=ml, mh=mh, vl=vl, vh=vh)
  }

# must be able to accept interval arguments

minmaxmedian <- function(min, max, median, name=''){
  p <- ii()       
  u <- rep(median, Pbox$steps)
  u[p <= 0.5] <- min
  p <- jjj()
  d <- rep(median, Pbox$steps)
  d[0.5 < p] <- max
  ml <- (min + median)/2
  mh <- (median + max)/2
  vl <- 0
  vh <- (max - min) * (max - min) / 4
  pbox(u=u, d=d, shape='{min, max, median}', name=name, ml=ml, mh=mh, vl=vl, vh=vh)
  }

# must be able to accept interval arguments

minmaxmedianismode <- function(min, max, m, name=''){
  p <- ii()
  u <- rep(m, Pbox$steps)
  u[p <= 0.5] <- p[p <= 0.5] * 2.0 * (m - min) + min
  p <- jjj()
  d <- rep(m, Pbox$steps)
  d[0.5 < p] <- (p[0.5 < p] - 0.5) * 2.0 * (max - m) + m
  ml <- (min + 3 * m) / 4
  mh <- (3 * m + max) / 4
  vl <- 0
  vh <- (max - min) * (max - min) / 4
  pbox(u=u, d=d, shape='{min, max, median=mode}', name=name, ml=ml, mh=mh, vl=vl, vh=vh)
  }

# must be able to accept interval arguments

minmaxpercentile <- function(min, max, fraction, percentile, name=''){
  p <- ii()
  u <- rep(percentile, Pbox$steps)
  u[p <= fraction] <- min
  p <- jjj()
  d <- rep(percentile, Pbox$steps)
  d[fraction < p] <- max
  ml <- fraction * min + (1 - fraction) * percentile
  mh <- fraction * percentile + (1 - fraction) * max
  vl <- 0
  vh <- (max - min) * (max - min) / 4
  pbox(u=u, d=d, shape='{min, max, percentile}', name=name, ml=ml, mh=mh, vl=vl, vh=vh)
  }

# must be able to accept interval arguments

symmeanstd <- function(mean, std, name=''){
  p <- iii()
  u <- rep(mean,Pbox$steps)
  u[p <= 0.5] <- mean - std / sqrt(2 * p[p <= 0.5])
  p <- jjj()
  d <- rep(mean, Pbox$steps)
  d[0.5 < p] <- mean + std / sqrt(2 * (1 - p[0.5 < p]))
  pbox(u=u, d=d, shape='{symmetric, mean, std}', name=name, ml=mean, mh=mean, vl=std^2, vh=std^2)
  }

# must be able to accept interval arguments

maxmean <- function(max, mean, name='') pbox(negate.pbox(minmean(-(max),-(mean))),name=name)

# this mmms function is not correct (it's missing the Berleant-Myers correction); it's also not intervalized...which should be env of (left mean, RIGHT sd) and (right mean, RIGHT sd)
mmms <- function(min,max,mean,std, name='') pbox(imp.pbox(meanstd(mean,std),imp.pbox(minmean(min,mean), maxmean(max,mean))),ml=left(mean), mh=right(mean), vl=left(std)^2,vh=right(std)^2,name=name)

mmmv <- function(a,b,c,d) return(mmms(a,b,c,sqrt(d)))

uniminmax <- function(min,max,mode, name='') pbox(mode + conv.pbox(U(0,1), minmax(min-mode,max-mode), '*'),name=name)

# the following version should be used when interval parameters are supported
unimmmv <- function(min, max, mean, var, mode, name='') {
  pbox(
    imp.pbox(
        mmms(min,max,mean,sqrt(bigger.interval(0, var))),
        mode+conv.pbox(
                        U(0,1),
                        mmms(min-mode,max-mode,2*(mean-mode),sqrt(bigger.interval(0, 3*var-(mean-mode)^2)))
                       )
        )
    ,name=name
    )
  }

unimmmv <- function(min, max, mean, var, mode, name='') {
  #cat('entering unimmmv\nmin:',min,'\nmax:',max,'\nmean:',mean,'\nvar:',var,'\nmode:',mode,'\n')
  #A <- mmms(min,max,mean,sqrt(var))
  #show(A)
  #B <- U(0,1)
  #E <- mmms(min-mode,max-mode,2*(mean-mode),sqrt(3*var-(mean-mode)^2)) 
  #cat('conv.pbox(B,E):')
  #F <- conv.pbox(B,E,'*')
  #show(F)
  #cat('mode+conv.pbox:')
  #G <- mode+F
  #show(G)
  #cat('imp.pbox:')
  #H <- imp.pbox(A,G)
  #show(H)
  pbox(
           imp.pbox(
                     mmms(min,max,mean,sqrt(var)),
                     mode+conv.pbox(
                                    U(0,1),
                                    mmms(min-mode,max-mode,2*(mean-mode),sqrt(3*var-(mean-mode)^2)) 
                                   , '*')
                   )
           ,name=name
          )
  }

# must be able to accept interval arguments

unimmms <- function(min, max, mean, std, mode, name='') unimmmv(min, max, mean, std^2, mode, name=name)

# must be able to accept interval arguments

unimodal <- unimodal.pbox <- function(pb,mode) imp.pbox(pb, unimmms(left(pb),right.pbox(pb),mean.pbox(pb),sd.pbox(pb), mode))

# must be able to accept interval arguments

minmaxmeanismedian <- function(min, max, m, name='') pbox(imp.pbox(minmaxmean(min,max,m), minmaxmedian(min,max,m)),name=name)

# must be able to accept interval arguments

minmaxmeanismode <- function(min, max, m, name='')  pbox(unimodal.pbox(minmaxmean(min,max,m),m),name=name)

whatIknow <- function(    # NOT BEST POSSIBLE 
  min=NULL,
  max=NULL,
  mean=NULL,
  median=NULL,
  mode=NULL,              # implies unimodality too
  std=NULL,
  var=NULL,
  cv=NULL,
  percentiles=NULL,       # array of pairs (percentage, percentile)
  coverages=NULL,
  shape=NULL,             # array of strings that might include 'unimodal', 'symmetric', 'positive' (implies min), 'nonnegative' (implies min), 'concave', 'convex', 'increasinghazard', 'decreasinghazard', 'discrete', 'integervalued', 'continuous', '', '', '', '', 'normal', 'lognormal', etc.
  data=NULL,              # data alone evokes KS confidence bands, data with a named shape evokes the Chen & Iles confidence bands 
  confidence=0.95,
  dependencies=NULL,
  #units=NULL,
  ..., debug=FALSE) {
  #pb <- pbox(dep=dependencies) #, units=units)
  pb <- pbox(dids=dependencies) #, units=units)
  plotting <- c(Pbox$plotting, Pbox$plotting.every)
  Pbox$plotting <- FALSE
  Pbox$plotting.every <- FALSE
  tryCatch( {
  if ('positive' %in% shape)     if (missing(min)) min <- 0 else min <- pmax(min,0)     # should be 0+ (zeroplus) if we had infinitessimals
  if ('nonnegative' %in% shape)  if (missing(min)) min <- 0 else min <- pmax(min,0)
  if (debug) cat('1 ')
  if (present(mode))             shape <- c(shape, 'unimodal')
  if (debug) cat('2 ')
  if (present(min,max))          try(pb <- c(pb,minmax(min,max)))
  if (debug) cat('3 ')
  if (present(min,mean))         try(pb <- c(pb,minmean(min,mean)))
  if (debug) cat('4 ')
  if (present(max,mean))         try(pb <- c(pb,negate(minmean(-max,-mean))))
  if (debug) cat('5 ')
  if (present(min,max,mean))     try(pb <- c(pb,minmaxmean(min,max,mean)))
  if (debug) cat('6 ')
  if (present(min,max,mode))     try(pb <- c(pb,minmaxmode(min,max,mode)))
  if (debug) cat('7 ')
  if (present(min,max,median))   try(pb <- c(pb,minmaxmedian(min,max,median)))
  if (debug) cat('8 ')
  if (present(min,mean,std))     try(pb <- c(pb,min+posmeanstd(mean-min,std)))
  if (debug) cat('9 ')
  if (present(max,mean,std))     try(pb <- c(pb,negate(-max+posmeanstd(-mean+max,std))))
  if (debug) cat('10 ')
  if (present(min,max,mean,std)) try(pb <- c(pb,mmms(min,max,mean,std)))
  if (debug) cat('11 ')
  if (present(min,max,mean,var)) try(pb <- c(pb,mmms(min,max,mean,sqrt(var))))
  if (debug) cat('12 ')
  if (present(mean,std))         try(pb <- c(pb,meanstd(mean,std)))
  if (debug) cat('13 ')
  if (present(mean,cv))          try(pb <- c(pb,meanstd(mean,mean*cv)))
  if (debug) cat('14 ')
  
  if (debug) cat('15 ')
  #minmaxpercentile
  if (debug) cat('16 ')
  
  if (debug) cat('17 ')
# if ('unimodal' %in% shape) { if (!present(mode)) mode <- interval()
  if (debug) cat('18 ')
  if ('unimodal' %in% shape && present(mode)) {
  if (debug) cat('19 ')
    if (present(min,max))          try(pb <- c(pb, uniminmax(min,max,mode)))
  if (debug) cat('20 ')
    if (present(min,max,mean,std)) try(pb <- c(pb, unimmms(min,max,mean,std,mode)))
  if (debug) cat('21 ')
    if (present(min,max,mean,var)) try(pb <- c(pb, unimmmv(min,max,mean,var,mode)))
  if (debug) cat('22 ')
    # can I not otherwise use naked unimodality?
  if (debug) cat('23 ')
    }
  if (debug) cat('24 ')

  if (debug) cat('25 ')
  if ('symmetric' %in% shape && present(mean,std))                     try(pb <- imp(pb,symmeanstd(mean,std)))
  if (debug) cat('26 ')

  if (debug) cat('27 ')
  if (present(min,max,mean,median) && isTRUE(all.equal(mean,median))) try(pb <- c(pb, minmaxmeanismedian(min,max,mean)))
  if (debug) cat('28 ')
  if (present(min,max,mean,mode)   && isTRUE(all.equal(mean,mode)))   try(pb <- c(pb, minmaxmeanismedian(min,max,mean)))
  if (debug) cat('29 ')
  if (present(min,max,median,mode) && isTRUE(all.equal(median,mode))) try(pb <- c(pb, minmaxmedianismode(min,max,median)))
  if (debug) cat('30 ')
  
  if (debug) cat('31 ')
  if ('normal' %in% shape)       try(pb <- c(pb,normal(mean=mean,std=std,var=var,cv=cv,median=median,mode=mode,...)))
  if (debug) cat('32 ')
  if ('lognormal' %in% shape)    try(pb <- c(pb,lognormal(mean=mean,std=std,var=var,cv=cv,median=median,mode=mode...)))
  if (debug) cat('33 ')
    }, finally = {
                  Pbox$plotting <- plotting[[1]]
                  Pbox$plotting.every <- plotting[[2]]
                  if (Pbox$plotting) try(plot.pboxlist(pb))
                 })
  pb
  }

known <- function(...) imp.pbox(whatIknow(...))

"""



##########################################################################
# Functions for elementary logical and Bayes' rule calculations
##########################################################################

def ratiokm(k,m) : return(1/(1+m/k)) # frequentist characterization of the binomial rate

def ratioKN(k,n) : return(k/n) # frequentist characterization of the binomial rate

def jeffkm(k,m) :
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given k successes and m failures 
  # randomly observed in k + m independent Bernoulli trials 
  if ((k < 0)  or (m < 0)) : stop('Improper arguments to function jeffkm')
  return(beta(k+0.5,m+0.5)) 

def jeffKN(k,n) : 
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given only k successes out 
  # of n randomly observed independent Bernoulli trials 
  if ((k < 0)  or (n < k)) : stop('Improper arguments to function jeffKN')
  return(beta(k+0.5,n-k+0.5))
  
def Bppv(p,s,t) : return(1/(1+((1/p-1)*(1-t))/s))

def Bnpv(p,s,t) : return(1/(1+(1-s)/(t*(1/p-1))))

# the mk argments below can be km, KN, jeffkm, jeffKN, ratiokm, or ratioKN

def ppv(pk,pm,sk,sm,tk,tm, mk=km) : return(Bppv(mk(pk,pm), mk(sk,sm), mk(tk,tm)))

def npv(pk,pm,sk,sm,tk,tm, mk=km) : return(Bnpv(mk(pk,pm), mk(sk,sm), mk(tk,tm)))

def ANDi(x,y) :
  nx = len(x)
  ny = len(y)
  if (nx==1) and (ny==1) : return(x*y)
  return(env(x[0:many]*y[0:many], x[many:nx]*y[many:ny]))

def ORi(x,y) :    # x+y-xy == 1-(1-x)*(1-y)
  nx = len(x)
  ny = len(y)
  if (nx==1) and (ny==1) : return(1-(1-x)*(1-y))
  return(1-(1-x[1:many] )*(1-y[1:many]), 1-(1-x[(many+1):nx] )*(1-y[(many+1):ny]))



"""
##########################################################################
# Bayesian inference constructors
##########################################################################

##########################################################################
# Additional compound and conjugate distribution constructors not already in R for use in Bayesian inference

# When there is no data, the posterior is the prior distribution.  
# As the sample size increases, the posterior tends to the data 
# manifested as the likelihood.  Typically, as the prior hyper-
# parmeters (a,b) increase, the posterior grows away from the 
# likelihood and towards the prior.  Generally, the default values 
# for the prior hyperparameters yield the uninformative case. 

# compound distributions already defined
# BB, betabinomial, gammaexponential, NB, negativebinomial, poissonbinomial

BCbernoulli <- function(x,a=0.5,b=0.5,only=TRUE) {
  # data x is an array of zeros and ones (failures and successes)
  # a and b are the parameters of the prior beta
  # see also the km( ) and KN( ) functions
  s = sum(x)
  n = len(x)
  lk = beta(s+1, n-s+1)
  pp = betabinomial(1, s+a, n-s+b)
  if (only) return(pp) else return(list(pr = beta(a, b), po = beta(s+a, n-s+b), pp = pp, lk = lk ))}
##############################################################################    
#par(mfrow=c(3,4))
#options(digits=3)
#x = c(0,0,0,1,1)
#  BCbernoulli(x,only=FALSE) #Jeffreys is the default
#  b = BCbernoulli(x,0,0)  # Haldane prior
#  b; title(paste('Haldane', mean(b)))
#  abline(h=1-sum(x)/len(x),col='green')
#  b = BCbernoulli(x)  # Jeffreys prior
#  b; title(paste('Jeffreys prior', mean(b)))
#  b = BCbernoulli(x,1,1)  # Bayes-Laplace prior
#  b; title(paste('Bayes-Laplace', mean(b)))
#  b = BCbernoulli(x,2,2)  # Walley prior
#  b; title(paste('Walley', mean(b)))
#x = c(1,1,rep(0,8))
#  BCbernoulli(x,only=FALSE) #Jeffreys is the default
#  b = BCbernoulli(x,0,0)  # Haldane prior
#  b; title(paste('Haldane', mean(b)))
#  abline(h=1-sum(x)/len(x),col='green')
#  b = BCbernoulli(x)  # Jeffreys prior
#  b; title(paste('Jeffreys prior', mean(b)))
#  b = BCbernoulli(x,1,1)  # Bayes-Laplace prior
#  b; title(paste('Bayes-Laplace', mean(b)))
#  b = BCbernoulli(x,2,2)  # Walley prior
#  b; title(paste('Walley', mean(b)))
#b = BCbernoulli(x,0,0,only=FALSE); Haldane = b$pr; pl(0,1, xlab='Haldane'); green(Haldane); title('prior')
#b = BCbernoulli(x,only=FALSE);      Jeffreys = b$pr;  pl(0,1, xlab='Jeffreys'); blue(Jeffreys); title('prior')
#b = BCbernoulli(x,1,1,only=FALSE); Laplace = b$pr;  pl(0,1, xlab='Laplace'); black(Laplace); title('prior')
#b = BCbernoulli(x,2,2,only=FALSE); Walley = b$pr;   pl(0,1, xlab='Walley'); gray(Walley); title('prior')

BCbinomial <- function(N, k,n,a=0.5,b=0.5,only=TRUE) {
  # data k is the count of successes, and n is a corresponding number of trials
  # both k and n may be arrays, but they must have the same length
  # a and b are the parameters of the prior beta
  # N is the number of trials to use for the posterior predictive distribution
  # see also the km( ) and KN( ) functions
  s = sum(k)
  sn = sum(n)
  lk = beta(s+1, sn-s+1)
  pp = betabinomial(N, s+a, sn-s+b)
  if (only) return(pp) else return(list(
    pr = beta(a, b),
    po = beta(s+a, sn-s+b),   
    pp = pp,
    lk = lk
    ))}
##############################################################################    
## https://stats.stackexchange.com/questions/512148/beta-binomial-vs-updating-a-prior-beta-distribution
## the BCbernoulli function is specified in terms of 1's (successes) and 0's (failures)
## the BCbinomial function is specified in terms of successes and number of TRIALS
## the beta and betabinomial distributions are specified in terms of numbers of successes and number of FAILURES 
#par(mfrow=c(1,1))
#b = BCbinomial(1,7,7+18,5,12,only=FALSE)
#bb = BCbernoulli(c(rep(1,7),rep(0,18)),5,12,only=FALSE)
#b
#bb
#pl(0,1)
#cyan(beta(5,12));  blue(b$pr);  cyan(bb$pr,lty='dotted')# should all be the same
#black(beta(7+5, 18+12));  gray(b$po);  black(bb$po,lty='dotted')# should all be the same
#khaki(betabinomial(1,12,30));  green(b$pp);  khaki(bb$pp,lty='dotted')# should all be the same
#pl(0,1)
#bc = BCbinomial(1,c(),c(),2,2,only=FALSE)
#blue(bc$pr)
#gray(bc$po)
#green(bc$pp)
##############################################################################    
#par(mfrow=c(3,5))
#options(digits=3)
#N = 10
#k = c(0,2)
#n = c(3,2)  # that is, the Bernoulli trials were {{0,0,0},{1,1}}
#2/5
#sum(k)/sum(n)
#  BCbinomial(N,k,n,only=FALSE) #Jeffreys is the default
#  b1 = BCbinomial(N,k,n,0,0)  # Haldane prior
#  pl(0,N); green(b1); title(paste('Haldane', mean(b1)))
#  b2 = BCbinomial(N,k,n)  # Jeffreys prior
#  pl(0,N); blue(b2); title(paste('Jeffreys prior', mean(b2)))
#  b3 = BCbinomial(N,k,n,1,1)  # Bayes-Laplace prior
#  pl(0,N); black(b3); title(paste('Bayes-Laplace', mean(b3)))
#  b4 = BCbinomial(N,k,n,2,2)  # Walley prior
#  pl(0,N); gray(b4); title(paste('Walley', mean(b4)))
#  pl(0,N); green(b1); blue(b2); black(b3); gray(b4)
## x = c(1,1,rep(0,8)) = {1,1,0,0,0,0,0,0,0,0}
#k = c(2,0,0,0)
#n = c(3,3,2,2)
#2/10
#sum(k)/sum(n)
#  BCbinomial(N,k,n,only=FALSE) #Jeffreys is the default
#  b1 = BCbinomial(N,k,n,0,0)  # Haldane prior
#  pl(0,N); green(b1); title(paste('Haldane', mean(b1)))
#  b2 = BCbinomial(N,k,n)  # Jeffreys prior
#  pl(0,N); blue(b2); title(paste('Jeffreys prior', mean(b2)))
#  b3 = BCbinomialN(k,n,1,1)  # Bayes-Laplace prior
#  pl(0,N); black(b3); title(paste('Bayes-Laplace', mean(b3)))
#  b4 = BCbinomial(N,k,n,2,2)  # Walley prior
#  pl(0,N); gray(b4); title(paste('Walley', mean(b4)))
#  pl(0,N); green(b1); blue(b2); black(b3); gray(b4)
#b1 = BCbinomial(N,k,n,0,0,only=FALSE); Haldane = b1$pr; pl(0,1, xlab='Haldane'); green(Haldane); title('prior')
#b2 = BCbinomial(N,k,n,only=FALSE);      Jeffreys = b2$pr;  pl(0,1, xlab='Jeffreys'); blue(Jeffreys); title('prior')
#b3 = BCbinomial(N,k,n,1,1,only=FALSE); Laplace = b3$pr;  pl(0,1, xlab='Laplace'); black(Laplace); title('prior')
#b4 = BCbinomial(N,k,n,2,2,only=FALSE); Walley = b4$pr;   pl(0,1, xlab='Walley'); gray(Walley); title('prior')
#pl(0,1); green(Haldane); blue(Jeffreys); black(Laplace); gray(Walley); title('prior')

BCpoisson = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default hyperparameters seem to be the uninformative case
  s = sum(x)
  n = len(x)
  lk = gamma(shape=s+1, rate=n)
  pr = gamma(shape=a, rate=b)   #prr = rgamma(many,shape=a,rate=b)
  po = gamma(shape = a + s, rate = b + n)  
  #por = rgamma(r_(r),shape = a + s, rate = b + n)
  #ppr = mc(rpois(MC$many,por))   
  pp = negativebinomial(a + s, 1-1/(1+b + n))  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
##############################################################################    
#doit = function(x,a,b) {
#  pl(0,18)
#  points(x,rep(-0.017,len(x)),col='red')
#  bc = BCpoisson(x,a,b,only=FALSE)
#  blue(bc$pr)
#  gray(bc$po)
#  green(bc$pp)
#  edf(bc$ppr)
#  }
#par(mfcol=c(4,2))
#a = 2; b = 4
#doit(rpois(3000,5),a,b)   # should be the same as poisson(5)
#doit(rpois(10,5),a,b)   
#doit(c(3,4,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior
#a = 1; b = 1/50
#doit(rpois(3000,5),a,b)   # should be the same as poisson(5)
#doit(rpois(10,5),a,b)   
#doit(c(3,4,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior

BCgeometric = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default values for the hyperparameters a and b are the uninformative case, although they will make the function crash if x=NULL
  s = sum(x)
  n = len(x)
  lk = beta(n+1, s+1)
  pr = beta(a,b)  
  po = beta(a + n, b + s)  
  por = rbeta(r_(r), a + n, b + s) 
  pp = mc(rgeom(MC$many,por))   # do we know an analytical formula?  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
##############################################################################    
#doit = function(x,a,b) {
#  pl(0,18)
#  points(x,rep(-0.017,len(x)),col='red')
#  bc = BCgeometric(x,a,b,only=FALSE)
#  blue(bc$pr)
#  gray(bc$po)
#  green(bc$pp)
#  }
#par(mfcol=c(4,2))
#a = 1/2; b = 1/2
#doit(rgeom(3000,0.5),a,b)   # should be the same as geometric(0.5)
#doit(rgeom(10,0.5),a,b)   
#doit(c(3,14,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior
#a = 1; b = 10
#doit(rgeom(3000,0.5),a,b)   # should be the same as geometric(0.5)
#doit(rgeom(10,0.5),a,b)   
#doit(c(3,14,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior

BCuniform.knownmin = BCuniform = function(x,A,a=A,b=A+1,r=runif(MC$many),only=TRUE) {  
  #  x_i ~ uniform(A,theta), that is, from a uniform distribution whose minimum is A and whose maximum needs to be established
  lk = pareto(max(x), len(x)+1)
  pr = pareto(a, b)
  po = pareto(max(a,max(x)), b+len(x)) # Masatoshi says max(x) is m
  por = qpareto(r_(r),max(a,max(x)), b+len(x)) # Masatoshi says max(x) is m
  pp = mc(runif(MC$many,A,por))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
#x = runif(25,5,13)
#bc = BCuniform(x,A=5,a=1,b=1,only=FALSE) 
#bc

BCnegativebinomial = function(x,R,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default hyperparameters are the uninformative case
  s = sum(x)
  n = len(x)
  lk = beta(R*n+1, s+1)
  pr = beta(a, b)
  po = beta(a+R*n, b+s)
  por = rbeta(r_(r),a+R*n, b+s)
  pp = mc(rnbinom(MC$many,R,por))   # do we know an analytical formula?  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

BCnormal.knownsigma = function(x,sigma,m0=mean(x),s0=sd(x),r=runif(MC$many),only=TRUE) {
  # the default hyperparameters seem to be uninformative (using s0=sd(x)/sqrt(len(x) makes the posterior differ from the likelihood more strongly)
  # increasing s0 makes the prior more uninformative, unlike the typical behaviour of (a,b) hyperparameters in other functions
  s = sum(x)
  n = len(x)
  pr = normal(m0,s0)  
  lk = normal(s/n,sigma/sqrt(n))
  mprime = (m0/s0^2+s/sigma^2)/(1/s0^2+n/sigma^2)
  sprime = 1/sqrt(1/s0^2+n/sigma^2)
  if (abs(s0<1e-20)) {mprime = m0; sprime = 0}
  po = normal(mprime, sprime)
  por = rnorm(r_(r),mprime, sprime)
  pp = normal(mprime, sqrt(sprime^2 + sigma^2))    # pp = mc(rnorm(MC$many,por,sigma)) 
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

BCnormal.knownmu = function(x,mu,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the hyperparameter defaults are the uninformative case (so the pr won't be defined as it would be improper)
  n = len(x)
  s = sum((x-mu)^2)
  pr = sqrt(inversegamma(shape=a,rate=b)) # we parameterize N with sd, not var
  po = sqrt(inversegamma(shape=a + n/2, rate=b + s/2))
  por = sqrt(1/rgamma(MC$many,shape=a + n/2, rate=b + s/2))
  pp = mc(rnorm(MC$many,mu,por))     
  #ppa = mean(x)+sd(x)*sqrt(1+1/n)*student(n-1)   
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}
#  
#x = rnorm(1000,15,2)
#bc = BCnormal.knownmu(x,15,5,5)
#bc
#edf(x)
#  
#x = rnorm(10,15,2)
#bc = BCnormal.knownmu(x,15,5,5)
#bc
#edf(x)
#  
#x = rnorm(1000,15,2)
#bc = BCnormal.knownsigma(x,2,10,2)
#bc
#edf(x)
#  
#x = rnorm(10,15,2)
#bc = BCnormal.knownsigma(x,2,10,2)
#bc
#edf(x)

rnormgamma <- function(n, mu, lambda, alpha, beta) {
  # normal-gamma deviates: (1) Sample tau from a gamma distribution with parameters alpha and beta, (2) Sample x from a normal distribution with mean mu and variance 1/(lambda * tau)
  # E(x) = mu;  E(tau) = alpha/beta
  if (len(n) > 1) n = len(n)
  tau <- rgamma(n, alpha, beta)
  x <- rnorm(n, mu, sqrt(1/(lambda*tau)))    # tau and x are NOT independent, as tau is used to compute x
  #cat(mean(x), mu,'  ', mean(tau), alpha/beta,'\n')
  data.frame(x = x, tau = tau)                              
}

seepr = function(m,l,a,b) {pr = rnormgamma(MC$many,m, l, a, b); prx = mc(pr$x); prt = mc(pr$t); pl(min(left(prx),left(prt)), max(right(prx),right(prt))); blue(prx); cyan(prt); title(paste(m,l,a,b)) }

BCnormal = function(x, mu0=mean(x), lambda0=1, alpha0=1, beta0=1, only=TRUE) {
  # x_i ~ N(mu, 1/tau)
  # The selection of the prior for a normal involves chosing values for 4 parameters.  The first argument 
  # mu0 is your guess about the mean of the normal data, and second is related to the dispersion of this 
  # estimate about the mean.  You can use the seepr( ) function to visualize the priors you select for 
  # BCnormal.  In practice, setting mu0 to mean(x) and the other values to one seems to often give 
  # reasonable results, so they have been specified as defaults, but I'm sure this violates some crucial
  # Bayesian stricture I'm not aware.
  pr = rnormgamma(MC$many,mu0, lambda0, alpha0, beta0)  # hmm...I guess pr and po and pp should all be correlated, right?  And they should have an intended correlation with ab outside caller via r_(r)
  n = len(x)
  xbar = mean(x)
  s = sd(x)
  po = rnormgamma(MC$many, (lambda0*mu0 + n*xbar)/(lambda0+n), lambda0+n, alpha0+n/2, beta0+(n*s+(lambda0*n*(xbar-mu0)^2)/(lambda0 + n))/2 )
  pp = mc(rnorm(MC$many,po$x,1/po$tau))
  pr = list(x=mc(pr$x), tau=mc(pr$tau))
  po = list(x=mc(po$x), tau=mc(po$tau))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}
#x = rnorm(20,15000,2)
#m = mean(x)
#l = 1
#a = 1
#b = 1
#seepr(m,l,a,b)
#bc = BCnormal(x,m,l,a,b)
#bc
#edf(x)
#mean(bc)
#sd(bc)

BCexponential = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the prior and posterior estimate the MEAN of the exponential
  # which is the reciprocal of its RATE parameter used by rexp()
  # the default hyperparameters are the uninformative case
  sm = 'exponential(theta)'
  s = sum(x)
  n = len(x)
  lk = gamma(shape=n+1,rate=s)         # reciprocated in the returned list
  pr = gamma(shape=a, rate=b)           # reciprocated in the returned list
  po = gamma(shape=a+n, rate=b+s)  # reciprocated in the returned list
  por = rgamma(r_(r), shape=a+n,  rate=b+s)
  pp = mc(rexp(MC$many, por))
  if (only) return(pp) else return(list(pr = 1/pr, po = 1/po, pp = pp, lk = 1/lk))}

BCpareto.knownmin = BCpareto = function(x,xm,a=0,b=0,r=runif(MC$many),only=TRUE) { 
  # default hyperparameters seem to be the uninformative case
  s = sum(log(x/xm))
  n = len(x)
  sm = 'pareto(xm,theta)'
  lk = gamma(shape=n, rate=s)  # this is just a guess; prolly wrong as po not always between pr and this lk
  #pr = gamma(a, b)
  #po = gamma(1/(1/a+s), b+n) # po = gamma(1/(a+s), b+len(n))
  #por = rgamma(r_(r), 1/(1/a+s), b+n)
  #pp = mc(qpareto(r_(r),xm,por))
  pr = gamma(shape=a, rate=b)
  po = gamma(shape=a+n, rate=b+s) 
  por = rgamma(r_(r), shape=a+n, rate=b+s)
  pp = mc(qpareto(r_(r),xm,por))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

#BCweibull <- function(x,a,b,c,d0,v=0) {
#  # x_i ~ weibull(shape=beta, scale=theta)
#  # https://www.johndcook.com/CompendiumOfConjugatePriors.pdf, page 33f
#  
#  n = len(x)
#  #sxb = sum(x^beta)
#  px = prod(x)
#  #lk = (beta/theta)^n * prod(x) ^ (beta-1) * exp(- sum(x^beta) / theta)
#  
#  #pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); ifelse((0<beta) and (0<theta), beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K,0)}
#  pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K}
#  #need to evaluate normalization factor K
#  
#  po = function(beta,theta) pr(beta,theta, a=a+n, b=b-log(px), c=c+n, d0=d0, v=n) 
#
#  # the pr and po are bivariate, but the pp is just univariate, so it seems like we should be able to smash [not marginalize, right?] all the pr and po values into a single array with which to create the pp
#  # or should pobeta and potheta just be the marginalizations?
#
#  pp = mc(qweibull(r_(r),shape=pobeta,scale=potheta))
#  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}

#x = qweibull(runif(30), shape=2, scale=3)
#x = c(0.9, 1.52, 1.10)
#a = 20.0
#b = 2.0
#c = 6.0
#d0 = 2.5
#v=0
#  n = len(x)
#  px = prod(x)
#pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K} #need to evaluate normalization factor K
#po = function(beta,theta) pr(beta,theta, a=a+n, b=b-log(px), c=c+n, d0=d0, v=n) 
  
# old naming scheme  
#bernoulli.uniform = function(x) {s = sum(x); beta(1+s,1+len(x)-s)}
#bernoulli.beta = function(x,a,b) {s = sum(x); beta(a+s, b+len(x)-s)
#binomial.beta = function(x,N,a,b) {s = sum(x); beta(a+s, b+sum(N)-s)
#negativebinomial.beta = function(x,r,a,b)  beta(a+r*len(x), b+sum(x))
#poisson.gamma = function(x,k,theta) gamma(k + sum(x), 1/(len(x)+1/theta))
#poisson.gamma = function(x,a,b) gamma(a + sum(x), b + len(x))
#hypergeometric.betabinomial = function(x,N,a,b) betabinomial(a+sum(x), b+sum(N) - sum(x))
#geometric.beta = function(x,a,b) beta(a+len(x), b + sum(x))
#normal.normal = function(s, mu, sigma) {v=1/(sigma^2) + len(x)/s; normal((mu/(sigma^2) + sum(x)/(s^2)/v, 1/sqrt(v))) }
#exponential.gamma = function(x,a,b) gamma(a+len(x),b+sum(x))



###############################################################################
# Maximum a posteriori 
###############################################################################


# Wikipedia [https://en.wikipedia.org/wiki/Maximum_likelihood_estimation] says
# From the perspective of Bayesian inference, MLE is generally equivalent to 
# maximum a posteriori (MAP) estimation with a prior distribution that is uniform 
# in the region of interest. In frequentist inference, MLE is a special case of 
# an extremum estimator, with the objective function being the likelihood.


# Piech and Sahami (2017) https://web.stanford.edu/class/archive/cs/cs109/cs109.1196/reader/11%20Parameter%20Estimation.pdf
# describe the Maximum A Posteriori (MAP) method for parameter estimation.
# The paradigm of MAP is that we should chose the value for our parameters
# that is the most likely given the data.  The MAP estimate is the mode of 
# the posterior distribution for the parameter.  Parameters for which conjugate
# priors are known include
# 
# Parameter           Prior distribution
# Bernoulli p         Beta
# Binomial p          Beta
# Poisson λ           Gamma
# Exponential λ       Gamma
# Multinomial pi      Dirichlet
# Normal µ            Normal
# Normal σ2           Inverse Gamma


"""



##########################################################################
# Maximum entropy distribution constructors 
##########################################################################

def MEminmax(min, max) : return(uniform(min,max))

def MEminmaxmean(min, max, mean) : return(sawinconrad(min,mean,max)) #http://mathoverflow.net/questions/116667/whats-the-maximum-entropy-probability-distribution-given-bounds-a-b-and-mean, http://www.math.uconn.edu/~kconrad/blurbs/analysis/entropypost.pdf for discussion of this solution.

def MEmeansd(mean, sd) : return(normal(mean, sd))

def MEminmean(min,mean) : return(min+exponential(mean-min))

#def MEmeangeomean(mean, geomean)

def MEdiscretemean(x,mu,steps=10,iterations=50) : # e.g., MEdiscretemean(1:10,2.3)
  x = np.array(x)
  def fixc(x,r) : return(1/np.sum(r**x))
  r = br = 1
  c = bc = fixc(x,r)
  d = bd = (mu - np.sum((c*r**x)*x))**2
  for j in range(steps) :
    step = 1/(j+1)
    for i in range(iterations) :
      r = np.abs(br + (np.random.uniform() - 0.5) * step)
      c = fixc(x,r)
      d = (mu - np.sum((c*r**x)*x))**2
      if d < bd :
        br = r
        bc = c
        bd = d
  w = bc*br**x
  w = w / np.sum(w) # needed? 
  
  z = np.array([])
  k = len(x)
  for i in range(k) : z = np.concatenate((z, np.repeat(x[i],w[i]*many)))
  if len(z)>=many :
      z = z[0:many] 
  else : 
      z = np.concatenate((z, np.random.choice(x, size=many-len(z), p=w)))
  np.random.shuffle(z)   # shuffles z in place
  return(z)

def MEquantiles(v,p) :
  if len(v) != len(p) : stop('Inconsistent array lengths for quantiles')
  if (min(p) < 0) or (1 < max(p)) : stop('Improper probability for quantiles') # ensure 0 <= p <= 1
  if not (min(p) == 0 and max(p)==1) : stop('Probabilities must start at zero and go to one for quantiles')
  if (np.any(np.diff(p)<0)) : stop('Probabilities must increase for quantiles') # ensure montone probabilities
  if (np.any(np.diff(v)<0)) : stop('Quantiles values must increase') # ensure montone quantiles
  x = np.repeat(np.inf,many)
  r = np.random.uniform(size=many)
  # np.where is Python's version of R's ifelse function
  for i in range(len(p)-1) : 
      x = np.where((p[i]<=r) & (r<p[i+1]), v[i]+(r-p[i])*(v[i+1]-v[i])/(p[i+1]-p[i]),  x)
  return(x)

def MEdiscreteminmax(min,max) : return(np.minimum(np.trunc(uniform(min,max+1)),max))

def MEmeanvar(mean, var) : return(MEmeansd(mean,np.sqrt(var)))

def MEminmaxmeansd(min, max, mean, sd) : return(beta1((mean - min) / (max - min),  sd/(max - min) ) * (max - min) + min)

def MEmmms(min, max, mean, sd) : return(beta1((mean - min) / (max - min),  sd/(max - min) ) * (max - min) + min)

def MEminmaxmeanvar(min, max, mean, var) : return(MEminmaxmeansd(min,max,mean,np.sqrt(var)))



###############################################################################
# Miscellaneous: PERT, Fermi methods, mean-normal range, KS, EDF, Antweiler
###############################################################################

def antweiler(x) : return(triangular(min=min(x), mode=3*np.mean(x)-max(x)-min(x), max=max(x))) # https://wernerantweiler.ca/blog.php?item=2019-06-05       #**

def betapert(min, max, mode) :  # N.B.  Not in numerical order!
  mu = (min + max + 4*mode)/6
  if (abs(mode-mu)<1e-8) :
      alpha1 = alpha2 = 3 
  else :
      alpha1 = (mu - min)*(2*mode - min - max)/((mode - mu)*(max - min))
      alpha2 = alpha1*(max - mu)/(mu - min) 
  return(min + (max - min) * beta(alpha1, alpha2))  

def mnr(n,many=10000) :
   xL = xU = np.random.normal(size=many)
   for i in range(n-1) : 
     xx = np.random.normal(size=many)
     xL = np.minimum(xL,xx)
     xU = np.maximum(xU,xx) 
   return(np.mean(xU - xL))

def fermilnorm(x1, x2, n=None, pr=0.9) :
   gm = np.sqrt(x1*x2)
   if n is None : gsd = np.sqrt(x2/x1) ** (1/sps.norm.ppf(pr)) # qnorm(pr)
   else : gsd = np.exp((np.log(x2) - np.log(x1)) / mnr(n))
   return(np.log((gm, gsd)))

def ferminorm(x1, x2, n=None, pr=0.9) :
   m = (x1 + x2) / 2
   if n is None : s = (x2 - x1) / (2 * sps.norm.ppf(pr))  # qnorm(pr) 
   else : s = (x2 - x1) / mnr(n)
   return(np.array((m, s)))
 
def approxksD95(n) :
    from scipy.interpolate import CubicSpline
    # approximations for the critical level for Kolmogorov-Smirnov statistic D,
    # for confidence level 0.95. Taken from Bickel & Doksum, table IX, p.483
    # and Lienert G.A.(1975) who attributes to Miller,L.H.(1956), JASA
    if n > 80 : return(1.358 /(np.sqrt(n) + .12 + .11/np.sqrt(n))) # Bickel&Doksum, table IX,p.483
    else :
        x = np.array((1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80)) # from Lienert
        y = np.array((.975,   .84189, .70760, .62394, .56328, # 1:5
                      .51926, .48342, .45427, .43001, .40925, # 6:10
                      .33760, .29408, .24170, .21012,         # 15,20,30,40
                      .18841, .17231, .15975, .14960))         # 50,60,70,80
        f = CubicSpline(x, y, bc_type='natural')
        return(f(n))

def ks(x, conf=0.95, min=None, max=None) :
  if conf != 0.95 : stop('Cannot currently handle confidence levels other than 95%')
  h = histogram(x)
  mn = np.min(x)
  mx = np.max(x)
  if min is None : min = mn - (mx-mn) / 2
  if max is None : max = mx + (mx-mn) / 2
  lots = int(approxksD95(len(x)) * many)
  Lfermi = np.concatenate((lots*[min],h))
  Rfermi = np.concatenate((h,lots*[max]))
  Lfermi.sort()
  Rfermi[::-1].sort()  
  return(np.concatenate((Lfermi[0:many],Rfermi[0:many]))) # should prolly shuffle
  
def ferminormconfband(x1, x2, n, pr=0.9, conf=0.95, bOt=0.001, tOp=0.999) :
  if conf != 0.95 : stop('Cannot handle confidence levels other than 95%')
  m,s = ferminorm(x1,x2,n,pr)
  lots = int(approxksD95(n) * many)
  BOT = sps.norm.ppf(bOt,m,s)
  TOP = sps.norm.ppf(tOp,m,s)
  Lfermi = np.concatenate((lots*[BOT],sps.norm.rvs(m,s,size=many)))
  Rfermi = np.concatenate((lots*[TOP],sps.norm.rvs(m,s,size=many)))
  Lfermi.sort()
  Rfermi[::-1].sort() 
  return(np.concatenate((Lfermi[0:many],Rfermi[0:many]))) # should prolly shuffle
  
def fermilnormconfband(x1, x2, n, pr=0.9, conf=0.95, bOt=0.001, tOp=0.999) :
  if conf != 0.95 : stop('Cannot handle confidence levels other than 95%')
  mlog,slog = fermilnorm(x1,x2,n,pr)
  d = lognormal2(mlog,slog)
  lots = int(approxksD95(n) * many)
  BOT = sps.lognorm.ppf(bOt,s=slog,scale=np.exp(mlog))
  TOP = sps.lognorm.ppf(tOp,s=slog,scale=np.exp(mlog))
  Lfermi = np.concatenate((lots*[BOT],d))
  Rfermi = np.concatenate((lots*[TOP],d))
  Lfermi.sort()
  Rfermi[::-1].sort() 
  return(np.concatenate((Lfermi[0:many],Rfermi[0:many]))) # should prolly shuffle



"""

OPi = function(x,y,op) { # op can be '+', '-', '*', '/', '^', 'pmin', or 'pmax'
  nx = len(x)
  ny = len(y) 
  if (op=='-') return(OPi(x, c(-y[(many+1):ny], -y[1:many]),'+'))
  if (op=='/') return(OPi(x, c(1/y[(many+1):ny], 1/y[1:many]),'*'))
  if ((nx==1) and (ny==1)) return(do.call(op,list(x,y)))
  if (nx==1) return(c(do.call(op,list(x,y[1:many])), do.call(op,list(x,y[(many+1):ny]))))
  if (ny==1) return(c(do.call(op,list(x[1:many],y)), do.call(op,list(x[(many+1):nx],y))))
  c(do.call(op,list(x[1:many],y[1:many])), do.call(op,list(x[(many+1):nx],y[(many+1):ny])))
  }

opi = function(x,y,op) {  # in case it is not obvious what OPi is doing
  nx = len(x)
  ny = len(y) 
  if ((nx==1) and (ny==1)) return(do.call(op,list(x,y)))
  if (nx==1) return(opi(rep(x,2*many),y))
  if (ny==1) return(opi(x,rep(y,2*many)))
  if (op=='+') return(c(x[1:many]+y[1:many], x[(many+1):nx]+y[(many+1):ny]))
  if (op=='-') return(opi(x, c(-y[(many+1):ny], -y[1:many]),'+'))
  if (op=='*') return(c(x[1:many]*y[1:many], x[(many+1):nx]*y[(many+1):ny]))
  if (op=='/') return(opi(x, c(1/y[(many+1):ny], 1/y[1:many]),'*'))
  if (op=='^') return(c(x[1:many]^y[1:many], x[(many+1):nx]^y[(many+1):ny]))
  if ((op=='min') || (op=='pmin')) return(c(pmin(x[1:many],y[1:many]), pmin(x[(many+1):nx],y[(many+1):ny])))
  if ((op=='max') || (op=='pmax')) return(c(pmax(x[1:many],y[1:many]), pmax(x[(many+1):nx],y[(many+1):ny])))
  stop('ERROR unknown operator in opi')
  }

plotbox <- function(b,new=TRUE,col='blue',lwd=2,xlim=range(b[is.finite(b)]),ylim=c(0,1),xlab='',ylab='Prob',...) {
  edf <- function (x, col, lwd, ...) {
      n = len(x)
      s <- sort(x)
      lines(c(s[[1]],s[[1]]),c(0,1/n),lwd=lwd,col=col,...)
      for (i in 2:n) lines(c(s[[i-1]],s[[i]],s[[i]]),c(i-1,i-1,i)/n,col=col,lwd=lwd,...)
      }
  b = ifelse(b==-Inf, xlim[1] - 10, b)
  b = ifelse(b==Inf, xlim[2] + 10, b)
  if (new) plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
  if (len(b) < many) edf(b,col,lwd) else
  edf(c(min(b),max(b),b[1:min(len(b),many)]),col,lwd)
  if (many < len(b)) edf(c(min(b),max(b),b[(many+1):len(b)]),col,lwd)
  }

"""





###############################################################################
# Testing and exercising the functions
###############################################################################

###############################################################################
# Maximum likelihood estimation for lognormal using sps.fit()

# So far, I'm not impressed with sps ML fitting. What am I doing wrong?
# It looks like I'm not doing anything wrong.  It just doesn't work so well 
# with really small data sets.

print('****** 1')

# N = 10
w = np.array([2.912,2.5565,2.9077,4.6462,3.5,2.2677,4.6362,3.017,3.9792,4.6102])
print(sps.lognorm.fit(w)) # (10.6565, 2.26770, 0.03194)
L = sps.lognorm.rvs(s=10.656, loc=2.2677, scale=0.0319,size=many)
LL = MLlognormal(w)
LLL = sMLlognormal(w)
pl((-100,1e3)); edf(L,'r'); edf(LL); edf(LLL,'g'); #edf(w)
edf(LLL,'g'); edf(w)
print(np.mean(L));print(np.mean(LL));print(np.mean(LLL));print(np.mean(w))
plt.show()

# N = 30
W = sps.lognorm.rvs(s=2,scale=np.exp(3),size=30)
print(sps.lognorm.fit(W)) 
L = sps.lognorm.rvs(*sps.lognorm.fit(W),size=many)
LL = MLlognormal(W)
LLL = sMLlognormal(W)
pl((-100,1e3)); edf(L,'r'); edf(LL); edf(LLL,'g'); edf(W,'k')
#edf(LLL,'g'); edf(W)
print(np.mean(L));print(np.mean(LL));print(np.mean(LLL));print(np.mean(W))
plt.show()

# N = 100
W = sps.lognorm.rvs(s=2,scale=np.exp(3),size=100)
print(sps.lognorm.fit(W)) 
L = sps.lognorm.rvs(*sps.lognorm.fit(W),size=many)
LL = MLlognormal(W)
LLL = sMLlognormal(W)
pl((-100,1e3)); edf(L,'r'); edf(LL); edf(LLL,'g'); edf(W,'k')
#edf(LLL,'g'); edf(W)
print(np.mean(L));print(np.mean(LL));print(np.mean(LLL));print(np.mean(W))
plt.show()

###############################################################################
# Miscellaneous PERT and Fermi estimates

print('****** 2')

x = np.random.normal(size=25)
edf(antweiler(x))
 
edf(betapert(5, 10, 6),'r')

ferminorm(12, 16)  # array([14., 1.560])
edf(normal(*ferminorm(12, 16)),'g') 

fermilnorm(16, 32)  # array([3.12, 0.27])
edf(lognormal2(*fermilnorm(16, 32)),'k') 

plt.show()

###############################################################################
# Fermi and KS confidence bands

print('****** 3')

bOt = 0.001
tOp = 0.999
m,s = ferminorm(2,10,100,.9)
BOT = sps.norm.ppf(bOt,m,s)
TOP = sps.norm.ppf(tOp,m,s)
n=normal(m,s)
f=ferminormconfband(2,10,100)
red(n); edf(f)

m,s = ferminorm(12,20,100,.9)
n=normal(m,s)
f=ferminormconfband(12,20,100,bOt=0.00001,tOp=0.99999)
red(n); edf(f)

mlog,slog = fermilnorm(22,27,50,.9)
n = lognormal2(mlog,slog)
f = fermilnormconfband(22,27,50)  # n = 50, lognormal
red(n); edf(f)

w = 36 + 2*np.random.normal(size=25)
k = ks(w)
edf(k); red(w)

plt.show()

###############################################################################
# beta distribution constructors should be able to handle funky parameters

print('****** 4')

x = beta(1,1)   # U(0,1)
print(x.mean()) # 0.5

x = beta(1,0)   # 1
print(x.mean()) # 1

x = beta(0,1)   # 0
print(x.mean()) # 0

x = beta(0,0)   # bernoulli(0.5), or should it be [0,1]?   
print(x.mean()) # 0.5, or the interval [0,1]

###############################################################################
# Pareto

print('****** 5')

# picture in Wikipedia https://en.wikipedia.org/wiki/Pareto_distribution
m=1;pl((0,5)); edf(pareto(m,1));edf(pareto(m,2));edf(pareto(m,3));edf(pareto(m,3000));plt.title('pareto');plt.show() 

m=2;pl((0,5)); edf(pareto(m,1));edf(pareto(m,2));edf(pareto(m,3));edf(pareto(m,3000));plt.title('pareto');plt.show()

###############################################################################
# power function distribution

print('****** 6')

pl((0,1)); pl((0,1)); pl((0,1))

edf(powerfunction(1,1)); plt.show()
edf(powerfunction(2,1)); plt.show()
edf(powerfunction(3,1)); plt.show()
edf(powerfunction(4,1)); plt.show()
    
edf(powerfunction(1,2)); plt.show()
edf(powerfunction(2,2)); plt.show()
edf(powerfunction(3,2)); plt.show()
edf(powerfunction(4,2)); plt.show()

edf(powerfunction(1,3)); plt.show()
edf(powerfunction(2,3)); plt.show()
edf(powerfunction(3,3)); plt.show()
edf(powerfunction(4,3)); plt.show()

###############################################################################
# Laplace

print('****** 7')

# picture on Wikipedia https://en.wikipedia.org/wiki/Laplace_distribution
m=0;pl((-10,10));edf(laplace(m,1));edf(laplace(m,2));edf(laplace(m,4));edf(laplace(-5,4));plt.title('laplace');plt.show()  
 
###############################################################################
# I don't think the fits to loguniform are correct

print('****** 8')

w = np.random.uniform(2,5,size=200); edf(w); edf(MMloguniform(w)); edf(MLloguniform(w)); plt.show()
w = np.random.uniform(2,5,size=20); edf(w); edf(MMloguniform(w)); edf(MLloguniform(w)); plt.show()

###############################################################################
# are we sure the asterisk operator works?  (it's the best thing about Python)
# the red and green fitted distributions should be the same, modulo MC error

print('****** 9')

def AMLgumbel(x) :
    loc, scale = sps.gumbel_r.fit(x)
    return(gumbel(loc,scale))

def BMLgumbel(x) : return(gumbel(*sps.gumbel_r.fit(w)))    

w = np.array([2.91247063, 2.55651104, 2.90768457, 4.64622234, 3.49995966,
       2.26770086, 4.63619271, 3.01703563, 3.97919485, 4.61017778,
       2.00292333, 3.13348299, 4.68998771, 2.30031397, 2.14102056,
       4.23825192, 2.56982047, 4.86396995, 3.79969706, 4.00203139])
g1 = AMLgumbel(w); g2 = BMLgumbel(w)
edf(w); edf(g1,'r'); edf(g2,'k'); plt.show()

w = np.array([2.91247063, 2.55651104, 2.90768457, 4.64622234, 3.49995966,
       2.26770086, 4.63619271, 3.01703563, 3.97919485, 4.61017778])
g1 = AMLgumbel(w); g2 = BMLgumbel(w)
edf(w); edf(g1,'r'); edf(g2,'k'); plt.show()

G1,G2 = sps.gumbel_r.fit(w)
print(G1,G2)
print(*sps.gumbel_r.fit(w))
print(sps.gumbel_r.fit(w))

###############################################################################
# parameterisations: Rayleigh, power function, gamma, inverse gamma, gammaexponential 

print('****** 10')

# parameterization of Rayliegh distributions does NOT match Risk Calc and pba.r

# raleigh() ***************************************


# parameterization for power function distributions seems to match Risk Calc and pba.r

edf(powerfunction(1,1));edf(powerfunction(2,1));edf(powerfunction(3,1));edf(powerfunction(4,1))
edf(powerfunction(1,2));edf(powerfunction(2,2));edf(powerfunction(3,2));edf(powerfunction(4,2))
edf(powerfunction(1,3));edf(powerfunction(2,3));edf(powerfunction(3,3));edf(powerfunction(4,3));plt.show()

# parameterizations for gamma and inversegamma distributions don't match Risk Calc
# but I've just updated pba.r [15 Nov 2024] so it now agrees with these conventions

# The parameterizations for gamma and inversegamma distributions 
# match with their Wikipedia articles (as of 15 November 2024).
# Note that the wrinkle is that one's scale is the other's rate.

######### examples from the Wikipedia pages

print('display Reyleigh examples')
# https://en.wikipedia.org/wiki/Rayleigh_distribution
#blue(rayleigh(0.5)); green(rayleigh(1)); red(rayleigh(2)); cyan(rayleigh(3)); purple(rayleigh(4)); plt.show()

# https://en.wikipedia.org/wiki/Gamma_distribution
red(gamma(shape=1, scale=2));orange(gamma(shape=2, scale=2));yellow(gamma(shape=3, scale=2));green(gamma(shape=5, scale=1)); black(gamma(shape=9, scale=0.5));blue(gamma(shape=7.5, scale=1));purple(gamma(shape=0.5, scale=1));plt.show()

# https://en.wikipedia.org/wiki/Inverse-gamma_distribution
red(inversegamma(shape=1, scale=1));green(inversegamma(shape=2, scale=1));blue(inversegamma(shape=3, scale=1));cyan(inversegamma(shape=3, scale=0.5)); plt.show()

# there is no Wikipedia article on gammaexponential, so instead we can check the
# picture created below with the picture made by R (immediately further down)
#
# Python
many = 100000
g = gammaexponential(shape=1, scale=1); pl(0,20); blue(g)
g = gammaexponential(shape=2, scale=1); pl(0,20); red(g)
g = gammaexponential(shape=1, scale=2); pl(0,20); black(g)
g = gammaexponential(shape=2, scale=2); pl(0,20); yellow(g)
g = gammaexponential(shape=1, scale=11); pl(0,20); cyan(g)
g = gammaexponential(shape=1, scale=0.1); pl(0,20); purple(g)

# # R
# source('pba BETTER.r')
# rbyc()
# pl(0,20)
# g = gammaexponential(shape=1, scale=1);   blue(g)
# g = gammaexponential(shape=2, scale=1);   red(g)
# g = gammaexponential(shape=1, scale=2);   black(g)
# g = gammaexponential(shape=2, scale=2);   yellow(g)
# g = gammaexponential(shape=1, scale=11);  cyan(g)
# g = gammaexponential(shape=1, scale=0.1); purple(g)

######### moments are correct now
shape = 4
scale = 6
rate = 1/scale
#rbyc(3,2)
tOp = 0.99999

print('display Reyleigh examples')
#r = rayleigh(scale)
#print(scale * np.sqrt(np.pi/2),   r.mean())           # mean
#print((4-np.pi)*scale**2/2,       r.var())            # var 

p = powerfunction(scale,shape)  # Risk Calc is the reference for the moments [scale=b,shape=c]
print(scale/(1+1/shape),   p.mean())                  # mean
print(scale**2/((1+2/shape)*(shape+1)**2), p.var())   # var # maybe not defined, unless truncated

g = gamma(shape=shape, scale=scale)
print(shape * scale,      g.mean())                   # mean
print(shape * scale**2,   g.var())                    # var 

ig = inversegamma(shape=shape, scale=scale)
print(scale / (shape - 1),               ig.mean())   # mean
print(scale**2 / ((shape-1)**2*(shape-2)), ig.var())    # var

######### reciprocation

def compare(x,y) : edf(x,'b',lw=7); edf(y,'y',lw=3); plt.show()

SCALE = 6
SHAPE = 4
RATE = 1/SCALE

compare(1/zbuff(gamma(shape=SHAPE, rate=RATE)), inversegamma(shape=SHAPE, rate=1/RATE))
compare(1/zbuff(gamma(shape=SHAPE, rate=RATE)), inversegamma(shape=SHAPE, scale=1/SCALE))
compare(1/zbuff(gamma(shape=SHAPE, rate=RATE)), inversegamma(shape=SHAPE, rate=SCALE))
compare(1/zbuff(gamma(shape=SHAPE, rate=RATE)), inversegamma(shape=SHAPE, scale=RATE))
compare(1/zbuff(gamma(shape=SHAPE, scale=SCALE)), inversegamma(shape=SHAPE, scale=1/SCALE))
compare(1/zbuff(gamma(shape=SHAPE, scale=SCALE)), inversegamma(shape=SHAPE, rate=1/RATE))
compare(1/zbuff(gamma(shape=SHAPE, scale=SCALE)), inversegamma(shape=SHAPE, scale=RATE))
compare(1/zbuff(gamma(shape=SHAPE, scale=SCALE)), inversegamma(shape=SHAPE, rate=SCALE))

# Using the Monte Carlo model for distributions simplifies things a little here.
# You can just reciprocate a gamma and you'll get the appropriate inverse gamma.
# You don't even need the zbuff function to protect against division by zero, as 
# the Monte Carlo deviates will never be exactly zero.  Having an inversegamma 
# constructor is handy in Risk Calc and pba.r because their distribution models 
# are discretizations rather than Monte Carlo assemblages.

#def gamma(shape,rate=1,scale=None) :
#    if scale is None : scale = 1/rate
#    rate = 1/scale
#    return(sps.gamma.rvs(a=shape,scale=1/rate,size=many))

def inversegammaRECIP(shape, scale=None, rate=None) : 
    if scale is None : scale = 1/rate
    return(1/gamma(shape=shape,scale=1/scale))

def inversegammaSCIPY(shape, scale=None, rate=None) : 
    if scale is None and not rate is None : scale = 1/rate
    return(sps.invgamma.rvs(a=shape,scale=scale,size=many))

a = inversegammaRECIP(scale=SCALE, shape=SHAPE)
b = inversegammaSCIPY(scale=SCALE, shape=SHAPE)
compare(a,b)

###############################################################################
# c-box and distribution-free p-box constructors

print('****** 11')

def sh(x,t,Data=None) : 
    if Data is None : Data = data
    plt.title(t)
    edf(Data,'y')
    edf(x)
    plt.show()

k = 22
m = 11
n = k + m
fdata = np.concatenate((m*[0],k*[1]))
bdata = np.random.uniform(size=25) > 0.35
idata = np.round(np.random.uniform(size=25) * 16)
data = np.random.uniform(size=25) * 30
x2 = 5 + np.random.uniform(size=25) * 30
error = np.random.normal(size=25)

x=km(k,m);                                    sh(x,'km',fdata)
x=KN(k,n);                                    sh(x,'KN',fdata)
x=FKN(k,n);                                   sh(x,'FKN',fdata)
x=CBbernoulli(bdata);                         sh(x,'CBbernoulli',bdata)
x=CBbernoulli_p(bdata);                       sh(x,'CB p',bdata)             
x=CBbinomial(n,idata);                        sh(x,'CBbinomial(n)',idata)
x=CBbinomial_p(n,idata);                      sh(x,'CB p',idata)
#x=CBbinomialnp(x);                           sh(x,'CB')
#x=CBbinomialnp_n(x);                         sh(x,'CB')
#x=CBbinomialnp_p(x);                         sh(x,'CB') 
x=CBpoisson(idata);                           sh(x,'CBpoisson',idata)
x=CBpoisson_mean(idata);                      sh(x,'CB mean',idata)
x=CBexponential(data);                        sh(x,'CBexponential')          
x=CBexponential_mean(data);                   sh(x,'CB mean')                
x=CBexponential_lambda(data);                 sh(x,'CB lambda')              
x=CBnormal(data);                             sh(x,'CBnormal')
x=CBnormal_mu(data);                          sh(x,'CB mu')
x=CBnormal_sigma(data);                       sh(x,'CB sigma')
x=CBlognormal(data);                          sh(x,'CBlognormal')
x=CBlognormal_mu(data);                       sh(x,'CB mu')
x=CBlognormal_sigma(data);                    sh(x,'CB sigma')
x=CBuniform(data);                            sh(x,'CBuniform')
x=CBuniform_midpoint(data);                   sh(x,'CB midpoint')
x=CBuniform_width(data);                      sh(x,'CB width')
x=CBuniform_minimum(data);                    sh(x,'CB minimum')
x=CBuniform_maximum(data);                    sh(x,'CB maximum')
x=CBnonparametric(data);                      sh(x,'CBnonparametric')
x=CBnormal_meandifference(data, x2);          sh(x,'CB normal mean difference, should be ~5')
#x=CBnonparametric_deconvolution(data, error); sh(x,'CB deconvolution')       ##

###############################################################################
# maxent constructors
"""
print('****** 12')

def sh(x,t) : 
    edf(x)
    plt.title(t)
    plt.show()

x=MEminmax(min=10, max=14);              sh(x,'MEminmax(10,14)')
x=MEminmaxmean(min=10, max=14, mean=11); sh(x,'MEminmaxmean(10,14,11)')      ##
x=MEmeansd(mean=20, sd=1);               sh(x,'MEmeansd(20,1)')
x=MEminmean(min=1,mean=2);               sh(x,'MEminmean(1,2)')
#x=MEmeangeomean(mean=12, geomean=10);    sh(x,'MEmeangeomean(12,10)')       ##
x=MEdiscretemean(x=[1,2,3,4,5,6],mu=2.3);sh(x,'discretemean([1,2,3,4,5,6],2.3)') # e.g., MEdiscretemean(1:10,2.3)
x=MEquantiles(v=np.array((0,1,2,3,5)),p=np.array((0,.03,.3,.36,1))); sh(x,'MEquantiles([0,1,2,3,5],p=[0,.03,.3,.36,1]')
x=MEdiscreteminmax(min=21,max=45);       sh(x,'MEdiscreteminmax(21,45)')
x=MEmeanvar(10,3);                       sh(x,'MEmeanvar(10,3)')
x=MEminmaxmeansd(10,20,13,1);            sh(x,'MEminmaxmeansd(10,20,13,1)')
x=MEmmms(min=10, max=20, mean=13, sd=2); sh(x,'MEmmms(10,20,13,2)')
x=MEminmaxmeanvar(0,1,0.8,0.1);          sh(x,'MEminmaxmeanvar(0,1,0.8,0.1)')

###############################################################################
# maximum likelihood constructors

print('****** 13')

def sh(x,t,Data=None) : 
    if Data is None : Data = data
    plt.title(t)
    edf(Data,'y')
    edf(x)
    plt.show()

data = np.random.uniform(size=25)
datat = 2 * np.random.normal(size=25)
data2 = 2 * data
data10 = np.round(data*10) 
datap = data + 1
data10p = np.round(data*10) + 1
data100 = 1+np.round(data * 100)
dataBB = np.array((0,0,0,0,0,0,0,0,0,0,0,0,3,24,104,286,670,1033,1343,1112,829,478,181,45,7))
N = 12

x=sMLbernoulli(data); sh(x,'sMLbernoulli(data)')
x=sMLnormal(data); sh(x,'sMLnormal(data)')
x=sMLgaussian(data); sh(x,'sMLgaussian(data)')
x=sMLexponential(data); sh(x,'sMLexponential(data)')
x=sMLpoisson(data); sh(x,'sMLpoisson(data)')
x=sMLgeometric(data); sh(x,'sMLgeometric(data)')
x=sMLgumbel(data); sh(x,'sMLgumbel(data)')
x=sMLpascal(data); sh(x,'sMLpascal(data)')
x=sMLuniform(data); sh(x,'sMLuniform(data)')
x=sMLrectangular(data); sh(x,'sMLrectangular(data)')
x=sMLpareto(data); sh(x,'sMLpareto(data)')
x=sMLlaplace(data); sh(x,'sMLlaplace(data)')
x=sMLdoubleexponential(data); sh(x,'sMLdoubleexponential(data)')
x=sMLlognormal2(data); sh(x,'sMLlognormal2(data)')
x=sMLlognormal(data); sh(x,'sMLlognormal(data)')
x=sMLloguniform(data); sh(x,'sMLloguniform(data)')
x=sMLweibull(data); sh(x,'sMLweibull(data)')
x=sMLgamma(data); sh(x,'sMLgamma(data)')

###############################################################################
# alternative maximum likelihood constructors (using sps)

print('****** 14')


x=MLbernoulli(data); sh(x,'MLbernoulli(data)')                           


data = beta(2,3)[0:25]; x=MLbeta(data); sh(x,'MLbeta(data)')
data = betabinomial(2,3)[0:25]; x=MLbetabinomial(data); sh(x,'MLbetabinomial(data)')
x=MLbinomial(data); sh(x,'MLbinomial(data)')
x=MLchisquared(data); sh(x,'MLchisquared(data)')
x=MLexponential(data); sh(x,'MLexponential(data)')
x=MLF(data); sh(x,'MLF(data)')
x=MLgamma(data); sh(x,'MLgamma(data)')
x=MLgammaexponential(data); sh(x,'MLgammaexponential(data)')
x=MLgeometric(data); sh(x,'MLgeometric(data)')
x=MLgumbel(data); sh(x,'MLgumbel(data)')
x=MLlaplace(data); sh(x,'MLlaplace(data)')
x=MLlogistic(data); sh(x,'MLlogistic(data)')
x=MLlognormal(data); sh(x,'MLlognormal(data)')
x=MLloguniform(data); sh(x,'MLloguniform(data)')
x=MLnegativebinomial(data); sh(x,'MLnegativebinomial(data)')
x=MLnormal(data); sh(x,'MLnormal(data)')
x=MLpareto(data); sh(x,'MLpareto(data)')
x=MLpoisson(data); sh(x,'MLpoisson(data)')
x=MLpowerfunction(data); sh(x,'MLpowerfunction(data)')
x=MLrayleigh(data); sh(x,'MLrayleigh(data)')
x=MLstudent(data); sh(x,'MLstudent(data)')
x=MLtriangular(data); sh(x,'MLtriangular(data)')
x=MLuniform(data); sh(x,'MLuniform(data)')

###############################################################################
# method of matching moments (MoM) constructors

print('****** 15')

def sh(x,t,Data=None) : 
    if Data is None : Data = data
    plt.title(t)
    edf(Data,'y')
    edf(x)
    plt.show()

data = np.random.uniform(size=25)
datat = 2 * np.random.normal(size=25)
data2 = 2 * data
data10 = np.round(data*10) 
datap = data + 1
data10p = np.round(data*10) + 1
data100 = 1+np.round(data * 100)
dataBB = np.array((0,0,0,0,0,0,0,0,0,0,0,0,3,24,104,286,670,1033,1343,1112,829,478,181,45,7))
N = 12

x=MMbernoulli(data);                sh(x,'MMbernoulli(data)')
x=MMbeta(data);                     sh(x,'MMbeta(data)')
x=MMbetabinomial(N,dataBB);         sh(x,'MMbetabinomial(12,dataBB)',dataBB)  
x=MMbinomial(data10p);              sh(x,'MMbinomial(data100)',data10) 
x=MMchisquared(data);               sh(x,'MMchisquared(data)')
x=MMexponential(data);              sh(x,'MMexponential(data)')
x=MMF(data2);                       sh(x,'MMF(data2)',data2)               
x=MMgamma(data10);                  sh(x,'MMgamma(data10)',data10)
x=MMgaussian(data);                 sh(x,'MMgaussian(data)')
x=MMgeometric(data10p);             sh(x,'MMgeometric(data10p)',data10p)
x=MMpascal(data10p);                sh(x,'MMpascal(data10p)',data10p)
x=MMgumbel(data);                   sh(x,'MMgumbel(data)')
x=MMextremevalue(data);             sh(x,'MMextremevalue(data)')
x=MMlognormal(data);                sh(x,'MMlognormal(data)')
x=MMlaplace(data);                  sh(x,'MMlaplace(data)')
x=MMdoubleexponential(data);        sh(x,'MMdoubleexponential(data)')
x=MMlogistic(data);                 sh(x,'MMlogistic(data)')
x=MMloguniform(data);               sh(x,'MMloguniform(data)')
x=MMnormal(data);                   sh(x,'MMlognormal(data)')
x=MMpareto(data);                   sh(x,'MMpareto(data)')
x=MMpoisson(data);                  sh(x,'MMpoisson(data)')
x=MMpowerfunction(data);            sh(x,'MMpowerfunction(data)')
x=MMt(datat);                       sh(x,'MMt(datat)',datat)                        
x=MMstudent(datat);                 sh(x,'MMstudent(datat)',datat)                  
x=MMuniform(data);                  sh(x,'MMuniform(data)')
x=MMrectangular(data);              sh(x,'MMrectangular(data)')
x=MMtriangular(data);               sh(x,'MMtriangular(data)')

###############################################################################
# maximum likelihood constructors

print('****** 13')




###############################################################################
# alternative maximum likelihood constructors

print('****** 14')


###############################################################################
# Bayes constructors

print('****** 16')


###############################################################################
# maximum a posteriori constructors

print('****** 17')

















###############################################################################
# bestiary of precise distributions

print('****** 18')

def sh(x,t) : 
    edf(x)
    plt.title(t)
    plt.show()

x=bernoulli(p=0.25);               sh(x,'bernoulli(p=0.25)') 
x=beta(a=2,b=3) ;                  sh(x,'beta(a=2,b=3)') 
x=betabinomial2(size=10,v=2,w=3);  sh(x,'betabinomial2(size=10,v=2,w=3)') 
x=betabinomial(size=10,v=2,w=3);   sh(x,'betabinomial(size=10,v=2,w=3)')  
x=binomial(12,0.4);                sh(x,'binomial(size=12,p=0.4)')
x=chisquared(v=6);                 sh(x,'chisquared(v=6)')
x=exponential(mean=2);             sh(x,'exponential(mean=2)') 
x=F(6,11);                         sh(x,'F(df1=6,df2=11)')
x=gamma(shape=4,rate=2);           sh(x,'gamma(shape=4,rate=2)')
x=gammaexponential(shape=4,rate=2);sh(x,'gammaexponential(shape=4,rate=2)')
x=geometric(m=0.3);                sh(x,'geometric(m=0.3)')
x=gumbel(2,4);                     sh(x,'gumbel(loc=2,scale=4)')
x=inversechisquared(14);           sh(x,'inversechisquared(df=14)') 
x=inversegamma(shape=2,scale=4);   sh(x,'inversegamma(shape=2,scale=4)')
x=laplace(a=4,b=5);                sh(x,'laplace(a=4,b=5)') 
x=logistic(2,3);                   sh(x,'logistic(loc=2,scale=3)')
x=lognormal(m=2,s=1);              sh(x,'lognormal(m=10,s=1)')
x=lognormal2(mlog=-2,slog=1);      sh(x,'lognormal2(mlog=-2,slog=1)')
x=loguniform(min=2, max=6);        sh(x,'loguniform(min=2, max=6)')
x=negativebinomial(size=10,prob=0.25); sh(x,'negativebinomial(size=10,prob=0.25)') 
x=normal(m=5,s=1) ;                sh(x,'normal(m=5,s=1)')
x=pareto(mode=3, c=2);             sh(x,'pareto(mode=3, c=2)')
x=poisson(m=4);                    sh(x,'poisson(m=4)')
#x=rayleigh(4,3);                   sh(x,'rayleigh(4,3)')
x=sawinconrad(2,4,9) ;             sh(x,'student(2,4,9)')
x=student(v=5) ;                   sh(x,'student(v=5)')
x=triangular(2,5,11);              sh(x,'triangular(2,5,11)')
x=uniform(a=2,b=4) ;               sh(x,'uniform(a=2,b=4)') 

print('****** 19')

x = np.array([ 1,  2,  3,  4,  5,  6,  7,  8,  9])
w = np.array([ 1,  5,  1,  1,  1,  1,  2,  3, 15])
q = mixture(x,w)
h = histogram(x)
edf(q)
edf(h)
plt.show()
print(np.mean(x),np.mean(h), 'but', np.mean(q), np.sum(x * (w/np.sum(w))))

print('****** 20')

"""
 
"""


###############################################################################
# Repeat some of the above calculations in R to compare and check the results
###############################################################################

###############################################################################
# power function distribution

rbyc(3,4)
powerfunction(1,1);powerfunction(2,1);powerfunction(3,1);powerfunction(4,1)
powerfunction(1,2);powerfunction(2,2);powerfunction(3,2);powerfunction(4,2)
powerfunction(1,3);powerfunction(2,3);powerfunction(3,3);powerfunction(4,3)

###############################################################################
# gamma-exponential compound distribution

#source('pba BETTER.r')
rbyc()
pl(0,20)
g = gammaexponential(shape=1, scale=1);   blue(g)
g = gammaexponential(shape=2, scale=1);   red(g)
g = gammaexponential(shape=1, scale=2);   black(g)
g = gammaexponential(shape=2, scale=2);   yellow(g)
g = gammaexponential(shape=1, scale=11);  cyan(g)
g = gammaexponential(shape=1, scale=0.1); purple(g)

###############################################################################
# maxent

cat('****** 12\n')

sh = function(x,t) { pl(x); edf(x,1); title(t) }
rbyc(3,4)
x=MEminmax(min=10, max=14);              sh(x,'MEminmax(10,14)')
x=MEminmaxmean(min=10, max=14, mean=11); sh(x,'MEminmaxmean(10,14,11)')      ##
x=MEmeansd(mean=20, sd=1);               sh(x,'MEmeansd(20,1)')
x=MEminmean(min=1,mean=2);               sh(x,'MEminmean(1,2)')
#x=MEmeangeomean(mean=12, geomean=10);    sh(x,'MEmeangeomean(12,10)')       ##
x=MEdiscretemean(x=c(1,2,3,4,5,6),mu=2.3);sh(x,'discretemean(c(1,2,3,4,5,6),2.3)') # e.g., MEdiscretemean(1:10,2.3)
x=MEquantiles(c(0,1,2,3,5),c(0,.03,.3,.36,1)); sh(x,'MEquantiles([0,1,2,3,5],p=[0,.03,.3,.36,1]')
x=MEdiscreteminmax(min=21,max=45);       sh(x,'MEdiscreteminmax(21,45)')
x=MEmeanvar(10,3);                       sh(x,'MEmeanvar(10,3)')
x=MEminmaxmeansd(10,20,13,1);            sh(x,'MEminmaxmeansd(10,20,13,1)')
x=MEmmms(min=10, max=20, mean=13, sd=2); sh(x,'MEmmms(10,20,13,2)')
x=MEminmaxmeanvar(0,1,0.8,0.1);          sh(x,'MEminmaxmeanvar(0,1,0.8,0.1)')

###############################################################################
# method of matching moments

cat('****** 13\n')

sh = function(x,t,Data=NULL) { 
  if (is.null(Data)) Data = data
  pl(x)
  title(t)
  edf(Data,col='green')
  edf(x)
  }

data = runif(25)
datat = 2 * rnorm(25)
data2 = 2 * data
data10 = round(data*10)
datap = data + 1
data10p = round(data*10) + 1
data100 = 1+round(data * 100)
N = int(max(data100))

rbyc(5,5)
x=MMbernoulli(data);                sh(x,'MMbernoulli(data)')
x=MMbeta(data);                     sh(x,'MMbeta(data)')
x=MMbetabinomial(N,data100);        sh(x,'MMbetabinomial(int(max(data)),data100)',data100)  
x=MMbinomial(data10p);              sh(x,'MMbinomial(data100)',data10) 
x=MMchisquared(data);               sh(x,'MMchisquared(data)')
x=MMexponential(data);              sh(x,'MMexponential(data)')
x=MMF(data2);                       sh(x,'MMF(data2)',data2)               
x=MMgamma(data10);                  sh(x,'MMgamma(data10)',data10)
x=MMgaussian(data);                 sh(x,'MMgaussian(data)')
x=MMgeometric(data10p);             sh(x,'MMgeometric(data10p)',data10p)
x=MMpascal(data10p);                sh(x,'MMpascal(data10p)',data10p)
x=MMgumbel(data);                   sh(x,'MMgumbel(data)')
x=MMextremevalue(data);             sh(x,'MMextremevalue(data)')
x=MMlognormal(data);                sh(x,'MMlognormal(data)')
x=MMlaplace(data);                  sh(x,'MMlaplace(data)')
x=MMdoubleexponential(data);        sh(x,'MMdoubleexponential(data)')
x=MMlogistic(data);                 sh(x,'MMlogistic(data)')
x=MMloguniform(data);               sh(x,'MMloguniform(data)')
x=MMnormal(data);                   sh(x,'MMlognormal(data)')
x=MMpareto(data);                   sh(x,'MMpareto(data)')
x=MMpoisson(data);                  sh(x,'MMpoisson(data)')
x=MMpowerfunction(data);            sh(x,'MMpowerfunction(data)')
#x=MMt(datat);                       sh(x,'MMt(datat)',datat)                        
x=MMstudent(datat);                 sh(x,'MMstudent(datat)',datat)                  
x=MMuniform(data);                  sh(x,'MMuniform(data)')
#x=MMrectangular(data);              sh(x,'MMrectangular(data)')
x=MMtriangular(data);               sh(x,'MMtriangular(data)')

###############################################################################
# bestiary

sh = function(x,t) {
    edf(x,new=TRUE)
    title(t)
    }
rbyc(5,5)
x=bernoulli(p=0.25);               sh(x,'bernoulli(p=0.25)') 
x=beta(2,3) ;                      sh(x,'beta(2,3)') 
x=betabinomial(n=10,2,3);          sh(x,'betabinomial(size=10,2,3)')  
x=binomial(12,0.4);                sh(x,'binomial(size=12,0.4)')
x=chisquared(6);                   sh(x,'chisquared(6)')
x=exponential(mean=2);             sh(x,'exponential(mean=2)') 
x=F(6,11);                         sh(x,'F(df1=6,df2=11)')
x=gamma(shape=4,rate=2);           sh(x,'gamma(shape=4,rate=2)')
x=gammaexponential(shape=4,rate=2);sh(x,'gammaexponential(shape=4,rate=2)')
x=geometric(prob=0.3);             sh(x,'geometric(prob=0.3)')
x=gumbel(2,4);                     sh(x,'gumbel(loc=2,scale=4)')
x=laplace(a=4,b=5);                sh(x,'laplace(a=4,b=5)') 
x=logistic(2,3);                   sh(x,'logistic(loc=2,scale=3)')
x=lognormal(m=2,s=1);              sh(x,'lognormal(m=10,s=1)')
x=SL(meanlog=-2,stdlog=1);         sh(x,'SL(meanlog=-2,stdlog=1)')  # lognormal should work here instead of SL
x=loguniform(min=2, max=6);        sh(x,'loguniform(min=2, max=6)')
x=negativebinomial(size=10,prob=0.25); sh(x,'negativebinomial(size=10,prob=0.25)') 
x=normal(m=5,s=1) ;                sh(x,'normal(m=5,s=1)')
x=pareto(mode=3, c=2);             sh(x,'pareto(mode=3, c=2)')
x=poisson(4);                      sh(x,'poisson(4)')
x=powerfunction(4,3);              sh(x,'powerfunction(4,3)')
x=rayleigh(4,3);                   sh(x,'rayleigh(4,3)')
x=sawinconrad(2,4,9) ;             sh(x,'student(2,4,9)')
x=student(df=5) ;                  sh(x,'student(df=5)')
x=triangular(2,5,11);              sh(x,'triangular(2,5,11)')
x=uniform(2,4) ;                   sh(x,'uniform(2,4)') 



# -----------------------------------------------------------------------------
# End of IOANNA6.PY LIBRARY
# -----------------------------------------------------------------------------
'''  









































''' CHAPTER

   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
   TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS     TESTS & CHECKS    
'''

# ----------------------------------------------------------------------
# Some tests, checks, examples and exercises to check the implementation
# ----------------------------------------------------------------------

"""
# If you want to compare results against what RAMAS Risk Calc gives you, be
# sure to set

PbO.steps = 100
PbO.bOt = 0.005
PbO.tOp = 0.995

# and remember that Risk Calc distinguishes + and |+|, as well as < and |<|.
"""   


"""   
x = I(2, 5)
x.left()           # 2
x.right()          # 5

y = I(4, 3)
y.left()           # 3
y.right()          # 4

z = I('7.5')
z.left()           # 7.45
z.right()          # 7.55

v = x + y 
print(v)           # [5, 9]
w = x * y * z
w                  # [44.7, 151.0]

u = U(0, 1)        # uniform(0,1) precise distribution
n = N(0, 1)        # normal(0,1) precise distribution
u                  # Pbox(range=[0.0, 1.0], mean=0.5)
print(n)           # Pbox(range=[-3.09, 3.09], mean=0.0)
beta(2, 3)         # Pbox(range=[0.0, 1.0], mean=0.4)  

N(I(5,6), 1)       # p-box of normal distributions with mean in [5,6]
KN(2, 10)          # c-box about a probability given 2 successes in 10 trials
MMM(0, 10, 1)      # p-box of all distributions with unit mean and range [0,10]

s = u + n + I(0,1) # sum of two distributions and an interval is a p-box
s                  # Pbox(range=[-3.09, 5.09], mean=[0.5, 1.5])
e = env(u, n)      # envelope of two distributions is a p-box
e                  # Pbox(range=[-3.09, 3.09], mean=[0.0, 0.5])
KN(2,10)*KN(2,100) # logical conjunction (AND) of two c-boxes

U(I(1,2), I(3,4))  # Pbox(range=[1.0, 4.0], mean=[2.0, 3.0])
U(I(1,3), I(2,4))  # Pbox(range=[1.0, 4.0], mean=[1.5, 3.5])    
B(0,0)             # Pbox(range=[0.0, 1.0], mean=[0.0, 1.0])
B(I(0.01,1), I(3,4))  # Pbox(range=[0.0, 1.0], mean=[0.0, 0.25])

# Plotting
plot(x)            # interval
plot(x,form='e')   # ellipse, 't' for triangle
plot(s)            # p-box
plot(e,fmt='b:')   # blue with dotted lines
plt.show()

def both(a,b): plot(b, fmt='b'); plot(a) 
    
both(beta(0,0),  beta(0.01,0.01))  # dunno
both(beta(0,10), beta(0.01,10))    # zero
both(beta(10,0), beta(10,0.01))    # one
both(km(0,0),    km(0.01,0.01))    # dunno
both(km(0,10),   km(0.01,10))
both(km(10,0),   km(10,0.01))
both(KN(0,0),    KN(0.01,0.01))    # dunno
both(KN(0,10),   KN(0.01,10))
both(FKN(0,0),   FKN(0.01,0.01))   # dunno
both(FKN(0,10),  FKN(0.01,10))     # zero

env(Interval(1,2),3,5,Interval(6,7),9)
env(Interval(1,2),3,5,Interval(6,7),3)
env(Interval(1,2),N(5,1),6)

Interval(1,2) + Interval(3,4)  # [4, 6]  unlike Risk Calc, only p-boxes are automatically plotted
Interval(1,2) * Interval(3,4)  # [3, 8]
Interval(1,2) - Interval(3,4)  # [-3, -1]
Interval(1,2) / Interval(3,4)  # [0.25, 0.6667]

12 + Interval(1,2)          # [13.0, 14.0]
Interval(1,2) + 12          # [13.0, 14.0]
2 * Interval(1,2)           # [2.0, 4.0]
Interval(1,2) * 2           # [2.0, 4.0]
12 - Interval(1,2)          # [10.0, 11.0]
Interval(1,2) - 12          # [-11.0, -10.0]
Interval(1,2) * -1          # [-2.0, -1.0]
-1 * Interval(1,2)          # [-2.0, -1.0]
-Interval(1,2)              # [-2.0, -1.0] 
2 / Interval(3,4)           # [0.5, 0.667]
Interval(3,4)/2             # [1.5, 2.0]

U(1,2) + U(3,4)             # Pbox(range=[4.0, 6.0], mean=5.0)
U(1,2) * U(3,4)             # Pbox(range=[3.0, 8.0], mean=5.0)
U(1,2) - U(3,4)             # Pbox(range=[-3.0, -1.0], mean=-2.0)
U(1,2) / U(3,4)             # Pbox(range=[0.25, 0.6667], mean=[1.787, 1.787])

2 + N(5,1)      # Pbox(range=[3.9, 10.09], mean=7.0)
N(5,1) + 2      # Pbox(range=[3.9, 10.09], mean=7.0)
2 * N(5,1)      # Pbox(range=[3.8, 16.18], mean=10.0)
N(5,1) * 2      # Pbox(range=[3.8, 16.18], mean=10.0)
2 - N(5,1)      # Pbox(range=[-6.09, 0.09], mean=-3.0)
N(5,1) - 2      # Pbox(range=[-0.09, 6.09], mean=3.0)
1/N(5,1)        # Pbox(range=[0.123, 0.5236], mean=[0.2081, 0.2102])
N(5,1)/10       # Pbox(range=[0.19, 0.809], mean=5.1)
N(5,1) * -1     # Pbox(range=[-8.09, -1.91], mean=-5.0)
-1 * N(5,1)     # Pbox(range=[-8.09, -1.91], mean=-5.0)
-N(5,1)         # Pbox(range=[-8.09, -1.91], mean=-5.0)
2  + -N(5,1)    # Pbox(range=[-6.09, 0.090], mean=-3.0)

N(5,1)+Interval(2,3)    # Pbox(range=[3.9, 11.09], mean=[7.0, 8.0])
Interval(2,3)+N(5,1)    # Pbox(range=[3.9, 11.09], mean=[7.0, 8.0])
N(5,1)*Interval(2,3)    # Pbox(range=[3.8, 24.27], mean=[7.0, 8.0])
Interval(2,3)*N(5,1)    # Pbox(range=[3.8, 24.27], mean=[7.0, 8.0])
N(5,1)-Interval(2,3)    # Pbox(range=[-1.09, 6.09], mean=[2.0, 3.0])
Interval(2,3)-N(5,1)    # Pbox(range=[-6.09, 1.09], mean=[-3.0, -2.0])
N(5,1)/Interval(3,4)    # Pbox(range=[0.477, 2.69], mean=[1.25, 1.67])
Interval(5,6)/N(5,1)    # Pbox(range=[0.618, 3.142], mean=[1.04, 1.26])

A = Interval(3,4)/N(5,1)
B = 3/N(5,1)
C = 4/N(5,1)
plot(A)
plot(A); plot(B,fmt='b')
plot(A); plot(B,fmt='b'); plot(C,fmt='g')





aa = [0.51, 0.49, -0.1, -0.51]
tab = '\t'
print('a',tab,'trunc(a)',tab,'round(a)',tab,'ceil(a)',tab,'floor(a)')
for a in aa: print(a, tab, trunc(a),tab,round(a),tab,tab,ceil(a),tab,tab,floor(a))

#      a    trunc(a)  round(a)   ceil(a)  floor(a)
#    0.51       0.0       1.0       1.0       0.0
#    0.49       0.0       0.0       1.0       0.0
#    -0.1      -0.0      -0.0      -0.0      -1.0
#   -0.51      -0.0      -1.0      -0.0      -1.0


# subplotting

fig, ax = init_splot(4, 5, sharex=True, sharey=True)

splot(N(5,1))                   # draw a p-box
lines([0,10],[0,1])             # draw a line
lines([5,10],[0,1],c='b')       # draw a colored line  

splot(U(0,10))                  # another p-box in a new subplot
title("This is a title")

splot(T(0,10,1))                # yet another p-box in a new subplot
lines(U(0,1),c='r')             # superimpose another colored p-box

splot(U(0,10),c='g')            # yet another p-box in a new subplot, in color

splot(U(0,10),c='xkcd:grey',main='THIS IS A TITLE')            # yet another p-box in a new subplot, in color



def Aabove(a,c):
    A = above(a,c)
    plot(a,c='xkcd:grey',lw=1)
    plot(A,lw=1)
    
a = N(5,1)
plot(a,c='xkcd:grey',lw=1)
Aabove(a,0)
Aabove(a,1)
Aabove(a,2)
Aabove(a,3)
Aabove(a,4)
Aabove(a,5)
Aabove(a,6)
Aabove(a,7)
Aabove(a,8)
#Aabove(a,9)
#Aabove(a,10)

def Abelow(a,c):
    A = below(a,c)
    plot(a,c='xkcd:grey',lw=1)
    plot(A,lw=1)
    
a = N(5,1)
plot(a,c='xkcd:grey',lw=1)
#Abelow(a,0)
#Abelow(a,1)
Abelow(a,2)
Abelow(a,3)
Abelow(a,4)
Abelow(a,5)
Abelow(a,6)
Abelow(a,7)
Abelow(a,8)
Abelow(a,9)
Abelow(a,10)

# panoply of mass reassigning functions

massreassignexamples()




a = N(5,1);     A = N(2,1)
b = I(0.5,3);   B = I(-2,3)
c = 1.2;        C = 0;        CC = -1

log(a)          # Pbox(range=[0.6469, 2.0906], mean=[1.5845, 1.5917])
log(A)          # Pbox(range=[-inf, 1.6273], mean=[-inf, -inf])
log(b)          # [-0.6931, 1.0986]
log(B)          # [-inf, 1.0986]
log(c)          # np.float64(0.1823)
log(C)          # np.float64(-inf)
log(CC)         # np.float64(nan)
log(N(-1,1))    # notice the ordinate axis is not [0,1]


log(math.e)         # 1
log(0.2)            # -1.6094379124341003
log(0)              # -inf
log(-1)             # nan
log([1,2,3])        # array([0.        , 0.69314718, 1.09861229])
log([0,1,2,3])      # array([      -inf, 0.        , 0.69314718, 1.09861229])
log([0,1,-2,3])     # array([      -inf, 0.        ,        nan, 1.09861229])
log(Interval(1,3))  # [0.0, 1.0986122886681098]
log(Interval(0,3))  # [-inf, 1.0986122886681098]
log(Interval(-1,3)) # [nan, 1.0986122886681098]
#log(Interval(-2,-1))# ValueError: Could not take log
log(N(5,1))         # Pbox(range=[0.6469, 2.0906], mean=[1.5845, 1.5917])
log(N(3,1))         # Pbox(range=[-inf, 1.8066], mean=[-inf, 1.0376])
log(N(2,1))         # Pbox(range=[-inf, 1.6273], mean=[-inf, 0.5776])
#log(-N(5,1))        # ValueError: Could not take log


sqrt(math.e)         # 1.6487
sqrt(0.2)            # 0.4472
sqrt(0)              # 0
sqrt(-1)             # nan
sqrt([1,2,3])        # array([1.        , 1.41421356, 1.73205081])
sqrt([0,1,2,3])      # array([0.        , 1.        , 1.41421356, 1.73205081])
sqrt([0,1,-2,3])     # array([0.        , 1.        ,        nan, 1.73205081])
sqrt(Interval(1,3))  # [1.0, 1.7320]
sqrt(Interval(0,3))  # [0.0, 1.732]
sqrt(Interval(-1,3)) # [0.0, 1.7320]
#sqrt(Interval(-2,-1))# ValueError: sqrt undefined for negative values
sqrt(N(5,1))         # Pbox(range=[1.3819, 2.8443], mean=[2.2207, 2.2281])
sqrt(N(3,1))         # Pbox(range=[0.0, 2.4678], mean=[1.6976, 1.7100])
sqrt(N(2,1))         # Pbox(range=[0.0, 2.2561], mean=[1.3721, 1.3948])
#sqrt(-N(5,1))        # ValueError: sqrt undefined for negative values



# cos ... make code to check answers in Risk Calc
import random
for j in range(20):
    x = round(Interval(10*random.random(), 10*random.random()),3)
    r = cos(x)
    R = str(r)
    print('\ncos(',x,' radians) ===',R)
    
    

# tangent function
import random
for j in range(20):
    x = round(Interval(random.random(), random.random())*2*np.pi-np.pi,3)
    r = tan(x)
    R = str(r)
    print('\ntan(',x,' radians) ===',R)
    a = seq(-np.pi, np.pi)
    plt.plot([left(x),-3,-3,right(x)],[left(r),left(r),right(r),right(r)],c='xkcd:grey')
    plt.plot(a, np.minimum(3,np.maximum(-3,np.tan(a))), 'b')
    plot(x)
    plt.title('tan('+str(x)+')='+R)
    plt.show()
for j in range(20):
    x = round(Interval(10*(random.random()-0.2), 10*(random.random()-0.2)),3)
    r = tan(x)
    R = str(r)
    print('\ntan(',x,' radians) ===',R)

 

#------------------------------------------------------------------------------
# test KS_band(x) against histogram(x)

# this routine is only for testing histogram
def KS_band(values, mn=None, mx=None, conf=0.95, two_sided=True):
    x, y = as_vectors(values)
    allscalars = np.all(x == y)
    if mn is None: mn = np.min(x) if y is None else min(np.min(x), np.min(y))
    if mx is None: mx = np.max(x) if y is None else max(np.max(x), np.max(y))
    x = np.sort(np.asarray(x))
    n = len(x)
    if n == 0: raise ValueError("KS_band: empty sample")
    d = KS_critical(n, conf=conf, two_sided=two_sided)
    F = (np.arange(1, n+1) / n) # empirical CDF values at order statistics F_n(x_k) = k/n for k=1..n
    pL,pU = PbO.ii(), PbO.jjj()
    QL, QU = np.empty_like(pL, dtype=float), np.empty_like(pU, dtype=float)
    # left quasi-inverse, for each p, we want the smallest x such that F_n(x) >= p - d
    for i, p in enumerate(pL):
        target = max(p - d, 0.0)
        if target <= 0.0: QL[i] = x[0]
        else:
            k = np.searchsorted(F, target, side='left')
            if k >= n: QL[i] = x[-1]
            else: QL[i] = x[k]
    # right quasi-inverse, for each p, we want the smallest x such that F_n(x) >= p + d
    for i, p in enumerate(pU):
        target = min(p + d, 1.0)
        if target <= 0.0: QU[i] = x[0]
        else:
            k = np.searchsorted(F, target, side='left')
            if k >= n: QU[i] = x[-1]
            else: QU[i] = x[k] 
    #print('np.all(x == y) is',np.all(x == y),'\nx=...',x[-5:],'\ny=...',y[-5:])
    if allscalars: return fatten(Pbox(QL, QU),0,mn,mx)
    else: return env(Pbox(QL, QU), KS_band(y,mn,mx,conf,two_sided))
   
def h(d):
    a = KS_band(d)
    b = histogram(d) 
    red(a)
    blue(b)

import random

def randi(m=0,M=1,prop=False):
    e = [random.random() * (M-m) + m, random.random() * (M-m) + m]
    if (prop): e.sort()
    return(Interval(e[0],e[1])) 
    
def randinti(m=0,M=100,prop=False):
    e = [random.random() * (M-m) + m, random.random() * (M-m) + m]
    e = [math.trunc(_) for _ in e]
    if (prop): e.sort()
    return(Interval(e[0],e[1])) 
    
h([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])
h([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
h([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
h([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21])
h([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
m = 50; h(np.array(range(m)))
h(np.array(range(2,m+2)))
d1 = np.array(range(m));    d2 = np.array(range(2,m+2))
ha1 = KS_band(d1);          ha2 = KS_band(d2);          
hb1 = histogram(d1);        hb2 = histogram(d2) 
red(env(ha1,ha2));          blue(env(hb1,hb2))

D = as_intervals(np.array(range(m)), np.array(range(2,m+2)))
Ha = KS_band(D)
Hb = histogram(D) 
red(Ha); blue(Hb)

m = 50; h(as_intervals(np.array(range(m)), np.array(range(2,m+2))))

h([randinti(i, i+2) for i in range(20)])
h([randinti() for i in range(40)])
h([randinti(i, i+2) for i in range(200)])

d = [2.64, 3.157, 4.566, 6.879, 7.346, 4.767, 5.395, 7.726, 5.855, 13.663, 0.381, 2.431, 03.171, 10.139, 0.715, 10.563, 11.058, 2.929, 12.151, 3, 4, Interval(2,5), (1,7), np.array([6,9])]; h(d)
d = [2.64, 3.157, 4.566, 6.879, 7.346, 4.767, 5.395, 7.726, 5.855, 13.663, 0.381, 2.431, 03.171, 10.139, 0.715, 10.563, 11.058, 2.929, 12.151, 3, 4, Interval(2,5), (1,7), np.array([6,9])]
x,y = as_vectors(d)
x = c( 2.64 ,  3.157,  4.566,  6.879,  7.346,  4.767,  5.395,  7.726, 5.855, 13.663,  0.381,  2.431,  3.171, 10.139,  0.715, 10.563, 11.058,  2.929, 12.151,  3.   ,  4.   ,  2.   ,  1.   ,  6.   )
y = c( 2.64 ,  3.157,  4.566,  6.879,  7.346,  4.767,  5.395,  7.726, 5.855, 13.663,  0.381,  2.431,  3.171, 10.139,  0.715, 10.563, 11.058,  2.929, 12.151,  3.   ,  4.   ,  5.   ,  7.   ,  9. )




#------------------------------------------------------------------------------
# test significant digit display and interpretation


def test_sigorder():
    assert sigorder('32125') == 1       
    assert sigorder('32125.0') == -1
    assert sigorder('32120.') == 1
    assert sigorder('321250') == 2      
    assert sigorder('3212000') == 4     
    assert sigorder('32120.001') == -3     #
    assert sigorder('32120.000') == -3
    assert sigorder('32121.010') == -3     #
    assert sigorder('32120100') == 3    
    assert sigorder('.32120100') == -8  
    assert sigorder('200') == 3            #
    assert sigorder('210') == 2         
    assert sigorder('201') == 1         
    assert sigorder('200.') == 1
    assert sigorder('200.0') == -1
    assert sigorder('200.00') == -2

def test_lastsigfig():
    assert lastsigfig('32125') == '5'
    assert lastsigfig('32125.0') == '0'
    assert lastsigfig('321250') == '5'
    assert lastsigfig('3212000') == '2'
    assert lastsigfig('32120.001') == '1'
    assert lastsigfig('32120.000') == '0'
    assert lastsigfig('32121.010') == '0'
    assert lastsigfig('32120100') == '1'

def test_about():
    lo, hi = about('200')
    assert abs(lo - 173.6951) < 1e-4
    assert abs(hi - 226.3049) < 1e-4

    lo, hi = about(v=200, r=3, f=0)
    assert abs(lo - 173.6951) < 1e-4
    assert abs(hi - 226.3049) < 1e-4

    lo, hi = about('200.')
    assert abs(lo - 193.82) < 0.01
    assert abs(hi - 206.18) < 0.01

    lo, hi = about(v=200, r=0, f=0)
    assert abs(lo - 197.0046) < 1e-4
    assert abs(hi - 202.9954) < 1e-4

    lo, hi = about(v=200, r=-1, f=0)
    assert abs(lo - 198.5481) < 1e-4
    assert abs(hi - 201.4519) < 1e-4

test_sigorder()
test_lastsigfig()
test_about()

def test_format_sigdigs():
    cases = [
        # (input, D, expected_output)
        (1000.3, 1, '1.e+3'),
        (1000.3, 2, '1.0e+3'),
        (1000.3, 3, '1.00e+3'),
        (1000.3, 4, '1000.'),
        (1000.3, 5, '1000.3'),
        (1000.3, 6, '1000.30'),
        (1000.003, 6, '1000.00'),
        (1000.007, 6, '1000.01'),
        (1234.5678, 9, '1234.56780'),
        (1234.5678, 8, '1234.5678'),
        (1234.5678, 7, '1234.568'),
        (1234.5678, 6, '1234.57'),
        (1234.5678, 5, '1234.6'),
        (1234.5678, 4, '1235'),
        (1234.5678, 3, '1230'),
        (1234.5678, 2, '1200'),
        (1234.5678, 1, '1000'),
        (6666.6666, 1, '7000'),
        (6666.6666, 2, '6700'),
        (6666.6666, 3, '6670'),
        (6666.6666, 4, '6667'),
        (6666.6666, 5, '6666.7'),
        (6666.6666, 6, '6666.67'),
        (6666.6666, 7, '6666.667'),
        (6666.6666, 8, '6666.6666'),
        (6666.6666, 9, '6666.66660'),
    ]

    for x, D, expected in cases:
        result = format_sigdigs(x, D)
        assert result == expected, f"format_sigdigs({x}, {D}) returned {result}, expected {expected}"





#fermi(normal,0.2, 2.3)     # N(1.25, 0.8193193533759979)
#fermi(lognormal,0.2, 2.3)    # lognormal3(0.6782329983125268, 2.5931847639688455) 
#fermi(normal,data=qnormal([0.839, 0.388, 0.899, 0.838, 0.547, 0.748, 0.416, 0.051, 0.035,
#       0.329, 0.792, 0.325, 0.063, 0.165, 0.102, 0.117, 0.612, 0.618,
#       0.407, 0.387, 0.297, 0.489, 0.583, 0.974, 0.838, 0.858, 0.848,
#       0.174, 0.343, 0.43 ],4)) # ~N(4,1)



"""









""" CHAPTER

  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   
  RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   RESEARCH   

#------------------------------------------------------------------------------
# Compare performances of recursive and nonrecursive qpbox() implementations
#------------------------------------------------------------------------------

# To repeat this simulation, uncomment the old qpbox() and RENAME the new, 
# nonrecursive qpbox() to qpboxN().

import time

def benchmark(label, func, repeats=100):
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = func()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    avg = sum(times) / repeats
    print(f" {avg:.6f} seconds: {label}")
    return result, times

ii = PbO.ii();  jj = PbO.jj()
def test_uniform_():  return  qpbox(qunif, ii, jj, I(0.5,2.5), Interval(0,1), Interval(2,3))
def test_uniform_N(): return qpboxN(qunif, ii, jj, I(0.5,2.5), Interval(0,1), Interval(2,3))
u, t = benchmark("Recursive uniform", test_uniform_)
uN,tN= benchmark("Nonrecursive uniform", test_uniform_N)
print('')

def test_uniform_():  return  qpbox(qunif, ii, jj, I(0.5,2.5), Interval(0,3), Interval(2,4))
def test_uniform_N(): return qpboxN(qunif, ii, jj, I(0.5,2.5), Interval(0,3), Interval(2,4))
u, t = benchmark("Recursive uniform overlap", test_uniform_)
uN,tN= benchmark("Nonrecursive uniform overlap", test_uniform_N)
print('')

iii = PbO.iii();  jjj = PbO.jjj()
def test_normal_():  return  qpbox(qnorm, iii, jjj, Interval(4,6), Interval(4,6), Interval(0.8,1.2))
def test_normal_N(): return qpboxN(qnorm, iii, jjj, Interval(4,6), Interval(4,6), Interval(0.8,1.2))
n, t = benchmark("Recursive N([4,6],[0.8,1.2])", test_normal_)
nN,tN= benchmark("Nonrecursive N([4,6],[0.8,1.2])", test_normal_N) 
print('')
    
def test_beta_():  return  qpbox(qbeta, ii, jj, I(0,1), Interval(1,3), Interval(2,4))
def test_beta_N(): return qpboxN(qbeta, ii, jj, I(0,1), Interval(1,3), Interval(2,4))
u, t = benchmark("Recursive beta", test_uniform_)
uN,tN= benchmark("Nonrecursive beta", test_uniform_N)
print('')

lo =I(0,1);  hi = I(8,9);  mu = I(2,3)
def test_mmm_():  return  qpbox([qMMML,qMMMR], ii, jj, mu, lo, hi, mu)
def test_mmm_N(): return qpboxN([qMMML,qMMMR], ii, jj, mu, lo, hi, mu)
u, t = benchmark("Recursive minmaxmean", test_mmm_)
uN,tN= benchmark("Nonrecursive minmaxmean", test_mmm_N)
print('')

# output from the comparisons show the new, nonrecursive algorithm to be better
 # ---------------------------------------------
 # 0.006313 seconds: Recursive uniform
 # 0.003042 seconds: Nonrecursive uniform

 # 0.003779 seconds: Recursive uniform overlap
 # 0.003526 seconds: Nonrecursive uniform overlap

 # 0.004656 seconds: Recursive N([4,6],[0.8,1.2])
 # 0.003420 seconds: Nonrecursive N([4,6],[0.8,1.2])

 # 0.003944 seconds: Recursive beta
 # 0.002932 seconds: Nonrecursive beta

 # 0.001708 seconds: Recursive minmaxmean
 # 0.000177 seconds: Nonrecursive minmaxmean
 # --------------------------------------------
 # 0.005087 seconds: Recursive uniform
 # 0.002878 seconds: Nonrecursive uniform

 # 0.003703 seconds: Recursive uniform overlap
 # 0.002748 seconds: Nonrecursive uniform overlap

 # 0.003961 seconds: Recursive N([4,6],[0.8,1.2])
 # 0.003152 seconds: Nonrecursive N([4,6],[0.8,1.2])

 # 0.004228 seconds: Recursive beta
 # 0.002921 seconds: Nonrecursive beta

 # 0.001667 seconds: Recursive minmaxmean
 # 0.000164 seconds: Nonrecursive minmaxmean
 # ---------------------------------------------
 

#------------------------------------------------------------------------------
# KS confidence bands
#------------------------------------------------------------------------------

# The R version of this library still uses the Miller table to for the KS bands

def KSDmax(n, conf=0.95):
    def MillerD(n, alpha, A):
        return math.sqrt(math.log(1.0 / alpha) / (2*n)) - 0.16693/n - A*(n**-1.5)
    ks80 = [1.0, 0.90000, 0.68377, 0.56481, 0.49265, 0.44698, 0.41037, 0.38148,
            0.35831, 0.33910, 0.32260, 0.30829, 0.29577, 0.28470, 0.27481,
            0.26588, 0.25778, 0.25039, 0.24360, 0.23735, 0.23156]
    ks90 = [1.0, 0.95000, 0.77639, 0.63604, 0.56522, 0.50945, 0.46799, 0.43607,
            0.40962, 0.38746, 0.36866, 0.35242, 0.33815, 0.32549, 0.31417,
            0.30397, 0.29472, 0.28627, 0.27851, 0.27136, 0.26473]
    ks95 = [1.0, 0.97500, 0.84189, 0.70760, 0.62394, 0.56328, 0.51926, 0.48342,
            0.45427, 0.43001, 0.40925, 0.39122, 0.37543, 0.36143, 0.34890,
            0.33760, 0.32733, 0.31796, 0.30936, 0.30143, 0.29408]
    ks98 = [1.0, 0.99000, 0.90000, 0.78456, 0.68887, 0.62718, 0.57741, 0.53844,
            0.50654, 0.47960, 0.45662, 0.43670, 0.41918, 0.40362, 0.38970,
            0.37713, 0.36571, 0.35528, 0.34569, 0.33685, 0.32866]
    ks99 = [1.0, 0.99500, 0.92929, 0.82900, 0.73424, 0.66853, 0.61661, 0.57581,
            0.54179, 0.51332, 0.48893, 0.46770, 0.44905, 0.43247, 0.41762,
            0.40420, 0.39201, 0.38086, 0.37062, 0.36117, 0.35241]
    def table_lookup(tbl, n):
        return tbl[n] if 0 <= n < len(tbl) else None
    if conf == 0: return 0.0
    if conf == 0.80:
        if n > 20: return MillerD(n, 0.10, 0.00256)
        return table_lookup(ks80, n)
    if conf == 0.90:
        if n > 20: return MillerD(n, 0.05, 0.0526)
        return table_lookup(ks90, n)
    if conf == 0.98:
        if n > 20: return MillerD(n, 0.01, 0.20562)
        return table_lookup(ks98, n)
    if conf == 0.99:
        if n > 20: return MillerD(n, 0.005, 0.28464)
        return table_lookup(ks99, n)
    if n > 20: return MillerD(n, 0.025, 0.11282)
    return table_lookup(ks95, n)

import math
import random

# Confidence levels supported by KSDmax (Miller tables + asymptotics)
SUPPORTED_CONFS = [0.80, 0.90, 0.95, 0.98, 0.99]
# n values to test: small, medium, large, plus some randoms
N_SMALL  = list(range(1, 21))          # 1..20 (table region)
N_MID    = [25, 30, 40, 50, 75, 100]   # mid-range
N_LARGE  = [150, 200, 300, 500, 1000]  # large n
N_RANDOM = [random.randint(21, 1000) for _ in range(50)]
N_VALUES = sorted(set(N_SMALL + N_MID + N_LARGE + N_RANDOM))
DECIMAL_PLACES = 4  # Tolerance: number of decimal places that must match
TOL = 10 ** (-DECIMAL_PLACES)

def compare_ks_values(n, conf, verbose_mismatch=True):
    v1 = KS_critical(n, conf=conf, two_sided=True)
    v2 = KSDmax(n, conf=conf)
    diff = abs(v1 - v2)
    match = diff <= TOL
    if (not match) and verbose_mismatch:
        print(f"Mismatch at n={n}, conf={conf:.2f}: "
              f"KS_critical={v1:.10f}, KSDmax={v2:.10f}, diff={diff:.10g}")
    return match, diff, v1, v2

def run_ks_consistency_tests():
    '''Compare KS algorithms for many n and conf values.'''
    total_tests = 0
    mismatches = 0
    max_diff = 0.0
    worst_case = None
    print(f"Testing KS_critical vs KSDmax with tolerance {TOL} "
          f"(~{DECIMAL_PLACES} decimal places)...")
    for conf in SUPPORTED_CONFS:
        for n in N_VALUES:
            total_tests += 1
            match, diff, v1, v2 = compare_ks_values(n, conf, verbose_mismatch=True)
            if not match:
                mismatches += 1
                if diff > max_diff:
                    max_diff = diff
                    worst_case = (n, conf, v1, v2, diff)
    print("\n------------------------------------------------------------")
    print(f"Total tests run: {total_tests}")
    print(f"Mismatches:      {mismatches}")
    if worst_case is not None:
        n, conf, v1, v2, diff = worst_case
        print(f"Worst mismatch at n={n}, conf={conf:.2f}:")
        print(f"  KS_critical = {v1:.10f}")
        print(f"  KSDmax      = {v2:.10f}")
        print(f"  abs diff    = {diff:.10g}")
    else:
        print(f"All values matched within {DECIMAL_PLACES} decimal places.")
    print("------------------------------------------------------------")

# The function run_ks_consistency_tests() compares the approximate Miller 
# against the exact MTW algorithm and finds remarkably good agreement:
# ------------------------------------------------------------
# Total tests run: 390
# Mismatches:      165
# Worst mismatch at n=20, conf=0.80:
#   KS_critical = 0.2315186231
#   KSDmax      = 0.2315600000
#   abs diff    = 4.137685868e-05
# ------------------------------------------------------------



# -----------------------------------------------------------------------------
# End of RESEARCH
# -----------------------------------------------------------------------------
"""  








""" CHAPTER

 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE
 WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE      WOODPILE

#------------------------------------------------------------------------------
# The stuff in the woodpile is obsolete or non-working code, kept for reference
#------------------------------------------------------------------------------

# the beta1() distribution (parameterized by its mean and sd)

# I could just make a qbeta1() function and use it in a call to qpbox()

# v and w are monotone over s
# v and w are simple, parabolic-looking functions of m

for m in np.linspace(0.001,0.999,100):
    s = np.linspace(0, abs(1-2)/2)
    tmp = m*(1-m)/(s*s) - 1   # tmp has repeated variables
    v = m * tmp
    w = (1-m) * tmp
    plt.plot(s,w)
    d = np.diff(w)
    if not all(d<0): print('nonmonotone',d)

for m in np.linspace(0.001,0.999,100):
    s = np.linspace(0, abs(1-2)/2)
    tmp = m*(1-m)/(s*s) - 1   # tmp has repeated variables
    v = m * tmp
    w = (1-m) * tmp
    plt.plot(s,v)
    d = np.diff(v)
    if not all(d<0): print('nonmonotone',d)

for s in np.linspace(0, abs(1-2)/2):
    m = np.linspace(0.001,0.999,100)
    tmp = m*(1-m)/(s*s) - 1   # tmp has repeated variables
    v = m * tmp
    w = (1-m) * tmp
    plt.plot(m,w)
    d = np.diff(w)
    if not all(d<0): print('nonmonotone',d)

for s in np.linspace(0, abs(1-2)/2):
    m = np.linspace(0.001,0.999,100)
    tmp = m*(1-m)/(s*s) - 1   # tmp has repeated variables
    v = m * tmp
    w = (1-m) * tmp
    plt.plot(m,v)
    d = np.diff(v)
    if not all(d<0): print('nonmonotone',d)



# struggles with domain-aware unary functions

# def exp(x):
#     if isinstance(x,Pbox): return Pbox(np.exp(x.u),np.exp(x.d))
#     if isinstance(x,Interval): return Interval(np.exp(left(x)),np.exp(right(x)))
#     return np.exp(x)    

# def handle_domain(lo, hi, bad_lo, bad_hi, f, name, why, error_flag):
#     if (np.any(bad_lo) or np.any(bad_hi)):
#         if error_flag: raise ValueError(f"{name} {why}")
#         flo = np.where(bad_lo, -np.inf, f(lo))
#         fhi = np.where(bad_hi, -np.inf, f(hi))
#         return flo, fhi
#     return f(lo), f(hi)

# def log(x, base=math.e):   # ln=log(x),  log10(x)=log((x,10),  log2(x)=log(x,2)
#     ctx = np.errstate(divide='ignore', invalid='ignore') if IvO.suppress_np_warnings else np.errstate()
#     f = lambda z: np.log(z)/np.log(base)
#     with ctx:
#         if isinstance(x, Pbox):
#             return Pbox(handle_domain(x.u, x.d, x.u<=0, x.d<=0, f, "log", IvO.why_nonpositive, IvO.error_log))
#         if isinstance(x, Interval):
#             return Interval(handle_domain(x.lo, x.hi, x.lo<=0, x.hi<=0, f, "log", IvO.why_nonpositive, IvO.error_log))
#         r = f(x) # scalars, arrays, etc.
#         return float(r) if np.isscalar(r) else r

# def domain_transform(x, f, in_domain, name, why):
#     # Scalars ---------------------------------------------------------------
#     if isscalar(x):
#         if not in_domain(x):
#             if IvO.error_log: raise ValueError(f"{name} {why}")
#             return float('nan')
#         return f(x)
#     # Intervals -------------------------------------------------------------
#     if isinstance(x, Interval):
#         lo, hi = left(x), right(x)
#         lo_ok, hi_ok = in_domain(lo), in_domain(hi)
#         if lo_ok and hi_ok: return Interval(f(lo), f(hi))
#         if IvO.error_log: raise ValueError(f"{name} {why}")
#         # interval repair: drop invalid parts
#         if not lo_ok and not hi_ok: return Interval(float('nan'))
#         if not lo_ok: return Interval(f(hi), float('inf') if f(hi)==float('inf') else f(hi))
#         if not hi_ok: return Interval(f(lo), float('inf') if f(lo)==float('inf') else f(lo))
#         #####
#         if not lo_ok: return Interval(f(hi), f(hi))
#         if not hi_ok: return Interval(f(lo), f(lo))
#         #####
#     # P-boxes ---------------------------------------------------------------
#     if isinstance(x, Pbox):
#         # truncate to domain
#         mask_u = in_domain(x.u)
#         mask_d = in_domain(x.d)
#         if not (mask_u.any() and mask_d.any()):
#             if IvO.error_log: raise ValueError(f"{name} {why}")
#             # no valid mass at all
#             return Pbox([float('nan')],[float('nan')])
#         # renormalize
#         u = x.u[mask_u]
#         d = x.d[mask_d]
#         u = (u - u.min())/(u.max()-u.min()) if len(u)>1 else u
#         d = (d - d.min())/(d.max()-d.min()) if len(d)>1 else d
#         # apply transform
#         return Pbox(f(u), f(d))
#     # Arrays ---------------------------------------------------------------
#     arr = np.asarray(x)
#     mask = in_domain(arr)
#     if not mask.all():
#         if IvO.error_log: raise ValueError(f"{name} {why}")
#         arr = arr.copy()
#         arr[~mask] = np.nan
#     return f(arr)

# def log(x, base=math.e):
#     f = lambda z: np.log(z)/np.log(base)
#     return domain_transform(x, f, lambda z: z>0, "log", IvO.why_nonpositive)

def slambertW(x):  # argument must be larger than -1/e
    prec = 1e-12
    w = rep(NA,len(x))
    w = ifelse(500 < x, log(x - 4.0) - (1.0 - 1.0/log(x)) * log(log(x)), w)
    lx1 = np.log(x + 1.0)
    w = ifelse((-1/exp(1) < x) & (x <= 500), 0.665 * (1 + 0.0195 * lx1) * lx1 + 0.04, w)
    i = 1
    diff = np.full(1,len(x))
    while (i < 100):
        i = i + 1
        wew = w * exp(w)
        wpew = (w+1) * exp(w)
        diff = abs((x-wew)/wpew)
        w = w-(wew-x)/(wpew-(w+2)*(wew-x)/(2*w+2))
    return w
    
def lambertW(x, k=0, tol=1e-8, *args, **kwargs):
    from scipy.special import lambertw
    f = lambda z: lambertw(z, k=0, tol=1e-8)
    return domaintonp(x, Interval(-1/np.e,Inf), 'undefined for such negative values', slambertW, *args, **kwargs)

"""

"""
# Copilot suggested the 'singledispatch' constructions, but we later agreed 
# that this was not actually better than straightforward 'def' functions.
 
from functools import singledispatch
from numbers import Number

# # foo           # Copilot says, yes, we need ALL FIVE of these registrations!
# @singledispatch
# def foo(x): return np.foo(x)
# @exp.register
# def _(x: Number): return np.foo(x)
# @exp.register
# def _(x: np.generic): return np.foo(x) 
# @exp.register
# def _(x: Interval): return Interval(np.foo(x.lo), np.foo(x.hi))
# @exp.register
# def _(x: Pbox): return Pbox(np.foo(x.leftside()), np.foo(x.rightside()))


# left
@singledispatch
def left(x): return min(x)# default: works for lists, tuples, arrays, iterables
@left.register
def _(x: Number): return x                                      # Python scalar
@left.register
def _(x: np.generic): return float(x) # numpy scalar np.float64, np.int64, etc.
@left.register
def _(x: np.ndarray):                                             # numpy array
    if x.size == 1: return float(x)      # scalar-like
    if x.size == 2: return float(x[0])   # interval-like
    return min(x)                        # iterable
@left.register
def _(x: Interval): return x.left()
@left.register
def _(x: Pbox): return x.left()

# right
@singledispatch
def right(x): return max(x)
@right.register
def _(x: Number): return x
@right.register
def _(x: np.generic): return float(x) 
@right.register
def _(x: np.ndarray):
    if x.size == 1: return float(x)      # scalar-like
    if x.size == 2: return float(x[1])   # interval-like
    return max(x)                        # iterable
@right.register
def _(x: Interval): return x.right()
@right.register
def _(x: Pbox): return x.right()

# leftside
@singledispatch
def leftside(x): return min(x)
@leftside.register
def _(x: Number): return x
@leftside.register
def _(x: Interval): return x.left()
@leftside.register
def _(x: Pbox): return x.leftside()
@leftside.register
def _(x: list): return [left(elem) for elem in x]   

# rightside
@singledispatch
def rightside(x): return max(x)
@rightside.register
def _(x: Number): return x
@rightside.register
def _(x: Interval): return x.right()
@rightside.register
def _(x: Pbox): return x.rightside()
@rightside.register
def _(x: list): return [right(elem) for elem in x]   

# steps
@singledispatch
def steps(x): raise TypeError(f"steps() not defined for type {type(x)}")
@steps.register
def _(x: Number): return 1
@steps.register
def _(x: Interval): return 1
@steps.register
def _(x: Pbox): return x.n

# mean
@singledispatch
def mean(x): raise TypeError(f"mean() not defined for type {type(x)}")
@mean.register
def _(x: Number): return x
@mean.register
def _(x: Interval): return x
@mean.register
def _(x: Pbox): return x.mean()

# sd (standard deviation)
@singledispatch
def sd(x,pop=True): raise TypeError(f"sd() not defined for type {type(x)}")
@sd.register
def _(x: Number,pop=True): return 0.0
@sd.register
def _(x: Interval,pop=True): return Interval(0.0, abs(x.hi - x.lo)/2)
@sd.register
def _(x: Pbox,pop=True):
    # Plug in your extremal-variance algorithm here.
    # Example structure:
    # sl, sh = sd_bounds_from_pbox(x.u, x.d)
    return (None, None)
    
# breadth
@singledispatch
def breadth(x): raise TypeError(f"breadth() not defined for type {type(x)}")
@breadth.register
def _(x: Number): return 0
@breadth.register
def _(x: Interval): return right(x) - left(x)
@breadth.register
def _(x: Pbox): return np.sum(x.d - x.u) / x.n

"""

"""
Before the change-at-a-distance crisis:    

    def __init__(self, u, d=None, ml=None, mh=None):
        if is_interval(u) and d is None: u,d = left(u), right(u)
        if isscalar(u) or isscalar(d): 
            many = max(long(u), long(d))
            #many = PbO.steps
            if isscalar(u): u = [u] * many
            if isscalar(d): d = [d] * many
        uu = list(u)
        dd = list(d if d is not None else u)
        if len(uu) != len(dd): raise ValueError("Left and right sides must have same length")
        for i in range(len(uu)):
            if is_missing(uu[i]): # NA or None
                if i==0 : uu[i] = float("-inf")
                else:     uu[i] = uu[i-1]           # outward rounding downward
            #print(type(uu[i]), uu[i])
        for i in reversed(range(len(dd))):
            if is_missing(dd[i]): # NA or None
                if i==len(dd)-1: dd[i] = float("inf")             
                else:            dd[i] = dd[i+1]      # outward rounding upward
        u = np.asarray(uu, dtype=float)
        d = np.asarray(dd, dtype=float)        
        if not is_monotone(u): raise ValueError("Left side nonmonotonic")
        if not is_monotone(d): raise ValueError("Right side nonmonotonic")
        self.u = u
        self.d = d
        self.n = len(u)
        self.ml = float(np.mean(u) if ml is None else ml)
        self.mh = float(np.mean(d) if mh is None else mh)

    def copy(self):
        return Pbox(self.u.copy(), self.d.copy(), self.ml, self.mh)
    
    def left(self):
        return float(self.u[0])

    def right(self):
        return float(self.d[-1])

    def leftside(self):
        return self.u

    def rightside(self):
        return self.d

    def steps(self):
        return self.n

    def mean(self):    # use ends(mean()) if you need to iterate both endpoints
        return Interval(self.ml, self.mh)
            
    def cut(self, p, tight=True):
        if p < 0 or p > 1: raise ValueError("Second argument for cut must be a probability between zero and one")
        n = self.n  
        if tight:   
            p_long = p * n
            fractional = (p_long % 1) == 0
            idx_u = min(n, (1 if fractional else 0) + math.ceil(p_long))
            idx_d = max(1, math.ceil(p_long))
            return Interval(self.u[idx_u-1], self.d[idx_d-1])
        if p == 1: lower = n
        else:
            if (p % (1/n)) == 0: lower = round(p * n)
            else: lower = math.ceil(p * n)
        if p == 0: upper = 1
        else:
            if (p % (1/n)) == 0: upper = round(p * n) + 1
            else: upper = math.floor(p * n) + 1
        return Interval(self.u[max(lower, 1)-1], self.d[min(upper, n)-1])
        
    def median(self):  # conservative w.r.t. discretization  
        #Returns median conservative with respect to discretization.
        return Interval(self.u[self.n // 2 - (1-self.n % 2)], self.d[self.n // 2])
    # # the conservative median=cut(0.5,False) rather than the optimistic one
    # setoption(steps=10)  # even number of discretization steps
    # a = N(5,1)    
    # m = median(a) # the optimistic median would be 5, rather than [4.75,5.25]
    # print(m)
    # plt.plot([left(m),left(m),right(m),right(m)],[0,1,1,0])
    # plt.plot([1,7],[0.5,0.5]);  plot(a)    
    # plt.show()   
    # setoption(steps=9)  # odd number of discretization steps   
    # a = N(5,1)    
    # m = median(a) # same as the optimistic median
    # print(m)
    # plt.plot([left(m),left(m),right(m),right(m)],[0,1,1,0])
    # plt.plot([1,7],[0.5,0.5]);  plot(a)    
    
    def iqr(self,tight):
        return Interval(left(self.cut(0.25,tight)), right(self.cut(0.75,tight)))
    
    def __repr__(self):
        return (f"Pbox(range=[{self.left()}, {self.right()}], "
                f"mean="+Interval(self.ml, self.mh,auto=False).__repr__()+")")
     
    def identical(self, other):
        if not isinstance(other, Pbox): other = Pbox(other)
        return all(same(x,y) for x,y in zip(self.u, other.u)) \
               and all(same(x,y) for x,y in zip(self.d, other.d)) \
               and same(self.ml, other.ml) and same(self.mh, other.mh)

    def __neg__(self):
        return Pbox(-self.d[::-1], -self.u[::-1],
                    ml=-self.mh, mh=-self.ml)

    def __abs__(self):   # doesn't look perfect; cf. abs(N(0,1)) with Risk Calc
        if 0 <= left(self): return(self)   
        if right(self) <= 0: return(-self) 
        u = np.sort(abs(self.u))
        d = np.sort(abs(self.d))
        return Pbox(u=u, d=d, ml=np.mean(u), mh=np.mean(d)) 

    def __add__(self, other):
        if isinstance(other, Pbox):
            return conv_pbox(self, other, op="+")
        if isinstance(other, Interval):
            return self + as_pbox(other)
        return Pbox(self.u + other, self.d + other,
                    ml=self.ml + other, mh=self.mh + other)

    __radd__ = __add__
  
    def __sub__(self, other):
        other = as_pbox(other)
        return conv_pbox(self, -other, op="+")

    def __rsub__(self, other): return as_pbox(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, Pbox):
            return conv_pbox(self, other, op="*")
        if isinstance(other, Interval):
            return self * as_pbox(other)
        if other >= 0:
            return Pbox(self.u * other, self.d * other,
                        ml=self.ml * other, mh=self.mh * other)
        return - (self * (-other))

    __rmul__ = __mul__
    
    def __truediv__(self, other): # self/other
        if is_pbox(other): return self * reciprocate(other)
        if is_interval(other) or is_scalar(other):
            return self * reciprocate(Pbox(other))
        return NotImplemented  
    
    def __rtruediv__(self, other): return reciprocate(self) * other # other/self

    def __eq__(self, other): # Equality is not generally meaningful for distributions or p-boxes
        raise ValueError("Equality comparisons are not meaningful for distributions or p-boxes")

    def __ne__(self, other): # Equality is not generally meaningful for distributions or p-boxes
        raise ValueError("Equality comparisons are not meaningful for distributions or p-boxes")
     
    def prob(self, s=0):
        return Interval( len(self.d[self.d<s])/self.n, len(self.u[self.u<=s])/self.n)

    def xprob(self, s=0):              # required by the inequality comparisons
        return Interval( len(self.d[self.d<=s])/self.n, len(self.u[self.u<s])/self.n)

    # A < B = Pr(A − B < 0)
    # A ≤ B = Pr(A − B ≤ 0)
    # A > B = Pr(B − A < 0)
    # A ≥ B = Pr(B − A ≤ 0)        
     
    def __lt__(self, other): return (self-other).prob()
    def __le__(self, other): return (self-other).xprob()   
    def __gt__(self, other): return (other-self).xprob()
    def __ge__(self, other): return (other-self).prob()   
    def __rlt__(self, other): return Pbox(other).__lt__(self)
    def __rle__(self, other): return Pbox(other).__le__(self)
    def __rgt__(self, other): return Pbox(other).__gt__(self)
    def __rge__(self, other): return Pbox(other).__ge__(self)
    
    def summary(self):
        return {     
            'name':       '',
            'units':      '',
            'shape':      '',
            'mean':       mean(self),
            'sd':         sd(self),
            'var':        var(self),
            'breadth':    breadth(self),
            'iqwidth':    self.iqr().width(),
            'iqr':        self.iqr(),
            'support':    support(self),
            'left':       left(self),
            'pc01':       cut(self, 0.01),
            'pc05':       cut(self, 0.05),
            'pc25':       cut(self, 0.25),
            'median':     cut(self, 0.50),
            'pc75':       cut(self, 0.75),
            'pc95':       cut(self, 0.95),
            'pc99':       cut(self, 0.99),
            'right':      right(self),
            'steps':      steps(self) }

def summarize(x: Pbox): 
    labs = {     
        'line0':      'Summary',
        'name':       '  Name:                ',
        'units':      '  Units:               ',
        'shape':      '  Shape:               ',
        'mean':       '  Average:             ',
        'sd':         '  Std dev:             ',
        'var':        '  Variance:            ',
        'breadth' :   '  Breadth:             ',
        'iqwidth':    '  Interquartile width: ',
        'iqr':        '  Interquartile range: ',
        'support':    '  Support range:       ',
        'line1':      '  Order statistics',
        'left':       '    Left (min) value:  ',
        'pc01':       '    1st percentile:    ',
        'pc05':       '    5th percentile:    ',
        'pc25':       '    25th percentile:   ',
        'median':     '    Median (50th%ile): ',
        'pc75':       '    75th percentile:   ',
        'pc95':       '    95th percentile:   ',
        'pc99':       '    99th percentile:   ',
        'right':      '    Right (max) value: ',
        'steps':      '  Discretizations:     ',         }
    xs = x.summary()
    for k,v in labs.items(): print(v, xs.get(k, '')) 
"""    

"""


# def interval_sin(I):
#     a, b = I.lo, I.hi
#     # endpoint values
#     sa = math.sin(a)
#     sb = math.sin(b)
#     lo = min(sa, sb)
#     hi = max(sa, sb)
#     # check if interval contains a maximum (+1)
#     # maxima at pi/2 + 2k*pi
#     k1 = math.ceil((a - math.pi/2) / (2*math.pi))
#     k2 = math.floor((b - math.pi/2) / (2*math.pi))
#     if k1 <= k2:
#         hi = 1.0
#     # check if interval contains a minimum (-1)
#     # minima at 3*pi/2 + 2k*pi
#     k1 = math.ceil((a - 3*math.pi/2) / (2*math.pi))
#     k2 = math.floor((b - 3*math.pi/2) / (2*math.pi))
#     if k1 <= k2:
#         lo = -1.0
#     return Interval(lo, hi)


# NOT USED
def interval_template(I, f_scalar, extrema_points, extrema_values, domain=None):
   if I not in domain:
        if IvO.quieterrors: 
            try: I = retaindomain(I,domain) 
            except ValueError:
                if isinstance(I,Interval): raise ValueError("Interval outside "+f_scalar.__name__+" domain")
    a,  b  = I.lo,         I.hi
    fa, fb = f_scalar(a),  f_scalar(b)
    lo, hi = min(fa, fb),  max(fa, fb)

    # evaluate extrema inside the interval
    for x0 in extrema_points(a, b):
        fx0 = extrema_values(x0)
        lo = min(lo, fx0)
        hi = max(hi, fx0)

    return Interval(lo, hi)


def sin_extrema_points(a, b):
    # maxima at pi/2 + 2k*pi
    # minima at 3*pi/2 + 2k*pi
    pts = []

    # maxima
    k1 = math.ceil((a - math.pi/2) / (2*math.pi))
    k2 = math.floor((b - math.pi/2) / (2*math.pi))
    pts += [math.pi/2 + 2*k*math.pi for k in range(k1, k2+1)]

    # minima
    k1 = math.ceil((a - 3*math.pi/2) / (2*math.pi))
    k2 = math.floor((b - 3*math.pi/2) / (2*math.pi))
    pts += [3*math.pi/2 + 2*k*math.pi for k in range(k1, k2+1)]

    return pts

def sin_extrema_values(x):
    # maxima = +1, minima = -1
    return 1.0 if abs((x - math.pi/2) % (2*math.pi)) < 1e-12 else -1.0

def interval_sin(I):
    return interval_template(I, math.sin, sin_extrema_points, sin_extrema_values)




#################################################
# MIN R implementation

op = 'pmin'
for (i in 1:n)
        {
          j <- i:n
          k <- n:i
          zd[[i]]  <- min(do.call(op,list(x@d[j],y@d[k])))  #zd[[i]] <- min(x@d[j] + y@d[k])
          j <- 1:i
          k <- i:1
          zu[[i]]  <- max(do.call(op,list(x@u[j],y@u[k])))  #zu[[i]] <- max(x@u[j] + y@u[k])
        }


#################################################
# MIN C++ implementation

    case RandomNbr::minimum:
      {
          for (i=0; i<n; i++)
          {
              outlier = toobig;
              for (j = i; j<n; j++)
              {
                  if (x.d[j] < y.d[i - j + n - 1])
                      here = x.d[j];
                  else
                      here = y.d[i - j + n - 1];
                  if (here<outlier) outlier = here;
              }
              z.d[i] = outlier;

              outlier = -toobig;
              for (j = 0; j<=i; j++)
              {
                  if (x.u[j] < y.u[i - j])
                      here = x.u[j];
                  else
                      here = y.u[i - j];
                  if (here>outlier) outlier = here;
              }
              z.u[i] = outlier;
          }
          if (x.right()<y.left())
          {
              z.mymean = x.mean();
              z.myvar = x.variance();
          }
          else if (y.right()<x.left())
          {
              z.mymean = y.mean();
              z.myvar = y.variance();
          }
          else
          {
              z.mymean = imp(env(intmin(x.left(),y.left()), intmin(x.mean(),y.mean())),
                      x.mean()+y.mean()-env(intmax(x.right(), y.right()), intmax(x.mean(),y.mean())));
              z.myvar = env(intmax(x.variance(),y.variance()),0.0);
          }
          // This causes the infamous min(3,N(5,1)) bug
          //z.mymean = imp(z.mymean, VKmeanminimum(x, y, RandomNbr::dw));
      }  //case
      break;

    

# class Interval:
#     def __init__(self, lo, hi=None, auto=True):
#         '''
#         Interval(12,22)                                   # [12.0, 22.0]
#         Interval(22,12)                                   # [12.0, 22.0]
#         Interval(12)                                      # 12.0
#         Interval(12,)                                     # 12.0
#         Interval(12,None)                                 # 12.0
#         Interval(12,NA)                                   # [12.0, inf]
#         Interval(12,np.inf)                               # [12.0, inf]
#         Interval(12,float('inf'))                         # [12.0, inf]
#         Interval(NA,NA)                                   # [-inf, inf]
#         Interval('12')                                    # [11.5, 12.5]
#         Interval(Interval(12,13), Interval(21,22))        # [12, 22]
#         Interval(Interval(12,50), Interval(9,22))         # [12, 22]
#         '''
#         def _scalarize(x):          # numPy scalar (np.float64, np.int64, etc.)
#             if isinstance(x, np.generic): return float(x)
#             if isinstance(x, np.ndarray) and x.size == 1: return float(x)
#             return x
#         lo = _scalarize(lo)
#         hi = _scalarize(hi)
#         if is_na(lo): lo = float("-inf")
#         if is_na(hi): hi = float("inf")
#         if isinstance(lo, str) and hi is None:     # significant-digit interval
#             lo_val, hi_val = sgnumber(lo)
#             self.lo = float(lo_val)
#             self.hi = float(hi_val)
#         elif hi is None and isinstance(lo, Interval):       # existing interval
#             self.lo = float(lo.lo)
#             self.hi = float(lo.hi)
#         elif hi is None and np.isscalar(lo):              # degenerate interval
#             self.lo = float(lo)
#             self.hi = float(lo)
#         elif hi is None and isinstance(lo, np.ndarray) and lo.size == 2:
#             self.lo = float(lo[0])
#             self.hi = float(lo[1])
#         elif hi is None and hasattr(lo,"__len__") and len(lo)==2: # tuple, list
#             self.lo = float(lo[0])
#             self.hi = float(lo[1])
#         elif hi is not None:                                  # explicit bounds
#             self.lo = float(left(lo))
#             self.hi = float(right(hi))
#         else: raise ValueError("Bad interval input")
#         if self.hi < self.lo and IvO.autocorrect and auto:
#             self.lo, self.hi = self.hi, self.lo     


#################################################
# Scott's early draft Python implementation
    
MAX

def pmax(x,y,op='max'):
    toobig = inf   
    n = steps(x)
    zu = np.full(n,0)
    zd = np.full(n,0)
    if   op == "max": 
        for i in range(n):       
            outlier = toobig
            for j in range(i,n):    # (j = i; j<n; j++)
                if (x.d[j] > y.d[i - j + n - 1]):  here = x.d[j]
                else: here = y.d[i - j + n - 1]
                if here < outlier: outlier = here
            zd[i] = outlier

            outlier = -toobig
            for j in range(i+1):    # (j = 0; j<=i; j++)
                if (x.u[j] > y.u[i - j]): here = x.u[j]
                else: here = y.u[i - j]
                if here > outlier: outlier = here
            zu[i] = outlier
        if right(y) < left(x): ml, mh = mean(x)
        elif right(x) < left(y): ml, mh = mean(y)
        else: ml, mh = imp(env(smax(mean(x),mean(y)), smax(right(x),right(y))), mean(x)+mean(y)-env(smin(left(x),left(y)), smin(mean(x),mean(y))))
        # ml, mh = imp(z.mymean, VKmeanmaximum(x, y, RandomNbr::dw))   # causes the max(3,N(5,1)) bug
    return Pbox(zu,zd,ml,mh)
  
# ..... it got better!

def frechet_min(x, y):
    n = len(x.d)
    assert len(y.d) == n == len(x.u) == len(y.u)
    zd = np.empty(n, dtype=float)
    zu = np.empty(n, dtype=float)
    for i in range(n):
        outlier = inf
        for j in range(i, n):
            k = i - j + n - 1
            here = min(x.d[j], y.d[k])
            if here < outlier: outlier = here
        zd[i] = outlier
        outlier = -inf
        for j in range(0, i + 1):
            k = i - j
            here = min(x.u[j], y.u[k])
            if here > outlier: outlier = here
        zu[i] = outlier
    return Pbox(zu, zd)    

def frechet_max(x, y):
    n = len(x.d)
    assert len(y.d) == n == len(x.u) == len(y.u)
    zd = np.empty(n, dtype=float)
    zu = np.empty(n, dtype=float)
    for i in range(n):
        outlier = inf
        for j in range(i, n):
            k = i - j + n - 1
            here = max(x.d[j], y.d[k])
            if here < outlier: outlier = here
        zd[i] = outlier
        outlier = -inf
        for j in range(0, i + 1):
            k = i - j
            here = max(x.u[j], y.u[k])
            if here > outlier: outlier = here
        zu[i] = outlier
    if right(y) < left(x): ml, mh = mean(x)
    elif right(x) < left(y): ml, mh = mean(y)
    else: ml, mh = imp(env(smax(mean(x),mean(y)), smax(right(x),right(y))), mean(x)+mean(y)-env(smin(left(x),left(y)), smin(mean(x),mean(y))))
    ml, mh = imp(Interval(np.mean(zu),np.mean(zd)),Interval(ml,mh))
    # ml, mh = imp(z.mymean, VKmeanmaximum(x, y, RandomNbr::dw))   # causes the max(3,N(5,1)) bug
    return Pbox(zu, zd, ml, mh)    
  
  
def frechetconvFAST(x, y, op='+'):  # uses a "blisterinly fast" Toeplitz‑like index matrix
    if op=='-': return frechetconv(x,(-y),'+')
    if op=='/': return frechetconv(x,reciprocate.pbox(y),'*')

    xd, xu = x.d, x.u     # OMIT THESE
    yd, yu = y.d, y.u
    n = len(xd)

    if   op == '+':        f = np.add
    elif op == '*':        f = np.multiply
    elif op == min:        f = np.minimum
    elif op == max:        f = np.maximum
    else:                  raise ValueError("Unsupported op")

    # LOWER ENVELOPE (d)
    # Build full pairwise matrix for lower side
    Md = f(xd[:,None], yd[::-1][None,:])   # shape (n,n)

    # Extract anti-diagonals: flip left-right, then take diagonals
    Md_flip = np.fliplr(Md)
    zd = np.array([Md_flip.diagonal(i) for i in range(n)])
    zd = np.min(zd, axis=1)   # min over each anti-diagonal

    # UPPER ENVELOPE (u)
    Mu = f(xu[:,None], yu[None,:])         # shape (n,n)

    # Extract forward diagonals
    zu = np.array([Mu.diagonal(i) for i in range(n)])
    zu = np.max(zu, axis=1)   # max over each diagonal

    # Mean/variance logic (unchanged)
    ml = -np.inf
    mh = np.inf
    if op == '+':        ml, mh = ends(x.mean() + y.mean())
    elif op == '*':        ml, mh = ends(x.mean() * y.mean())
    elif op == min:
        if right(y) < left(x): ml, mh = mean(y)
        elif right(x) < left(y): ml, mh = mean(x)
        else: ml, mh = imp(env(smin(mean(x),mean(y)), smin(right(x),right(y))),
                mean(x)+mean(y)-env(smax(left(x),left(y)), smax(mean(x),mean(y))))
        ml, mh = imp(Interval(np.mean(zu),np.mean(zd)), Interval(ml,mh))
    elif op == max:
        if right(y) < left(x): ml, mh = mean(x)
        elif right(x) < left(y): ml, mh = mean(y)
        else: ml, mh = imp(env(smax(mean(x),mean(y)), smax(right(x),right(y))),
                mean(x)+mean(y)-env(smin(left(x),left(y)), smin(mean(x),mean(y))))
        ml, mh = imp(Interval(np.mean(zu),np.mean(zd)), Interval(ml,mh))
    return Pbox(u=zu, d=zd, ml=ml, mh=mh)    
 



# Assumed available in your environment:
#   class Pbox(u, d, ml=None, mh=None)
#   class Interval(lo, hi)
#   PbO.steps
#   left(x), right(x), mean(x), ends(x), smin(a,b), smax(a,b)
#   imp(a, b), env(a, b)
#   straddles(x), is_zero(x)
#   reciprocate.pbox(y)

def _frechet_pairwise_envelopes(xd, xu, yd, yu, op):
    '''
    Mixed-resolution Fréchet envelopes for one of: '+', '*', min, max.
    xd, xu, yd, yu are 1D arrays (possibly different lengths).
    Returns zd_raw, zu_raw with length nx + ny - 1.
    '''
    xd = np.asarray(xd, float)
    xu = np.asarray(xu, float)    # OMIT THESE
    yd = np.asarray(yd, float)
    yu = np.asarray(yu, float)

    nx, ny = len(xd), len(yd)
    n_out = nx + ny - 1

    # Choose pairwise operator
    if op == '+':        f = np.add
    elif op == '*':        f = np.multiply
    elif op == min:        f = np.minimum
    elif op == max:        f = np.maximum
    else:        raise ValueError(f"Unsupported op {op}")

    # LOWER ENVELOPE (d): anti-diagonals of f(xd, reversed yd)
    Md = f(xd[:, None], yd[::-1][None, :])  # shape (nx, ny)
    Md_flip = np.fliplr(Md)

    zd_raw = np.empty(n_out, dtype=float)
    # anti-diagonals: offsets from -(nx-1) to +(ny-1)
    for idx, k in enumerate(range(-(nx - 1), ny)):
        diag = Md_flip.diagonal(k)
        if op in ('+', '*', min):
            # Fréchet min / sum / product: min over pairwise values
            zd_raw[idx] = np.min(diag)
        elif op is max:
            # Fréchet max: min over max-values
            zd_raw[idx] = np.min(diag)

    # UPPER ENVELOPE (u): forward diagonals of f(xu, yu)
    Mu = f(xu[:, None], yu[None, :])        # shape (nx, ny)

    zu_raw = np.empty(n_out, dtype=float)
    for idx, k in enumerate(range(-(nx - 1), ny)):
        diag = Mu.diagonal(k)
        # For all four ops, upper envelope is max over pairwise values
        zu_raw[idx] = np.max(diag)

    return zd_raw, zu_raw


def _resample_envelopes(zd_raw, zu_raw, n_target):
    '''
    Resample raw envelopes (length m) to length n_target using
    monotone-preserving interpolation in index space.
    '''
    m = len(zd_raw)
    if m == n_target: return zd_raw.copy(), zu_raw.copy()
    old_idx = np.linspace(0.0, 1.0, m)
    new_idx = np.linspace(0.0, 1.0, n_target)
    zd = np.interp(new_idx, old_idx, zd_raw)
    zu = np.interp(new_idx, old_idx, zu_raw)
    return zd, zu


def frechetconv(x, y, op='+'):
    '''
    Mixed-resolution Fréchet convolution for +, -, *, /, min, max.
    Result is resampled to PbO.steps discretization levels.
    '''
    # Reduce to +, *, min, max
    if op == '-':        return frechetconv(x, -y, '+')
    if op == '/':        return frechetconv(x, reciprocate.pbox(y), '*')

    # Sign logic for multiplication
    if op == '*':
        if straddles(x) or straddles(y):        raise NotImplementedError("Fréchet product straddling zero not implemented")
        if is_zero(x) or is_zero(y):            return 0  # or Pbox.zero(...) in your style
        if (right(x) <= 0) and (right(y) <= 0): return frechetconv(-x, -y, '*')
        if right(x) <= 0:                       return -frechetconv(-x, y, '*')
        if right(y) <= 0:                       return -frechetconv(x, -y, '*')

    # Core envelopes (mixed resolution)
    zd_raw, zu_raw = _frechet_pairwise_envelopes(x.d, x.u, y.d, y.u, op)

    # Resample to PbO.steps
    n_target = PbO.steps
    zd, zu = _resample_envelopes(zd_raw, zu_raw, n_target)

    # Mean / variance logic
    ml, mh = -np.inf, np.inf
    if op == '+':          ml, mh = ends(mean(x) + mean(y))
    elif op == '*':          ml, mh = ends(mean(x) * mean(y))
    elif op is min:
        if right(y) < left(x):            ml, mh = mean(y)
        elif right(x) < left(y):            ml, mh = mean(x)
        else: ml, mh = imp(env(smin(mean(x), mean(y)), smin(right(x), right(y))),
                mean(x) + mean(y) - env(smax(left(x), left(y)), smax(mean(x), mean(y))))
        ml, mh = imp(Interval(np.mean(zu), np.mean(zd)), Interval(ml, mh))
    elif op is max:
        if right(y) < left(x):            ml, mh = mean(x)
        elif right(x) < left(y):            ml, mh = mean(y)
        else: ml, mh = imp(env(smax(mean(x), mean(y)), smax(right(x), right(y))),
                mean(x) + mean(y) - env(smin(left(x), left(y)), smin(mean(x), mean(y))))
        ml, mh = imp(Interval(np.mean(zu), np.mean(zd)), Interval(ml, mh))
    else: raise ValueError(f"Unsupported operation {op}")
    print(zu)
    print(zd)
    return Pbox(u=zu, d=zd, ml=ml, mh=mh)






###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Significant digit displays and interpretations


def sformat_sigdigs(a, D):               # format_sigdigs(100/3,3) yields '33.3'
    fmt = f"{{:.{D}g}}"
    return f"{fmt.format(a)}"

def cformat_sigdigs(x, D):
    '''Format x with exactly D significant digits, preserving significance rules.'''
    sign = "-" if x < 0 else ""
    x = abs(x)
    if x == 0: return sign + "0." + "0"*(D-1)    # zero is a special case

    
    fmt = f"{{:.{D}g}}"
    s = fmt.format(x)                # get correct rounded value using g-format

    # If scientific notation appears, we must normalize ourselves
    if "e" in s or "E" in s:
        mantissa, exp = s.lower().split("e")
        exp = int(exp)

        # Normalize mantissa to D significant digits with trailing zeros
        if "." in mantissa:
            digits = mantissa.replace(".", "")
            digits = digits + "0"*(D - len(digits))
            mantissa = digits[0] + "." + digits[1:]
        else:
            # integer mantissa
            digits = mantissa + "0"*(D - len(mantissa))
            mantissa = digits[0] + "." + digits[1:]

        return f"{sign}{mantissa}e{exp:+d}"

    # Step 2: If s already contains a decimal point, enforce trailing zeros
    if "." in s:
        # Remove exponent (already handled)
        whole, frac = s.split(".")
        digits = whole + frac

        # Pad or trim to exactly D significant digits
        digits = digits + "0"*(D - len(digits))
        digits = digits[:D]

        # Reinsert decimal point in the correct place
        if len(whole) >= D:
            # No fractional digits needed
            return sign + whole[:D] + "."
        else:
            return sign + whole + "." + digits[len(whole):]

    # Step 3: No decimal point → integer-like string
    # Determine how many digits we have
    digits = s
    digits = digits + "0"*(D - len(digits))
    digits = digits[:D]

    # If all digits are to the left of decimal, append decimal point
    if len(digits)<=len(s): return sign + digits + "."
    return sign + digits[:len(s)] + "." + digits[len(s):]      # else add point

def test_format_sigdigs(x,D,s): 
    print('format_sigdigs('+str(x)+','+str(D)+")  # '"+s+"'"+"     '"+sformat_sigdigs(x,D)+"'"+"     '"+cformat_sigdigs(x,D)+"'")


#The string displayed after the # is what your algorithm yields.


test_format_sigdigs(1000.3,1, '1.e+3') #   '1.e+3'   '1.e+3'
test_format_sigdigs(1000.3,2, '1.0e+3') #   '1.0e+3'   '1.0e+3' *
test_format_sigdigs(1000.3,3, '1.00e+3') #     *
test_format_sigdigs(1000.3,4, '1000.') #     *
test_format_sigdigs(1000.3,5, '1000.3') #
test_format_sigdigs(1000.3,6, '1000.30') #     *
test_format_sigdigs(1234.5678,9, '1234.56780') #     *
test_format_sigdigs(1234.5678,8, '1234.5678') #
test_format_sigdigs(1234.5678,7, '1234.568') #
test_format_sigdigs(1234.5678,6, '1234.57') #
test_format_sigdigs(1234.5678,5, '1234.6') #
test_format_sigdigs(1234.5678,4, '1235') #    @ '1235.'    
test_format_sigdigs(1234.5678,3, '1230') #   @'1.23e+3'        *
test_format_sigdigs(1234.5678,2, '1200') # @'1.2e+3'        *
test_format_sigdigs(1234.5678,1, '1000') #  @'1.e+3'        *
test_format_sigdigs(6666.6666,1, '7000') #   @ '7.e+3'         *
test_format_sigdigs(6666.6666,2, '6700') #  @ '6.7e+3'       *
test_format_sigdigs(6666.6666,3, '6670') # @ '6.67e+3'      *
test_format_sigdigs(6666.6666,4, '6667') # @ '6667.'   
test_format_sigdigs(6666.6666,5, '6666.7') #
test_format_sigdigs(6666.6666,6, '6666.67') #
test_format_sigdigs(6666.6666,7, '6666.667') #
test_format_sigdigs(6666.6666,8, '6666.6666') #
test_format_sigdigs(6666.6666,9, '6666.66660') #     *

test_format_sigdigs(1000.003,6, '1000.00')
test_format_sigdigs(1000.007,6, '1000.01')


Below is a side-by-side comparison of the outputs of the original function format_sigdigs() and your revised function with what I think is the correct output.

The first string after the '#' is what I think the answer SHOULD BE.  The second string is the result returned by my original, very simple function, and the third string is the result given by your latest version of the code.  I've annotated the end of each line with an S if the second string disagrees with the first, and with a C if the third string disagrees with the first.  We can see that there are several more S than C marks, indicating that your latest version of the function is a big improvement.

The first and last lines with a C for the inputs (1234.5678,4) and (6666.6666,4) are not really errors, although we'd slightly perfer not to have the decimal point in the output if it's not actually needed in that instance.  But we much prefer to have an unneeded decimal point here rather than to not have it where it's needed, as for the output for the inputs (1000.3,4).

The other lines with the C marker indicate a more serious complaint, which applies to both the original function and your revised version of the function.  In these cases, scientific notation is used, which of course is not wrong at all, but the better answers in these cases are the simple decimal numbers.  Definitely a rule of this algorithm should be that, if an explicit decimal expansion characterizes the value and its precision correctly, it is prefered over scientific notation.  Code that obeys this rule is much more humane.  Can you please update your revised code to reflect this rule?

               Input           Correct     Scott       Copilot  
format_sigdigs(1000.3,1)     # '1.e+3'     '1e+03'     '1.e+3'
format_sigdigs(1000.3,2)     # '1.0e+3'     '1e+03'     '1.0e+3'        S
format_sigdigs(1000.3,3)     # '1.00e+3'     '1e+03'     '1.00e+3'      S
format_sigdigs(1000.3,4)     # '1000.'     '1000'     '1000.'           S
format_sigdigs(1000.3,5)     # '1000.3'     '1000.3'     '1000.3'
format_sigdigs(1000.3,6)     # '1000.30'     '1000.3'     '1000.30'      S
format_sigdigs(1000.003,6)   # '1000.00'     '1000'     '1000.00'         S
format_sigdigs(1000.007,6)   # '1000.01'     '1000.01'     '1000.01'
format_sigdigs(1234.5678,9)  # '1234.56780'     '1234.5678'     '1234.56780'    S
format_sigdigs(1234.5678,8)  # '1234.5678'     '1234.5678'     '1234.5678'
format_sigdigs(1234.5678,7)  # '1234.568'     '1234.568'     '1234.568'
format_sigdigs(1234.5678,6)  # '1234.57'     '1234.57'     '1234.57'
format_sigdigs(1234.5678,5)  # '1234.6'     '1234.6'     '1234.6'
format_sigdigs(1234.5678,4)  # '1235'     '1235'     '1235.'              C
format_sigdigs(1234.5678,3)  # '1230'     '1.23e+03'     '1.23e+3'       S C
format_sigdigs(1234.5678,2)  # '1200'     '1.2e+03'     '1.2e+3'        S C
format_sigdigs(1234.5678,1)  # '1000'     '1e+03'     '1.e+3'           S C
format_sigdigs(6666.6666,1)  # '7000'     '7e+03'     '7.e+3'           S C
format_sigdigs(6666.6666,2)  # '6700'     '6.7e+03'     '6.7e+3'          S C
format_sigdigs(6666.6666,3)  # '6670'     '6.67e+03'     '6.67e+3'          S C
format_sigdigs(6666.6666,4)  # '6667'     '6667'     '6667.'               C
format_sigdigs(6666.6666,5)  # '6666.7'     '6666.7'     '6666.7'
format_sigdigs(6666.6666,6)  # '6666.67'     '6666.67'     '6666.67'
format_sigdigs(6666.6666,7)  # '6666.667'     '6666.667'     '6666.667'
format_sigdigs(6666.6666,8)  # '6666.6666'     '6666.6666'     '6666.6666'
format_sigdigs(6666.6666,9)  # '6666.66660'     '6666.6666'     '6666.66660'      S

I'd hoped that some of my complaints about your latest algorithm would not appear in the second strings created by the old, simple algorithm.  But, apart from the minor issue with nonessential decimal points, mostly the C complaints are a subset of the S complaints.  Nevertheless, my programmer's instinct is to ask whether an efficient solution is perhaps to create the output as before (or maybe in two ways) and then check whether it has the right number of significant digits using sgnumber() before investing in a post-hoc emergency correction.

def sgnumber(user_input: str):     # number ± its significant-digit imprecision
    user_input = user_input.strip().lower()
    tens = '0'
    if 'e' in user_input: mantissa, tens = user_input.split('e', 1)
    else: mantissa = user_input
    if '.' in mantissa: j = len(mantissa.split('.')[1])
    #else: j = len(mantissa.split('0', 1)) - len(mantissa) + 1           
    else: j = len(mantissa.rstrip('0')) - len(mantissa) 
    pm = 10**(-j) * 10**int(tens) / 2
    #print('input:',user_input,', mantissa:',mantissa, ', j:',j, ', tens:',tens, ', pm:',pm)
    return([float(user_input)-pm, float(user_input)+pm])

       



###############################################################################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Fermi stuff







def fermi(shape, lower=None, upper=None, pr=0.9, data=None):
    def fermi_N(x1, x2, n=None, pr=0.9):
        m = (x1 + x2) / 2
        if n is None: s = (x2 - x1) / (2 * qnormal(pr)) 
        else: s = (x2 - x1) / mean_normal_range(n)
        return (m, s)
    def fermi_L(x1, x2, n=None, pr=0.9):
        gm = np.sqrt(x1*x2)
        if n is None: gsd = np.sqrt(x2/x1) ** (1/qnormal(pr)) 
        else: gsd = np.exp((np.log(x2) - np.log(x1)) / mean_normal_range(n))
        print(gm,gsd)
        return np.log((gm, gsd))
    def fermi_W(x1, x2, n=None, pr=0.9):
        p1, p2 = (1 - pr) / 2, (1 + pr) / 2  # symmetric quantiles
        if n is None: r = x2 / x1  # adjust range if sample size is given
        else: r = np.exp((np.log(x2) - np.log(x1)) / mean_normal_range(n))
        A = -np.log(1 - p1)
        B = -np.log(1 - p2)
        k = np.log(B/A) / np.log(r)  # solve for shape k 
        lam = x1 / (A ** (1/k)) # scale λ
        return k, lam
    def fermi_G(x1, x2, n=None, pr=0.9):
        m = (x1 + x2) / 2       # Fermi-style mean and sd
        if n is None: s = (x2 - x1) / (2 * qnormal(pr))
        else: s = (x2 - x1) / mean_normal_range(n)
        alpha = (m / s)**2   # convert mean/sd to gamma parameters
        theta = s**2 / m
        return alpha, theta
    def fermi_logistic(x1, x2, n=None, pr=0.9):
        p1 = (1 - pr) / 2
        p2 = (1 + pr) / 2
        L = np.log(p2/(1-p2)) - np.log(p1/(1-p1))
        if n is None: s = (x2 - x1) / L
        else: # adjust effective range
            R = (x2 - x1) / mean_normal_range(n)   
            s = R / L
        mu = (x1 + x2) / 2
        return mu, s
    if data is None: n = None
    else: lower, upper, n = *support(data), len(data)
    if shape in ('normal','N', normal): return normal(*fermi_N(lower,upper,n,pr))
    if shape in ('lognormal','L', lognormal): return lognormal2(*fermi_L(lower,upper,n,pr))
    if shape in ('weibull','W', weibull): return weibull(*fermi_W(lower, upper, n, pr))
    if shape in ('gamma','G', gamma): return gamma(*fermi_G(lower, upper, n, pr))
    if shape in ('logistic','S', logistic): return logistic(*fermi_logistic(lower, upper, n, pr))
    raise ValueError('unknown distribution')
    
        
from scipy.optimize import root

def fermi_from_quantiles(
    ppf,                  # function ppf(p, *theta)
    x1, x2,               # elicited lower, upper
    n=None, pr=0.9,       # sample size, central probability
    theta0=None,          # initial guess for parameters
    transform=None,       # optional: (to_unconstrained, from_unconstrained)
):
    '''
    Generic Fermi-style parameter inference from two symmetric quantiles.
    '''

    # symmetric quantiles
    p1 = (1 - pr) / 2
    p2 = (1 + pr) / 2

    # optional sample-range adjustment (normal-based)
    if n is not None:
        # effective 'central' range that would give same expected sample range
        R = (x2 - x1) / mean_normal_range(n)
        m = (x1 + x2) / 2
        x1_eff = m - R/2
        x2_eff = m + R/2
    else:
        x1_eff, x2_eff = x1, x2

    # default transform: identity
    if transform is None:
        def to_u(theta):   return np.array(theta, dtype=float)
        def from_u(u):     return u
    else:
        to_u, from_u = transform

    # initial guess
    if theta0 is None:
        # crude: center and scale from elicited interval
        m = (x1_eff + x2_eff) / 2
        s = (x2_eff - x1_eff) / (2 * qnormal(pr))
        theta0 = (m, s)

    u0 = to_u(theta0)

    def residual(u):
        theta = from_u(u)
        q1 = ppf(p1, *theta)
        q2 = ppf(p2, *theta)
        return np.array([q1 - x1_eff, q2 - x2_eff])

    sol = root(residual, u0)
    if not sol.success:
        raise RuntimeError("Fermi quantile inversion failed: " + sol.message)

    theta_hat = from_u(sol.x)
    return tuple(theta_hat)    
   









# ============================================================
# Quantile wrappers (using your conventions)
# ============================================================

def qnormal(p, mu=0, sigma=1):
    return sps.norm(loc=mu, scale=sigma).ppf(p)

# def qlognormal(p, m, s):
#     # mean/sd -> mu/sigma
#     sigma = np.sqrt(np.log(1 + (s*s)/(m*m)))
#     mu = np.log(m) - 0.5*sigma*sigma
#     return sps.lognorm.ppf(p, s=sigma, scale=np.exp(mu))

# def qweibull(p, k, lam):
#     return sps.weibull_min.ppf(p, c=k, scale=lam)

# def qgamma(p, shape, scale=1, rate=None):
#     if rate is not None:
#         scale = 1/rate
#     return sps.gamma.ppf(p, shape, scale=scale)

# def qlogistic(p, mu, s):
#     return sps.logistic.ppf(p, loc=mu, scale=s)

# ============================================================
# Constructors (adapt these to your qpbox/PbO setup)
# ============================================================

# def normal(mu, sigma):
#     # replace with your qpbox-based constructor
#     return (mu, sigma)

# def lognormal(m, s):
#     # replace with your qpbox-based constructor
#     return (m, s)

# def weibull(k, lam):
#     # replace with your qpbox-based constructor
#     return (k, lam)

# def gamma(shape, scale=1, rate=None):
#     if rate is not None:
#         scale = 1/rate
#     # replace with your qpbox-based constructor
#     return (shape, scale)

# def logistic(mu, s):
#     # replace with your qpbox-based constructor
#     return (mu, s)

# ============================================================
# mean_range with memoization (for scale families only)
# ============================================================

_mean_range_cache = {}

def mean_range(ppf, n, sims=200000):
    '''Monte Carlo estimate of E[max - min] for n IID samples.'''
    if n <= 1:
        return 1.0  # effectively "no correction"
    key = (id(ppf), n)
    if key in _mean_range_cache:
        return _mean_range_cache[key]
    u = np.random.rand(sims, n)
    x = ppf(u)
    mr = np.mean(x.max(axis=1) - x.min(axis=1))
    _mean_range_cache[key] = mr
    return mr

# ============================================================
# Transforms (constrained -> unconstrained)
# ============================================================

def normal_transform():
    def to_u(theta):
        mu, sigma = theta
        return np.array([mu, np.log(sigma)])
    def from_u(u):
        mu, log_sigma = u
        return (mu, np.exp(log_sigma))
    return to_u, from_u

def lognormal_transform():
    # work in mean/sd space, enforce positivity
    def to_u(theta):
        m, s = theta
        return np.log([m, s])
    def from_u(u):
        m, s = np.exp(u)
        return (m, s)
    return to_u, from_u

def weibull_transform():
    def to_u(theta):
        k, lam = theta
        return np.log([k, lam])
    def from_u(u):
        k, lam = np.exp(u)
        return (k, lam)
    return to_u, from_u

def gamma_transform():
    def to_u(theta):
        a, t = theta
        return np.log([a, t])
    def from_u(u):
        a, t = np.exp(u)
        return (a, t)
    return to_u, from_u

def logistic_transform():
    def to_u(theta):
        mu, s = theta
        return np.array([mu, np.log(s)])
    def from_u(u):
        mu, log_s = u
        return (mu, np.exp(log_s))
    return to_u, from_u

# ============================================================
# Core quantile inversion
# ============================================================

def fermi_from_quantiles(ppf, x1, x2, pr=0.9, n=None,
                         theta0=None, transform=None,
                         range_ppf=None):
    '''
    Solve for parameters theta such that:
        ppf(p1, *theta) = x1
        ppf(p2, *theta) = x2

    If range_ppf is provided and n>1, treat [x1,x2] as a sample range
    and shrink/expand it using the distribution-specific mean range
    of the *scale family* defined by range_ppf.
    '''

    p1 = (1 - pr) / 2.0
    p2 = (1 + pr) / 2.0

    # Transform
    if transform is None:
        def to_u(theta): return np.array(theta, float)
        def from_u(u):   return tuple(u)
    else:
        to_u, from_u = transform

    # Provisional initial guess in parameter space (normal-based)
    if theta0 is None:
        m = (x1 + x2) / 2.0
        s = (x2 - x1) / (2.0 * sps.norm.ppf((1+pr)/2.0))
        theta0 = (m, s)

    # Range correction for scale families only
    if (n is not None) and (n > 1) and (range_ppf is not None):
        mr = mean_range(range_ppf, n)
        R = (x2 - x1) / mr
        m = (x1 + x2) / 2.0
        x1_eff = m - R/2.0
        x2_eff = m + R/2.0
    else:
        x1_eff, x2_eff = x1, x2

    u0 = to_u(theta0)

    def residual(u):
        theta = from_u(u)
        q1 = ppf(p1, *theta)
        q2 = ppf(p2, *theta)
        return np.array([q1 - x1_eff, q2 - x2_eff])

    sol = root(residual, u0)
    if not sol.success:
        raise RuntimeError("Fermi quantile inversion failed: " + sol.message)

    return from_u(sol.x)

# ============================================================
# Unified Fermi dispatcher (constructor-only)
# ============================================================

def fermi(shape, lower, upper, n=None, pr=0.9):
    '''
    Fermi elicitation:
      - shape is a constructor: normal, lognormal, weibull, gamma, logistic
      - [lower, upper] is interpreted as central pr interval (quantiles)
        unless n>1 and shape is a scale family (normal, logistic),
        in which case it's treated as a sample range and corrected.
    '''

    # Identify distribution by constructor
    if shape is normal:
        ppf, transform = qnormal, normal_transform()
        # unit-scale ppf for range correction
        range_ppf = lambda u: qnormal(u, 0.0, 1.0)

    elif shape is lognormal:
        ppf, transform = qlognormal, lognormal_transform()
        range_ppf = None  # not a pure scale family

    elif shape is weibull:
        ppf, transform = qweibull, weibull_transform()
        range_ppf = None  # shape-dependent, not pure scale

    elif shape is gamma:
        ppf, transform = qgamma, gamma_transform()
        range_ppf = None  # shape-dependent, not pure scale

    elif shape is logistic:
        ppf, transform = qlogistic, logistic_transform()
        range_ppf = lambda u: qlogistic(u, 0.0, 1.0)

    else:
        raise ValueError(f"Unsupported distribution constructor: {shape}")

    theta = fermi_from_quantiles(
            ppf,
            lower, upper,
            pr=pr,
            n=n,
            theta0=None,
            transform=transform,
            range_ppf=range_ppf
        )
    print(theta)
    return shape(*theta)



# def fermi(shape, lower, upper, n=None, pr=0.9):
#     if   shape is normal:    ppf, transform = qnormal, normal_transform()
#     elif shape is lognormal: ppf, transform = qlognormal, lognormal_transform()
#     elif shape is weibull:   ppf, transform = qweibull, weibull_transform()
#     elif shape is gamma:     ppf, transform = qgamma, gamma_transform()
#     elif shape is logistic:  ppf, transform = qlogistic, logistic_transform()
#     else: raise ValueError(f"Unsupported distribution constructor: {shape}")
#     theta = fermi_from_quantiles(ppf,lower,upper,pr=pr,n=n,theta0=None,transform=transform)
#     print('theta',theta)
#     return shape(*theta)  # build distribution using the constructor







#boundary following but NOT RIGOROUS

import numpy as np
import scipy.stats as sps

# ------------------------------------------------------------
# Lognormal PPF in mean/sd parameterization
# ------------------------------------------------------------

def qlognormal(p, m, s):
    sigma = np.sqrt(np.log(1 + (s*s)/(m*m)))
    mu    = np.log(m) - 0.5*sigma*sigma
    return sps.lognorm.ppf(p, s=sigma, scale=np.exp(mu))

def mean_range_lognormal(m, s, n, sims=50000):
    if n <= 1:
        return 0.0
    u = np.random.rand(sims, n)
    x = qlognormal(u, m, s)
    return np.mean(x.max(axis=1) - x.min(axis=1))

# ------------------------------------------------------------
# Main routine: boundary-based p-box construction
# ------------------------------------------------------------

def fermi_pbox_lognormal_boundary(lower, upper, n,
                                  tol_mean=0.05,
                                  tol_range=0.10,
                                  m_grid=40, s_grid=40,
                                  sims_range=50000):
    '''
    Returns:
      QL(p)  -- lower quantile envelope on PbO.ii()
      QU(p)  -- upper quantile envelope on PbO.jjj()
      m_interval = (m_min, m_max)
      s_interval = (s_min, s_max)
      boundary   -- list of boundary (m,s) pairs
    '''

    # --------------------------------------------------------
    # Targets
    # --------------------------------------------------------
    m_target = (lower + upper) / 2.0
    R_target = upper - lower

    # crude normal-based scale guess
    if n > 1 and R_target > 0:
        mr_norm = 3.077  # approx for n≈10; refine if desired
        s0 = R_target / mr_norm
    else:
        s0 = max(R_target, 1.0)

    # --------------------------------------------------------
    # Parameter grid
    # --------------------------------------------------------
    m_vals = np.linspace(0.5*m_target, 1.5*m_target, m_grid)
    s_vals = np.linspace(0.3*s0, 3.0*s0, s_grid)

    admissible = []

    # --------------------------------------------------------
    # Filter admissible (m,s)
    # --------------------------------------------------------
    for m in m_vals:
        for s in s_vals:
            if m <= 0 or s <= 0:
                continue

            # mean constraint
            if abs(m - m_target)/m_target > tol_mean:
                continue

            # range constraint
            R_hat = mean_range_lognormal(m, s, n, sims=sims_range)
            if abs(R_hat - R_target)/R_target > tol_range:
                continue

            admissible.append((m, s))

    if not admissible:
        raise RuntimeError("No admissible lognormal parameters found.")

    # --------------------------------------------------------
    # Extract boundary of admissible region
    # --------------------------------------------------------
    admissible = np.array(admissible)
    m_list = admissible[:,0]
    s_list = admissible[:,1]

    m_interval = (m_list.min(), m_list.max())
    s_interval = (s_list.min(), s_list.max())

    boundary = []

    # For each m-bin, keep min and max s
    m_bins = np.linspace(m_interval[0], m_interval[1], 50)
    for i in range(len(m_bins)-1):
        lo, hi = m_bins[i], m_bins[i+1]
        mask = (m_list >= lo) & (m_list < hi)
        if np.any(mask):
            s_sub = s_list[mask]
            boundary.append((m_list[mask][np.argmin(s_sub)], s_sub.min()))
            boundary.append((m_list[mask][np.argmax(s_sub)], s_sub.max()))

    # For each s-bin, keep min and max m
    s_bins = np.linspace(s_interval[0], s_interval[1], 50)
    for i in range(len(s_bins)-1):
        lo, hi = s_bins[i], s_bins[i+1]
        mask = (s_list >= lo) & (s_list < hi)
        if np.any(mask):
            m_sub = m_list[mask]
            boundary.append((m_sub.min(), s_list[mask][np.argmin(m_sub)]))
            boundary.append((m_sub.max(), s_list[mask][np.argmax(m_sub)]))

    boundary = np.array(boundary)

    # --------------------------------------------------------
    # Quantile envelopes on your canonical grids
    # --------------------------------------------------------
    pL = PbO.ii()   # left grid
    pU = PbO.jjj()  # right grid

    QL = np.full_like(pL, np.inf, dtype=float)
    QU = np.full_like(pU, -np.inf, dtype=float)

    for (m, s) in boundary:
        QL = np.minimum(QL, qlognormal(pL, m, s))
        QU = np.maximum(QU, qlognormal(pU, m, s))

    return QL, QU, m_interval, s_interval, boundary









# PR-MODE QUANTILE-FITTING

from scipy.integrate import quad
from scipy.optimize import root

# ============================================================
# Mean functions for each family
# ============================================================

def mean_normal(theta):
    mu, sigma = theta
    return mu

def mean_lognormal(theta):
    m, s = theta
    # by construction, m is the mean
    return m

def mean_weibull(theta):
    k, lam = theta
    return sps.weibull_min.mean(c=k, scale=lam)

def mean_gamma(theta):
    shape, scale = theta
    return shape * scale

def mean_logistic(theta):
    mu, s = theta
    return mu

# ============================================================
# Transforms (constrained -> unconstrained) for root-finding
# ============================================================

def normal_transform():
    def to_u(theta):
        mu, sigma = theta
        return np.array([mu, np.log(sigma)])
    def from_u(u):
        mu, log_sigma = u
        return (mu, np.exp(log_sigma))
    return to_u, from_u

def lognormal_transform():
    def to_u(theta):
        m, s = theta
        return np.log([m, s])
    def from_u(u):
        m, s = np.exp(u)
        return (m, s)
    return to_u, from_u

def weibull_transform():
    def to_u(theta):
        k, lam = theta
        return np.log([k, lam])
    def from_u(u):
        k, lam = np.exp(u)
        return (k, lam)
    return to_u, from_u

def gamma_transform():
    def to_u(theta):
        a, t = theta
        return np.log([a, t])
    def from_u(u):
        a, t = np.exp(u)
        return (a, t)
    return to_u, from_u

def logistic_transform():
    def to_u(theta):
        mu, s = theta
        return np.array([mu, np.log(s)])
    def from_u(u):
        mu, log_s = u
        return (mu, np.exp(log_s))
    return to_u, from_u

# ============================================================
# Core: range-mode Fermi for a given family
# ============================================================
def fermi_range_family(ppf, mean_fun, transform, lower, upper, n,
                       theta0=None):
    '''
    Option A: range mode.
    Inputs:
      - ppf: quantile function ppf(p, *theta)
      - mean_fun(theta): E_theta[X]
      - transform: (to_u, from_u) for constrained params
      - [lower, upper]: sample min/max from size n
      - n: sample size
    Solves:
      mean_fun(theta)      = (lower + upper)/2
      mean_range_ppf(ppf, theta, n) = upper - lower
    '''
    m_target = (lower + upper) / 2.0
    R_target = upper - lower

    to_u, from_u = transform

    # crude initial guess: center at m_target, scale from naive normal trick
    if theta0 is None:
        # use normal-based range heuristic for scale-ish start
        if R_target > 0 and n > 1:
            mr_norm = mean_range_ppf(lambda u, mu, sigma: qnormal(u, 0.0, 1.0),
                                     (0.0, 1.0), n)
            s0 = R_target / mr_norm
        else:
            s0 = R_target if R_target > 0 else 1.0
        theta0 = (m_target, s0)

    u0 = to_u(theta0)

    def residual(u):
        theta = from_u(u)
        m = mean_fun(theta)
        R = mean_range_ppf(ppf, theta, n)
        return np.array([m - m_target, R - R_target])

    sol = root(residual, u0)
    if not sol.success:
        raise RuntimeError("Fermi range inversion failed: " + sol.message)

    return from_u(sol.x)

# ============================================================
# Unified Fermi (Option A only: n required, pr ignored)
# ============================================================

def fermi(shape, lower, upper, n):
    '''
    Range-mode Fermi elicitation (Option A):
      - shape is a constructor: normal, lognormal, weibull, gamma, logistic
      - [lower, upper] are sample min/max from a sample of size n
      - pr is ignored; we match mean and expected sample range.
    '''

    if n is None or n <= 1:
        raise ValueError("Range-mode Fermi requires n > 1")

    if shape is normal:
        # closed form: mean-range trick
        mr = mean_range_ppf(lambda u, mu, sigma: qnormal(u, 0.0, 1.0),
                            (0.0, 1.0), n)
        mu = (lower + upper) / 2.0
        sigma = (upper - lower) / mr
        return normal(mu, sigma)

    elif shape is lognormal:
        theta = fermi_range_family(
            ppf=lambda p, m, s: qlognormal(p, m, s),
            mean_fun=mean_lognormal,
            transform=lognormal_transform(),
            lower=lower,
            upper=upper,
            n=n
        )
        return lognormal(*theta)

    elif shape is weibull:
        theta = fermi_range_family(
            ppf=lambda p, k, lam: qweibull(p, k, lam),
            mean_fun=mean_weibull,
            transform=weibull_transform(),
            lower=lower,
            upper=upper,
            n=n
        )
        return weibull(*theta)

    elif shape is gamma:
        theta = fermi_range_family(
            ppf=lambda p, a, t: qgamma(p, a, t),
            mean_fun=mean_gamma,
            transform=gamma_transform(),
            lower=lower,
            upper=upper,
            n=n
        )
        return gamma(*theta)

    elif shape is logistic:
        theta = fermi_range_family(
            ppf=lambda p, mu, s: qlogistic(p, mu, s),
            mean_fun=mean_logistic,
            transform=logistic_transform(),
            lower=lower,
            upper=upper,
            n=n
        )
        return logistic(*theta)

    else:
        raise ValueError(f"Unsupported distribution constructor: {shape}")








def fermi_from_quantiles(ppf, x1, x2, pr=0.9, n=None, theta0=None, transform=None):
    p1 = (1 - pr) / 2
    p2 = (1 + pr) / 2
    if transform is None: # default transform
        def to_u(theta): return np.array(theta, float)
        def from_u(u): return tuple(u)
    else: to_u, from_u = transform
    if theta0 is None: # provisional initial guess in parameter space
        m = (x1 + x2) / 2
        s = (x2 - x1) / (2 * sps.norm.ppf((1+pr)/2))
        theta0 = (m, s)
    if n is not None: # distribution-specific range correction using theta0
        R = (x2 - x1) / mean_range(lambda u: ppf(u, *theta0), n)
        m = (x1 + x2) / 2
        x1_eff = m - R/2
        x2_eff = m + R/2
    else: x1_eff, x2_eff = x1, x2
    u0 = to_u(theta0)
    def residual(u):
        theta = from_u(u)
        q1 = ppf(p1, *theta)
        q2 = ppf(p2, *theta)
        return np.array([q1 - x1_eff, q2 - x2_eff])
    sol = root(residual, u0)
    if not sol.success: raise RuntimeError("Fermi quantile inversion failed: " + sol.message)
    return from_u(sol.x)



# def fermi(shape, lower, upper, n=None, pr=0.9):
#     if   shape is normal:    ppf, transform = qnormal, normal_transform()
#     elif shape is lognormal: ppf, transform = qlognormal, lognormal_transform()
#     elif shape is weibull:   ppf, transform = qweibull, weibull_transform()
#     elif shape is gamma:     ppf, transform = qgamma, gamma_transform()
#     elif shape is logistic:  ppf, transform = qlogistic, logistic_transform()
#     else: raise ValueError(f"Unsupported distribution constructor: {shape}")
#     theta = fermi_from_quantiles(ppf,lower,upper,pr=pr,n=n,theta0=None,transform=transform)
#     print('theta',theta)
#     return shape(*theta)  # build distribution using the constructor




def fermiF(f,a,b,n=None): #,pr=0.9):
    
    fab = fermi(f,a,b,n) #,pr=pr)
    plot(fab)
    #plt.title(str(n))
    plt.plot([a,a,b,b],[1,0,0,1],'xkcd:grey',ls=':')
    A,B = ends(support(fab))
    
    pr = 0.9
    
    alpha = (1-pr)/2
    plt.plot([A,B,B,A],[1-alpha,1-alpha,alpha,alpha],'r',ls=':',lw=1)
    #plt.show()
  
a = N((52+86)/2, (86-52)/3.077) 
b = N(69, 11.06)
c = fermi(normal,52,86,10)     
red(a); blue(b); cyan(c)
    

fermiF(normal, 52, 86, n=10)   # ≈ N(69, 11.06)


# Lognormal, range-mode:
fermiF(lognormal, 1.0, 5.0, n=20)

# Gamma, range-mode:
fermiF(gamma, 0.2, 3.0, n=15)



# memoizable mean range curve
# ppf = qnormal
# nn = [n for n in range(2,100)] + [125, 150, 200, 250, 300, 400, 600, 900, 1000]
# mr = []
# for n in nn : mr.append( mean_range(ppf, n) )
# plt.plot(nn,mr,'-o')

fermiF(normal, a=0.2, b=3) 
fermiF(lognormal, a=0.2, b=3) 
fermiF(weibull, a=0.2, b=3) 
fermiF(gamma, a=0.2, b=3) 
fermiF(logistic, a=0.2, b=3) 
fermiF(normal, 0.2, 3, pr=0.8)   
"""
    














"""
# EXTENDED FERMI METHODS



import numpy as np
import scipy.stats as sps
from math import sqrt, log
from scipy.integrate import quad

# --------------------------------------------------------------------
# Helpers you already have
# --------------------------------------------------------------------

def mean_normal_range(n):
    '''Expected range of n independent N(0,1) samples.'''
    phi = sps.norm.pdf
    Phi = sps.norm.cdf
    def log_f(x): return np.log(n) + np.log(phi(x)) + (n-1)*np.log(Phi(x))
    xs = np.linspace(-8, 8, 2001)
    m = np.max(log_f(xs))
    def integrand(x): return x * np.exp(log_f(x) - m)
    I = quad(integrand, -8, 8, limit=200)[0]
    return 2 * I * np.exp(m)

# lognormal PPF in mean/sd parameterization
def qlognormal(p, m, s):
    sigma = np.sqrt(np.log(1 + (s*s)/(m*m)))
    mu    = np.log(m) - 0.5*sigma*sigma
    return sps.lognorm.ppf(p, s=sigma, scale=np.exp(mu))

def mean_range_lognormal(m, s, n, sims=100000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = qlognormal(u, m, s)
    return np.mean(x.max(axis=1) - x.min(axis=1))

# --------------------------------------------------------------------
# Generic 2-parameter Fermi p-box engine (rigorous boundary version)
# --------------------------------------------------------------------

def fermi_pbox_generic_2d(lower, upper, n,
                          ppf,          # ppf(p, t1, t2)
                          mean_fun,     # mean(t1, t2)
                          mean_range_fun,  # mean_range(t1, t2, n)
                          t1_vals, t2_vals,
                          rigor=True,
                          tol_mean=0.05,
                          tol_range=0.10,
                          sims_range=None):
    '''
    Generic rigorous Fermi p-box constructor for 2-parameter families.
    Returns a Pbox(QL, QU, m_min, m_max), where m_* are mean bounds.
    '''

    m_target = (lower + upper) / 2.0
    R_target = upper - lower

    admissible = []
    exterior   = []

    # classify grid points
    for t1 in t1_vals:
        for t2 in t2_vals:
            theta = (t1, t2)
            m = mean_fun(*theta)
            mean_ok = abs(m - m_target)/m_target <= tol_mean

            R_hat = mean_range_fun(*theta, n=n) if sims_range is None \
                    else mean_range_fun(*theta, n=n, sims=sims_range)
            range_ok = abs(R_hat - R_target)/R_target <= tol_range

            if mean_ok and range_ok:
                admissible.append(theta)
            else:
                exterior.append(theta)

    if not admissible:
        raise RuntimeError("no admissible parameters found")

    admissible = np.array(admissible)
    exterior_set = set((round(t1,12), round(t2,12)) for (t1,t2) in exterior)

    # rigorous boundary: admissible + first exterior neighbors
    boundary = set()
    t1_step = t1_vals[1] - t1_vals[0] if len(t1_vals) > 1 else 0.0
    t2_step = t2_vals[1] - t2_vals[0] if len(t2_vals) > 1 else 0.0

    for (t1, t2) in admissible:
        boundary.add((t1, t2))
        if rigor:
            for dt1 in [-t1_step, 0, t1_step]:
                for dt2 in [-t2_step, 0, t2_step]:
                    if dt1 == 0 and dt2 == 0:
                        continue
                    t1b = t1 + dt1
                    t2b = t2 + dt2
                    key = (round(t1b,12), round(t2b,12))
                    if key in exterior_set:
                        boundary.add((t1b, t2b))

    boundary = np.array(list(boundary))
    t1_list = boundary[:,0]
    t2_list = boundary[:,1]

    # mean interval (epistemic)
    means = np.array([mean_fun(t1, t2) for (t1, t2) in boundary])
    m_interval = (means.min(), means.max())

    # quantile envelopes on your canonical grids
    pL = PbO.ii()
    pU = PbO.jjj()
    QL = np.full_like(pL, np.inf, dtype=float)
    QU = np.full_like(pU, -np.inf, dtype=float)

    for (t1, t2) in boundary:
        QL = np.minimum(QL, ppf(pL, t1, t2))
        QU = np.maximum(QU, ppf(pU, t1, t2))

    return Pbox(QL, QU, *m_interval)

# --------------------------------------------------------------------
# 1. Normal (p-boxified)
# --------------------------------------------------------------------

def ppf_normal(p, mu, sigma):
    return sps.norm(loc=mu, scale=sigma).ppf(p)

def mean_normal(theta1, theta2):
    mu, sigma = theta1, theta2
    return mu

def mean_range_normal(mu, sigma, n, sims=None):
    # exact: sigma * E[range of N(0,1)]
    return sigma * mean_normal_range(n)

def fermi_normal_pbox(lower, upper, n,
                      rigor=True, tol_mean=0.05, tol_range=0.10,
                      mu_grid=40, s_grid=40):
    m_target = (lower + upper)/2
    R_target = upper - lower
    mr = mean_normal_range(n)
    s0 = R_target / mr if mr > 0 else max(R_target, 1.0)
    mu_vals = np.linspace(0.5*m_target, 1.5*m_target, mu_grid)
    s_vals  = np.linspace(0.3*s0, 3.0*s0, s_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_normal,
                                 mean_normal,
                                 mean_range_normal,
                                 mu_vals, s_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range)

# original crisp Fermi normal (for reference)
def fermi_normal_crisp(lower, upper, n):
    m = (lower + upper)/2
    R = upper - lower
    sigma = R / mean_normal_range(n)
    return normal(m, sigma)

# --------------------------------------------------------------------
# 2. Lognormal (your version, generalized through the engine)
# --------------------------------------------------------------------

def ppf_lognormal(p, m, s):
    return qlognormal(p, m, s)

def mean_lognormal(m, s):
    return m

def mean_range_lognormal_theta(m, s, n, sims=100000):
    return mean_range_lognormal(m, s, n, sims=sims)

def fermi_lognormal(lower, upper, n,
                    rigor=True, tol_mean=0.05, tol_range=0.10,
                    m_grid=40, s_grid=40, sims_range=50000):
    m_target = (lower + upper)/2
    R_target = upper - lower
    if n > 1 and R_target > 0:
        mr_norm = 3.077
        s0 = R_target / mr_norm
    else:
        s0 = max(R_target, 1.0)
    m_vals = np.linspace(0.5*m_target, 1.5*m_target, m_grid)
    s_vals = np.linspace(0.3*s0, 3.0*s0, s_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_lognormal,
                                 mean_lognormal,
                                 mean_range_lognormal_theta,
                                 m_vals, s_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range,
                                 sims_range=sims_range)

# --------------------------------------------------------------------
# 3. Uniform (true scale)
# --------------------------------------------------------------------

def ppf_uniform(p, m, h):
    return m + h*(2*p - 1)

def mean_uniform(m, h):
    return m

def mean_range_uniform(m, h, n, sims=None):
    return 2*h*(n-1)/(n+1)

def fermi_uniform(lower, upper, n,
                  rigor=True, tol_mean=0.05, tol_range=0.10,
                  m_grid=40, h_grid=40):
    m_target = (lower + upper)/2
    R_target = upper - lower
    h0 = R_target*(n+1)/(2*(n-1)) if n > 1 else R_target/2
    m_vals = np.linspace(0.5*m_target, 1.5*m_target, m_grid)
    h_vals = np.linspace(0.3*h0, 3.0*h0, h_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_uniform,
                                 mean_uniform,
                                 mean_range_uniform,
                                 m_vals, h_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range)

# --------------------------------------------------------------------
# 4. Laplace (true scale)
# --------------------------------------------------------------------

def ppf_laplace(p, mu, b):
    return sps.laplace(loc=mu, scale=b).ppf(p)

def mean_laplace(mu, b):
    return mu

def mean_range_laplace(mu, b, n, sims=None):
    H = np.sum(1/np.arange(1, n)) if n > 1 else 0.0
    return 2*b*H

def fermi_laplace(lower, upper, n,
                  rigor=True, tol_mean=0.05, tol_range=0.10,
                  mu_grid=40, b_grid=40):
    m_target = (lower + upper)/2
    R_target = upper - lower
    H = np.sum(1/np.arange(1, n)) if n > 1 else 1.0
    b0 = R_target/(2*H)
    mu_vals = np.linspace(0.5*m_target, 1.5*m_target, mu_grid)
    b_vals  = np.linspace(0.3*b0, 3.0*b0, b_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_laplace,
                                 mean_laplace,
                                 mean_range_laplace,
                                 mu_vals, b_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range)

# --------------------------------------------------------------------
# 5. Logistic (true scale)
# --------------------------------------------------------------------

def ppf_logistic(p, mu, s):
    return sps.logistic(loc=mu, scale=s).ppf(p)

def mean_logistic(mu, s):
    return mu

def mean_range_logistic(mu, s, n, sims=50000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = ppf_logistic(u, mu, s)
    return np.mean(x.max(axis=1) - x.min(axis=1))

def fermi_logistic(lower, upper, n,
                   rigor=True, tol_mean=0.05, tol_range=0.10,
                   mu_grid=40, s_grid=40, sims_range=50000):
    m_target = (lower + upper)/2
    R_target = upper - lower
    s0 = R_target/4.0 if n > 1 else max(R_target, 1.0)
    mu_vals = np.linspace(0.5*m_target, 1.5*m_target, mu_grid)
    s_vals  = np.linspace(0.3*s0, 3.0*s0, s_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_logistic,
                                 mean_logistic,
                                 mean_range_logistic,
                                 mu_vals, s_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range,
                                 sims_range=sims_range)

# --------------------------------------------------------------------
# 6. Weibull (shape family)
# --------------------------------------------------------------------

def ppf_weibull(p, k, lam):
    return sps.weibull_min(c=k, scale=lam).ppf(p)

def mean_weibull(k, lam):
    return sps.weibull_min(c=k, scale=lam).mean()

def mean_range_weibull(k, lam, n, sims=50000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = ppf_weibull(u, k, lam)
    return np.mean(x.max(axis=1) - x.min(axis=1))

def fermi_weibull(lower, upper, n,
                  rigor=True, tol_mean=0.05, tol_range=0.10,
                  k_grid=40, lam_grid=40, sims_range=50000):
    R_target = upper - lower
    k_vals   = np.linspace(0.3, 5.0, k_grid)
    lam_vals = np.linspace(0.1*R_target, 5.0*R_target, lam_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_weibull,
                                 mean_weibull,
                                 mean_range_weibull,
                                 k_vals, lam_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range,
                                 sims_range=sims_range)

# --------------------------------------------------------------------
# 7. Gamma (shape family)
# --------------------------------------------------------------------

def ppf_gamma(p, a, t):
    return sps.gamma(a, scale=t).ppf(p)

def mean_gamma(a, t):
    return a*t

def mean_range_gamma(a, t, n, sims=50000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = ppf_gamma(u, a, t)
    return np.mean(x.max(axis=1) - x.min(axis=1))

def fermi_gamma(lower, upper, n,
                rigor=True, tol_mean=0.05, tol_range=0.10,
                a_grid=40, t_grid=40, sims_range=50000):
    R_target = upper - lower
    a_vals = np.linspace(0.3, 10.0, a_grid)
    t_vals = np.linspace(0.1*R_target, 5.0*R_target, t_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_gamma,
                                 mean_gamma,
                                 mean_range_gamma,
                                 a_vals, t_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range,
                                 sims_range=sims_range)

# --------------------------------------------------------------------
# 8. Pareto (shape family)
# --------------------------------------------------------------------

def ppf_pareto(p, xm, alpha):
    return sps.pareto(b=alpha, scale=xm).ppf(p)

def mean_pareto(xm, alpha):
    if alpha <= 1:
        return np.inf
    return xm*alpha/(alpha-1)

def mean_range_pareto(xm, alpha, n, sims=50000):
    if n <= 1: return 0.0
    u = np.random.rand(sims, n)
    x = ppf_pareto(u, xm, alpha)
    return np.mean(x.max(axis=1) - x.min(axis=1))

def fermi_pareto(lower, upper, n,
                 rigor=True, tol_mean=0.05, tol_range=0.10,
                 xm_grid=40, a_grid=40, sims_range=50000):
    R_target = upper - lower
    xm_vals = np.linspace(0.1*lower, lower, xm_grid)
    a_vals  = np.linspace(1.1, 10.0, a_grid)
    return fermi_pbox_generic_2d(lower, upper, n,
                                 ppf_pareto,
                                 mean_pareto,
                                 mean_range_pareto,
                                 xm_vals, a_vals,
                                 rigor=rigor,
                                 tol_mean=tol_mean,
                                 tol_range=tol_range,
                                 sims_range=sims_range)


approx.ksD95 <- function(n) {
    # approximations for the critical level for Kolmogorov-Smirnov statistic D,
    # for confidence level 0.95. Taken from Bickel & Doksum, table IX, p.483
    # and Lienert G.A.(1975) who attributes to Miller,L.H.(1956), JASA
    ifelse(n > 80,
           1.358 /( sqrt(n) + .12 + .11/sqrt(n)),##Bickel&Doksum, table IX,p.483

           splinefun(c(1:9, 10, 15, 10 * 2:8),# from Lienert
                     c(.975,   .84189, .70760, .62394, .56328,# 1:5
                       .51926, .48342, .45427, .43001, .40925,# 6:10
                       .33760, .29408, .24170, .21012,# 15,20,30,40
                       .18841, .17231, .15975, .14960)) (n))
  }

# def histogram(x, mn=None, mx=None, conf=0.95):
#     x, y = as_vectors(x)
#     x = np.asarray(x)
#     n = len(x)
#     if mn is None: mn = np.min(x) if y is None else min(np.min(x), np.min(y))
#     if mx is None: mx = np.max(x) if y is None else max(np.max(x), np.max(y))
#     # sort indices
#     sx = np.argsort(x)
#     # canonical grid size
#     steps = PbO.steps
#     # R indexing trick:
#     # u <- x[sx[1+0:(steps-1)/((steps)/n)]]
#     idx_u = sx[np.arange(steps) * n // steps]
#     u = x[idx_u]
#     # d <- x[sx[1+1:steps/((steps)/(n))]]
#     idx_d = sx[np.arange(steps) * n // steps] # this indexiing is a bit better than the R version's
#     d = x[idx_d]
#     #d[-1] = np.max(x)
#     if steps % n == 0: d = u.copy()        # if steps divisible by n, set d = u
#     if np.all(x==y): A = Pbox(u, d)  #, shape='histogram') # WON'T WORK FOR INTERVAL DATA ARRAYS
#     else: A = env(Pbox(u, d), histogram(y))
#     # fatten by Kolmogorov bound
#     dm = KSDmax(n, conf)
#     return fatten(A, dm=dm, leftbound=mn, rightbound=mx)







#             THIS IS THE WOODPILE  (append above this line)
#------------------------------------------------------------------------------
# End of WOODPILE
#------------------------------------------------------------------------------
          
"""

