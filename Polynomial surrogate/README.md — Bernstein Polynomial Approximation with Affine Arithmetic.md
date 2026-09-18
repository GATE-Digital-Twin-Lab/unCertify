# Bernstein Polynomial Approximation with Affine Arithmetic

## Overview

This project implements a **Bernstein polynomial approximation framework combined with Affine Arithmetic (AA)** for rigorous range bounding of nonlinear functions.

The main objective is to obtain tighter bounds for functions for which direct Interval Arithmetic (IA) or Affine Arithmetic can produce overly conservative enclosures because of dependency, wrapping, and repeated-variable effects.

The approach combines:

1. **Domain subintervalization**
2. **Polynomial approximation on each subdomain**
3. **Bernstein coefficient computation**
4. **Affine Arithmetic evaluation of the polynomial**
5. **Rigorous residual bounding**
6. **Global hull computation over all subdomains**

The computational domain used by the current framework is normalised to

$[0,1]^d.$

A separate polynomial approximation is constructed on each subdomain, and the resulting local bounds are combined to obtain a global enclosure.

---

## Important clarification: interpolation nodes vs. polynomial basis

A key point in this implementation is that the parameter `node_type` controls **only the interpolation/sample points used to compute the Bernstein coefficients**. It does **not** change the polynomial basis.

For example, when

```python
node_type = "chebyshev"
```

the approximation is still a **Bernstein polynomial**. Chebyshev-Lobatto points are used only as the interpolation nodes at which the function is sampled before solving for the Bernstein coefficients.

Thus:

| `node_type` | Interpolation nodes | Polynomial basis |
|---|---|---|
| `"bernstein"` | Bernstein nodes | Bernstein |
| `"chebyshev"` | Chebyshev-Lobatto nodes | Bernstein |
| `"legendre"` | Legendre-Lobatto nodes | Bernstein |

In all three cases, the resulting approximation has the form

$\[
p_B(t)=
\sum_{i_1=0}^{n_1}\cdots\sum_{i_d=0}^{n_d}
C_{i_1,\ldots,i_d}
\prod_{k=1}^{d}
B_{i_k,n_k}(t_k),
\]$

where

$\[
B_{k,n}(t) =
\binom{n}{k}t^k(1-t)^{n-k}
\]$

are Bernstein basis functions.

Therefore, `node_type` should be interpreted as a choice of **coefficient-generation/interpolation grid**, rather than a choice of polynomial family.

For example,

```python
nvec = [3, 3]
node_type = "chebyshev"
```

means that the implementation:

1. generates a tensor-product grid of Chebyshev-Lobatto points,
2. evaluates the function at those points,
3. solves the interpolation systems to obtain the Bernstein coefficients,
4. represents the resulting approximation in the Bernstein basis,
5. evaluates the Bernstein polynomial using Affine Arithmetic.

This distinction is particularly important for the **Bernstein convex-hull property**. Regardless of whether the coefficients were obtained using Bernstein, Chebyshev-Lobatto, or Legendre-Lobatto interpolation points, the final polynomial is represented in the Bernstein basis. Consequently,

$\[
\min_i C_i
\leq
p_B(t)
\leq
\max_i C_i,
\qquad t\in[0,1]^d,
\]$

providing a coefficient-based enclosure of the polynomial.

---

# Methodology

The overall procedure is illustrated below:

```text
Original domain [0,1]^d
          │
          ▼
   Subintervalization
          │
          ▼
  Local subdomain X_i
          │
          ▼
   Map X_i → [0,1]^d
          │
          ▼
 Select interpolation nodes
(Bernstein / Chebyshev / Legendre)
          │
          ▼
 Evaluate f on interpolation grid
          │
          ▼
 Compute Bernstein coefficients
          │
          ▼
 Construct Bernstein polynomial p_B
          │
          ├───────────────┐
          ▼               ▼
      AA evaluation   Coefficient
       of p_B          convex hull
          │               │
          ▼               ▼
       p_AA          p_B interval
          │               │
          └───────┬───────┘
                  ▼
           Residual bound
             r = f - p_B
                  │
                  ▼
       Hull over all subdomains
                  │
                  ▼
        Global rigorous bound
```

---

# 1. Domain subintervalization

The normalised domain is

$\[
X=[0,1]^d.
\]$

The domain can be divided into a user-specified number of subintervals in each dimension.

For example,

```python
split_interval([4, 4])
```

divides a two-dimensional domain into

$\[
4\times4=16
\]$

subdomains.

Each subdomain has the form

$\[
X_i =
[a_1,b_1]\times\cdots\times[a_d,b_d].
\]$

The polynomial approximation is constructed independently on each subdomain.

Subintervalization is useful because it reduces the size of the region over which the approximation and range evaluation are performed. It can also reduce dependency and wrapping effects in the Affine Arithmetic evaluation.

The global enclosure is subsequently obtained by taking the hull of the local enclosures.

---

# 2. Mapping each subdomain to $`[0,1]^d`$

For a subdomain

$\[
X_i=
\prod_{k=1}^{d}[a_k,b_k],
\]$

the physical coordinates are mapped to normalised coordinates using

$\[
t_k =
\frac{x_k-a_k}{b_k-a_k}.
\]$

Therefore,

$\[
t_k\in[0,1].
\]$

The Bernstein polynomial is constructed in this normalised coordinate system.

The inverse mapping is

$\[
x_k=a_k+(b_k-a_k)t_k.
\]$

This mapping allows the same Bernstein basis definition to be used on every subdomain.

---

# 3. Bernstein basis

The one-dimensional Bernstein basis of degree \(n\) is

$\[
B_{k,n}(t) =
\binom{n}{k}
t^k(1-t)^{n-k},
\qquad
k=0,\ldots,n.
\]$

The basis satisfies two important properties:

### Non-negativity

For

$\[
t\in[0,1],
\]$

we have

$\[
B_{k,n}(t)\geq0.
\]$

### Partition of unity

The basis functions satisfy

$\[
\sum_{k=0}^{n}B_{k,n}(t)=1.
\]$

These properties lead directly to the Bernstein convex-hull property used for polynomial range bounding.

---

# 4. Bernstein basis recurrence

The implementation does not need to explicitly construct each basis function using the expression

$\[
\binom{n}{k}t^k(1-t)^{n-k}.
\]$

Instead, `AA_evaluationBr.py` constructs the basis recursively.

The recurrence starts with

$\[
B_{0,0}(t)=1.
\]$

For increasing degree,

$\[
B_{0,m}(t) =
(1-t)B_{0,m-1}(t),
\]$

$\[
B_{k,m}(t) =
(1-t)B_{k,m-1}(t)
+
tB_{k-1,m-1}(t),
\]$

and

$\[
B_{m,m}(t) =
tB_{m-1,m-1}(t).
\]$

This recurrence is particularly convenient for Affine Arithmetic because the basis functions can be constructed using the overloaded AA operations.

---

# 5. Tensor-product Bernstein polynomial

For a $\(d\)-$ dimensional problem with degree vector

$\[
\mathbf n=(n_1,\ldots,n_d),
\]$

the tensor-product Bernstein approximation is

$\[
p_B(t) =
\sum_{i_1=0}^{n_1}
\cdots
\sum_{i_d=0}^{n_d}
C_{i_1,\ldots,i_d}
\prod_{k=1}^{d}
B_{i_k,n_k}(t_k).
\]$

Here,

$\[
C_{i_1,\ldots,i_d}
\]$

are the Bernstein coefficients.

For example, with

```python
nvec = [3, 3]
```

the polynomial has

$\[
(3+1)(3+1)=16
\]$

Bernstein coefficients.

---

# 6. Interpolation nodes

The code supports three choices of interpolation nodes:

### Bernstein nodes

For degree $\(n\)$,

$\[
t_k=\frac{k}{n},
\qquad
k=0,\ldots,n.
\]$

### Chebyshev-Lobatto nodes

The Chebyshev-Lobatto nodes are generated on the requested physical interval and are then mapped to the normalized Bernstein domain.

### Legendre-Lobatto nodes

The Legendre-Lobatto nodes are generated using the Legendre polynomial formulation and similarly mapped to the normalized domain.

The available option is selected through

```python
node_type
```

with

```python
node_type = "bernstein"
```

```python
node_type = "chebyshev"
```

or

```python
node_type = "legendre"
```

Again, this parameter determines the **interpolation nodes**, not the polynomial basis.

---

# 7. Computing the Bernstein coefficients

The function values are first evaluated on the tensor-product interpolation grid.

For each dimension, an interpolation matrix is constructed from the Bernstein basis:

$\[
A_{jk} =
B_{k,n}(t_j).
\]$

The corresponding linear system is solved to obtain the Bernstein coefficients.

In one dimension,

$\[
A C = F,
\]$

where

- \(A\) is the Bernstein interpolation matrix,
- \(C\) is the vector of Bernstein coefficients,
- \(F\) contains the function values at the interpolation nodes.

For multiple dimensions, the interpolation problem is solved dimension-by-dimension using tensor-product structure.

The implementation uses

```python
np.linalg.solve
```

rather than explicitly computing matrix inverses.

---

# 8. Affine Arithmetic evaluation

Once the Bernstein coefficients have been computed, the polynomial can be evaluated using the Affine Arithmetic implementation.

The physical subdomain is first mapped to normalized coordinates:

$\[
t_k =
\frac{x_k-a_k}{b_k-a_k}.
\]$

The Bernstein basis functions are then constructed using AA arithmetic.

The polynomial is evaluated as

$\[
p_B(t) =
\sum_i C_i B_i(t).
\]$

Since the coefficients \(C_i\) are numerical values and the basis functions \(B_i(t)\) are represented using Affine Arithmetic, the resulting polynomial evaluation is itself an affine form.

The resulting enclosure is stored as an `AffineScalar`.

---

# 9. Bernstein convex-hull property

One of the main advantages of the Bernstein representation is its convex-hull property.

Because

$\[
B_i(t)\geq0
\]$

and

$\[
\sum_iB_i(t)=1,
\]$

the polynomial

$\[
p_B(t)=\sum_i C_iB_i(t)
\]$

is a convex combination of its Bernstein coefficients.

Therefore,

$\[
\min_i C_i
\leq
p_B(t)
\leq
\max_i C_i.
\]$

Consequently, if the minimum and maximum Bernstein coefficients are

$\[
C_{\min} =
\min_i C_i
\]$

and

$\[
C_{\max} =
\max_i C_i,
\]$

then

$\[
p_B(X_i)
\subseteq
[C_{\min},C_{\max}].
\]$

This provides a simple coefficient-based enclosure of the polynomial without requiring direct interval evaluation of the polynomial expression.

---

# 10. Residual computation

The polynomial surrogate is used to separate the original function into

$\[
f(x)=p_B(x)+r(x),
\]$

where

$\[
r(x)=f(x)-p_B(x)
\]$

is the approximation residual.

The residual is important because the final rigorous range bound can be constructed from the polynomial approximation and the residual enclosure.

Two residual calculations are implemented.

---

## 10.1 Direct Affine Arithmetic residual

The polynomial is evaluated using Affine Arithmetic:

$\[
p_{\mathrm{AA}}(X_i).
\]$

The residual is then computed directly as

$\[
r_{\mathrm{AA}}(X_i) =
f_{\mathrm{AA}}(X_i) -
p_{\mathrm{AA}}(X_i).
\]$

This approach preserves the dependency information available through the AA representation.

The implementation returns this quantity as

```python
residual_AA
```

---

## 10.2 Bernstein coefficient-based residual

A second residual enclosure uses the Bernstein convex-hull property.

Suppose the function has the AA enclosure

$\[
f(X_i)
\subseteq
[f_i^-,f_i^+]
\]$

and the Bernstein coefficients provide

$\[
p_B(X_i)
\subseteq
[p_i^-,p_i^+].
\]$

Then

$\[
r(X_i) =
f(X_i)-p_B(X_i)
\]$

is bounded by

$\[
r(X_i)
\subseteq
\left[
f_i^- - p_i^+,
\;
f_i^+ - p_i^-
\right].
\]$

In the implementation,

$\[
p_i^-=\min_j C_j,
\qquad
p_i^+=\max_j C_j.
\]$

Therefore,

$\[
r_B(X_i)
\subseteq
\left[
f_i^- - C_{\max},
\;
f_i^+ - C_{\min}
\right].
\]$

This is stored as

```python
residual_Bernstein
```

and provides a residual enclosure based on the Bernstein coefficient range.

---

# 11. Global bound through subintervalization

For each subdomain \(X_i\), the method produces a local residual enclosure

$\[
r(X_i)
\subseteq
[r_i^-,r_i^+].
\]$

The global residual enclosure is obtained by taking the hull of all local intervals:

$\[
r(X)
\subseteq
\mathrm{hull}
\left(
\bigcup_i r(X_i)
\right).
\]$

Similarly, global hulls can be computed for

- \(f(X)\),
- \(p_B(X)\),
- \(p_{\mathrm{AA}}(X)\),
- \(r_{\mathrm{AA}}(X)\),
- \(r_B(X)\).

The function

```python
hull_intervals(...)
```

performs this global aggregation.

---

# 12. Why subintervalization is used

Applying the polynomial approximation over the entire domain may produce a relatively large approximation error.

Instead, the domain is divided into smaller regions and a separate polynomial is constructed on each region.

For example,

```text
1 × 1
2 × 2
4 × 4
8 × 8
16 × 16
...
```

Each subdomain has a local polynomial approximation.

As the subdomains become smaller, the polynomial only needs to represent the function over a smaller region. This can reduce approximation error and may also reduce overestimation associated with nonlinear AA operations.

The global result is then obtained by taking the hull over all local results.

The framework does not assume that every individual interval width must decrease monotonically with increasing subdivision. The actual behavior depends on the function, polynomial degree, interpolation nodes, and AA evaluation.

---

# 13. Project structure

The main files are organized as follows:

```text
Bernstein/
│
├── BernsteinFiles.py
├── AA_evaluationBr.py
├── test_functions.py
├── subintervalization.py
└── run_Br_poly.py
```

### `BernsteinFiles.py`

Contains the numerical Bernstein approximation framework:

- interpolation node generation,
- tensor-product grid generation,
- function evaluation on the grid,
- Bernstein basis functions,
- Bernstein interpolation matrices,
- tensor-product Bernstein coefficient computation,
- Bernstein polynomial evaluation.

Main functions include:

```python
get_nodes()
bernstein_grid_nd()
evaluate_function_on_grid()
reshape_grid_values()
bernstein_basis()
bernstein_matrix()
bernstein_coefficients()
evaluate_bernstein_1d()
evaluate_bernstein_nd()
bernstein_approximation()
```

---

### `AA_evaluationBr.py`

Contains the Affine Arithmetic implementation of the Bernstein polynomial evaluation.

Main functions include:

```python
polylist_Bernstein()
map_to_bernstein()
bernstein_basis_AA()
evaluate_bernstein_AA()
residual_eval_Bernstein_AA()
```

This module connects the numerical Bernstein approximation with the Affine Arithmetic framework.

---

### `test_functions.py`

Contains the benchmark functions in both numerical and Affine Arithmetic forms.

Current benchmark functions include:

- Branin
- Ackley
- Eggholder

For example:

```python
Branin_numpy()
test_BraninAA()

Ackley_numpy()
test_AckleyAA()

Eggholder_numpy()
test_EggholderAA()
```

---

### `subintervalization.py`

Provides utilities for domain subdivision and interval aggregation.

Main functions include:

```python
split_interval()
interval_width()
hull_intervals()
```

---

### `run_Br_poly.py`

Provides the experiment driver.

It controls:

- benchmark function,
- polynomial degree,
- interpolation node type,
- Affine Arithmetic options,
- subdivision levels,
- result reporting.

---

# 14. Example configuration

A typical experiment can be configured as follows:

```python
from functools import partial

from test_functions import (
    Ackley_numpy,
    test_AckleyAA
)

fun_numpy = partial(
    Ackley_numpy,
    d=2
)

fun_AA = partial(
    test_AckleyAA,
    d=2
)

nvec = [3, 3]

node_type = "chebyshev"

cheb = False

split_values = [
    1,
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256
]
```

Here,

```python
nvec = [3, 3]
```

specifies a degree-3 Bernstein polynomial in each dimension.

The choice

```python
node_type = "chebyshev"
```

means that Chebyshev-Lobatto points are used for interpolation, while the resulting polynomial is still represented in the Bernstein basis.

---

# 15. Running the experiments

The experiment driver evaluates the approximation for different levels of subintervalization.

For example:

```python
run_experiment(
    fun_numpy,
    fun_AA,
    nvec,
    node_type,
    cheb,
    split_values
)
```

The results can then be reported for each subdivision level.

Typical quantities include:

```text
Global f hull
Global Bernstein coefficient hull
Global p_AA hull
Global residual AA hull
Global residual Bernstein hull
```

as well as their corresponding interval widths.

---

# 16. Interpretation of the reported intervals

The framework reports several different enclosures because they represent different stages of the approximation and bounding process.

### Function hull

$\[
f(X)
\]$

represents the global Affine Arithmetic enclosure of the original function.

### Bernstein coefficient hull

$\[
p_{\mathrm{Bern}}(X)
\]$

is obtained from the minimum and maximum Bernstein coefficients.

### AA polynomial hull

$\[
p_{\mathrm{AA}}(X)
\]$

is obtained by directly evaluating the Bernstein polynomial using Affine Arithmetic.

### Direct residual AA hull

$\[
r_{\mathrm{AA}}(X) =
f_{\mathrm{AA}}(X)-p_{\mathrm{AA}}(X).
\]$

### Bernstein residual hull

$\[
r_{\mathrm{Bern}}(X)
\]$

is obtained by combining the AA enclosure of the function with the Bernstein coefficient-based enclosure of the polynomial.

These bounds are intentionally reported separately so that the effect of the polynomial representation and the AA evaluation can be examined independently.

---

# 17. Benchmark functions

The current framework supports several benchmark functions.

## Branin

The normalised input variables are mapped to the standard Branin domain before evaluating the function.

The implementation contains both a NumPy version and an Affine Arithmetic version.

---

## Ackley

The Ackley function is implemented for configurable dimensions:

```python
Ackley_numpy(x, d=2)
```

and

```python
test_AckleyAA(x, d=2, cheb=False)
```

The normalised variables in \([0,1]^d\) are mapped to the standard Ackley input range.

---

## Eggholder

The Eggholder function is also available in numerical and Affine Arithmetic forms.

---

# 18. Dependencies

The implementation requires:

- Python 3
- NumPy
- SciPy

The main numerical functionality uses:

```python
import numpy as np
```

and SciPy routines such as:

```python
from scipy.special import comb
from scipy.special import roots_legendre
```

The Affine Arithmetic implementation additionally requires the project-specific

```python
Affine_ArithmeticClassV3.py
```

module.

---

# 19. Role of Affine Arithmetic

Affine Arithmetic is used to obtain conservative enclosures while maintaining correlation information between repeated occurrences of uncertain variables.

The polynomial basis functions are constructed directly using Affine Arithmetic operations.

For example, operations such as

```python
(1 - x) * basis[k]
```

and

```python
x * basis[k]
```

are evaluated using the overloaded AA arithmetic.

This allows the Bernstein polynomial to be evaluated as an affine form rather than treating every occurrence of a variable as an independent interval quantity.

The resulting AA enclosure can therefore be compared with the simpler coefficient-based Bernstein enclosure.

---

# 20. Role of the Bernstein representation

The Bernstein representation is particularly useful for range bounding because of its basis properties.

For

$\[
t\in[0,1],
\]$

the basis functions satisfy

$\[
B_{k,n}(t)\geq0
\]$

and

$\[
\sum_{k=0}^{n}B_{k,n}(t)=1.
\]$

Consequently, the polynomial value is always contained in the convex hull of its Bernstein coefficients:

$\[
p_B(t)
\in
\left[
\min_k C_k,
\max_k C_k
\right].
\]$

This property allows the polynomial approximation to be bounded directly from its coefficients.

The combination of this property with subintervalization provides a convenient framework for constructing local polynomial surrogates and bounding their approximation residuals.

---

# 21. Overall algorithm

For each subdomain \(X_i\):

### Step 1: Define the local domain

$\[
X_i=
\prod_{k=1}^{d}[a_k,b_k].
\]$

### Step 2: Generate interpolation nodes

Select one of:

- Bernstein nodes,
- Chebyshev-Lobatto nodes,
- Legendre-Lobatto nodes.

### Step 3: Evaluate the function

Compute

$\[
F_j=f(x_j)
\]$

at the tensor-product interpolation points.

### Step 4: Compute Bernstein coefficients

Solve the interpolation systems to obtain

$\[
C_{i_1,\ldots,i_d}.
\]$

### Step 5: Construct the Bernstein polynomial

$\[
p_B(t) =
\sum_i C_iB_i(t).
\]$

### Step 6: Evaluate using Affine Arithmetic

Construct

$\[
p_{\mathrm{AA}}(X_i).
\]$

### Step 7: Obtain a Bernstein coefficient enclosure

Use

$\[
p_B(X_i)
\subseteq
[C_{\min},C_{\max}].
\]$

### Step 8: Compute the residual enclosure

Using

$\[
r=f-p_B,
\]$

compute

$\[
r(X_i)
\subseteq
\left[
f_i^- - C_{\max},
f_i^+ - C_{\min}
\right].
\]$

### Step 9: Aggregate all subdomains

Finally,

$\[
r(X)
\subseteq
\mathrm{hull}
\left(
\bigcup_i r(X_i)
\right).
\]$

The same procedure can be used to obtain global enclosures for the function and polynomial.

---

# 22. Summary

The implemented framework combines **Bernstein polynomial approximation**, **Affine Arithmetic**, and **domain subintervalization** to obtain conservative range bounds for nonlinear functions.

The key features are:

- Tensor-product Bernstein polynomial approximation
- Multiple interpolation-node choices
- Chebyshev-Lobatto and Legendre-Lobatto sampling options
- Recursive Bernstein basis construction
- Tensor-product coefficient computation
- Affine Arithmetic polynomial evaluation
- Bernstein convex-hull based polynomial bounds
- Direct AA residual evaluation
- Bernstein coefficient-based residual evaluation
- Domain subintervalization
- Global hull aggregation
- Benchmark functions including Branin, Ackley, and Eggholder

Most importantly, the implementation distinguishes between **how the Bernstein coefficients are obtained** and **how the polynomial itself is represented**:

> Choosing `node_type="chebyshev"` does **not** construct a Chebyshev polynomial. It uses Chebyshev-Lobatto points as interpolation nodes to compute the coefficients of a polynomial that is still represented in the **Bernstein basis**.

This distinction is central when interpreting the coefficient convex-hull bound and the resulting residual enclosure.
