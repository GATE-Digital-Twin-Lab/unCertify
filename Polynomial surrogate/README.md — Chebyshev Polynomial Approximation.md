# Chebyshev Polynomial Approximation with Affine Arithmetic

This folder contains the implementation of a **Chebyshev polynomial approximation framework combined with affine arithmetic (AA)** for rigorous range bounding of nonlinear functions.

The main objective is to approximate a nonlinear function locally by a polynomial and then use affine arithmetic to obtain a conservative bound on the approximation error (residual). The computational domain can be subdivided into multiple subdomains to reduce the overestimation associated with direct interval or affine-arithmetic evaluation.

---

## 1. Overview

For a function

$\[
f(x), \qquad x\in[0,1]^d,
\]
$

the method constructs a local Chebyshev polynomial approximation

$
\[
p(x) \approx f(x)
\]
$
on each subdomain.

The residual is then defined as

$
\[
r(x)=f(x)-p(x).
\]
$
Instead of relying only on the approximation error observed at the interpolation nodes, the residual is evaluated using affine arithmetic:
$
\[
r(X)=f(X)-p(X),
\]
$
which provides a conservative enclosure of the residual over the complete subdomain.

The resulting global enclosure is obtained by taking the hull of the local interval enclosures over all subdomains.

The overall procedure is therefore
$
\[
\boxed{
f
\;\longrightarrow\;
p
\;\longrightarrow\;
r=f-p
\;\longrightarrow\;
\text{AA bound on }r
}
\]
$
combined with domain subintervalization.

---

## 2. Main Components

The implementation is organized into the following files.

### `ChebyshevFiles.py`

Contains the numerical construction of the Chebyshev polynomial approximation.

The main steps are:

1. Generate Chebyshev-Lobatto nodes.
2. Construct the multidimensional tensor-product grid.
3. Evaluate the original function at the grid points.
4. Reshape the function values into an N-dimensional tensor.
5. Compute the Chebyshev coefficients using a multidimensional DCT-I.
6. Return the resulting coefficient tensor.

The main function is:

```python
chebyshev_approximation(func, nvec, combo)
```

which returns

```text
nodes
points
F_grid
C
```

where `C` contains the Chebyshev coefficients.

---

### `AA_evaluationCh.py`

Contains the affine-arithmetic evaluation of the Chebyshev polynomial and the residual.

The main functions are:

```python
map_to_chebyshev(...)
polylist_Ch(...)
chebyshev_basis_AA(...)
evaluate_chebyshev_AA(...)
residual_eval_AA(...)
```

The physical subdomain
$
\[
[a_k,b_k]
\]
$
is mapped to the standard Chebyshev interval
$
\[
[-1,1]
\]
$
using
$
\[
\xi_k =
\frac{2x_k-(a_k+b_k)}
     {b_k-a_k}.
\]
$
The Chebyshev basis is then generated using the three-term recurrence
$
\[
T_0(\xi)=1,
\]
$
$
\[
T_1(\xi)=\xi,
\]
$
and, for $\(k\geq1\),$
$
\[
T_{k+1}(\xi)
=
2\xi T_k(\xi)-T_{k-1}(\xi).
\]
$
The multidimensional polynomial is represented as a tensor-product expansion
$
\[
p(\xi_1,\ldots,\xi_d)
=
\sum_{i_1=0}^{n_1}
\cdots
\sum_{i_d=0}^{n_d}
C_{i_1,\ldots,i_d}
\prod_{k=1}^{d}
T_{i_k}(\xi_k).
\]
$
The polynomial is evaluated using affine arithmetic.

---

### `subintervalization.py`

Provides utilities for domain subdivision and global interval hull computation.

The original normalized domain is
$
\[
[0,1]^d.
\]
$
The function

```python
split_interval(splits_per_dim)
```

divides each dimension into the requested number of equal subintervals.

For example,

```python
split_interval([4, 4])
```

creates
$
\[
4\times4=16
\]
$
two-dimensional subdomains.

The function

```python
hull_intervals(intervals)
```

computes the global interval hull from all local interval bounds.

---

### `test_functions.py`

Contains NumPy and affine-arithmetic implementations of the benchmark functions.

Currently included:

- Branin
- Ackley
- Eggholder

Each benchmark generally has two implementations:

```python
Branin_numpy(...)
test_BraninAA(...)

Ackley_numpy(...)
test_AckleyAA(...)

Eggholder_numpy(...)
test_EggholderAA(...)
```

The NumPy implementation is used to construct the polynomial approximation, while the AA implementation is used for rigorous evaluation of the original function.

---

## 3. Computational Workflow

For every subdomain, the following procedure is performed.

### Step 1 — Generate subdomains

For a specified number of splits per dimension,

```python
subdomains = split_interval(splits_per_dim)
```

the normalized domain `[0,1]^d` is partitioned into smaller regions.

---

### Step 2 — Construct an affine representation

For every subdomain,

```python
X_sub = AffineArray.from_intervals(combo)
```

creates an affine-arithmetic representation of the uncertain variables.

---

### Step 3 — Construct the local Chebyshev approximation

The function

```python
chebyshev_approximation(
    func_numpy,
    nvec,
    combo
)
```

generates the Chebyshev-Lobatto grid and evaluates the numerical function at the grid points.

The Chebyshev coefficients are then computed using a multidimensional DCT-I.

The resulting approximation is
$
\[
p(x)=
\sum_{\mathbf{k}}
C_{\mathbf{k}}
T_{\mathbf{k}}(x).
\]
$
---

### Step 4 — Evaluate the original function using AA

The affine-arithmetic implementation of the original function is evaluated over the entire subdomain:

```python
f_AA = func_AA(X_sub, ...)
```

This produces an enclosure

$
\[
f(X_{\text{sub}})
\subseteq
[f_{\mathrm{lo}},f_{\mathrm{hi}}].
\]
$
---


### Step 5 — Evaluate the polynomial using AA

The same affine variables are used to evaluate the Chebyshev polynomial:

```python
p_AA = evaluate_chebyshev_AA(...)
```

This gives an enclosure

$
\[
p(X_{\text{sub}})
\subseteq
[p_{\mathrm{lo}},p_{\mathrm{hi}}].
\]
$
---


### Step 6 — Bound the residual

The residual is calculated as

$
\[
r(X_{\text{sub}})
=
f(X_{\text{sub}})
-
p(X_{\text{sub}}).
\]
$

This is implemented by

```python
residual = f - p
```

and provides a conservative affine-arithmetic enclosure

$
\[
r(X_{\text{sub}})
\subseteq
[r_{\mathrm{lo}},r_{\mathrm{hi}}].
\]
$
---


## 4. Global Bounds

After every subdomain has been evaluated, the local bounds are collected:

```python
f_intervals
p_intervals
residual_intervals
```

The global hull is then calculated:

```python
f_hull = hull_intervals(f_intervals)
p_hull = hull_intervals(p_intervals)
residual_hull = hull_intervals(residual_intervals)
```

Thus, the global residual enclosure is
$
\[
r(X)
\subseteq
\operatorname{hull}
\left(
\bigcup_i r(X_i)
\right),
\]
$
where
$
\[
X=\bigcup_i X_i.
\]
$
The widths of the resulting global bounds are also reported:

```python
Width(f)
Width(p)
Width(residual)
```

These quantities can be used to study the effect of polynomial approximation and domain refinement on the resulting enclosure.

---

## 5. Subintervalization

The main experiment varies the number of subdivisions per dimension:

```python
split_values = [
    1, 2, 4, 8, 16,
    32, 64, 128, 256
]
```

For a two-dimensional problem, this corresponds to:

| Splits | Number of subdomains |
|---:|---:|
| 1 × 1 | 1 |
| 2 × 2 | 4 |
| 4 × 4 | 16 |
| 8 × 8 | 64 |
| 16 × 16 | 256 |
| 32 × 32 | 1024 |
| 64 × 64 | 4096 |
| 128 × 128 | 16384 |
| 256 × 256 | 65536 |

Increasing the number of subdomains makes the polynomial approximation local to smaller regions of the domain. The resulting residual bounds can therefore be studied as a function of the level of domain refinement.

---

## 6. Polynomial Degree

The polynomial degree is specified independently in each dimension:

```python
nvec = [3, 3]
```

For a two-dimensional problem, this produces a polynomial of degree 3 in each coordinate.

More generally,

```python
nvec = [n1, n2, ..., nd]
```

specifies the degree in each dimension.

The resulting tensor-product polynomial contains
$
\[
\prod_{k=1}^{d}(n_k+1)
\]
$
Chebyshev basis terms.

For example,

```python
nvec = [3, 3]
```

requires
$
\[
(3+1)(3+1)=16
\]
$
basis terms.

---

## 7. Benchmark Selection

The benchmark function can be selected in the configuration section.

For example, for the two-dimensional Ackley function:

```python
from functools import partial

fun_numpy = partial(Ackley_numpy, d=2)
fun_AA = partial(test_AckleyAA, d=2)
```

The `partial` function allows the dimension to be specified while keeping the remaining framework independent of the particular benchmark.

For example, a different Ackley dimension can be selected using:

```python
fun_numpy = partial(Ackley_numpy, d=5)
fun_AA = partial(test_AckleyAA, d=5)
```

---

## 8. Affine Arithmetic

The framework uses the `AffineArray` and `AffineScalar` classes from:

```text
Affine_ArithmeticClassV3.py
```

The affine representation separates affine uncertainty coefficients from non-affine/remainder errors.

For an affine variable, the enclosure is represented using a center and uncertainty terms. This allows correlations between repeated occurrences of the same uncertain variable to be retained during arithmetic operations.

Nonlinear operations such as

```python
sqrt()
exp()
log()
sin()
cos()
pow()
```

are implemented by the affine-arithmetic class.

The `cheb` parameter can be passed to these operations when Chebyshev-aware affine arithmetic is required.

---

## 9. Important Distinction: Approximation vs. Rigorous Bounding

The Chebyshev coefficients are obtained from floating-point evaluations of the original function at Chebyshev-Lobatto nodes.

Therefore, the polynomial construction itself is a numerical approximation step.

The rigorous bounding step is performed separately by evaluating
$
\[
f(X)-p(X)
\]
$
using affine arithmetic over each subdomain.

Consequently, the framework separates:

1. **Numerical polynomial construction**
2. **Affine-arithmetic evaluation**
3. **Residual enclosure**
4. **Domain subdivision**
5. **Global hull construction**

This distinction is important because agreement between the polynomial and the function at the interpolation nodes alone does not constitute a rigorous global error bound.

---

## 10. Running the Experiment

After configuring the benchmark function and polynomial degree, run the main experiment:

```bash
python run_Ch_poly.py
```

The script reports, for each subdivision level:

```text
Splits = ...
Number of subdomains = ...
Global f hull        = ...
Global p hull        = ...
Global residual hull = ...
Width(f)             = ...
Width(p)             = ...
Width(residual)      = ...
```

For example:

```text
Splits = 8 x 8
Number of subdomains = 64
Global f hull        = ...
Global p hull        = ...
Global residual hull = ...
Width(f)             = ...
Width(p)             = ...
Width(residual)      = ...
```

These results can then be used to compare the effect of subintervalization on the function enclosure, polynomial enclosure, and residual enclosure.

---

## 11. File Structure

A typical directory structure is:

```text
Chebyshev/
│
├── run_Ch_poly.py
├── ChebyshevFiles.py
├── AA_evaluationCh.py
├── subintervalization.py
├── test_functions.py
└── Affine_ArithmeticClassV3.py
```

### File responsibilities

| File | Purpose |
|---|---|
| `run_Ch_poly.py` | Main experiment and subintervalization study |
| `ChebyshevFiles.py` | Chebyshev nodes, grids, coefficients, and numerical approximation |
| `AA_evaluationCh.py` | AA evaluation of Chebyshev polynomial and residual |
| `subintervalization.py` | Domain subdivision and global interval hull |
| `test_functions.py` | NumPy and AA benchmark functions |
| `Affine_ArithmeticClassV3.py` | Affine arithmetic implementation |

---

## 12. Summary

The implemented approach combines **local Chebyshev polynomial approximation**, **affine-arithmetic evaluation**, and **domain subintervalization**.

For each subdomain \(X_i\),
$
\[
f(X_i)
\]
$
is evaluated directly using affine arithmetic, while a local polynomial
$
\[
p_i(x)
\]
$
is constructed numerically and evaluated using affine arithmetic. The residual

$
\[
r_i(X_i)=f(X_i)-p_i(X_i)
\]
$

is then enclosed using affine arithmetic.

The final global bounds are obtained by taking the hull of the local enclosures.

The framework therefore provides a systematic way to investigate whether replacing direct evaluation of a nonlinear function with a polynomial surrogate plus a rigorously bounded residual can reduce overestimation while maintaining conservative bounds.
