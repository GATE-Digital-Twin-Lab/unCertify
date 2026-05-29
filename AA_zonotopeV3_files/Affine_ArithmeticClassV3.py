import numpy as np


class AffineScalar:

    def __init__(self, parent=None, idx=None, x0=None, xi=None):

        # view mode
        if parent is not None:

            self.parent = parent
            self.idx = idx

            self._standalone = False

        # standalone mode
        else:

            self._x0 = float(x0)

            self._xi = np.asarray(xi, dtype=float)

            self._standalone = True
    @property
    def x0(self):

        if self._standalone:

            return self._x0

        return self.parent.x0[self.idx]
    @property
    def xi(self):

        if self._standalone:

            return self._xi

        return self.parent.E[self.idx]

    @property
    def interval(self):

        r = np.sum(np.abs(self.xi))

        return (
            self.x0 - r,
            self.x0 + r
        )
    
    @staticmethod
    def _align(xi1, xi2):

        n1 = len(xi1)

        n2 = len(xi2)

        if n1 < n2:

            xi1 = np.pad(
                xi1,
                (0, n2 - n1)
            )

        elif n2 < n1:

            xi2 = np.pad(
                xi2,
                (0, n1 - n2)
            )

        return xi1, xi2

    def __add__(self, other):

        # affine + affine
        if isinstance(other, self.__class__):

            xi1, xi2 = self._align(
                self.xi,
                other.xi
            )

            return AffineScalar(
                x0 = self.x0 + other.x0,
                xi = xi1 + xi2
            )

        # affine + scalar
        elif isinstance(other, (int, float)):

            return AffineScalar(
                x0 = self.x0 + other,
                xi = self.xi.copy()
            )

        raise TypeError(
            "other must be AffineScalar, int, or float"
        )
    
    def __mul__(self, other):

        # affine * affine
        if isinstance(other, self.__class__):

            # ---------- align dimensions ----------

            xi1, xi2 = self._align(
                self.xi,
                other.xi
            )

            # ---------- affine part ----------

            x0_new = self.x0 * other.x0

            xi_aff = (
                self.x0 * xi2
                +
                other.x0 * xi1
            )

            # ---------- improved error term (26) ----------

            v = xi1 * xi2

            v_pos = np.maximum(v, 0.0)

            v_neg = np.maximum(-v, 0.0)

            diag_term = max(
                np.sum(v_pos),
                np.sum(v_neg)
            )

            offdiag = 0.0

            n = len(xi1)

            for i in range(n):

                for j in range(i + 1, n):

                    offdiag += abs(
                        xi1[i] * xi2[j]
                        +
                        xi1[j] * xi2[i]
                    )

            e = diag_term + offdiag

            # ---------- append fresh noise symbol ----------

            xi_new = np.append(
                xi_aff,
                e
            )

            return AffineScalar(
                x0 = x0_new,
                xi = xi_new
            )

        # affine * scalar
        elif isinstance(other, (int, float)):

            return AffineScalar(
                x0 = self.x0 * other,
                xi = self.xi * other
            )

        raise TypeError(
            "other must be AffineScalar, int, or float"
        )
    
    def __rmul__(self, other):

        return self * other
    
    def __neg__(self):

        return AffineScalar(
            x0 = -self.x0,
            xi = -self.xi
        )
    
    def __sub__(self, other):

        return self + (-other)
    
    def __rsub__(self, other):

        return other + (-self)

    def __radd__(self, other):

        return self + other

    def __repr__(self):

        return (
            f"AffineScalar("
            f"x0={self.x0}, "
            f"xi={self.xi})"
        )


class AffineArray:

    def __init__(self, x0, E):

        self.x0 = np.asarray(x0, dtype=float)

        self.E = np.asarray(E, dtype=float)

    @classmethod
    def from_intervals(cls, intervals):

        intervals = np.asarray(intervals, dtype=float)

        lb = intervals[:, 0]
        ub = intervals[:, 1]

        x0 = 0.5 * (lb + ub)

        r = 0.5 * (ub - lb)

        E = np.diag(r)

        return cls(x0, E)

    @property
    def interval(self):

        rad = np.sum(np.abs(self.E), axis=1)

        lb = self.x0 - rad
        ub = self.x0 + rad

        return np.column_stack((lb, ub))

    def __getitem__(self, idx):

        return AffineScalar(self, idx)

    def __repr__(self):

        return (
            f"AffineArray(\n"
            f"x0={self.x0},\n"
            f"E=\n{self.E}\n)"
        )