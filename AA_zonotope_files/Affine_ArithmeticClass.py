# %%
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  23, 2026

@author: Premjit Saha
"""

import numpy as np
import sympy as sp
import cvxpy as cp
import itertools
from scipy.optimize import root
from scipy.stats import qmc

import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class AAZonotopeModel:
    """
    N-D affine approximation and Chebyshev minimax model.
    """

    def __init__(self, vars, f_target, V_base):
        """
        Parameters
        ----------
        vars : list[sympy.Symbol]
            Independent variables
        f_target : sympy.Expr
            Target function f(x)
        V_base : ndarray, shape (N, 2)
            Interval bounds for each variable
        """

        # ---- store symbolic data ----
        self.vars = vars
        self.f_target = f_target

        # ---- interval geometry ----
        self.V_base = np.asarray(V_base, dtype=float)
        self.dim = len(vars)

        self.x0 = 0.5 * (self.V_base[:, 0] + self.V_base[:, 1])
        self.x1 = 0.5 * (self.V_base[:, 1] - self.V_base[:, 0])

        # ---- lambdified functions ----
        self.f_fun = sp.lambdify(self.vars, self.f_target, modules="numpy")

        self.grad = [sp.diff(self.f_target, v) for v in self.vars]
        self.grad_fun = sp.lambdify(self.vars, self.grad, modules="numpy")

        self.bounds =(-1.0, 1.0)

    def _grad_vec(self, z):
        return np.asarray(self.grad_fun(*z), dtype=float).ravel()
    
    
    def chebyshev_fit(
        self,
        n_starts=None,
        M=None
    ):
        """
        Chebyshev minimax affine fit.
        """
        
        xmin, xmax = self.bounds

        # ---- multi-start stationary points ----
        initial_guesses = np.random.uniform(
            xmin, xmax, size=(n_starts, self.dim)
        )

        solutions = []
        for x0 in initial_guesses:
            sol = root(self._grad_vec, x0, method="hybr")
            if sol.success:
                solutions.append(sol.x)

        solutions = np.array(solutions)

        if solutions.size == 0:
            solutions = solutions.reshape(0, self.dim)

        # ---- boundary vertices ----
        eps_corners = np.array(
            list(itertools.product([xmin, xmax], repeat=self.dim))
        )


        all_points = np.vstack((eps_corners, solutions))
        all_points = np.unique(np.round(all_points, 18), axis=0)

        # ---- Latin Hypercube samples ----
        sampler = qmc.LatinHypercube(d=self.dim)
        samples = qmc.scale(
            sampler.random(n=M),
            xmin * np.ones(self.dim),
            xmax * np.ones(self.dim)
        )

        eps_all = np.vstack((all_points, samples))
        eps_all = np.unique(np.round(eps_all, 18), axis=0)

        # ---- physical samples ----
        X = np.array([
            self.x0 + np.diag(self.x1) @ e for e in eps_all
        ])

        # ---- target values ----
        y = np.asarray(self.f_fun(*eps_all.T), dtype=float).ravel()

        # ---- build V matrix ----
        V = np.column_stack([np.ones(len(X)), X])

        # ---- Chebyshev LP ----
        c = cp.Variable(self.dim + 1)
        E = cp.Variable(nonneg=True)

        constraints = [
            V @ c - y <= E,
            V @ c - y >= -E
        ]

        problem = cp.Problem(cp.Minimize(E), constraints)
        problem.solve()

        self.c_cheb = c.value
        self.E_cheb = E.value

        return self.c_cheb, self.E_cheb
    

    def affine_min_rangeLP(self, M=None, solver=None):
        """
        Minimum-range affine approximation via LP.
        """

        # ---- gradient bounds ----
        xmin, xmax = self.bounds

        g_low = np.asarray(
        self.grad_fun(*(xmin * np.ones(self.dim))),
        dtype=float
        )

        g_high = np.asarray(
        self.grad_fun(*(xmax * np.ones(self.dim))),
        dtype=float
        )

        d_min = np.minimum(g_low, g_high)
        d_max = np.maximum(g_low, g_high)

        alpha = cp.Variable(self.dim)
        z = cp.Variable(self.dim)

        objective = cp.Minimize(cp.sum(z))

        constraints = [
            z >=  cp.multiply(self.x1, alpha),
            z >= -cp.multiply(self.x1, alpha),
            z >= 0,
            alpha >= d_min,
            alpha <= d_max
        ]

        prob = cp.Problem(objective, constraints)
        prob.solve(solver=solver)

        alpha_star = alpha.value

        # -----------------------------
        # Sample domain (corners + LHS)
        # -----------------------------
        
        sampler = qmc.LatinHypercube(d=self.dim)
        samples = sampler.random(n=M)

        lower = np.ones(self.dim) * xmin
        upper = np.ones(self.dim) * xmax

        samples_scaled = qmc.scale(samples, lower, upper)

        eps_corners = np.array(list(itertools.product([-1, 1], repeat=self.dim)))

        eps_all = np.vstack((samples_scaled, eps_corners))

        X = np.array([
            self.x0 + np.diag(self.x1) @ e for e in eps_all
        ])

        f_vals = self.f_fun(*eps_all.T)
        g_vals = f_vals - X @ alpha_star

        gamma = 0.5 * (g_vals.max() + g_vals.min())
        delta = 0.5 * (g_vals.max() - g_vals.min())

        self.c_lp = np.hstack((gamma, alpha_star))
        self.delta_lp = delta

        return self.c_lp, self.delta_lp
    

    def zonotope(self, c, E):
        """
        Construct zonotope from affine coefficients.

        Returns
        -------
        center : ndarray, shape (dim+1,)
        G      : ndarray, shape (dim+1, dim+1)
        vertices : ndarray, shape (2^(dim+1), dim+1)
        f_range : list [f_min, f_max]
        """

        c0 = c[0]
        c_lin = c[1:]

        y0 = c0 + c_lin @ self.x0
        center = np.hstack([self.x0, y0])

        G = np.zeros((self.dim + 1, self.dim + 1))
        for i in range(self.dim):
            G[i, i] = self.x1[i]
            G[-1, i] = c_lin[i] * self.x1[i]

        G[-1, -1] = E

        eps = np.array(
            list(itertools.product([-1, 1], repeat=self.dim + 1))
        )

        vertices = np.array([center + G @ e for e in eps])

        f_range = [float(vertices[:, -1].min()), float(vertices[:, -1].max())]

        return center, G, vertices, f_range
    

    def plot_zonotope_3d(
        self,
        vertices,
        var_names,
        f_name,
        grid_points=None,
        surface_alpha=None,
        poly_alpha=None
    ):
        """
        Plot a 3D zonotope/polytope with the true function surface.

        Parameters
        ----------
        vertices : ndarray, shape (N, 3)
            Zonotope/polytope vertices [x, y, f]
        var_names : list[str]
            Names of the independent variables (length 2)
        f_name : str
            Name/label of the function output
        grid_points : int
            Resolution of surface plot
        surface_alpha : float
            Transparency of surface
        poly_alpha : float
            Transparency of polytope
        """

        # --------------------------------------------------
        # Dimension check
        # --------------------------------------------------

        if vertices.shape[1] != 3:
            raise ValueError(
                f"Plotting only supported for 3D objects. "
                f"Received vertices with shape {vertices.shape}."
            )

        if len(var_names) != 2:
            raise ValueError("var_names must have length 2.")
        
        # --------------------------------------------------
        # Build surface grid
        # --------------------------------------------------
        x_vals = np.linspace(vertices[:, 0].min(), vertices[:, 0].max(), grid_points)
        y_vals = np.linspace(vertices[:, 1].min(), vertices[:, 1].max(), grid_points)

        eps_vals = np.linspace(-1, 1, grid_points)

        X, Y = np.meshgrid(x_vals, y_vals)
        X1, Y1 = np.meshgrid(eps_vals, eps_vals)

        Z = self.f_fun (X1, Y1)

        # --------------------------------------------------
        # Convex hull of zonotope
        # --------------------------------------------------

        hull = ConvexHull(vertices)
        faces = [vertices[simplex] for simplex in hull.simplices]

        # --------------------------------------------------
        # Plot
        # --------------------------------------------------
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        poly = Poly3DCollection(
        faces,
        alpha=poly_alpha,
        edgecolor="black"
        )
        ax.add_collection3d(poly)

        # Vertices
        ax.scatter(
            vertices[:, 0],
            vertices[:, 1],
            vertices[:, 2],
            color="red",
            s=20
        )

        # Surface
        ax.plot_surface(
            X, Y, Z,
            cmap="viridis",
            alpha=surface_alpha
        )

        # Axis limits
        mins = vertices.min(axis=0)
        maxs = vertices.max(axis=0)

        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])

        # Labels
        ax.set_xlabel(var_names[0])
        ax.set_ylabel(var_names[1])
        ax.set_zlabel(f_name)

        plt.tight_layout()
        plt.show()

    def plot_zonotope_2d(
        self,
        cvec,
        vertices,
        xlabel="x",
        ylabel="f(x)",
        title="2D Zonotope",
        grid_points=None,
        fill_alpha=None
    ):
        """
        Plot a 2D zonotope / polytope.
        Parameters
        ----------
        cvec : optimal parameters from zonotope approximation

        vertices : ndarray, shape (N, 2)
            Zonotope vertices [x, f]
        xlabel, ylabel : str
            Axis labels
        title : str
            Plot title
        grid_points : float
            Resolution of surface plot
        fill_alpha : float
            Transparency of zonotope fill
        """

        vertices = np.asarray(vertices)

        if vertices.shape[1] != 2:
            raise ValueError(
                f"2D plot requires vertices of shape (N,2), "
                f"got {vertices.shape}"
            )

        # ----------------------------------
        # Order vertices counter-clockwise
        # ----------------------------------
        centroid = vertices.mean(axis=0)

        angles = np.arctan2(
            vertices[:, 1] - centroid[1],
            vertices[:, 0] - centroid[0]
        )

        order = np.argsort(angles)
        vertices_ord = vertices[order]

        # close polygon
        vertices_closed = np.vstack([vertices_ord, vertices_ord[0]])

        # ----------------------------------
        # Plot
        # ----------------------------------
        plt.figure(figsize=(6, 6))

        # True function 
        x_vals = np.linspace(vertices[:, 0].min(), vertices[:, 0].max(), grid_points)
        eps_vals = np.linspace(-1, 1, grid_points)
        y_vals = self.f_fun (eps_vals)
        plt.plot(x_vals, y_vals, "r-", label="True function")

        # Affine approximation
        y_aff = cvec[0] + cvec[1]*x_vals
        plt.plot(x_vals, y_aff, "b--", label="Affine approx")

        # Zonotope boundary
        plt.plot(
            vertices_closed[:, 0],
            vertices_closed[:, 1],
            color="black",
            linewidth=2,
            label="Zonotope boundary"
        )

        # Vertices
        plt.scatter(
            vertices_ord[:, 0],
            vertices_ord[:, 1],
            color="gray",
            zorder=3
        )

        # Fill zonotope
        plt.fill(
            vertices_ord[:, 0],
            vertices_ord[:, 1],
            alpha=fill_alpha
        )

        # Center
    
        plt.scatter(
            centroid[0],
            centroid[1],
            color="black",
            marker="x",
            s=100,
            label="Center"
        )

        plt.axis("equal")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.show()



   
       

  





