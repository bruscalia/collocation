import numpy as np
from scipy.optimize import root
from scipy.special import gamma, hyp2f1
import math


class OrthogonalCollocation:
    
    def __init__(self, fun, bc, n, a, x0=0.0, x1=1.0, vectorized=True):
        """
        Class for creating Orthogonal Collocation problems.
        
        When describing the shape of arguments, m is the number of decision variables and n is the number of elements along x coordinates.
        
        For details see Villadsen and Stewart (1967).
        
        Villadsen, J. V., and Warren E. Stewart. "Solution of boundary-value problems by orthogonal collocation." Chemical Engineering Science 22.11 (1967): 1483-1501.

        Parameters
        ----------
        
        fun : callable
            Function in the format ``f(x, y, dy, d2y, *args)`` that returns zeros in internal collocation points.
            Here x is a either a scalar or a (n,) shape vector, and y is a either a vector with shape (m,) or (m, n), in which each line corresponds to a single variable.
            The function must return array_like with y shape. These choices are related to the 'vectorized' parameter.
            
        bc : callable
            Function in the format ``bc(x, y, dy, d2y, *args)`` that returns a vector of zeros, shape (m,) in the boundary.
        
        n : int
            Number of internal collocation points.
        
        a : int
            Geometry: 
            
                - 1 for slabs
                - 2 for cylinders
                - 3 for spheres.
        
        x0 : float, optional
            Starting point in independent variable. Must be the symmetry point. Defaults to 0.0.
        
        x1 : float, optional
            Boundary point. Defaults to 1.0.
            
        vectorized : bool, optional
            Either or nor fun is implemented in vectorized fashion.
        """
        
        self.fun = FunctionCaller(fun, vectorized=vectorized)
        self.bc = bc
        self.n = n
        self.a = a
        self.x0 = x0
        self.x1 = x1
        self.scale_factor = (x1 - x0) ** (a - 1)
        
        points = OrthogonalCollocation._collocation_points(x0, x1, n, a)
        self.points = points["points"]
        self.x = points["x"]
        self.poly = points["poly"]
        self.Px2 = points["Px2"]
        
        matrices = OrthogonalCollocation._collocation_matrices(self.points, a)
        self.A = matrices["A"]
        self.B = matrices["B"]
        self.Q = matrices["Q"]
        self.R = matrices["R"]
        self.T = matrices["T"]
        self.W = matrices["W"]
    
    def collocate(self, y0, args=(), kwargs={}, **options):
        """
        Method to find y values in collocation points.

        Parameters
        ----------
        
        y0 : 2d array
            Initial estimates for collocation dependent variables. Shape (m, n).
            
        args : tuple, optional
            Additional arguments passed to functions. Defaults to ().
            
        kwargs : dict, optional
            Additional arguments passed to functions ignore for now. Defaults to {}.
            
        **options:
            keyword arguments passed to scipy.optimize.root.
        """
        
        y0 = np.array(y0).flatten()
        
        sol_root = root(self._obj_collocation, y0, args=args, **options)
        
        if not sol_root.success:
            print("Warning: solution did not converge.")
            self._root = sol_root
        
        self.y = sol_root.x.reshape([-1, self.n + 1])
        self._root = sol_root
    
    def effectiveness(self, fun, args=(), kwargs={}):
        """
        Integrate a function over the volume and obtain a coefficient compared to the surface value.

        Parameters
        ----------
        
        fun : callable
            Function to evaluate in the format ``f(x, y, *args, **kwargs)``
            
        args : tuple, optional
            Additional arguments passed to functions. Defaults to ().
            
        kwargs : dict, optional
            Additional arguments passed to functions ignore for now. Defaults to {}.

        Returns
        -------
        
        1d array
            Effectiveness factors for function outputs.
        """

        rates = fun(self.x, self.y, *args, **kwargs)
        int_func = rates.dot(self.W)
        int_surf = rates[..., -1] * self.W.sum()
        return int_func / int_surf
    
    def interpolate(self, x):
        """
        Calculates interpolated y values based on collocation points.

        Parameters
        ----------
        
        x : 1d array like
            Reference values of shape (n,).

        Returns
        -------
        
        2d array
            y values in shape (m, n)
        """
        x = np.array(x)
        
        x = ((x - self.x0) / (self.x1 - self.x0))
        X = np.atleast_1d(x).reshape([-1, 1]) ** (np.arange(self.n + 1) * 2)
        y = X.dot(self.theta).T
        return y
    
    @property
    def theta(self):
        """Interpolation polynomial coefficients."""
        
        return np.linalg.solve(self.Q, self.y.T)
    
    def Pi(self, i, x):
        """
        Evaluates the problem Jacobi polynomial of order i in x^2, Pi(x^2).

        Parameters
        ----------
        
        i : int
            Polynomial order.
            
        x : float or array like
            Independet variable.

        Returns
        -------
        
        float or array like
            Values of Pi(x^2)
        """
        x = np.array(x)
        
        return hyp2f1(-i, i + self.a/2 + 1, self.a/2, x**2)
    
    @property
    def ak(self):
        """Returns values of coefficients from the exapansion form as in Villadsen and Stewart (1967). \
            For vectorized implementations, theta is preferred.

        Returns
        -------
        
        array like
            Values of polynomial coefficients from original exapanded form.
        """
        
        return np.array([1 / OrthogonalCollocation.calc_const(i, self.a) *\
            (self.W * (self.y - self.y[:, [-1]])).dot(self.Pi(i, self.points))\
                for i in range(self.n)]).T
    
    def expansion(self, x):
        """
        Evaluates function in x using original expansion form from Villadsen and Stewart (1967), equation (6).

        Parameters
        ----------
        
        x : float or array like
            Value of independent variable.

        Returns
        -------
        
        float or array like
            Dependent variables evaluated in x.
        """
        x = np.array(x)
        
        return self.y[:, [-1]] + (1 - x**2) * np.sum([self.ak[:, [i]] * np.atleast_2d(self.Pi(i, x))\
            for i in range(self.n)], axis=0)
    
    @staticmethod
    def _collocation_points(x0, x1, n, a):
        
        poly = np.zeros(2 * n + 1)
        poly[0] = 1.
        polyx2 = np.zeros(n + 1)
        polyx2[0] = 1.
        
        for i in range(1, n + 1):
            
            _range_i = np.arange(i)
            p = np.prod(-n + _range_i) * np.prod(n + a/2 + 1 + _range_i)\
                / (np.prod(a/2 + _range_i) * math.factorial(i))

            poly[2 * i] = p
            polyx2[i] = p
            
        poly = np.flip(poly)
        all_roots = np.roots(poly)
        positive_roots = np.sort(all_roots[all_roots > 0])
        points = np.append(positive_roots, [1])
        x_points = points * (x1 - x0) + x0
        
        return {'points':points, 'x':x_points, "poly":poly, "Px2":polyx2}
    
    @staticmethod
    def _collocation_matrices(points, a):

        points = np.array(points)
        _len = len(points)

        base_mat = np.column_stack((points,) * _len)
        _range_n = np.arange(_len)
        
        Q = base_mat.copy() ** (2 * _range_n)
        Qinv = np.linalg.inv(Q)
        
        R = Q.copy() / base_mat.copy() * (2 * _range_n)
        
        T_base = R.copy() / base_mat.copy()
        T = T_base * (np.maximum(2 * _range_n - 1, 0) + (a - 1))
        
        A = np.matmul(R, Qinv)
        B = np.matmul(T, Qinv)
        
        WF = (1 ** (2 * _range_n + a)) / (2 * _range_n + a)\
            - (0 ** (2 * _range_n + a)) / (2 * _range_n + a)

        W = np.matmul(WF, Qinv)
        
        return {'A':A, 'B':B, 'Q':Q, 'R':R, 'T':T, 'W':W}
    
    @staticmethod
    def calc_const(i, a):
        
        C = gamma(a/2)**2 * gamma(i + 1) * gamma(i + 2) / (\
            (4*i + a + 2) * gamma(i + a/2) * gamma(i + a/2 + 1))
        
        return C
    
    def _evaluate(self, y, *args, **kwargs):
        
        dy_ = y.dot(self.A.T)
        d2y_ = y.dot(self.B.T) / self.scale_factor
        
        res = self.fun(self.x, y, dy_, d2y_, *args, **kwargs)
        res[..., -1] = self.bc(self.x[-1], y[..., -1], dy_[..., -1], d2y_[..., -1], *args, **kwargs)
        
        return res
    
    def _obj_collocation(self, y, *args, **kwargs):
        
        y = y.reshape([-1, self.n + 1])
        res = self._evaluate(y, *args, **kwargs).flatten()
        return res


class FunctionCaller:
    
    def __init__(self, fun, vectorized=True):
        
        self.fun = fun
        self.vectorized = vectorized
        self.caller = self._vectorized_call if vectorized else self._elementwise_call
    
    def __call__(self, x, y, dy, d2y, *args, **kwargs):
        return self.caller(x, y, dy, d2y, *args, **kwargs)
    
    def _vectorized_call(self, x, y, dy, d2y, *args, **kwargs):
        return self.fun(x, y, dy, d2y, *args, **kwargs)

    def _elementwise_call(self, x, y, dy, d2y, *args, **kwargs):
        return  np.column_stack([self.fun(xj, y[:, j], dy[:, j], d2y[:, j], *args, **kwargs) \
            for j, xj in enumerate(x)])
    
