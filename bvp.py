import numpy as np
from scipy.optimize import root
from scipy.special import gamma
import math


class OrthogonalCollocation:
    
    def __init__(self, fun, bc, n, a, x0=0.0, x1=1.0):
        
        self.fun = fun
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
        
        y0 = y0.flatten()
        
        sol_root = root(self._obj_collocation, y0, args=args, **options)
        
        if not sol_root.success:
            print("Warning: solution did not converge.")
            self._root = sol_root
        
        self.y = sol_root.x.reshape([-1, self.n + 1])
        self._root = sol_root
    
    def effectiveness(self, fun, args=(), kwargs={}):

        rates = fun(self.x, self.y, *args, **kwargs)
        int_func = rates.dot(self.W)
        int_surf = rates[..., -1] * self.W.sum()
        return int_func / int_surf
    
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
    
    def extrapolate(self, x):
        
        x = ((x - self.x0) / (self.x1 - self.x0))
        X = np.atleast_1d(x).reshape([-1, 1]) ** (np.arange(self.n + 1) * 2)
        y = X.dot(self.theta).T
        return y
    
    @property
    def theta(self):
        
        return np.linalg.solve(self.Q, self.y.T)
    
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