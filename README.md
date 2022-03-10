# collocation
A repository with a python implementation of Orthogonal Collocation by Villadsen and Stewart (1967) for solving second order boundary value problems.

See two simple examples and a real world problem in the notebook example_collocation.

¬pip install -e git+https://github.com/bruscalia/collocation#egg=collocation

## Usage

```
import numpy as np
from collocation.bvp import OrthogonalCollocation
```

```
#Internal function
def fun_1(x, y, dy, d2y, k):
    
    return d2y[0] + k * y[0] + 1

#Boundary function
def bc_1(x, y, dy, d2y, k):
    
    return dy[0] - 1

k = 1.0
```

```
#Number of collocatioin points
n_points = 6

#Create problem
problem_1 = OrthogonalCollocation(fun_1, bc_1, n_points, 1, x0=0.0, x1=1.0, vectorized=True)

#Initial estimation
y01 = np.zeros([1, n_points + 1])

#Collocation using scipy.optimize.root in backend
problem_1.collocate(y01, args=k, method="hybr", tol=1e-6)
```

This was developed as a part of the modeling in the published article: "Simulation and optimization of axial-flow and radial-flow reactors for dehydrogenation of ethylbenzene into styrene based on a heterogeneous kinetic model," Chem. Eng. Sci., vol. 244, p. 116805, 2021. doi:10.1016/j.ces.2021.116805.

The code is uploaded on ResearchGate and can be cited linked to doi:10.13140/RG.2.2.17223.78240.
