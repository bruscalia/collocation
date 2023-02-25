# collocation
A repository with a Python implementation of Orthogonal Collocation by Villadsen and Stewart (1967) for solving second order boundary value problems with symmetry.

See two simple examples and a real world problem in the notebook [example_collocation](https://github.com/bruscalia/collocation/blob/main/notebooks/example_collocation.ipynb).

## Install
As it is a very short package, we have not made it available via PyPi. So the user must either clone the repository using git (see code below) or download the files and use in a corresponding directory.

```
pip install git+https://github.com/bruscalia/collocation
```

## Usage

Consider the following system:

$$
\displaystyle
\begin{align}
  & \frac{d^2 y_1}{dx^2} + k_1 y_2 + 1 = 0\\
  & \frac{d^2 y_2}{dx^2} + k_2 log(1 + y_1) = 0\\
  & y_1 = 0, & \text{at } x = 1\\
  & y_2 - 1 = 0, & \text{at } x = 1\\
  & \frac{dy_1}{dx} = \frac{dy_2}{dx} = 0, & \text{at } x = 0
\end{align}
$$

First, import the necessary modules to solve it.

```python
import numpy as np
from collocation import OrthogonalCollocation
```

The user must define a function that returns zeros in internal points and another that returns zeros in the surface boundary.

```python
# Internal function
def fun(x, y, dy, d2y, k):
    
    return np.array([d2y[0] + k1 * y[1] + 1, d2y[1] + k2 * np.log(1 + y[0])])

# Boundary function
def bc(x, y, dy, d2y, k):
    
    return np.array([y[0], y[1] - 1])

k1 = 1.0
k2 = -1.0
```

Then must instantiate a problem using the **OrthogonalCollocation** class, define initial estimations and collocate points. The points are available in the *y* property of the problem and the method *interpolate* might provide values given *x* coordinates.

```python
# Create problem
problem = OrthogonalCollocation(fun, bc, 6, 1, x0=0.0, x1=1.0, vectorized=True)

# Initial estimation
y0 = np.zeros([2, n_points + 1])

# Collocation using scipy.optimize.root in backend
problem.collocate(y0, args=(k1, k2), method="hybr", tol=1e-6)
```

### Visualization
<p align="center">
  <img src="https://github.com/bruscalia/collocation/raw/main/images/example.png" alt="example"/>
</p>

## Citation
This was developed as a part of the modeling in the styrene reactor simulation project:

[Leite, B., Costa, A. O. S. & Costa Junior, E. F., 2021. Simulation and optimization of axial-flow and radial-flow reactors for dehydrogenation of ethylbenzene into styrene based on a heterogeneous kinetic model. Chem. Eng. Sci., Volume 244, p. 116805. doi:10.1016/j.ces.2021.116805.](https://doi.org/10.1016/j.ces.2021.116805)

The code is uploaded on ResearchGate and can be cited linked to doi:10.13140/RG.2.2.17223.78240.

## Original Article
Villadsen, J. & Stewart, W. E., 1967. Solution of boundary-value problems by orthogonal collocation. Chem. Eng. Sci., 22(11), pp. 1483-1501.

## Contact
e-mail: bruscalia12@gmail.com
