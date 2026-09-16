[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

===============================================================================
PROJECT: Laplace Limit Computation Engine
===============================================================================

OVERVIEW:
Calculates the Laplace Limit (L ≈ 0.6627434193491815...) to arbitrary 
precision (up to 1,000,000+ digits). The Laplace Limit defines the maximum 
eccentricity for which Kepler's equation for planetary orbits can be solved using 
a converging power series.

ALGORITHM & IMPLEMENTATION:
- Newton-Raphson Root Finding: Solves the transcendental equation:
    x * exp(sqrt(1 + x^2)) = 1 + sqrt(1 + x^2)
- Dynamic Precision Doubling: Starts at low precision and doubles the working 
  precision at each iteration to match the quadratic convergence rate of 
  Newton's method, minimizing CPU time spent on initial approximations.

## Usage

```bash
python "Laplace Limit.py" --help
```
