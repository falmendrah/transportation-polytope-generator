# Transportation polytope generator

This repository contains Python implementations of functions described in the paper "All Linear and Integer Programs are Slim 3-way Transportation Programs" which shows any rational convex polytope is polynomial-time representable as a three-way line-
sum transportation polytope of “slim” $(r, c, 3)$ format.

__Reference__

[De Loera and Onn, 2006] J. A. De Loera and S. Onn. All linear and integer programs are slim 3-way transportation
programs. SIAM Journal on Optimization, 17(3):806–821, 2006.

## Files

- `preprocessing.py`: Contains function to perform a coefficient reduction process on convex polytopes in standard form with "large" values in their defining matrix. This corresponds to stage 1 (section 3) in [De Loera and Onn, 2006].
- `plane_sum.py`: Contains the `plane_sum_entry_forbidden` class and related functions. This corresponds to stage 2 (section 3) in [De Loera and Onn, 2006].
- `slim_line_sum.py`: Contains the `slim_line_sum` class and related functions. This corresponds to stage 3 (section 3) in [De Loera and Onn, 2006].
- `embedding.py`: Contains the `slim_line_sum_representation()` function to represent convex polytopes as slim transportation polytopes as well as functions `embed_in_plane_sum(), embed_in_line_sum()` to map an integer point from a convex polytope to their image in the transportation polytope acoording to the linear isomorphism provided in the proof of the main result of [De Loera and Onn, 2006].
- `example_usage.ipynb`: Jupyter notebook demonstrating usage with examples.

## Installation

To use this code, clone this repository and install the required dependencies.

To import the functions to a Python file make sure to modify `sys.path` as needed to include the parent directory of `trans_polytope_repr`.