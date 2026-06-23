# Hybrid Root-Finding Algorithm

## Overview

This project implements a hybrid root-finding algorithm in MATLAB that combines Inverse Quadratic Interpolation (IQI) with the Bisection Method to efficiently and reliably solve nonlinear equations.

The algorithm is based on the idea behind Brent's Method, using fast interpolation steps when possible and falling back to bisection when convergence becomes unreliable.

## Features

- Inverse Quadratic Interpolation (IQI)
- Bisection fallback strategy
- Adaptive convergence criteria
- Function and root tolerance controls
- Interval validation and error handling

## Motivation

Traditional root-finding methods often involve a tradeoff between speed and reliability.

- Newton's Method converges quickly but may fail without a good initial guess.
- Bisection is reliable but converges slowly.

This implementation combines interpolation-based acceleration with the robustness of interval-based methods.

## Example

```matlab
f = @(x) x.^3 - x - 2;

Int.a = 1;
Int.b = 2;

params.root_tol = 1e-8;
params.func_tol = 1e-8;

[root, info] = modifiedzeroin3040422362(f, Int, params);
```

## Results

The algorithm successfully locates roots while maintaining numerical stability through interval validation and fallback logic.

## Skills Demonstrated

- Numerical Analysis
- Algorithm Design
- MATLAB Programming
- Scientific Computing
- Computational Mathematics
