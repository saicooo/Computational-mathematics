from math import sqrt, exp
from scipy.integrate import simpson
import numpy as np


def coordinate_transform(t_val, lower_bound, upper_bound):
    """Map from [-1, 1] to [lower_bound, upper_bound]"""
    return 0.5 * (upper_bound - lower_bound) * t_val + 0.5 * (lower_bound + upper_bound)


def gaussian_quadrature_8point(func, lower_bound, upper_bound):
    """
    Computes definite integral using 8-point Gaussian quadrature
    between lower_bound and upper_bound
    """
    # Precomputed nodes and weights for 8-point Gauss-Legendre quadrature
    quadrature_nodes = [-0.96028986, -0.79666648, -0.52553242, -0.18343464,
                        0.18343464, 0.52553242, 0.79666648, 0.96028986]

    quadrature_weights = [0.10122854, 0.22238103, 0.31370664, 0.36268378,
                          0.36268378, 0.31370664, 0.22238103, 0.10122854]

    integral_approximation = 0.0
    for node, weight in zip(quadrature_nodes, quadrature_weights):
        mapped_point = coordinate_transform(node, lower_bound, upper_bound)
        integral_approximation += weight * func(mapped_point)

    integral_approximation *= 0.5 * (upper_bound - lower_bound)
    return integral_approximation


def simpson_integration(func, lower_bound, upper_bound, n=1000):
    """
    Computes definite integral using Simpson's rule
    between lower_bound and upper_bound with n intervals
    """
    x = np.linspace(lower_bound, upper_bound, n + 1)
    y = [func(xi) for xi in x]
    return simpson(y=y, x=x)


def integrand_function(x_val):
    """
    Function to be integrated: (1/sqrt(1+x)) * exp(-x)
    """
    return (1 / sqrt(1 + x_val)) * exp(-x_val)


if __name__ == "__main__":
    integration_start = 0.0  # Lower integration limit
    integration_end = 1.0  # Upper integration limit

    # 8-point Gaussian quadrature
    custom_result = gaussian_quadrature_8point(integrand_function,
                                               integration_start,
                                               integration_end)

    # Simpson's method
    simpson_result = simpson_integration(integrand_function,
                                         integration_start,
                                         integration_end)

    print(f"Numerical integration of (1/sqrt(1+x))*exp(-x) from {integration_start} to {integration_end}:")
    print(f"8-point Gauss quadrature result: {custom_result:.10f}")
    print(f"Scipy simpson's rule result (n=1000): {simpson_result:.10f}")
    print("\nAbsolute difference:")
    print(f"Gauss vs scipy.simpson: {abs(custom_result - simpson_result):.2e}")
