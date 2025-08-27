import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange


def newton_interpolation(x_data, y_data, x):
    n = len(x_data)
    # Create divided differences table
    diff_table = [[0] * n for _ in range(n)]
    for i in range(n):
        diff_table[i][0] = y_data[i]

    # Fill the divided differences table
    for j in range(1, n):
        for i in range(n - j):
            diff_table[i][j] = (diff_table[i + 1][j - 1] - diff_table[i][j - 1]) / (x_data[i + j] - x_data[i])

    # Compute Newton interpolation polynomial
    result = diff_table[0][0]
    product_term = 1.0

    for j in range(1, n):
        product_term *= (x - x_data[j - 1])
        result += diff_table[0][j] * product_term

    return result


def read_data(filename):
    x_data = []
    y_data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                x, y = map(float, line.strip().split())
                x_data.append(x)
                y_data.append(y)
    return np.array(x_data), np.array(y_data)


try:
    # Read data from file
    filename = 'data.txt'
    x_data, y_data = read_data(filename)

    if len(x_data) == 0:
        raise ValueError("File contains no data")

    # Point for interpolation
    x_interp = 2.7344
    y_interp = newton_interpolation(x_data, y_data, x_interp)

    print(f"Interpolated value at x = {x_interp}: y = {y_interp:.4f}")

    # Compare with library function
    poly = lagrange(x_data, y_data)
    print(f'Library function result: {poly(x_interp)}')

    # Create points for smooth interpolation plot
    x_smooth = np.linspace(min(x_data), max(x_data), 500)
    y_smooth = [newton_interpolation(x_data, y_data, x) for x in x_smooth]

    # Plotting - vertical arrangement
    plt.figure(figsize=(10, 10))  # Adjusted figure size for vertical layout

    # 1. Plot of points connected by straight lines (top plot)
    plt.subplot(2, 1, 1)
    plt.plot(x_data, y_data, 'o-')
    plt.scatter([x_interp], [y_interp], color='red', label=f'x = {x_interp}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Original function points')
    plt.legend()
    plt.grid()

    # 2. Plot of Newton interpolation polynomial (bottom plot)
    plt.subplot(2, 1, 2)
    plt.plot(x_smooth, y_smooth, label='Newton polynomial')
    plt.scatter(x_data, y_data, color='blue', label='Original points')
    plt.scatter([x_interp], [y_interp], color='red', label=f'x = {x_interp}, y = {y_interp:.4f}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Newton polynomial approximation')
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()

except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
except Exception as e:
    print(f"An error occurred: {str(e)}")