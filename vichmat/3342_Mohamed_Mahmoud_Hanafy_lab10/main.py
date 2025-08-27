import numpy as np
import matplotlib.pyplot as plt


def target_function(x):
    return np.exp(x) - np.arccos(np.sqrt(x))


def get_chebyshev_nodes(start, end, num_nodes):
    k_indices = np.arange(1, num_nodes + 1)
    base_nodes = np.cos((2 * k_indices - 1) * np.pi / (2 * num_nodes))
    scaled_nodes = 0.5 * (start + end) + 0.5 * (end - start) * base_nodes
    return scaled_nodes


def compute_divided_differences(x_values, y_values):
    node_count = len(x_values)
    difference_table = np.zeros([node_count, node_count])
    difference_table[:, 0] = y_values

    for column in range(1, node_count):
        for row in range(node_count - column):
            numerator = difference_table[row + 1][column - 1] - difference_table[row][column - 1]
            denominator = x_values[row + column] - x_values[row]
            difference_table[row][column] = numerator / denominator

    return difference_table[0, :]


def evaluate_newton_polynomial(x_nodes, y_nodes, evaluation_points):
    node_count = len(x_nodes)
    coefficients = compute_divided_differences(x_nodes, y_nodes)
    result_values = np.zeros_like(evaluation_points, dtype=np.float64)

    for term_idx in range(node_count):
        term_value = np.ones_like(evaluation_points, dtype=np.float64) * coefficients[term_idx]
        for j in range(term_idx):
            term_value *= (evaluation_points - x_nodes[j])
        result_values += term_value

    return result_values


def handle_zoom_event(event):
    """Graph zoom event handler"""
    current_axes = event.inaxes
    if current_axes is None:
        return

    zoom_factor = 1.2 if event.button == 'up' else 1 / 1.2

    x_limits = current_axes.get_xlim()
    y_limits = current_axes.get_ylim()

    cursor_x, cursor_y = event.xdata, event.ydata
    if cursor_x is None or cursor_y is None:
        return

    new_x_limits = [
        cursor_x - (cursor_x - x_limits[0]) * zoom_factor,
        cursor_x + (x_limits[1] - cursor_x) * zoom_factor
    ]
    new_y_limits = [
        cursor_y - (cursor_y - y_limits[0]) * zoom_factor,
        cursor_y + (y_limits[1] - cursor_y) * zoom_factor
    ]

    current_axes.set_xlim(new_x_limits)
    current_axes.set_ylim(new_y_limits)
    plt.draw()


def visualize_interpolation(start, end, node_count):
    # Uniformly distributed nodes
    uniform_x = np.linspace(start, end, node_count)
    uniform_y = target_function(uniform_x)

    # Chebyshev nodes
    cheb_x = get_chebyshev_nodes(start, end, node_count)
    cheb_y = target_function(cheb_x)

    # Evaluation points
    plot_x = np.linspace(start, end, 10000)
    true_y = target_function(plot_x)

    # Interpolation polynomials
    interp_uniform = evaluate_newton_polynomial(uniform_x, uniform_y, plot_x)
    interp_cheb = evaluate_newton_polynomial(cheb_x, cheb_y, plot_x)

    # Create figure
    figure, axis = plt.subplots(figsize=(12, 8))
    figure.canvas.mpl_connect('scroll_event', handle_zoom_event)

    # Plot data
    axis.plot(plot_x, true_y, color='#000000', label='Original function', linewidth=2)
    axis.plot(plot_x, interp_uniform, '--', color='#32CD32',
             label=f'Uniform node interpolation (n={node_count})')
    axis.plot(plot_x, interp_cheb, '-.', color='#8A2BE2',
             label=f'Chebyshev node interpolation (n={node_count})')

    axis.scatter(uniform_x, uniform_y, color='#32CD32', label='Uniform nodes')
    axis.scatter(cheb_x, cheb_y, color='#8A2BE2', label='Chebyshev nodes')

    axis.set_title(f'Interpolation of eˣ - arccos(√x) on [{start}, {end}]')
    axis.set_xlabel('x')
    axis.set_ylabel('f(x)')
    axis.legend()
    axis.grid(True)

    # Coordinate axes
    axis.axhline(0, color='black', linewidth=1.0, linestyle='-', alpha=0.7)
    axis.axvline(0, color='black', linewidth=1.0, linestyle='-', alpha=0.7)

    plt.show()


def get_valid_interval():
    while True:
        start = float(input("Enter interval start (0 ≤ a ≤ 1): "))
        end = float(input("Enter interval end (0 ≤ b ≤ 1): "))
        if 0 <= start <= 1 and 0 <= end <= 1:
            return start, end
        print("Error: interval must be within [0, 1]")


if __name__ == "__main__":
    print("Function interpolation: eˣ - arccos(√x)")
    interval_start, interval_end = get_valid_interval()
    nodes = int(input("Enter number of interpolation nodes: "))
    visualize_interpolation(interval_start, interval_end, nodes)