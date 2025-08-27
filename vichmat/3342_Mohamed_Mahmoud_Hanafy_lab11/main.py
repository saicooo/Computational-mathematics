from math import fabs, sin, log, pi


def calculate_integrand(x):
    if x < LOWER_LIMIT:
        print(f"Boo<{x}")
    if x > UPPER_LIMIT:
        print(f"Foo>{x}")
    return log(sin(x))


def trapezoid_inner(h, n, x0, func):
    def x_point(i):
        return x0 + i * h

    total = 0
    for i in range(n):
        total += func(x_point(i)) + func(x_point(i + 1))
    total *= h / 
    return total


def trapezoid(a, b, epsilon, func):
    n = 1
    k = 2
    e0 = epsilon * ((1 << k) - 1)
    difference = e0 + 1

    step = (b - a) / n
    integral_h2 = trapezoid_inner(step, n, a, func)

    while difference >= e0:
        n *= 2
        integral_h = integral_h2
        step = (b - a) / n
        integral_h2 = trapezoid_inner(step, n, a, func)
        difference = fabs(integral_h2 - integral_h)

    return (integral_h2 - (integral_h2 - integral_h) / 3), n


def rectangle_inner(h, n, x0, func):
    x = x0
    total = 0
    for _ in range(n):
        total += func(x + h / 2)
        x += h
    total *= h
    return total


def rectangle(a, b, epsilon, func):
    n = 1
    k = 2
    e0 = epsilon * ((1 << k) - 1)
    difference = e0 + 1
    
    step = (b - a) / n
    integral_h2 = rectangle_inner(step, n, a, func)
    
    while difference >= e0:
        n *= 2
        integral_h = integral_h2
        step = (b - a) / n
        integral_h2 = rectangle_inner(step, n, a, func)
        difference = fabs(integral_h2 - integral_h)
    
    return (integral_h2 + (integral_h2 - integral_h) / 3), n


def simpson_inner(h, n, x0, func):
    def x_point(i):
        return x0 + i * h
    
    m = n // 2
    total = 0
    for i in range(m):
        total += (func(x_point(2 * i)) + 
                  4 * func(x_point(2 * i + 1)) + 
                  func(x_point(2 * i + 2)))
    total *= h / 3
    return total


def simpson(a, b, epsilon, func, echo=False):
    n = 1
    k = 4
    e0 = epsilon * ((1 << k) - 1)
    difference = e0 + 1
    
    step = (b - a) / n
    integral_h2 = simpson_inner(step, n, a, func)
    
    while difference >= e0:
        n *= 2
        integral_h = integral_h2
        step = (b - a) / n
        integral_h2 = simpson_inner(step, n, a, func)
        difference = fabs(integral_h2 - integral_h)
    
    return (integral_h2 - (integral_h2 - integral_h) / 15), n


if __name__ == "__main__":
    LOWER_LIMIT = pi/4
    UPPER_LIMIT = pi/2
    
    print("Epsilon\t\tRectangle\tTrapezoid\tSimpson\t\tRect_N\tTrap_N\tSimp_N")
    epsilon = 0.1
    while epsilon > 0.000001:
        rect_result, rect_n = rectangle(LOWER_LIMIT, UPPER_LIMIT, epsilon, calculate_integrand)
        trap_result, trap_n = trapezoid(LOWER_LIMIT, UPPER_LIMIT, epsilon, calculate_integrand)
        simp_result, simp_n = simpson(LOWER_LIMIT, UPPER_LIMIT, epsilon, calculate_integrand)
        print(f"{epsilon:.6f}\t{rect_result:.6f}\t{trap_result:.6f}\t{simp_result:.6f}\t"
              f"{rect_n}\t{trap_n}\t{simp_n}")
        epsilon /= 10