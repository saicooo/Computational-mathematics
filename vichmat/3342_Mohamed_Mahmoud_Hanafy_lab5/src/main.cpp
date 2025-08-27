#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

// Rounding function: rounds 'x' to the nearest multiple of 'delta'
double round(double x, double delta)
{
    if (delta <= 1E-9) // Check for invalid precision
    {
        puts("Invalid rounding precision specified\n");
        exit(1);
    }

    if (x > 0.0) // Round positive numbers
    {
        return delta * long(x / delta + 0.5);
    }
    else // Round negative numbers
    {
        return delta * long(x / delta - 0.5);
    }
}

// Function f(x) = e^x - arccos(sqrt(x))
double f(double x, double delta)
{
    if (x < 0.0 || x > 1.0) // Check if x is within the domain [0, 1]
    {
        return NAN;
    }

    double term1 = exp(x);                // e^x
    double term2 = acos(sqrt(x));         // arccos(sqrt(x))
    double result = term1 - term2;        // f(x) = e^x - arccos(sqrt(x))

    return round(result, delta); // Round the result to the specified precision
}

// First derivative of f(x): f'(x) = e^x + 1 / (2 * sqrt(x) * sqrt(1 - x))
double f1(double x)
{
    if (x <= 0.0 || x >= 1.0) // Check if x is within the domain (0, 1)
    {
        return NAN;
    }

    double term1 = exp(x);                              // e^x
    double term2 = 1.0 / (2.0 * sqrt(x) * sqrt(1.0 - x)); // 1 / (2 * sqrt(x) * sqrt(1 - x))
    double result = term1 + term2;                       // f'(x) = e^x + term2

    return result;
}

// Newton's method to find the root of f(x)
double newton(double x, double eps, int &k, double delta)
{
    double y, y1, dx, eps0;
    k = 0; // Initialize iteration counter

    double m1 = 2.44; // Lower bound for |f'(x)|
    double M2 = 6.3;  // Upper bound for |f''(x)|

    eps0 = sqrt(2 * m1 * eps / M2); // Tolerance for stopping condition

    do
    {
        y = f(x, delta); // Evaluate f(x)

        if (fabs(y) <= eps) // If f(x) is close enough to zero
        {
            return x;
        }

        y1 = f1(x); // Evaluate f'(x)

        if (fabs(y1) <= eps) // If the derivative is too small
        {
            puts("The derivative becomes zero\n");
            exit(1);
        }

        dx = y / y1; // Calculate the Newton step
        x -= dx;     // Update x
        k++;         // Increment iteration counter
    } while (fabs(dx) > eps0); // Continue until the step size is small enough

    return x; // Return the approximate root
}

// Display help message
void print_help()
{
    cout << "Usage: program [options]\n"
         << "Options:\n"
         << "  -i          Display the number of iterations\n"
         << "  -d delta    Set the value of delta (default 0.0001)\n"
         << "  -e EPS      Set the value of eps (default 0.0001)\n"
         << "  -t          Test table with various eps and delta\n"
         << "  -h          Show this message and exit\n";
}

// Generate a test table with various eps and delta values
void test_table()
{
    int k;
    double a, b, eps, x_N, x_0;

    a = 0.1;       // Start of interval
    b = 0.25;      // End of interval
    x_0 = 0.1;     // Initial guess for Newton's method

    double x = 0.154406; // Known root for comparison

    printf("eps\t\tdelta\t\ta\t\tb\t\tx\t\tn_m\t\tn\t\tk\tC\n");

    for (double delta = 0.1; delta >= 0.000001; delta /= 10) // Vary delta
    {
        for (eps = 0.1; eps >= 0.000001; eps /= 10) // Vary eps
        {
            x_N = newton(x_0, eps, k, delta); // Find the root using Newton's method
            printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
                   eps, delta, a, b, x_N, abs(x - x_N), k, eps >= abs(x - x_N));
        }
    }
}

int main(int argc, char *argv[])
{
    double delta = 0.0001; // Default delta value
    double eps = 0.0001;   // Default eps value
    bool additional_information = false; // Flag to display additional info

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i)
    {
        string arg = argv[i];
        if (arg == "-t") // Generate test table
        {
            test_table();
            return 0;
        }
        else if (arg == "-i") // Display number of iterations
        {
            additional_information = true;
        }
        else if (arg == "-d") // Set delta value
        {
            if (i + 1 < argc)
            {
                delta = atof(argv[++i]);
            }
            else
            {
                cerr << "Error: -d requires a delta value\n";
                return 1;
            }
        }
        else if (arg == "-e") // Set eps value
        {
            if (i + 1 < argc)
            {
                eps = atof(argv[++i]);
            }
            else
            {
                cerr << "Error: -e requires an eps value\n";
                return 1;
            }
        }
        else if (arg == "-h") // Display help message
        {
            print_help();
            return 0;
        }
        else // Handle unknown arguments
        {
            cerr << "Unknown argument: " << arg << "\n";
            print_help();
            return 1;
        }
    }

    double a = 0.1;     // Start of interval
    double b = 0.25;    // End of interval
    double x_0 = 0.1;   // Initial guess for Newton's method
    int k;              // Number of iterations

    // Find the root using Newton's method
    double result = newton(x_0, eps, k, delta);

    // Output the result
    cout << fixed << setprecision(6);
    cout << "x = " << result << "\n";
    if (additional_information) // Display additional info if requested
    {
        cout << "eps: " << eps << "\n";
        cout << "delta: " << delta << "\n";
        cout << "number of iterations: " << k << "\n";
        double m1 = 2.44;
        double M2 = 6.3;

        double eps0 = sqrt(2 * m1 * eps / M2); // Calculate eps0
        cout << "eps_0: " << eps0 << "\n";
    }

    return 0;
}