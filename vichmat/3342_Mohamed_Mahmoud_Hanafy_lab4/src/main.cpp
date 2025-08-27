#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

// Rounding function
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

// Chord method (secant method) to find the root of f(x) in the interval [a, b]
double horda(double a, double b, double eps, int &k, double delta)
{
    double e = fabs(eps) * 2.0; // Tolerance for the interval
    double fa = f(a, delta);    // f(a)
    double fb = f(b, delta);    // f(b)
    double x, fx;

    if (fa * fb > 0.0) // Check if the interval is valid (f(a) and f(b) must have opposite signs)
    {
        puts("Invalid interval specified\n");
        exit(1);
    }

    if (eps <= 0.0) // Check for invalid precision
    {
        puts("Invalid precision specified\n");
        exit(1);
    }

    k = 0; // Initialize iteration counter

    if (fa == 0.0) // If f(a) is already the root
    {
        return a;
    }

    if (fb == 0.0) // If f(b) is already the root
    {
        return b;
    }

    do
    {
        x = a - (b - a) * fa / (fb - fa); // Calculate the next approximation using the chord method
        fx = f(x, delta);                 // Evaluate f(x)

        if (fabs(fx) == 0) // If the root is found
        {
            return x;
        }

        if (fx * fa < 0.0) // If the root is in the left subinterval
        {
            b = x;
            fb = fx;
        }
        else // If the root is in the right subinterval
        {
            a = x;
            fa = fx;
        }

        k++; // Increment iteration counter
    } while (fabs(fx) >= eps); // Continue until the desired precision is achieved

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
    double a, b, eps, x;

    a = 0.0;
    b = 1.0; // Interval [0, 1], since arccos(sqrt(x)) is defined for x ∈ [0, 1]
    printf("eps\t\tdelta\t\ta\t\tb\t\tx\t\tn_m\t\tn\t\tk\tC\n");

    for (double delta = 0.1; delta >= 0.000001; delta /= 10) // Vary delta
    {
        for (eps = 0.1; eps >= 0.000001; eps /= 10) // Vary eps
        {
            x = horda(a, b, eps, k, delta); // Find the root using the chord method
            double n_m = eps / delta;         // Ratio of eps to delta
            double n = 1 / f1(x);             // Inverse of the derivative f'(x)
            int C = (eps / delta > n) ? 1 : 0; // Condition C

            printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
                   eps, delta, a, b, x, n_m, n, k, C);
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

    double a = 0.0; // Start of interval
    double b = 1.0; // End of interval
    int k;          // Number of iterations

    // Find the root using the chord method
    double result = horda(a, b, eps, k, delta);

    // Output the result
    cout << fixed << setprecision(6);
    cout << "x = " << result << "\n";
    if (additional_information) // Display additional info if requested
    {
        cout << "eps: " << eps << "\n";
        cout << "delta: " << delta << "\n";
        cout << "number of iterations: " << k << "\n";
    }

    return 0;
}