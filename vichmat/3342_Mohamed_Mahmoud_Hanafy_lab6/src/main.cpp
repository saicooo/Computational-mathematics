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

// Iteration function phi(x) = x - alpha * f(x)
double phi(double x, double delta)
{
    double alpha = 0.25; // Relaxation parameter
    return round(x - alpha * (f(x, delta)), delta);
}

// Simple iteration method to find the root of f(x)
double iter(double x0, double eps, int &k, double delta)
{
    double x = x0; // Initial guess
    double x_new;  // New approximation
    k = 0;        // Iteration counter

    do
    {
        x_new = phi(x, delta); // Compute the next approximation
        k++;                  // Increment iteration counter
        if (abs(x_new - x) < eps) // Check for convergence
        {
            break;
        }
        x = x_new; // Update the current approximation
    } while (true);

    return x_new; // Return the approximate root
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
    x_0 = 0.25;    // Initial guess for the iteration method

    double x = 0.154406; // Known root for comparison

    cout << "eps\t\tdelta\t\ta\t\tb\t\tx\t\tdx\t\tk\tC\n";

    for (double delta = 0.1; delta >= 0.000001; delta /= 10) // Vary delta
    {
        for (eps = 0.1; eps >= 0.000001; eps /= 10) // Vary eps
        {
            double result = iter(x_0, eps, k, delta); // Find the root using the iteration method
            double dx = abs(x - result);              // Difference between the result and the known root
            int C = (eps >= dx) ? 1 : 0;              // Condition C: 1 if eps >= dx, else 0

            printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
                   eps, delta, a, b, result, dx, k, C);
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
    int k;              // Number of iterations

    // Find the root using the iteration method
    double result = iter(b, eps, k, delta);

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