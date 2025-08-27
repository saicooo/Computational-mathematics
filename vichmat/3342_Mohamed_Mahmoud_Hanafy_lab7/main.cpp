#include <random>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>

#define MAX_SIZE 10000

const float MIN_NUMBER = -1000.0f;
const float MAX_NUMBER = 1000.0f;


// Main function implementing Gaussian elimination with partial pivoting
int gaussMethod(std::vector<std::vector<float>>& coeffs, std::vector<float>& constants, std::vector<float>& answers, int n) {
    const int min_size = 2;
    
    std::cout << "Gauss method is running..." << std::endl;

    // Validate matrix size
    if (n < min_size) return 2;
    if (n > MAX_SIZE) return 3;

    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Forward elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Find pivot row with maximum element in current column
        float current_max = 0.0f;
        int pivot_row = -1;
        
        for (int k = i; k < n; ++k) {
            float val = fabs(coeffs[k][i]);
            if (val > current_max) {
                current_max = val;
                pivot_row = k;
            }
        }
        
        if (pivot_row < 0) return 1;  // Singular matrix
            
        if (pivot_row != i) {
            // Swap rows for partial pivoting
            coeffs[i].swap(coeffs[pivot_row]);
            std::iter_swap(constants.begin() + i, constants.begin() + pivot_row);
        }
        
        // Eliminate current column below diagonal
        for (int k = i + 1; k < n; ++k) {
            float ratio = coeffs[k][i] / coeffs[i][i];
            for (int j = i; j < n; ++j) {
                coeffs[k][j] -= coeffs[i][j] * ratio;
            }
            constants[k] -= constants[i] * ratio;
        }
    }
    
    // Back substitution
    answers[n - 1] = constants[n - 1] / coeffs[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        float accum = constants[i];
        for (int j = n - 1; j > i; --j) {
            accum -= coeffs[i][j] * answers[j];
        }
        answers[i] = accum / coeffs[i][i];
    }

    // Calculate and display execution time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Gauss method finished working" << std::endl;
    if (duration.count() > 1000) {
        std::cout << "Total execution time: " << duration.count()/1000.0 << " seconds" << std::endl;
    } else {
        std::cout << "Total execution time: " << duration.count() << " ms" << std::endl;
    }
    
    return 0;
}

// Generates random matrix and vector within specified bounds
void randomMatrixGenerator(std::vector<std::vector<float>>& A, std::vector<float>& b, int n) {
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd() ^ static_cast<unsigned>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    ));
    
    std::uniform_real_distribution<float> distrib(MIN_NUMBER, MAX_NUMBER);
    
    std::cout << "Generating matrix..." << std::endl;
    
    // Generate random matrix
    for (int i = 0; i < n; i++) {
        std::vector<float> row;
        for (int j = 0; j < n; j++) {
            row.push_back(distrib(gen));
        }
        A.push_back(row);
    }
    
    // Generate random vector
    for (int j = 0; j < n; j++) {
        b.push_back(distrib(gen));
    }
    
    std::cout << "Generation complete!" << std::endl;
}

// Computes matrix condition number using SVD
float computeConditionNumber(const Eigen::MatrixXf& eigenMatrix) {
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(eigenMatrix);
    
    float maxSingularValue = svd.singularValues()(0);
    float minSingularValue = svd.singularValues()(svd.singularValues().size() - 1);

    if (minSingularValue == 0) {
        throw std::runtime_error("Matrix is singular");
    }
    return maxSingularValue / minSingularValue;
}

int main() {
    std::ofstream outputFile("../result.txt");

    // Get matrix size from user
    int n;
    std::cout << "Input matrix size: ";
    std::cin >> n;
    
    std::vector<std::vector<float>> A;
    std::vector<float> b;
    std::vector<float> X(n);

    randomMatrixGenerator(A, b, n);
    
    // Write matrix to file
    outputFile << "A = [ ";
    for (auto row: A) {
        for (auto x: row) {
            outputFile << x << " ";
        }
        outputFile << ";\n";
    }
    outputFile << "];\n";

    // Write vector to file
    outputFile << "b = [ ";
    for (auto x: b) {
        outputFile << x << "; ";
    }
    outputFile << "];\n";

    // Solve using Gaussian elimination
    int result = gaussMethod(A, b, X, n);

    // Prepare data for Eigen solver
    Eigen::MatrixXf A_Eigen(n, n);
    Eigen::VectorXf b_Eigen(n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_Eigen(i, j) = A[i][j];
        }
        b_Eigen(i) = b[i];
    }

    // Solve using Eigen library for comparison
    std::cout << "Eigen solver is running..." << std::endl;
    Eigen::VectorXf result_Eigen = A_Eigen.colPivHouseholderQr().solve(b_Eigen);
    std::cout << "Eigen solver finished working" << std::endl;

    // Compare results
    if (result == 0) {
        std::cout << "\nResult comparison (Gauss_method | Eigen):" << std::endl;
        for (int i = 0; i < std::min(10, static_cast<int>(X.size())); i++) {
            std::cout << "x[" << i << "]: " << X[i] << " | " << result_Eigen[i] << std::endl;
        }
        if (X.size() > 10) {
            std::cout << "... (showing first 10 of " << X.size() << " elements), results were written in file" << std::endl;
        }
        
        // Calculate difference between solutions
        Eigen::VectorXf diff = Eigen::Map<Eigen::VectorXf>(X.data(), X.size()) - result_Eigen;
        std::cout << "Norm of difference: " << diff.norm() << std::endl;
        
        // Compute and display condition number
        std::cout << "Computing condition number..." << std::endl;
        std::cout << "Condition number: " << computeConditionNumber(A_Eigen) << std::endl;
    }
    else {
        std::cout << "Input error. Code: " << result << std::endl;
    }

    outputFile.close();
    return 0;
}