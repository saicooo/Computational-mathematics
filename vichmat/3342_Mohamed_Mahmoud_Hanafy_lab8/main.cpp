#include <random>
#include <fstream>
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <algorithm>

const float MIN_NUMBER = -1000.0f;
const float MAX_NUMBER = 1000.0f;
const int N_MAX = 1000000;  // Changed from #define to const

// Implements Jacobi iterative method for solving linear systems
int jacobiMethod(Eigen::MatrixXf& C, Eigen::VectorXf& d, Eigen::VectorXf& X, int n, int precision) {
    // Validate input parameters
    if (C.norm() >= 1) {
        std::cerr << "Error: Matrix norm must be less than 1 for convergence" << std::endl;
        return 1;
    }
    else if (n < 2) {
        std::cerr << "Error: Matrix size too small (minimum 2)" << std::endl;
        return 2;
    }
    else if (n > N_MAX) {
        std::cerr << "Error: Matrix size exceeds maximum allowed (" << N_MAX << ")" << std::endl;
        return 3;
    }
    
    std::cout << "Jacobi method started..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    Eigen::VectorXf X_0(n);
    const float eps = 1.0f / pow(10, precision);
    std::cout << "Target precision: " << eps << std::endl;
    
    // Initialize with starting guess
    X_0 = X;
    Eigen::VectorXf X_1 = C * X_0 + d;

    int iterations = 1;
    const float convergence_factor = C.norm() / (1 - C.norm());
    
    // Main iteration loop
    while (convergence_factor * (X_1 - X_0).norm() > eps) {
        X_0 = X_1;
        X_1 = C * X_0 + d;
        iterations++;
        
        }
    

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Jacobi method completed successfully" << std::endl;
    std::cout << "Total iterations: " << iterations << std::endl;
    if (duration.count() > 1000) {
        std::cout << "Total execution time: " << duration.count()/1000.0 << " seconds" << std::endl;
    } else {
        std::cout << "Total execution time: " << duration.count() << " ms" << std::endl;
    }

    X = X_1;
    return 0;
}

// Generates diagonally dominant random matrix and vectors
void generateRandomMatrix(std::vector<std::vector<float>>& A, std::vector<float>& b, std::vector<float>& X, int n) {
    std::cout << "Generating random matrix of size " << n << "..." << std::endl;
    
    std::random_device rd;
    std::mt19937 gen(rd() ^ static_cast<unsigned>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()
    ));
    
    std::uniform_real_distribution<float> distrib(-1, 1);
    
    for (int i = 0; i < n; i++) {
        std::vector<float> row;
        float sum = 0;
        for (int j = 0; j < n; j++) {
            row.push_back(distrib(gen));
            sum += fabs(row[j]);
        }
        A.push_back(row);
        // Ensure diagonal dominance
        A[i][i] = sum * 1.5f;
    }

    for (int j = 0; j < n; j++) {
        b[j] = distrib(gen);
        X[j] = distrib(gen);
    }
    
    std::cout << "Matrix generation complete" << std::endl;
}

// Computes condition number using SVD
float computeConditionNumber(const Eigen::MatrixXf& eigenMatrix) {
    try {
        Eigen::JacobiSVD<Eigen::MatrixXf> svd(eigenMatrix);
        
        float maxSingularValue = svd.singularValues()(0);
        float minSingularValue = svd.singularValues()(svd.singularValues().size() - 1);

        if (minSingularValue == 0) {
            throw std::runtime_error("Matrix is singular (condition number is infinite)");
        }
        return maxSingularValue / minSingularValue;
    } catch (const std::exception& e) {
        std::cerr << "Error in condition number computation: " << e.what() << std::endl;
        return std::numeric_limits<float>::infinity();
    }
}

// Solves linear system using Eigen's built-in solver
Eigen::VectorXf solveLinearSystem(const Eigen::MatrixXf& eigenA, const Eigen::VectorXf& eigenB) {
    std::cout << "Solving system using Eigen..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    Eigen::VectorXf x = eigenA.colPivHouseholderQr().solve(eigenB);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "Eigen solver completed in " << duration.count() << " ms" << std::endl;
    return x;
}

int main() {
    std::ofstream outFile("output.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file" << std::endl;
        return 1;
    }

    // Get input parameters
    int n, precision; 
    std::cout << "Enter matrix size (n) and precision: ";
    std::cin >> n >> precision;
    
    std::vector<std::vector<float>> A;
    std::vector<float> b(n);
    std::vector<float> X(n);

    // Generate and prepare matrices
    generateRandomMatrix(A, b, X, n);

    Eigen::MatrixXf eigenA(n, n);
    Eigen::MatrixXf eigenC(n, n);
    Eigen::VectorXf eigenB(n);
    Eigen::VectorXf eigenD(n);
    Eigen::VectorXf eigenX(n);

    // Convert to Eigen format and prepare Jacobi matrices
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            eigenA(i, j) = A[i][j];
            eigenC(i, j) = (i != j) ? -A[i][j] / A[i][i] : 0;
        }
        eigenB(i) = b[i];
        eigenD(i) = b[i] / A[i][i];
        eigenX(i) = X[i];
    }
    
    // Solve using Jacobi method
    int result = jacobiMethod(eigenC, eigenD, eigenX, n, precision);
    
    if (result == 0) {
        // Compare with direct solution
        Eigen::VectorXf results_real = solveLinearSystem(eigenA, eigenB);
        
        std::cout << "\nSolution comparison (Jacobi | Eigen):" << std::endl;
        for (int i = 0; i < std::min(10, n); i++) {
            std::cout << "x[" << i << "]: " << eigenX(i) << " | " << results_real[i] << std::endl;
        }
        if (n > 10) {
            std::cout << "... (showing first 10 of " << n << " elements)" << std::endl;
        }
        
        Eigen::VectorXf diff = eigenX - results_real;
        std::cout << "\nNorm of difference: " << diff.norm() << std::endl;
        std::cout << "Computing condition number..." << std::endl;
        std::cout << "Condition number: " << computeConditionNumber(eigenA) << std::endl;
    }
    else {
        std::cerr << "Error in Jacobi method (code: " << result << ")" << std::endl;
    }

    outFile.close();
    return 0;
}