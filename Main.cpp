#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

// Function to calculate mean
double mean(const std::vector<double>& returns) {
    double sum = 0.0;
    for (double r : returns) {
        sum += r;
    }
    return sum / returns.size();
}

// Function to calculate standard deviation
double standardDeviation(const std::vector<double>& returns, double mean) {
    double sum = 0.0;
    for (double r : returns) {
        sum += (r - mean) * (r - mean);
    }
    return std::sqrt(sum / (returns.size() - 1));
}

// Function to perform Monte Carlo simulation
std::vector<double> monteCarloSimulation(double mean, double stdDev, int numSimulations) {
    std::vector<double> simulatedReturns(numSimulations);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(mean, stdDev);

    for (int i = 0; i < numSimulations; ++i) {
        simulatedReturns[i] = d(gen);
    }

    return simulatedReturns;
}

// Function to calculate VaR
double calculateVaR(const std::vector<double>& simulatedReturns, double confidenceLevel) {
    std::vector<double> sortedReturns = simulatedReturns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    int index = static_cast<int>((1.0 - confidenceLevel) * sortedReturns.size());
    return -sortedReturns[index];
}

// Function to calculate CVaR
double calculateCVaR(const std::vector<double>& simulatedReturns, double confidenceLevel) {
    std::vector<double> sortedReturns = simulatedReturns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    int index = static_cast<int>((1.0 - confidenceLevel) * sortedReturns.size());
    double VaR = -sortedReturns[index];

    // Calculate the average of all losses exceeding VaR
    double sumExceedingVaR = 0.0;
    int countExceedingVaR = 0;
    for (int i = 0; i < sortedReturns.size(); ++i) {
        if (-sortedReturns[i] > VaR) {
            sumExceedingVaR += -sortedReturns[i];
            ++countExceedingVaR;
        }
    }

    return sumExceedingVaR / countExceedingVaR;
}

int main() {
    // Historical returns data (example)
    std::vector<double> historicalReturns = {0.01, -0.02, 0.015, -0.005, 0.02, -0.01, 0.03, -0.015, 0.005, 0.01};

    // Parameters
    double confidenceLevel = 0.99;
    int numSimulations = 10000;

    // Calculate mean and standard deviation of historical returns
    double meanReturn = mean(historicalReturns);
    double stdDevReturn = standardDeviation(historicalReturns, meanReturn);

    // Perform Monte Carlo simulation
    std::vector<double> simulatedReturns = monteCarloSimulation(meanReturn, stdDevReturn, numSimulations);

    // Calculate VaR
    double VaR = calculateVaR(simulatedReturns, confidenceLevel);

    // Calculate CVaR
    double CVaR = calculateCVaR(simulatedReturns, confidenceLevel);

    // Output the result
    std::cout << "Value at Risk (VaR) at " << confidenceLevel * 100 << "% confidence level: " << VaR << std::endl;
    std::cout << "Conditional Value at Risk (CVaR) at " << confidenceLevel * 100 << "% confidence level: " << CVaR << std::endl;

    return 0;
}
