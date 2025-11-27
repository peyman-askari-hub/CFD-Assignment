#include "Headers.h"
#include "ErrorCalculator.h"
#include <chrono>

ResultSeries Laasonen::run() {
    ResultSeries results;

    double r = D * dt / (dx * dx);

    // Initial condition
    std::vector<double> T(N + 1, Tin);
    T[0] = T[N] = Tsur;

    // Tridiagonal coefficients (constant)
    std::vector<double> a(N - 1, -r);
    std::vector<double> b(N - 1, 1 + 2 * r);
    std::vector<double> c(N - 1, -r);

    double nextOutput = 0.0;
    double outputInterval = 0.1;

    // Save initial temperature
    results[0.0] = T;

    // Time stepping (start from 1!)
    for (int n = 1; n <= timeSteps; n++) {
        double t = n * dt;

        // Build RHS vector d
        std::vector<double> d(N - 1);
        for (int i = 1; i < N; i++)
            d[i - 1] = T[i];

        // Apply BCs
        d[0]     += r * Tsur;  // left boundary
        d[N - 2] += r * Tsur;  // right boundary

        // Solve tridiagonal system
        std::vector<double> T_inner = ThomasSolver::solve(a, b, c, d);

        // Update full T vector
        for (int i = 1; i < N; i++)
            T[i] = T_inner[i - 1];

        // Output at intervals
        if (t >= nextOutput - 1e-8) {
            results[t] = T;
            nextOutput += outputInterval;
        }
    }

    return results;
}

void Laasonen::computeErrorForDt(const std::vector<double>& dtList,
                              std::vector<double>& rmsList,
                              std::vector<double>& maxList,
                              std::vector<double>& runtimeList)
{
    rmsList.clear();
    maxList.clear();
    runtimeList.clear();

    for (double dtTest : dtList)
    {
        // Temporary solvers with new dt
        SolverContext ctxDt{L, D, Tin, Tsur, dx, totalTime, dtTest};
        Laasonen solver(ctxDt);
        Analytical analytical(ctxDt);

        auto start = std::chrono::high_resolution_clock::now();
        ResultSeries num = solver.run();
        ResultSeries ana = analytical.run();
        auto end = std::chrono::high_resolution_clock::now();

        double runtime = std::chrono::duration<double,std::milli>(end - start).count();

        double rms, maxErr;
        ErrorCalculator::compute(num, ana, rms, maxErr);

        rmsList.push_back(rms);
        maxList.push_back(maxErr);
        runtimeList.push_back(runtime);
    }
}
