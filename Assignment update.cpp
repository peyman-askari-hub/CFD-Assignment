#include <iostream>       // For standard input/output
#include <vector>         // For std::vector
#include <cmath>          // For math functions like sin, exp, pow, fabs
#include <iomanip>        // For controlling output formatting (precision, setw)
#include <fstream>        // For file input/output (CSV files)
#include <chrono>         // For measuring runtime of code
using namespace std;
using namespace std::chrono; // To easily use high_resolution_clock and duration

// Class to solve 1D heat conduction problem using multiple methods
class HeatEquation {
private:
    const double pi = 3.14159265358979323846; // Constant pi
    double L, D, Tin, Tsur, dx, dt, totalTime; // Geometry, material, initial/boundary conditions
    int N, timeSteps;                          // Number of spatial points and time steps
    vector<double> x;                           // Spatial grid
    int M_max;                                  // Maximum number of terms for analytical series
    double outputInterval;                      // Interval to print/save results

    // Thomas algorithm to solve tridiagonal matrix Ax = d
    vector<double> thomasSolver(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
        int n = b.size();
        vector<double> c_star(n, 0.0), d_star(n, 0.0), T(n, 0.0); // Temporary arrays for elimination

        c_star[0] = c[0] / b[0];  // Forward sweep
        d_star[0] = d[0] / b[0];  // Forward sweep

        for (int i = 1; i < n; i++) {
            double m = b[i] - a[i] * c_star[i - 1]; // Denominator
            c_star[i] = c[i] / m;                   // Modified upper diagonal
            d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m; // Modified RHS
        }

        T[n - 1] = d_star[n - 1]; // Back substitution
        for (int i = n - 2; i >= 0; i--)
            T[i] = d_star[i] - c_star[i] * T[i + 1];

        return T; // Return solution vector
    }

    // Calculate RMS and Maximum absolute error
    vector<double> calculateErrorVec(const vector<double>& analytical, const vector<double>& numerical) {
        double meanSquareError = 0.0; // RMS (L2 norm)
        double maxAbsoluteError = 0.0; // Max absolute error (L-infinity norm)

        for (int i = 0; i < analytical.size(); i++) {
            double error1 = fabs(analytical[i] - numerical[i]); // absolute error at point i
            meanSquareError += error1 * error1;
            maxAbsoluteError = max(maxAbsoluteError, error1);
        }
        meanSquareError = sqrt(meanSquareError / analytical.size()); // RMS normalization
        return { meanSquareError, maxAbsoluteError };
    }

    // Save solution to CSV file
    void saveToCSV(const vector<double>& T, double t, const string& methodName) {
        string filename = methodName + ".csv";
        ofstream fout(filename, ios::app); // append mode
        if (!fout.is_open()) {
            cerr << "Error opening " << filename << endl;
            return;
        }
        for (int i = 0; i <= N; i++)
            fout << t << "," << x[i] << "," << T[i] << "\n"; // write time, x, T
        fout.close();
    }

    // Print results to console
    void printResult(const vector<double>& T, double t, const string& methodName, int stepInterval = 10) {
        cout << "===== " << methodName << " | t = " << t << " hr =====\n";
        for (int i = 0; i <= N; i += stepInterval) { // Print every stepInterval point
            cout << "x=" << fixed << setprecision(2) << x[i]
                << " cm: T=" << fixed << setprecision(2) << T[i] << "C\n";
        }
        cout << "========================================\n\n";
    }

public:
    // Constructor: initialize all parameters
    HeatEquation(double L_, double D_, double Tin_, double Tsur_, double dx_, double dt_, double totalTime_, int M_max_ = 200, double outputInterval_ = 0.1)
        : L(L_), D(D_), Tin(Tin_), Tsur(Tsur_), dx(dx_), dt(dt_), totalTime(totalTime_), M_max(M_max_), outputInterval(outputInterval_) {
        N = static_cast<int>(L / dx); // number of spatial points
        timeSteps = static_cast<int>(totalTime / dt); // number of time steps
        x.resize(N + 1);
        for (int i = 0; i <= N; i++) x[i] = i * dx; // spatial grid
    }

    // Analytical solution using Fourier series
    vector<double> analyticalSolution(double t) {
        vector<double> T(N + 1, 0.0);
        for (int i = 0; i <= N; i++) {
            double sum = 0.0;
            for (int m = 1; m <= M_max; m++) {
                double term = (1 - pow(-1.0, m)) / (m * pi) * sin(m * pi * x[i] / L) * exp(-D * pow(m * pi / L, 2) * t);
                sum += term;
            }
            T[i] = Tsur + 2 * (Tin - Tsur) * sum;
        }
        T[0] = T[N] = Tsur; // enforce boundary conditions
        return T;
    }

    // Print analytical solution for all time steps
    void printAnalytical() {
        double nextOutput = 0.0;
        ofstream("analytical.csv", ios::trunc).close(); // clear file before writing
        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) { // output interval check
                vector<double> T = analyticalSolution(t);
                printResult(T, t, "Analytical");
                saveToCSV(T, t, "analytical");
                nextOutput += outputInterval;
            }
        }
    }

    // DuFort-Frankel explicit method
    void duFortFrankel() {
        double r = D * dt / (dx * dx);
        vector<double> T_prev(N + 1, Tin), T_curr(N + 1, Tin), T_next(N + 1, Tin);
        T_prev[0] = T_prev[N] = T_curr[0] = T_curr[N] = Tsur; // BCs
        double nextOutput = 0.0;
        ofstream("dufort.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T_curr, t, "DuFort-Frankel");
                saveToCSV(T_curr, t, "dufort");
                nextOutput += outputInterval;
            }

            if (n == 0) { // first step: forward Euler
                for (int i = 1; i < N; i++)
                    T_next[i] = T_curr[i] + r * (T_curr[i + 1] - 2 * T_curr[i] + T_curr[i - 1]);
            }
            else { // subsequent steps: DuFort-Frankel formula
                for (int i = 1; i < N; i++)
                    T_next[i] = ((1 - 2 * r) * T_prev[i] + 2 * r * (T_curr[i + 1] + T_curr[i - 1])) / (1 + 2 * r);
            }

            T_prev = T_curr;
            T_curr = T_next;
        }
    }

    // Richardson method (explicit, central difference in time)
    void richardsonMethod() {
        double r = D * dt / (dx * dx);
        vector<double> T_prev(N + 1, Tin), T_curr(N + 1, Tin);
        T_prev[0] = T_prev[N] = T_curr[0] = T_curr[N] = Tsur; // BCs
        double nextOutput = 0.0;
        ofstream("richardson.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T_curr, t, "Richardson");
                saveToCSV(T_curr, t, "richardson");
                nextOutput += outputInterval;
            }

            if (n == 0) { // first step: forward Euler
                for (int i = 1; i < N; i++)
                    T_curr[i] = T_prev[i] + r * (T_prev[i + 1] - 2 * T_prev[i] + T_prev[i - 1]);
            }
            else {
                vector<double> T_next(N + 1, Tin);
                T_next[0] = T_next[N] = Tsur; // BCs
                for (int i = 1; i < N; i++)
                    T_next[i] = T_prev[i] + 2 * r * (T_curr[i + 1] - 2 * T_curr[i] + T_curr[i - 1]);
                T_prev = T_curr;
                T_curr = T_next;
            }
        }
    }

    // Laasonen method (implicit, backward Euler)
    void laasonenMethod() {
        double r = D * dt / (dx * dx);
        vector<double> T(N + 1, Tin);
        T[0] = T[N] = Tsur; // BCs
        vector<double> a(N - 1, -r), b(N - 1, 1 + 2 * r), c(N - 1, -r); // tridiagonal matrix
        double nextOutput = 0.0;
        ofstream("laasonen.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T, t, "Laasonen");
                saveToCSV(T, t, "laasonen");
                nextOutput += outputInterval;
            }

            vector<double> d(N - 1, 0.0);
            for (int i = 1; i < N; i++) d[i - 1] = T[i];
            d[0] += r * Tsur; d[N - 2] += r * Tsur; // BC effect
            vector<double> T_inner = thomasSolver(a, b, c, d); // solve tridiagonal
            for (int i = 1; i < N; i++) T[i] = T_inner[i - 1];
        }
    }

    // Crank-Nicholson method (implicit, trapezoidal in time)
    void crankNicholsonMethod() {
        double r = D * dt / (2 * dx * dx);
        vector<double> T(N + 1, Tin);
        T[0] = T[N] = Tsur; // BCs
        vector<double> a(N - 1, -r), b(N - 1, 1 + 2 * r), c(N - 1, -r);
        double nextOutput = 0.0;
        ofstream("cranknicholson.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T, t, "Crank-Nicholson");
                saveToCSV(T, t, "cranknicholson");
                nextOutput += outputInterval;
            }

            vector<double> d(N - 1, 0.0);
            for (int i = 1; i < N; i++)
                d[i - 1] = r * T[i - 1] + (1 - 2 * r) * T[i] + r * T[i + 1]; // RHS
            d[0] += r * Tsur; d[N - 2] += r * Tsur; // BCs
            vector<double> T_inner = thomasSolver(a, b, c, d); // solve tridiagonal
            for (int i = 1; i < N; i++) T[i] = T_inner[i - 1];
        }
    }

    // Compute error for Laasonen method for different dt
    void LaasonenErrorMethod() {
        cout << "\n------------------------------------------------------------\n";
        cout << "Effect of Time Step (Δt) on Laasonen Accuracy&Time (dx=0.05)\n";
        cout << "Δt,RMS_Error,Max_Error,Runtime_ms\n";
        cout << "------------------------------------------------------------\n";

        double dt_values[] = { 0.01, 0.025, 0.05, 0.1 };
        ofstream errorFile("LaasonenError.csv", ios::trunc); // open CSV
        errorFile << "dt,RMS_Error,Max_Error,Runtime_ms\n";

        for (double current_dt : dt_values) {
            double r = D * current_dt / (dx * dx);  // compute r for current dt
            vector<double> T(N + 1, Tin);
            T[0] = T[N] = Tsur; // BCs
            int currentTimeSteps = static_cast<int>(totalTime / current_dt); // steps for this dt

            vector<double> a(N - 1, -r), b(N - 1, 1 + 2 * r), c(N - 1, -r); // tridiagonal

            auto start = high_resolution_clock::now(); // start timer

            for (int n = 0; n <= currentTimeSteps; n++) {
                vector<double> d(N - 1, 0.0);
                for (int i = 1; i < N; i++) d[i - 1] = T[i];
                d[0] += r * Tsur; d[N - 2] += r * Tsur;
                vector<double> T_inner = thomasSolver(a, b, c, d);
                for (int i = 1; i < N; i++) T[i] = T_inner[i - 1];
            }

            auto end = high_resolution_clock::now(); // end timer
            double ms = duration<double, std::milli>(end - start).count();

            vector<double> analyticalT = analyticalSolution(totalTime); // reference
            vector<double> errors = calculateErrorVec(analyticalT, T);  // RMS and max error

            cout << fixed << setprecision(3);
            cout << left << setw(10) << current_dt << setw(8) << errors[0] << setw(10) << errors[1] << setw(14) << ms << "\n";

            errorFile << current_dt << "," << errors[0] << "," << errors[1] << "," << ms << "\n"; // write to CSV
        }
        errorFile.close(); // close file
    }
};

int main() {
    HeatEquation wall(31.0, 93.0, 38.0, 149.0, 0.05, 0.01, 0.5); // initialize heat problem
    cout << "Starting calculations..." << endl;

    wall.printAnalytical();      // analytical solution
    wall.duFortFrankel();        // DuFort-Frankel method
    wall.richardsonMethod();     // Richardson method
    wall.laasonenMethod();       // Laasonen method
    wall.crankNicholsonMethod(); // Crank-Nicholson method
    wall.LaasonenErrorMethod();  // error analysis for Laasonen

    cout << "All methods completed. CSV files created for MATLAB plotting." << endl;
    return 0;
}
