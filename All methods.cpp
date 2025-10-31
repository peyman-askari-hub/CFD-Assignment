#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

class HeatEquation {
private:
    const double pi = 3.14159265358979323846;
    double L, D, Tin, Tsur, dx, dt, totalTime;
    int N, timeSteps;
    vector<double> x;
    int M_max;
    double outputInterval;

    vector<double> thomasSolver(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
        int n = b.size();
        vector<double> c_star(n, 0.0), d_star(n, 0.0), T(n, 0.0);
        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];
        for (int i = 1; i < n; i++) {
            double m = b[i] - a[i] * c_star[i - 1];
            c_star[i] = c[i] / m;
            d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
        }
        T[n - 1] = d_star[n - 1];
        for (int i = n - 2; i >= 0; i--)
            T[i] = d_star[i] - c_star[i] * T[i + 1];
        return T;
    }

    void saveToCSV(const vector<double>& T, double t, const string& methodName) {
        string filename = methodName + ".csv";
        ofstream fout(filename, ios::app);
        if (!fout.is_open()) {
            cerr << "Error opening " << filename << endl;
            return;
        }
        for (int i = 0; i <= N; i++)
            fout << t << "," << x[i] << "," << T[i] << "\n";
        fout.close();
    }

    void printResult(const vector<double>& T, double t, const string& methodName, int stepInterval = 10) {
        cout << "===== " << methodName << " | t = " << t << " hr =====\n";
        for (int i = 0; i <= N; i += stepInterval) {
            cout << "x=" << fixed << setprecision(2) << x[i]
                << " cm: T=" << fixed << setprecision(2) << T[i] << " °C\n";
        }
        cout << "========================================\n\n";
    }

public:
    HeatEquation(double L_, double D_, double Tin_, double Tsur_, double dx_, double dt_, double totalTime_, int M_max_ = 200, double outputInterval_ = 0.1)
        : L(L_), D(D_), Tin(Tin_), Tsur(Tsur_), dx(dx_), dt(dt_), totalTime(totalTime_), M_max(M_max_), outputInterval(outputInterval_) {
        N = static_cast<int>(L / dx);
        timeSteps = static_cast<int>(totalTime / dt);
        x.resize(N + 1);
        for (int i = 0; i <= N; i++) x[i] = i * dx;
    }

    // Analytical solution
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
        T[0] = T[N] = Tsur;
        return T;
    }

    void printAnalytical() {
        double nextOutput = 0.0;
        ofstream("analytical.csv", ios::trunc).close(); // clear file
        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                vector<double> T = analyticalSolution(t);
                printResult(T, t, "Analytical");
                saveToCSV(T, t, "analytical");
                nextOutput += outputInterval;
            }
        }
    }

    void duFortFrankel() {
        double r = D * dt / (dx * dx);
        vector<double> T_prev(N + 1, Tin), T_curr(N + 1, Tin), T_next(N + 1, Tin);
        T_prev[0] = T_prev[N] = T_curr[0] = T_curr[N] = Tsur;
        double nextOutput = 0.0;
        ofstream("dufort.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T_curr, t, "DuFort-Frankel");
                saveToCSV(T_curr, t, "dufort");
                nextOutput += outputInterval;
            }

            if (n == 0) {
                for (int i = 1; i < N; i++)
                    T_next[i] = T_curr[i] + r * (T_curr[i + 1] - 2 * T_curr[i] + T_curr[i - 1]);
            }
            else {
                for (int i = 1; i < N; i++)
                    T_next[i] = ((1 - 2 * r) * T_prev[i] + 2 * r * (T_curr[i + 1] + T_curr[i - 1])) / (1 + 2 * r);
            }

            T_prev = T_curr;
            T_curr = T_next;
        }
    }

    void richardsonMethod() {
        double r = D * dt / (dx * dx);
        vector<double> T_prev(N + 1, Tin), T_curr(N + 1, Tin);
        T_prev[0] = T_prev[N] = T_curr[0] = T_curr[N] = Tsur;
        double nextOutput = 0.0;
        ofstream("richardson.csv", ios::trunc).close();

        for (int n = 0; n <= timeSteps; n++) {
            double t = n * dt;
            if (t >= nextOutput - 1e-8) {
                printResult(T_curr, t, "Richardson");
                saveToCSV(T_curr, t, "richardson");
                nextOutput += outputInterval;
            }

            if (n == 0) {
                for (int i = 1; i < N; i++)
                    T_curr[i] = T_prev[i] + r * (T_prev[i + 1] - 2 * T_prev[i] + T_prev[i - 1]);
            }
            else {
                vector<double> T_next(N + 1, Tin);
                T_next[0] = T_next[N] = Tsur;
                for (int i = 1; i < N; i++)
                    T_next[i] = T_prev[i] + 2 * r * (T_curr[i + 1] - 2 * T_curr[i] + T_curr[i - 1]);
                T_prev = T_curr;
                T_curr = T_next;
            }
        }
    }

    void laasonenMethod() {
        double r = D * dt / (dx * dx);
        vector<double> T(N + 1, Tin);
        T[0] = T[N] = Tsur;
        vector<double> a(N - 1, -r), b(N - 1, 1 + 2 * r), c(N - 1, -r);
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
            d[0] += r * Tsur; d[N - 2] += r * Tsur;
            vector<double> T_inner = thomasSolver(a, b, c, d);
            for (int i = 1; i < N; i++) T[i] = T_inner[i - 1];
        }
    }

    void crankNicholsonMethod() {
        double r = D * dt / (2 * dx * dx);
        vector<double> T(N + 1, Tin);
        T[0] = T[N] = Tsur;
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
                d[i - 1] = r * T[i - 1] + (1 - 2 * r) * T[i] + r * T[i + 1];
            d[0] += r * Tsur; d[N - 2] += r * Tsur;
            vector<double> T_inner = thomasSolver(a, b, c, d);
            for (int i = 1; i < N; i++) T[i] = T_inner[i - 1];
        }
    }
};

int main() {
    HeatEquation wall(31.0, 93.0, 38.0, 149.0, 0.05, 0.01, 0.5);
    cout << "Starting calculations..." << endl;

    wall.printAnalytical();
    wall.duFortFrankel();
    wall.richardsonMethod();
    wall.laasonenMethod();
    wall.crankNicholsonMethod();

    cout << "All methods completed. CSV files created for MATLAB plotting." << endl;
    return 0;
}
