#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>


double F(double S, double t) {
    return -S + t + 1.0;
}

double solucao_exata(double t) {
    return std::exp(-t) + t;
}

double runge_kutta4(double S, double t, double dt) {
    double F1 = F(S, t);
    double S2 = S + (dt / 2.0) * F1;
    double F2 = F(S2, t + dt / 2.0);
    double S3 = S + (dt / 2.0) * F2;
    double F3 = F(S3, t + dt / 2.0);
    double S4 = S + dt * F3;
    double F4 = F(S4, t + dt);
    return S + (dt / 6.0) * (F1 + 2.0*F2 + 2.0*F3 + F4);
}

// adams-bashforth-moulton ordem 4
void adams_bashforth_moulton4(
    double S0, double t0, double dt, int n_steps,
    double tol = 1e-10, int max_iter = 50
) {
    std::vector<double> S(n_steps + 1, 0.0);
    std::vector<double> Fv(n_steps + 1, 0.0);
    std::vector<double> t(n_steps + 1, 0.0);

    S[0] = S0;
    t[0] = t0;
    Fv[0] = F(S[0], t[0]);

    for (int i = 0; i < 3; ++i) {
        t[i+1] = t[i] + dt;
        S[i+1] = runge_kutta4(S[i], t[i], dt);
        Fv[i+1] = F(S[i+1], t[i+1]);
    }

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n=== Metodo Preditor-Corretor de 4a Ordem (Adams-Bashforth-Moulton) ===\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::setw(6)  << "Passo"
              << std::setw(12) << "t_i"
              << std::setw(16) << "S_i (ABM)"
              << std::setw(16) << "S_i (exata)"
              << std::setw(14) << "Erro abs."
              << "\n";
    std::cout << std::string(80, '-') << "\n";

    for (int i = 0; i <= 3; ++i) {
        double exata = solucao_exata(t[i]);
        double erro = std::abs(S[i] - exata);
        std::cout << std::setw(4) << i << " (RK4)"
                  << std::setw(12) << t[i]
                  << std::setw(16) << S[i]
                  << std::setw(16) << exata
                  << std::setw(14) << erro << "\n";
    }

    for (int i = 3; i <= n_steps - 1; ++i) {
        t[i+1] = t[i] + dt;

        // predicao
        double S_pred = S[i] + (dt / 24.0) * (
            -9.0  * Fv[i-3]
            + 37.0 * Fv[i-2]
            - 59.0 * Fv[i-1]
            + 55.0 * Fv[i]
        );

        // correcao iterativa
        double S_cor = S_pred;
        double S_prev;
        int iter = 0;
        do {
            S_prev = S_cor;
            double F_pred = F(S_prev, t[i+1]);
            S_cor = S[i] + (dt / 24.0) * (
                  Fv[i-2]
                - 5.0  * Fv[i-1]
                + 19.0 * Fv[i]
                +  9.0 * F_pred
            );
            iter++;
        } while (std::abs(S_cor - S_prev) / std::abs(S_cor) >= tol && iter < max_iter);

        S[i+1]  = S_cor;
        Fv[i+1] = F(S[i+1], t[i+1]);

        double exata = solucao_exata(t[i+1]);
        double erro  = std::abs(S[i+1] - exata);
        std::cout << std::setw(4) << (i+1) << " (PC) "
                  << std::setw(12) << t[i+1]
                  << std::setw(16) << S[i+1]
                  << std::setw(16) << exata
                  << std::setw(14) << erro << "\n";
    }

    std::cout << std::string(80, '-') << "\n";
}

int main() {
    double S0 = 1.0;
    double t0 = 0.0;
    double tf = 2.0;
    double dt = 0.1;

    // Calculo seguro de n_steps
    int n_steps = static_cast<int>(std::round((tf - t0) / dt));

    std::cout << "\nPVI: dS/dt = -S + t + 1,   S(0) = 1\n";
    std::cout << "Solucao exata: S(t) = exp(-t) + t\n";
    std::cout << "Intervalo: [" << t0 << ", " << tf << "],  dt = " << dt
              << ",  N = " << n_steps << "\n";

    adams_bashforth_moulton4(S0, t0, dt, n_steps);

    return 0;
}