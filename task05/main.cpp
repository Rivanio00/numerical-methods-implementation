// Gauss-Legendre (2, 3, 4 pontos) com partição adaptativa
// Compilação: g++ -std=c++17 -O2 -o quadratura main.cpp

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>

struct Quadratura {
    std::vector<double> nodes;
    std::vector<double> weights;
};

struct Resultado {
    double valor;
    int    N_final;
    int    iteracoes;
    double erro_estimado;
};

Quadratura getGaussLegendre(int n) {
    Quadratura q;
    switch (n) {
        case 2:
            q.nodes   = { -0.5773502691896257,  0.5773502691896257 };
            q.weights = {  1.0,                  1.0                };
            break;
        case 3:
            q.nodes   = { -0.7745966692414834,  0.0,                0.7745966692414834 };
            q.weights = {  0.5555555555555556,  0.8888888888888889, 0.5555555555555556 };
            break;
        case 4:
            q.nodes   = { -0.8611363115940526, -0.3399810435848563,
                           0.3399810435848563,  0.8611363115940526 };
            q.weights = {  0.3478548451374538,  0.6521451548625461,
                           0.6521451548625461,  0.3478548451374538 };
            break;
        default:
            throw std::invalid_argument("n deve ser 2, 3 ou 4");
    }
    return q;
}

// Transforma [-1,1] → [xi,xf]: x(α) = (xi+xf)/2 + (xf-xi)/2·α
double gaussLegendreSimples(const std::function<double(double)>& f,
                            double xi, double xf, int n) {
    const Quadratura& q = getGaussLegendre(n);
    const double meio   = 0.5 * (xi + xf);
    const double metade = 0.5 * (xf - xi);
    double soma = 0.0;
    for (std::size_t k = 0; k < q.nodes.size(); ++k)
        soma += q.weights[k] * f(meio + metade * q.nodes[k]);
    return metade * soma;
}

double gaussLegendreParticionado(const std::function<double(double)>& f,
                                 double xi, double xf, int n, int N) {
    const double h = (xf - xi) / N;
    double total = 0.0;
    for (int i = 0; i < N; ++i)
        total += gaussLegendreSimples(f, xi + i * h, xi + (i + 1) * h, n);
    return total;
}

// Dobra N até |I_novo - I_anterior| < tol
Resultado integrarAdaptativo(const std::function<double(double)>& f,
                             double xi, double xf, int n,
                             double tol = 1e-6, int max_iter = 30) {
    int    N       = 1;
    double I_atual = gaussLegendreParticionado(f, xi, xf, n, N);
    double erro    = 0.0;
    int    iter    = 0;

    for (iter = 1; iter <= max_iter; ++iter) {
        double I_anterior = I_atual;
        N      *= 2;
        I_atual = gaussLegendreParticionado(f, xi, xf, n, N);
        erro    = std::abs(I_atual - I_anterior);
        if (erro < tol) break;
    }

    return { I_atual, N, iter, erro };
}

double integrando(double x) {
    double g = std::sin(2.0 * x) + 4.0 * x * x + 3.0 * x;
    return g * g;
}

int main() {
    const double xi  = 0.0;
    const double xf  = 1.0;
    const double tol = 1e-6;
    const double ref = 17.8764703;

    auto f = [](double x) { return integrando(x); };

    std::cout << "Problema: I = integral de 0 a 1 de (sin(2x) + 4x^2 + 3x)^2 dx\n";
    std::cout << "Valor de referencia: "
              << std::fixed << std::setprecision(7) << ref << "\n";
    std::cout << "Tolerancia: " << std::scientific << tol << "\n\n";

    for (int n : {2, 3, 4}) {
        Resultado r = integrarAdaptativo(f, xi, xf, n, tol);
        double erro_abs = std::abs(r.valor - ref);

        std::cout << "n = " << n << " pontos:\n";
        std::cout << "  Iteracoes:         " << r.iteracoes << "\n";
        std::cout << "  N (subintervalos): " << r.N_final << "\n";
        std::cout << "  Valor calculado:   "
                  << std::fixed << std::setprecision(7) << r.valor << "\n";
        std::cout << "  Erro absoluto:     "
                  << std::scientific << std::setprecision(2) << erro_abs << "\n\n";
    }

    std::cout << "--------------------------------------------------------------\n";
    std::cout << std::left
              << std::setw(6)  << "n"
              << std::setw(12) << "Iteracoes"
              << std::setw(10) << "N final"
              << std::setw(16) << "Valor"
              << std::setw(14) << "Erro abs"
              << "\n";
    std::cout << "--------------------------------------------------------------\n";

    for (int n : {2, 3, 4}) {
        Resultado r = integrarAdaptativo(f, xi, xf, n, tol);
        double erro_abs = std::abs(r.valor - ref);

        std::cout << std::left
                  << std::setw(6)  << n
                  << std::setw(12) << r.iteracoes
                  << std::setw(10) << r.N_final
                  << std::fixed    << std::setprecision(7)
                  << std::setw(16) << r.valor
                  << std::scientific << std::setprecision(2)
                  << erro_abs << "\n";
    }
    std::cout << "--------------------------------------------------------------\n";

    return 0;
}
