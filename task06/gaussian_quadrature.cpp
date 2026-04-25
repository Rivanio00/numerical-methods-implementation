#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <functional>

using namespace std;

struct Quadratura {
    vector<double> x;
    vector<double> w;
};

double f(double x)
{
    return x*x + 2*x + 1;
}

double integrar(const Quadratura& q, function<double(double)> f) {
    double soma = 0.0;

    for (size_t i = 0; i < q.x.size(); ++i) {
        soma += q.w[i] * f(q.x[i]);
    }

    return soma;
}

Quadratura hermite2 = {
    { -1.0/sqrt(2), 1.0/sqrt(2) },
    { sqrt(M_PI)/2.0, sqrt(M_PI)/2.0 }
};

Quadratura hermite3 = {
    { -sqrt(3.0/2.0), 0, sqrt(3.0/2.0) },
    { sqrt(M_PI)/6.0, 2*sqrt(M_PI)/3.0, sqrt(M_PI)/6.0 }
};

Quadratura hermite4 = {
    { -(sqrt((3+sqrt(6))/2.0)), -(sqrt((3-sqrt(6))/2.0)), (sqrt((3-sqrt(6))/2.0)), (sqrt((3+sqrt(6))/2.0)) },
    { sqrt(M_PI)/(4*(3+sqrt(6))), sqrt(M_PI)/(4*(3-sqrt(6))), sqrt(M_PI)/(4*(3-sqrt(6))), sqrt(M_PI)/(4*(3+sqrt(6))) }
};

Quadratura laguerre2 = {
    { 2 - sqrt(2), 2 + sqrt(2) },
    { (1.0/4.0)*(2 + sqrt(2)), (1.0/4.0)*(2 - sqrt(2)) }
};

Quadratura laguerre3 = {
    { 0.4157745568, 2.2942803603, 6.2899450829 },
    { 0.7110930099, 0.2785177336, 0.0103892565 }
};

Quadratura laguerre4 = {
    { 0.322547, 1.74576, 4.536620, 9.395070 },
    { 0.603154, 0.357419, 0.0388879, 0.000539295 }
};

Quadratura chebyshev2 = {
    { -1.0/sqrt(2), 1.0/sqrt(2) },
    { (M_PI)/2.0, (M_PI)/2.0 }
};

Quadratura chebyshev3 = {
    { -sqrt(3)/2.0, 0, sqrt(3)/2.0 },
    { (M_PI)/3.0, (M_PI)/3.0, (M_PI)/3.0 }
};

Quadratura chebyshev4 = {
    { 0.923879, 0.38683, -0.38683, -0.923879 },
    { (M_PI)/4.0, (M_PI)/4.0, (M_PI)/4.0, (M_PI)/4.0 }
};

int main()
{
    cout << fixed << setprecision(6);

    cout << "GAUSS-HERMITE" << endl;
    cout << "n=2: " << integrar(hermite2, f) << endl;
    cout << "n=3: " << integrar(hermite3, f) << endl;
    cout << "n=4: " << integrar(hermite4, f) << endl;

    cout << "\nGAUSS-LAGUERRE" << endl;
    cout << "n=2: " << integrar(laguerre2, f) << endl;
    cout << "n=3: " << integrar(laguerre3, f) << endl;
    cout << "n=4: " << integrar(laguerre4, f) << endl;

    cout << "\nGAUSS-CHEBYSHEV" << endl;
    cout << "n=2: " << integrar(chebyshev2, f) << endl;
    cout << "n=3: " << integrar(chebyshev3, f) << endl;
    cout << "n=4: " << integrar(chebyshev4, f) << endl;

    return 0;
}