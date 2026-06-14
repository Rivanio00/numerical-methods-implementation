#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double g = 10.0; // gravity
const double k = 0.25; // kg/s
const double m = 2.0;  // masss

// F(S, t) -> dS/dt
void F(double v, double y, double* dv, double* dy) {
    *dv = -g - (k/m) * v;
    *dy = v;
}

void rk3(double dt) {
    double v = 5.0;     // v0 speed
    double y = 200.0;   // height
    double t = 0.0;     // t0

    double ymax = y, tmax = 0.0;

    cout << "\n=== Delta_t = " << dt << " ===" << endl;

    while (y > 0.0) {
        // F1
        double dv1, dy1;
        F(v, y, &dv1, &dy1);

        // S_meio = Si + (dt/2)*F1
        double v_meio = v + (dt/2.0) * dv1;
        double y_meio = y + (dt/2.0) * dy1;

        // F2
        double dv2, dy2;
        F(v_meio, y_meio, &dv2, &dy2);

        // S_barra = Si + dt*(-F1 + 2*F2)
        double v_barra = v + dt * (-dv1 + 2.0*dv2);
        double y_barra = y + dt * (-dy1 + 2.0*dy2);

        // F3
        double dv3, dy3;
        F(v_barra, y_barra, &dv3, &dy3);

        // S_{i+1}
        double v_new = v + dt * (1.0/6.0*dv1 + 4.0/6.0*dv2 + 1.0/6.0*dv3);
        double y_new = y + dt * (1.0/6.0*dy1 + 4.0/6.0*dy2 + 1.0/6.0*dy3);

        t += dt;

        // altura maxima
        if (y_new > ymax) {
            ymax = y_new;
            tmax = t;
        }

        v = v_new;
        y = y_new;
    }

    cout << fixed << setprecision(6);
    cout << "y_max    = " << ymax   << " m"   << endl;
    cout << "t_max    = " << tmax   << " s"   << endl;
    cout << "t_total  = " << t      << " s"   << endl;
    cout << "v_impact = " << v      << " m/s" << endl;
}

int main() {
    cout << scientific << setprecision(6);

    double deltas[] = {1.0, 0.5, 0.1, 0.01, 0.001};
    int n = sizeof(deltas) / sizeof(deltas[0]);

    for (int i = 0; i < n; i++) {
        rk3(deltas[i]);
    }

    return 0;
}