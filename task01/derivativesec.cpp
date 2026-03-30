#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x)
{
    return std::sqrt(std::exp(3 * x) + 4 * x * x);
}

double fd2(double h, double x)
{
    return (-f(x + 2 * h) 
        + 16 * f(x + h) 
        - 30 * f(x) 
        + 16 * f(x - h) 
        - f(x - 2 * h)) 
        / (12 * h * h);
}

int main()
{
    double x = 2.0;
    double fx = 0;      // f(x)
    double fx2 = 0;     // f''(x)
    double fx2_pre = 0; // previous iteration

    double finalError = 0.00001;
    double error = 1;

    int step = 1;
    double stepSize = 0.5;

    cout << scientific << setprecision(10);

    while (error >= finalError) {
        fx = f(x);
        fx2 = fd2(stepSize, x);
        if (step == 1) {
            error = 1;
        } else {
            error = fabs((fx2 - fx2_pre) / fx2);
        }

        cout << "==============================" << endl;
        cout << "Step:        " << step << endl;
        cout << "Step Size:   " << setw(10) << stepSize << endl;
        cout << "f(x):        " << setw(15) << fx << endl;
        cout << "f''(x):      " << setw(15) << fx2 << endl;
        cout << "Error:       " << setw(15) << error << endl;

        fx2_pre = fx2;
        stepSize /= 2;
        step++;
    }
    cout << "==============================" << endl;
    cout << "Finished: Error < " << finalError << endl;
    return 0;
}