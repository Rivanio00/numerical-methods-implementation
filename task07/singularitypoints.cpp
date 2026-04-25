#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// I = 2
double f1simple(double s)
{
    double xs   = 0.5 + 0.5 * tanh(s);
    double dxds = 0.5 / (cosh(s) * cosh(s));
    return (1.0 / sqrt(xs)) * dxds;
}

double f1double(double s)
{
    double arg  = (M_PI / 2.0) * sinh(s);
    double xs   = 0.5 + 0.5 * tanh(arg);
    double dxds = 0.5 * (M_PI / 2.0) * cosh(s) / (cosh(arg) * cosh(arg));
    return (1.0 / sqrt(xs)) * dxds;
}

// I = 1.5708
double f2simple(double s)
{
    double xs   = -1.0 + tanh(s);
    double dxds = 1.0 / (cosh(s) * cosh(s));
    return (1.0 / sqrt(4.0 - xs * xs)) * dxds;
}

double f2double(double s)
{
    double arg  = (M_PI / 2.0) * sinh(s);
    double xs   = -1.0 + tanh(arg);
    double dxds = (M_PI / 2.0) * cosh(s) / (cosh(arg) * cosh(arg));
    return (1.0 / sqrt(4.0 - xs * xs)) * dxds;
}

double closedNewtonCotes(double h, double x, double (*f)(double))
{
    h = h/4; // xi ... xi + 4h
    return (2*h/45) * (7*f(x) + 32*f(x+h) + 12*f(x+ (2*h)) + 32*f(x + (h*3)) + 7*f(x+ (4*h)));
}

double openNewtonCotes(double h, double x, double (*f)(double))
{
    h = h/6; // xi+h ... xi + 5h
    return (h/10) * (33*f(x+h) - 42*f(x+ (2*h)) + 78*f(x +(h*3)) - 42*f(x+(4*h)) + 33*f(x+ (5*h)));
}

double integral(double a, double b, double (*f)(double), int method){
    int step = 1; // iteration
    double finalError = 0.000001; // tolerable error
    double error = 1; // inicial error
    double h = (b-a); // interval size
    int partitions = 1; // number of problem partitions
    double result = 0.0;
    double partialresult = 0.0;
    while (error > finalError){
        partialresult = 0.0;
        for (int i=0; i<partitions; i++){ // calculate for each partition
            if(method == 0){
                // h/partitions -> interval size in this partition
                // h/partitions)*i -> first point for each partition
                partialresult += closedNewtonCotes((h/partitions), a+((h/partitions)*i), f);
            }else{
                partialresult += openNewtonCotes((h/partitions), a+((h/partitions)*i), f);
            }
        }
        error = fabs((partialresult - result)/partialresult);
        result = partialresult;
        step++;
        partitions = partitions * 2;
    }
    return result;
}

int main()
{
    int method;
    cout << "Singularity at the limits of integration" << endl;
    cout << "Welcome! Please choose the integral method: " << endl;
    cout << "0 - Closed Newton-Cotes formulas\n";
    cout << "1 - Open Newton-Cotes formulas\n";
    cout << "-1 - Exit\n";
    cin >> method;
    if (method == -1){
        return 0;
    }

    double finalError = 0.000001; // tolerable error
    double error = 1; // inicial error
    int step = 1; // iteration

    double (*funcs[])(double) = {f1simple, f1double, f2simple, f2double};
    string names[] = {"f1 simple", "f1 double", "f2 simple", "f2 double"};

    cout << scientific << setprecision(10);

    for (int i = 0; i < 4; i++) {
        cout << "\n=== " << names[i] << " ===" << endl;
        
        double c = 1.0;
        double result = 0.0, error = 1.0;
        int iter = 0;
        
        while (error > finalError) {
            double newResult = integral(-c, c, funcs[i], method);
            error = fabs((newResult - result) / newResult);
            
            //cout << "c=" << c << " | I=" << newResult << " | erro=" << error << endl;
            
            result = newResult;
            c += 0.1;
            iter++;
        }
        cout << "c=" << c-0.1 << " | I=" << result << " | erro=" << error << " | iter=" << iter << endl;
    }
    
    cout << "Finished: Error < " << finalError << endl;
    return 0;
}