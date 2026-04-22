#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x)
{
    return pow((sin(2*x) + 4*x*x + 3*x), 2);
}

double closedNewtonCotes(double h, double x)
{
    h = h/4; // xi ... xi + 4h
    return (2*h/45) * (7*f(x) + 32*f(x+h) + 12*f(x+ (2*h)) + 32*f(x + (h*3)) + 7*f(x+ (4*h)));
}

double openNewtonCotes(double h, double x)
{
    h = h/6; // xi+h ... xi + 5h
    return (h/10) * (33*f(x+h) - 42*f(x+ (2*h)) + 78*f(x +(h*3)) - 42*f(x+(4*h)) + 33*f(x+ (5*h)));
}

int main()
{
    int method;
    cout << "Degree 4 substitution polynomials" << endl;
    cout << "Welcome! Please choose the integral method: " << endl;
    cout << "0 - Closed Newton-Cotes formulas\n";
    cout << "1 - Open Newton-Cotes formulas\n";
    cout << "-1 - Exit\n";
    cin >> method;

    double finalError = 0.000001; // tolerable error
    double error = 1; // inicial error
    int step = 1; // iteration


    if (method == 0){
        cout << "Set the range [a,b]" << endl;
    } else if (method == 1) {
        cout << "Set the range (a,b)" << endl;
    }else{
        return 0;
    }

    double a = 0.0;
    double b = 0.0;
    
    cout << "Enter the value of a:" << endl;
    cin >> a;
    cout << "Enter the value of b:" << endl;
    cin >> b;

    double h = (b-a); // interval size
    int partitions = 1; // number of problem partitions

    cout << scientific << setprecision(10);

    double result = 0.0;
    double partialresult = 0.0;

    while (error > finalError){
        partialresult = 0.0;
        for (int i=0; i<partitions; i++){ // calculate for each partition
            if(method == 0){
                // h/partitions -> interval size in this partition
                // h/partitions)*i -> first point for each partition
                partialresult += closedNewtonCotes((h/partitions), a+((h/partitions)*i));
            }else{
                partialresult += openNewtonCotes((h/partitions), a+((h/partitions)*i));
            }
        }
        cout << "Step: " << step << endl;
        cout << "Partial Result: " << partialresult << endl;
        error = fabs((partialresult - result)/partialresult);
        cout << "Error: " << error << endl;
        cout << "==============================" << endl;
        result = partialresult;
        step++;
        partitions = partitions * 2;
    }
    
    cout << "Finished: Error < " << finalError << endl;
    return 0;
}