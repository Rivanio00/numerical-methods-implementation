#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip> 

using namespace std;

double f(double x);
double forwardDifference(double x, double h, int d, int e);
double backwardDifference(double x, double h, int d, int e);
double centralDifference(double x, double h, int d, int e);

double f(double x)
{
    return pow(x, 4) + 2*pow(x, 3) - pow(x, 2) + 5*x + sin(x);
}

double undefined() {
    return NAN;
}

double forwardDifference(double x, double h, int d, int e)
{
    if (d == 1){
        switch(e){
            case 1: return (f(x+h)-f(x))/h;
            case 2: return (-3*f(x)+4*f(x+h)-f(x+2*h))/(2*h);
            case 3: return (-11*f(x)+18*f(x+h)-9*f(x+2*h)+2*f(x+3*h))/(6*h);
            case 4: return (-25*f(x)+48*f(x+h)-36*f(x+2*h)+16*f(x+3*h)-3*f(x+4*h))/(12*h);
        }
    }
    else if (d == 2){
        switch(e){
            case 1: return (f(x)-2*f(x+h)+f(x+2*h))/(h*h);
            case 2: return (2*f(x)-5*f(x+h)+4*f(x+2*h)-f(x+3*h))/(h*h);
            case 3: return (35*f(x)-104*f(x+h)+114*f(x+2*h)-56*f(x+3*h)+11*f(x+4*h))/(12*h*h);
            case 4: return (45*f(x)-154*f(x+h)+214*f(x+2*h)-156*f(x+3*h)+61*f(x+4*h)-10*f(x+5*h))/(12*h*h);
        }
    }
    else if (d == 3){
        switch(e){
            case 1: return (-f(x)+3*f(x+h)-3*f(x+2*h)+f(x+3*h))/(h*h*h);
            case 2: return (-5*f(x)+18*f(x+h)-24*f(x+2*h)+14*f(x+3*h)-3*f(x+4*h))/(2*h*h*h);
            case 3: return (-17*f(x)+71*f(x+h)-118*f(x+2*h)+98*f(x+3*h)-41*f(x+4*h)+7*f(x+5*h))/(4*h*h*h);
            case 4: return undefined();
        }
    }
    return 0;
}

double backwardDifference(double x, double h, int d, int e)
{
    if (d == 1){
        switch(e){
            case 1: return (f(x)-f(x-h))/h;
            case 2: return (3*f(x)-4*f(x-h)+f(x-2*h))/(2*h);
            case 3: return (11*f(x)-18*f(x-h)+9*f(x-2*h)-2*f(x-3*h))/(6*h);
            case 4: return (25*f(x)-48*f(x-h)+36*f(x-2*h)-16*f(x-3*h)+3*f(x-4*h))/(12*h);
        }
    }
    else if (d == 2){
        switch(e){
            case 1: return (f(x)-2*f(x-h)+f(x-2*h))/(h*h);
            case 2: return (2*f(x)-5*f(x-h)+4*f(x-2*h)-f(x-3*h))/(h*h);
            case 3: return (35*f(x)-104*f(x-h)+114*f(x-2*h)-56*f(x-3*h)+11*f(x-4*h))/(12*h*h);
            case 4: return (45*f(x)-154*f(x-h)+214*f(x-2*h)-156*f(x-3*h)+61*f(x-4*h)-10*f(x-5*h))/(12*h*h);
        }
    }
    else if (d == 3){
        switch(e){
            case 1: return (f(x)-3*f(x-h)+3*f(x-2*h)-f(x-3*h))/(h*h*h);
            case 2: return (5*f(x)-18*f(x-h)+24*f(x-2*h)-14*f(x-3*h)+3*f(x-4*h))/(2*h*h*h);
            case 3: return (17*f(x)-71*f(x-h)+118*f(x-2*h)-98*f(x-3*h)+41*f(x-4*h)-7*f(x-5*h))/(4*h*h*h);
            case 4: return undefined();
        }
    }
    return 0;
}

double centralDifference(double x, double h, int d, int e)
{
    if (d == 1){
        switch(e){
            case 1: return undefined();
            case 2: return (f(x+h)-f(x-h))/(2*h);
            case 3: return undefined();
            case 4: return (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
        }
    }
    else if (d == 2){
        switch(e){
            case 1: return undefined();
            case 2: return (f(x+h)-2*f(x)+f(x-h))/(h*h);
            case 3: return undefined();
            case 4: return (-f(x+2*h)+16*f(x+h)-30*f(x)+16*f(x-h)-f(x-2*h))/(12*h*h);
        }
    }
    else if (d == 3){
        switch(e){
            case 1: return undefined();
            case 2: return (-f(x+2*h)+2*f(x+h)-2*f(x-h)+f(x-2*h))/(2*h*h*h);
            case 3: return undefined();
            case 4: return (f(x+3*h)-8*f(x+2*h)+13*f(x+h)-13*f(x-h)+8*f(x-2*h)-f(x-3*h))/(8*h*h*h);
        }
    }
    return 0;
}

int main()
{
    cout << fixed << setprecision(8);
    int method;
    cout << "Welcome! Please choose the differentiation method: " << endl;
    cout << "0 - Forward Difference\n";
    cout << "1 - Backward Difference\n";
    cout << "2 - Central Difference\n";
    cout << "3 - Complete Test\n";
    cout << "-1 - Exit\n";
    cin >> method;
    if (method == -1) 
    {
        return 0;
    }
    if (method < 0 || method > 3)
    {
        cout << "Invalid method selected." << endl;
        return 1;
    }
    double x;
    double stepSize;
    int derivativeOrder;
    int errorOrder;
    cout << "Enter the value of x:" << endl;
    cin >> x;

    if (method == 3){
        double h_values[] = {1e-4, 1e-3, 5e-3};

        cout << "\n=============================================\n";
        cout << " Complete Test - Numerical Differentiation\n";
        cout << " f(x) at x = " << x << "\n";
        cout << "=============================================\n\n";

        for (int m = 0; m < 3; m++)
        {
            string methodName;
            if (m == 0) methodName = "FORWARD";
            if (m == 1) methodName = "BACKWARD";
            if (m == 2) methodName = "CENTRAL";

            cout << "===== " << methodName << " =====\n";

            for (int d = 1; d <= 3; d++)
            {
                double h = h_values[d-1];

                for (int e = 1; e <= 4; e++)
                {
                    double result;

                    switch(m)
                    {
                        case 0: result = forwardDifference(x, h, d, e); break;
                        case 1: result = backwardDifference(x, h, d, e); break;
                        case 2: result = centralDifference(x, h, d, e); break;
                    }

                    cout << "Method: " << methodName
                        << " | d=" << d
                        << " | O(h^" << e << ")"
                        << " | h=" << h;

                    if (isnan(result))
                        cout << " | Undefined formula";
                    else
                        cout << " | Result=" << result;

                    cout << endl;
                }
                cout << endl;
            }
        }

        return 0;
    }
    
    cout << "Enter the step size (h):" << endl;
    cin >> stepSize;
    if (stepSize <= 0)
    {
        cout << "Step size must be positive." << endl;
        return 1;
    }

    cout << "Enter the derivative order:" << endl;
    cin >> derivativeOrder;
    if (derivativeOrder < 0 || derivativeOrder > 3){
        cout << "Derivative order must be 1, 2 or 3" << endl;
        return 1;
    }

    cout << "Enter the error order:" << endl;
    cin >> errorOrder;
    if (errorOrder < 0 || errorOrder > 4){
        cout << "Derivative order must be 1, 2, 3 or 4" << endl;
        return 1;
    }


    cout << "Computing..." << endl;
    double result;
    if (derivativeOrder == 0){
        result = f(x);
    }else{
        switch(method)
        {
            case 0: 
                result = forwardDifference(x, stepSize, derivativeOrder, errorOrder);
                break;
            case 1: 
                result = backwardDifference(x, stepSize, derivativeOrder, errorOrder);
                break;
            case 2:
                result = centralDifference(x, stepSize, derivativeOrder, errorOrder);
                break;
        }
    }
    if (isnan(result)) {
        cout << "Method: " << method
            << " | d=" << derivativeOrder
            << " | O(h^" << errorOrder << ")"
            << " | h=" << stepSize
            << " | Undefined formula for the selected error\n";
    } else {
        cout << "Method: " << method
            << " | d=" << derivativeOrder
            << " | O(h^" << errorOrder << ")"
            << " | h=" << stepSize
            << " | Result=" << result << endl;
    }
    return 0;
}