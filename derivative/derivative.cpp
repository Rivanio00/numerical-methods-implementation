#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip> 

using namespace std;

double f(double x);
double forwardDifference(double x, double stepSize, int derivativeOrder);
double backwardDifference(double x, double stepSize, int derivativeOrder);
double centralDifference(double x, double stepSize, int derivativeOrder);

double f(double x)
{
    return pow(x, 4) + 2*pow(x, 3) - pow(x, 2) + 5*x + sin(x);
}

double forwardDifference(double x, double stepSize, int derivativeOrder)
{
    if (derivativeOrder == 1)
    {
        return (f(x+stepSize)-f(x))/stepSize;
    }
    return (forwardDifference(x+stepSize, stepSize, derivativeOrder-1) - forwardDifference(x, stepSize, derivativeOrder-1))/stepSize;
}

double backwardDifference(double x, double stepSize, int derivativeOrder)
{
    if (derivativeOrder == 1)
    {
        return (f(x)-f(x-stepSize))/stepSize;
    }
    return (backwardDifference(x, stepSize, derivativeOrder-1) - backwardDifference(x-stepSize, stepSize, derivativeOrder-1))/stepSize;
}

double centralDifference(double x, double stepSize, int derivativeOrder)
{
    if (derivativeOrder == 1)
    {
        return (f(x+stepSize)-f(x-stepSize))/(2*stepSize);
    }
    return (centralDifference(x+stepSize, stepSize, derivativeOrder-1) - centralDifference(x-stepSize, stepSize, derivativeOrder-1))/(2*stepSize);
}

int main()
{
    cout << fixed << setprecision(8);
    int method;
    cout << "Welcome! Please choose the differentiation method: " << endl;
    cout << "0 - Forward Difference\n";
    cout << "1 - Backward Difference\n";
    cout << "2 - Central Difference\n";
    cout << "-1 - Exit\n";
    cin >> method;
    if (method == -1) 
    {
        return 0;
    }
    if (method < 0 || method > 2)
    {
        cout << "Invalid method selected." << endl;
        return 1;
    }
    double x;
    double stepSize;
    int derivativeOrder;
    cout << "Enter the value of x:" << endl;
    cin >> x;

    cout << "Enter the step size (h):" << endl;
    cin >> stepSize;
    if (stepSize <= 0)
    {
        cout << "Step size must be positive." << endl;
        return 1;
    }

    cout << "Enter the derivative order:" << endl;
    cin >> derivativeOrder;
    if (derivativeOrder < 0){
        cout << "Derivative order must be non-negative." << endl;
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
                result = forwardDifference(x, stepSize, derivativeOrder);
                break;
            case 1: 
                result = backwardDifference(x, stepSize, derivativeOrder);
                break;
            case 2:
                result = centralDifference(x, stepSize, derivativeOrder);
                break;
        }
    }
    cout << "The result is: " << result << endl;
    return 0;
}