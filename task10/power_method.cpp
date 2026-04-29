#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void normalize(int n, double* v){
    double norm = 0.0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) v[i] /= norm;
}

void matVec(int n, double* M, double* v, double* result)
{
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += M[i*n + j] * v[j];
        }
    }
}

// method
void converge(int n, double* M)
{
    // initial vector v = (1,1,1,...)
    double v_old[n];
    double v_new[n];
    for (int i = 0; i < n; i++) v_new[i] = 1.0;

    double error = 1;
    double finalError = 0.000001;

    double eigenvalue_old = 0.0;
    double eigenvalue_new = 0.0;


    do {
        eigenvalue_old = eigenvalue_new;
        
        //  v_old = v_new;
        for (int i = 0; i < n; i++) v_old[i] = v_new[i];
        
        normalize(n, v_old);

        matVec(n, M, v_old, v_new);

        eigenvalue_new = 0.0;
        for (int i = 0; i < n; i++) eigenvalue_new += v_old[i] * v_new[i];

        error = fabs((eigenvalue_new - eigenvalue_old) / eigenvalue_new);
    }while (error > finalError);
    
    cout << "Result for problem " << n << "x" << n << " =" << endl;
    cout << "Final Eigenvalue: " << eigenvalue_new << endl;
    cout << "Final Error: " << error << endl;
    cout << "Final vector: " << endl << "|";
    for (int i = 0; i < n; i++) cout << v_old[i] << "|";
    cout << endl << "==========================" << endl;
}


int main()
{
    // Problem 1
    int n1 = 3;
    double M1[3][3] = {
        {5, 2, 1},
        {2, 3, 1},
        {1, 1, 2}
    };

    // Problem 2
    int n2 = 5;
    double M2[5][5] = {
        {40, 8,  4,  2, 1},
        { 8, 30, 12, 6, 2},
        { 4, 12, 20, 1, 2},
        { 2,  6,  1, 25, 4},
        { 1,  2,  2,  4, 5}
    };

    cout << scientific << setprecision(10);

    converge(n1, (double*)M1);
    converge(n2, (double*)M2);

    return 0;
}