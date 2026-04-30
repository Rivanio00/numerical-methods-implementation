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

void solveLU(int n, double* L, double* U, double* v, double* result)
{
    double y[n];

    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < i; j++) {
            sum += L[i*n + j] * y[j];
        }

        y[i] = (v[i] - sum) / L[i*n + i];
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += U[i*n + j] * result[j];
        }

        result[i] = (y[i] - sum) / U[i*n + i];
    }
}

void decompLU(int n, double* L, double* U){
    for (int j = 0; j<n; j++){
        if(U[j*n + j] == 0){
            cout << "The problem " << n << "x" << n << " havent LU decomposition" << endl;
            return;
        }
        for(int i=j+1; i<n; i++){
            L[i*n + j] = (U[i*n + j]/U[j*n + j]);
            U[i*n + j] = 0;
            for(int k=j+1; k<n; k++){
                U[i*n + k] = U[i*n + k] - (L[i*n + j] * U[j*n + k]);
            }
        }
    }
    cout << "LU decomposition:" << endl;
    cout << "L = " << endl;
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            cout << " | " << L[i*n + j];
        }
        cout << " | " << endl;
    }
    cout << "U = " << endl;
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            cout << " | " << U[i*n + j];
        }
        cout << " | " << endl;
    }
}

// method
double convergeInverse(int n, double* M)
{
    // initial vector v = (1,1,1,...)
    double v_old[n];
    double v_new[n];
    for (int i = 0; i < n; i++) v_new[i] = 1.0;

    double error = 1;
    double finalError = 0.000001;

    double eigenvalue_old = 0.0;
    double eigenvalue_new = 0.0;

    // U começa como cópia de M
    double U[n*n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            U[i*n + j] = ((double(*)[n])M)[i][j];
        }
    }

    // Identity matrix
    double L[n*n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            L[i*n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    decompLU(n, L, U);


    do {
        eigenvalue_old = eigenvalue_new;
        
        //  v_old = v_new;
        for (int i = 0; i < n; i++) v_old[i] = v_new[i];
        
        normalize(n, v_old);

        solveLU(n, L, U, v_old, v_new);

        eigenvalue_new = 0.0;
        for (int i = 0; i < n; i++) eigenvalue_new += v_old[i] * v_new[i];

        error = fabs((eigenvalue_new - eigenvalue_old) / eigenvalue_new);
    }while (error > finalError);
    
    double eigenvalue_final = (1/eigenvalue_new);


    cout << "Result for problem " << n << "x" << n << " =" << endl;
    cout << "Final Eigenvalue: " << eigenvalue_final << endl;
    cout << "Final Error: " << error << endl;
    cout << "Final vector: " << endl << "|";
    for (int i = 0; i < n; i++) cout << v_old[i] << "|";
    cout << endl << "=======================================================" << endl;
    return eigenvalue_final;
}

void convergeShifted(int n, double* M)
{
    double mu = 1.0;
    cout << "Enter shift (mu): ";
    cin >> mu;


    double A_shifted[n*n];

    // A' = A - mu * I
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j)
                A_shifted[i*n + j] = M[i*n + j] - mu;
            else
                A_shifted[i*n + j] = M[i*n + j];
        }
    }

    // Inverse Power
    double lambda_inverse = convergeInverse(n, A_shifted);

    // Final eigenvalue
    double lambda_final = lambda_inverse + mu;

    cout << "Shifted Power Result:" << endl;
    cout << "Lambda: " << lambda_final << endl << endl << endl;

}


int main()
{
    int method;
    cout << "Power Method" << endl;
    cout << "Welcome! Please choose the method: " << endl;
    cout << "0 - Inverse Power Method\n";
    cout << "1 - Shifted Inverse Power Method\n";
    cout << "-1 - Exit\n";
    cin >> method;
    if (method == -1) return 0;

    // Problem 1
    int n1 = 3;
    double M1[3][3] = {
        {5, 2, 1},
        {2, 3, 1},
        {1, 1, 2}
    };

    // Problem 2
    int n2 = 3;
    double M2[3][3] = {
        {-14, 1, -2},
        {1, -1, 1},
        {-2, 1, -11}
    };

    // Problem 3
    int n3 = 5;
    double M3[5][5] = {
        {40, 8,  4,  2, 1},
        { 8, 30, 12, 6, 2},
        { 4, 12, 20, 1, 2},
        { 2,  6,  1, 25, 4},
        { 1,  2,  2,  4, 5}
    };

    cout << scientific << setprecision(10);

    if (method == 0){
        convergeInverse(n1, (double*)M1);
        convergeInverse(n2, (double*)M2);
        convergeInverse(n3, (double*)M3);
    }else{
        convergeShifted(n1, (double*)M1);
        convergeShifted(n2, (double*)M2);
        convergeShifted(n3, (double*)M3);
    }

    return 0;
}