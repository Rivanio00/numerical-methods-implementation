#ifndef EIGEN_METHODS_H
#define EIGEN_METHODS_H

#include <cmath>
#include <iostream>
using namespace std;

void normalize(int n, double* v){
    double norm = 0.0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) v[i] /= norm;
}

void multiMatrix(int n, double* H, double* Hi, double* result){
    double C[n*n];
    for (int i = 0; i < (n*n); i++) C[i] = 0.0;

    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            for (int k = 0; k<n; k++){
                C[i*n + j] = C[i*n + j] + H[i*n + k] * Hi[k*n + j];
            }
            result[i*n + j] = C[i*n + j];
        }
    }
}

void subHouseHolder(double* H, double* A, int step, int n){
    double v[n];
    for (int i = 0; i < n; i++) v[i] = 0.0;

    for (int i = step+1; i < n; i++){
        v[i] = A[i*n + step];
    }

    double norma = 0.0;
    for (int i = step+1; i < n; i++) norma += v[i] * v[i];
    norma = sqrt(norma);

    double v_line[n];
    for (int i = 0; i < n; i++) v_line[i] = 0.0;
    v_line[step+1] = norma;

    double w[n];
    for (int i = 0; i < n; i++) w[i] = v[i] - v_line[i];

    double u[n];
    for (int i = 0; i < n; i++) u[i] = w[i];
    normalize(n, u);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double identity = (i == j) ? 1.0 : 0.0;
            H[i*n + j] = identity - 2.0 * u[i] * u[j];
        }
    }
}

void similarityTransformation(double* A, double* Hi, double* Ai, int n){
    double temp[n*n];
    multiMatrix(n, Ai, Hi, temp);

    double transposed[n*n];
    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            transposed[i*n + j] = Hi[j*n + i];

    double result[n*n];
    multiMatrix(n, transposed, temp, result);

    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            A[i*n + j] = result[i*n + j];
}

void houseHolderMethod(int n, double* A, double* A_final, double* H_out){
    double Hi[n*n];
    double Ai[n*n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Ai[i*n + j] = ((double(*)[n])A)[i][j];
            Hi[i*n + j] = 0.0;
        }
    }

    double H[n*n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            H[i*n + j] = (i == j) ? 1.0 : 0.0;

    for (int i = 0; i<(n-2); i++){
        subHouseHolder(Hi, Ai, i, n);
        similarityTransformation((double*)A, Hi, Ai, n);

        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                Ai[k*n + j] = ((double(*)[n])A)[k][j];

        double temp[n*n];
        multiMatrix(n, H, Hi, temp);
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                H[k*n + j] = temp[k*n + j];
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++){
            A_final[i*n + j] = A[i*n + j];
            H_out[i*n + j] = H[i*n + j];
        }
}

void matrizJacobiBaseadaNoElemento_ij_DeRvelha(int n, int i, int j, double* A, double* Jij){
    for (int k = 0; k < n*n; k++) Jij[k] = 0.0;
    for (int k = 0; k < n; k++) Jij[k*n + k] = 1.0;

    double diff = 0.0;
    double error = 0.000001;

    if (abs(A[i*n+j]) <= error) return;
    if (abs(A[j*n+j]) <= error){
        diff = (A[i*n+j] < 0) ? 3.14159265/2 : -3.14159265/2;
    } else {
        diff = atan2(-A[i*n + j], A[j*n + j]);
    }

    Jij[i*n + i] = cos(diff);
    Jij[j*n + j] = cos(diff);
    Jij[i*n + j] = sin(diff);
    Jij[j*n + i] = -sin(diff);
}

void qrDecomposition(int n, double* A, double* Q, double* R_new){
    double QT[n*n] = {};
    for(int i = 0; i<n; i++) QT[i*n + i] = 1;

    double R_old[n*n];
    for (int i = 0; i < n * n; i++) R_old[i] = A[i];

    double Jij[n*n] = {};

    for (int j = 0; j<(n-1); j++){
        for (int i = j+1; i<n; i++){
           matrizJacobiBaseadaNoElemento_ij_DeRvelha(n, i, j, R_old, Jij);
           multiMatrix(n, Jij, R_old, R_new);
           for (int i = 0; i < n * n; i++) R_old[i] = R_new[i];

            double temp[n*n];
            multiMatrix(n, Jij, QT, temp);
            for (int i = 0; i < n * n; i++) QT[i] = temp[i];
        }
    }

    for(int i = 0; i<n; i++)
        for(int j = 0; j<n; j++)
            Q[i*n + j] = QT[j*n + i];
}

double calcVal(int n, double* A_new){
    double result = 0.0;
    for (int i = 0; i<n; i++)
        for (int j = 0; j<n; j++)
            if (i > j) result += (A_new[i*n + j] * A_new[i*n + j]);
    return result;
}

void qrMethod(int n, double* A, double* A_new, double* P, double* lamb){
    double finalError = 0.000001;
    double val = 100.0;

    for(int i = 0; i<n; i++) P[i*n + i] = 1;

    double Q[n*n] = {};
    double R[n*n] = {};

    while (val > finalError){
        qrDecomposition(n, A, Q, R);
        multiMatrix(n, R, Q, A_new);

        for (int i = 0; i < n * n; i++) A[i] = A_new[i];

        double temp[n*n];
        multiMatrix(n, P, Q, temp);
        for (int i = 0; i < n * n; i++) P[i] = temp[i];

        val = calcVal(n, A_new);
    }
    for(int i = 0; i<n; i++) lamb[i] = A_new[i*n+i];
}

#endif
