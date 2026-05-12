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

    // v(step+1 : n) <- column step of A, from row step+1 to n
    for (int i = step+1; i < n; i++){
        v[i] = A[i*n + step];
    }

    // norma = ||v||
    double norma = 0.0;
    for (int i = step+1; i < n; i++) norma += v[i] * v[i];
    norma = sqrt(norma);

    // v_line = zero(n); v_line(step+1) = norma
    double v_line[n];
    for (int i = 0; i < n; i++) v_line[i] = 0.0;
    v_line[step+1] = norma;

    // w = v - v_line
    double w[n];
    for (int i = 0; i < n; i++) w[i] = v[i] - v_line[i];

    // u = w / ||w||
    double u[n];
    for (int i = 0; i < n; i++) u[i] = w[i];
    normalize(n, u);

    // H = I - 2 * u * u^T
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

    double transpose[n*n];
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            transpose[i*n + j] = Hi[j*n + i];
        }
    }

    double result[n*n];
    multiMatrix(n, transpose, temp, result);

    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            A[i*n + j] = result[i*n + j];
        }
    }
}

void houseHolderMethod(int n, double* A, double* A_final, double* H_out)
{
    // Ai <- A
    double Hi[n*n];
    double Ai[n*n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Ai[i*n + j] = ((double(*)[n])A)[i][j];
            Hi[i*n + j] = 0.0;
        }
    }


    // Identity Matrix
    double H[n*n];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            H[i*n + j] = (i == j) ? 1.0 : 0.0;
        }
    }


    for (int i = 0; i<(n-2); i++){
        // Hi <- subHouseHolder(Ai, i)
        subHouseHolder(Hi, Ai, i, n);

        // A <- Hi^t * Ai * Hi
        similarityTransformation((double*)A, Hi, Ai, n);

        // Ai <- A
        for (int k = 0; k < n; k++){
            for (int j = 0; j < n; j++){
                Ai[k*n + j] = ((double(*)[n])A)[k][j];
            }
        }

        // H <- H * Hi
        double temp[n*n];
        multiMatrix(n, H, Hi, temp);
        for (int k = 0; k < n; k++){
            for (int j = 0; j < n; j++){
                H[k*n + j] = temp[k*n + j];
            }
        }
    }
    


    cout << "Result for problem " << n << "x" << n << " =" << endl;
    cout << "Final Matrix A = " << endl;
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++){
            cout << A[i*n + j] << "|";
            A_final[i*n + j] = A[i*n + j];
        }
        cout << endl;
    }
    for (int i = 0; i < n*n; i++) H_out[i] = H[i];
}

void matrizJacobiBaseadaNoElemento_ij_DeRvelha(int n, int i, int j, double* A, double* Jij)
{
    // Jij <- I
    for (int k = 0; k < n*n; k++) Jij[k] = 0.0;
    for (int k = 0; k < n; k++) Jij[k*n + k] = 1.0;

    double diff = 0.0;

    double error = 0.000001; // tolerable error

    if (abs(A[i*n+j]) <= error){
        // considera Aij = 0, retorna identidade
        return;
    }
    if (abs(A[j*n+j]) <= error){
        if (A[i*n+j] < 0){
            diff = 3.14159265 / 2;
        }else{
            diff = (-1) * 3.14159265 / 2;
        }
    }else{
        // arctan()
        diff = atan2(-A[i*n + j], A[j*n + j]);
    }

    Jij[i*n + i] = cos(diff);
    Jij[j*n + j] = cos(diff);
    Jij[i*n + j] = sin(diff);
    Jij[j*n + i] = (-1) * sin(diff);
}

void qrDecomposition(int n, double* A, double* Q, double* R_new)
{
    // altera as matrizes Q e R
    double QT[n*n] = {};
    // QT <- I
    for(int i = 0; i<n; i++) QT[i*n + i] = 1;

    // R_old <- A
    double R_old[n*n];
    for (int i = 0; i < n * n; i++) {
        R_old[i] = A[i];
    }

    double Jij[n*n] = {};

    for (int j = 0; j<(n-1); j++){
        for (int i = j+1; i<n; i++){
            // construção matriz jacobi Jij
           matrizJacobiBaseadaNoElemento_ij_DeRvelha(n, i, j, R_old, Jij);

           multiMatrix(n, Jij, R_old, R_new);

           // R_old <- R_new
           for (int i = 0; i < n * n; i++) {
                R_old[i] = R_new[i];
            }

            // QT <- Jij * QT (acumula)
            double temp[n*n];
            multiMatrix(n, Jij, QT, temp);

            // QT <- temp
            for (int i = 0; i < n * n; i++) {
                QT[i] = temp[i];
            }
        }
    }
    
    // Q <- QT^t (transposta)
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            Q[i*n + j] = QT[j*n + i];
        }
    }
}

double calcVal(int n, double* A_new)
{
    double result = 0.0;
    for (int i =0; i<n; i++){
        for (int j =0; j<n; j++){
            if (i > j){
                result += (A_new[i*n + j] * A_new[i*n + j]);
            }
        }
    }
    return result;
}

void qrMethod(int n, double* A, double* A_new, double* P, double* lamb)
{
    double finalError = 0.000001; // tolerable error
    double val=100.0;

    // P <- I
    for(int i = 0; i<n; i++) P[i*n + i] = 1;

    double Q[n*n] = {}; // ortogonal
    double R[n*n] = {}; // triangular superior

    while (val > finalError){
        qrDecomposition(n, A, Q, R);

        // A_new = R*Q
        multiMatrix(n, R, Q, A_new);

        cout << "A_new Matrix (Partial Result) = " << endl;
        for (int i = 0; i < n; i++){
            cout << "|";
            for (int j = 0; j < n; j++){
                cout << A_new[i*n + j] << "|";
            }
            cout << endl;
        }
        cout << endl;

        // A <- A_new
        for (int i = 0; i < n * n; i++) {
            A[i] = A_new[i];
        }

        // P <- P * Q (acumula)
        double temp[n*n];
        multiMatrix(n, P, Q, temp);

        // P <- temp
        for (int i = 0; i < n * n; i++) {
            P[i] = temp[i];
        }

        // Verificar se a matriz A_new já é diagonal
        val = calcVal(n, A_new);
    }
    // lamb <- recebe os elementos da diagonal
    for(int i = 0; i<n; i++) lamb[i] = A_new[i*n+i];
}

int main()
{
    int n = 5;
    double A[5][5] = {
        {40, 8, 4, 2, 1},
        {8, 30, 12, 6, 2},
        {4, 12, 20, 1, 2},
        {2, 6, 1, 25, 4},
        {1, 2, 2, 4, 5}
    };
    cout << scientific << setprecision(10);

    double A_householder[n*n];
    double H[n*n];
    for(int i=0;i<n*n;i++) H[i]=0;

    for(int i=0;i<n;i++)
        H[i*n+i]=1;
    houseHolderMethod(n, (double*)A, A_householder, H);

    double A_final[n*n] = {}; // matriz digonal final
    double P[n*n] = {}; // P = Q1*Q2*...
    double lamb[n]; // autovalores de A
    qrMethod(n, A_householder, A_final, P, lamb);

    cout << "QR Matrix Final = " << endl;
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++){
            cout << A_final[i*n + j] << "|";
        }
        cout << endl;
    }
    cout << endl << endl;

    cout << "A eigenvalues = " << endl;
    cout << "|";
    for (int i = 0; i < n; i++){
        cout << lamb[i] << "|";
    }
    cout << endl;


    cout << endl << "A eigenvalues = " << endl;
    cout << "|";
    for (int i = 0; i < n; i++){
        cout << A_final[i*n+i] << "|";
    }
    cout << endl;

    cout << endl << "Eigenvectors of tridiagonal matrix = " << endl;
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++){
            cout << P[i*n + j] << "|";
        }
        cout << endl;
    }

    double HP[n*n];
    multiMatrix(n, H, P, HP);
    cout << endl << "Matrix H * P = " << endl;
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++){
            cout << HP[i*n + j] << "|";
        }
        cout << endl;
    }
    cout << endl;


    return 0;
}