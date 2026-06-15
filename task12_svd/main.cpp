#include <iostream>
#include <cmath>
#include <algorithm>
#include "eigen_methods.h"
using namespace std;

void printMatrix(const char* name, double* M, int rows, int cols) {
    cout << name << " =" << endl;
    for (int i = 0; i < rows; i++) {
        cout << "|";
        for (int j = 0; j < cols; j++)
            cout << M[i * cols + j] << "|";
        cout << endl;
    }
}

void multiMatrixRect(double* A, int rA, int cA,
                     double* B, int rB, int cB,
                     double* C)
{
    for (int i = 0; i < rA; i++) {
        for (int j = 0; j < cB; j++) {
            C[i * cB + j] = 0.0;
            for (int k = 0; k < cA; k++) {
                C[i * cB + j] += A[i * cA + k] * B[k * cB + j];
            }
        }
    }
}

void transpose(double* M, int rows, int cols, double* Mt) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            Mt[j * rows + i] = M[i * cols + j];
}

int main() {
    // PASSO 1: LEITURA DE DADOS
    int m, n;
    cout << "Digite m linhas e n colunas: ";
    cin >> m >> n;

    double* A = new double[m * n];

    cout << "Digite a matriz A (" << m << "x" << n << "):" << endl;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            cin >> A[i * n + j];
    cout << endl;
    printMatrix("A", A, m, n);


    // PASSO 2: FORMAR A MATRIZ SIMÉTRICA A_bar
    double* At = new double[n * m];
    transpose(A, m, n, At);

    int p;
    double* A_bar;

    if (m <= n) {
        p = m;
        A_bar = new double[p * p];
        multiMatrixRect(A, m, n, At, n, m, A_bar);
    } else {
        p = n;
        A_bar = new double[p * p];
        multiMatrixRect(At, n, m, A, m, n, A_bar);
    }

    
    // PASSO 3: AUTOVALORES E AUTOVETORES DE A_bar
    double* A_tri = new double[p * p];
    double* H     = new double[p * p];

    for (int i = 0; i < p * p; i++) H[i] = 0.0;
    for (int i = 0; i < p; i++) H[i * p + i] = 1.0;

    houseHolderMethod(p, A_bar, A_tri, H);

    double* A_diag = new double[p * p];
    double* P      = new double[p * p];
    double* lamb   = new double[p];

    for (int i = 0; i < p * p; i++) P[i] = 0.0;

    qrMethod(p, A_tri, A_diag, P, lamb);

    double* HP = new double[p * p];
    multiMatrix(p, H, P, HP);


    // PASSO 4: VALORES SINGULARES
    double* sigma = new double[p];
    for (int i = 0; i < p; i++) {
        if (lamb[i] > 0)
            sigma[i] = sqrt(lamb[i]);
        else
            sigma[i] = 0.0;
    }


    // PASSO 5: POSTO DA MATRIZ
    double tol = 1e-10;
    int rank = 0;
    for (int i = 0; i < p; i++) {
        if (sigma[i] > tol) rank++;
    }
    cout << endl << "========================================" << endl;
    cout << "Posto da matriz A: " << rank << endl;

    
    // PASSO 6: MONTAR U, SIGMA E V

    // ordenar em ordem decrescente
    int* idx = new int[p];
    for (int i = 0; i < p; i++) idx[i] = i;

    sort(idx, idx + p, [&sigma](int a, int b) {
        return sigma[a] > sigma[b];
    });

    double* sigma_sorted = new double[p];
    double* lamb_sorted  = new double[p];
    double* HP_sorted    = new double[p * p];
    for (int i = 0; i < p; i++) {
        sigma_sorted[i] = sigma[idx[i]];
        lamb_sorted[i]  = lamb[idx[i]];
    }
    for (int j = 0; j < p; j++)
        for (int i = 0; i < p; i++)
            HP_sorted[i * p + j] = HP[i * p + idx[j]];

    for (int i = 0; i < p; i++) {
        sigma[i] = sigma_sorted[i];
        lamb[i]  = lamb_sorted[i];
    }
    for (int i = 0; i < p * p; i++) HP[i] = HP_sorted[i];

    delete[] sigma_sorted;
    delete[] lamb_sorted;
    delete[] HP_sorted;
    delete[] idx;

    cout << endl << "========================================" << endl;
    cout << "Valores singulares (ordenados):" << endl;
    for (int i = 0; i < p; i++)
        cout << "  sigma[" << i << "] = " << sigma[i] << endl;

    // montar Sigma (m x n)
    double* Sigma = new double[m * n];
    for (int i = 0; i < m * n; i++) Sigma[i] = 0.0;
    int minMN = (m < n) ? m : n;
    for (int i = 0; i < minMN; i++)
        Sigma[i * n + i] = sigma[i];
    cout << endl << "========================================" << endl;
    printMatrix("Sigma", Sigma, m, n);

    // montar U (m x m) e V (n x n)
    double* U = new double[m * m];
    double* V = new double[n * n];
    for (int i = 0; i < m * m; i++) U[i] = 0.0;
    for (int i = 0; i < n * n; i++) V[i] = 0.0;

    if (m <= n) {
        // HP contém colunas de U
        for (int i = 0; i < m * m; i++) U[i] = HP[i];

        // V(:,j) = (1/sigma_j) * At * U(:,j)
        for (int j = 0; j < p; j++) {
            if (sigma[j] > tol) {
                double* u_col = new double[m];
                for (int k = 0; k < m; k++) u_col[k] = U[k * m + j];

                double* v_col = new double[n];
                for (int k = 0; k < n; k++) {
                    v_col[k] = 0.0;
                    for (int l = 0; l < m; l++)
                        v_col[k] += At[k * m + l] * u_col[l];
                    v_col[k] /= sigma[j];
                }
                for (int k = 0; k < n; k++) V[k * n + j] = v_col[k];

                delete[] u_col;
                delete[] v_col;
            }
        }
    } else {
        // HP contém colunas de V
        for (int i = 0; i < n * n; i++) V[i] = HP[i];

        // U(:,j) = (1/sigma_j) * A * V(:,j)
        for (int j = 0; j < p; j++) {
            if (sigma[j] > tol) {
                double* v_col = new double[n];
                for (int k = 0; k < n; k++) v_col[k] = V[k * n + j];

                double* u_col = new double[m];
                for (int k = 0; k < m; k++) {
                    u_col[k] = 0.0;
                    for (int l = 0; l < n; l++)
                        u_col[k] += A[k * n + l] * v_col[l];
                    u_col[k] /= sigma[j];
                }
                for (int k = 0; k < m; k++) U[k * m + j] = u_col[k];

                delete[] v_col;
                delete[] u_col;
            }
        }
    }

    cout << endl << "========================================" << endl;
    printMatrix("U", U, m, m);
    cout << endl;
    printMatrix("V", V, n, n);


    // PASSO 7: VERIFICAÇÃO U * Sigma * Vt = A
    double* Vt = new double[n * n];
    transpose(V, n, n, Vt);

    double* SigmaVt = new double[m * n];
    multiMatrixRect(Sigma, m, n, Vt, n, n, SigmaVt);

    double* result = new double[m * n];
    multiMatrixRect(U, m, m, SigmaVt, m, n, result);

    cout << endl << "========================================" << endl;
    printMatrix("U * Sigma * Vt", result, m, n);
    cout << endl;

    double frobenius = 0.0;
    for (int i = 0; i < m * n; i++) {
        double diff = A[i] - result[i];
        frobenius += diff * diff;
    }
    frobenius = sqrt(frobenius);
    cout << "||A - U*Sigma*Vt||_F = " << frobenius << endl;

    if (frobenius < 1e-6)
        cout << "Verificacao OK! A decomposição está correta." << endl;
    else
        cout << "Erro grande demais. Verifique a implementação!" << endl;


    // limpeza de memória
    delete[] A;
    delete[] At;
    delete[] A_bar;
    delete[] A_tri;
    delete[] H;
    delete[] A_diag;
    delete[] P;
    delete[] lamb;
    delete[] HP;
    delete[] sigma;
    delete[] Sigma;
    delete[] U;
    delete[] V;
    delete[] Vt;
    delete[] SigmaVt;
    delete[] result;
    return 0;
}