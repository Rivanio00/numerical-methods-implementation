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

void multiMatrix(double* H, double* Hi, int n, double* result){
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

    // v(step+1 : n) <- coluna step de A, da linha step+1 até n
    for (int i = step+1; i < n; i++){
        v[i] = A[i*n + step];
    }

    // norma = ||v||
    double norma = 0.0;
    for (int i = step+1; i < n; i++) norma += v[i] * v[i];
    norma = sqrt(norma);

    // v_linha = zeros(n); v_linha(step+1) = norma
    double v_linha[n];
    for (int i = 0; i < n; i++) v_linha[i] = 0.0;
    v_linha[step+1] = norma;

    // w = v - v_linha
    double w[n];
    for (int i = 0; i < n; i++) w[i] = v[i] - v_linha[i];

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
    multiMatrix(Ai, Hi, n, temp);

    double transpose[n*n];
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            transpose[i*n + j] = Hi[j*n + i];
        }
    }

    double result[n*n];
    multiMatrix(transpose, temp, n, result);

    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            A[i*n + j] = result[i*n + j];
        }
    }
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
void converge(int n, double* M, double* result)
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
    
    cout << endl << "Result for problem " << n << "x" << n << " =" << endl;
    cout << "Final Eigenvalue: " << eigenvalue_new << endl;
    cout << "Final Error: " << error << endl;
    cout << "Final vector: " << endl << "|";
    for (int i = 0; i < n; i++) cout << v_old[i] << "|";
    cout << endl << "==========================" << endl;
    for (int i = 0; i < n; i++) result[i] = v_old[i];
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
}


// method
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
        multiMatrix(H, Hi, n, temp);
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

    cout << "Final Matrix H = " << endl;
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++){
            cout <<  H[i*n + j] << "|";
            H_out[i*n + j] = H[i*n + j];
        }
        cout << endl;
    }
}

double convergeInverse(int n, double* M, double* eigenvector)
{
    double v_old[n];
    double v_new[n];
    for (int i = 0; i < n; i++) v_new[i] = 1.0;

    double error = 1;
    double finalError = 0.000001;

    double eigenvalue_old = 0.0;
    double eigenvalue_new = 0.0;

    double U[n*n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U[i*n + j] = ((double(*)[n])M)[i][j];

    double L[n*n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            L[i*n + j] = (i == j) ? 1.0 : 0.0;

    decompLU(n, L, U);

    do {
        eigenvalue_old = eigenvalue_new;
        for (int i = 0; i < n; i++) v_old[i] = v_new[i];
        normalize(n, v_old);
        solveLU(n, L, U, v_old, v_new);
        eigenvalue_new = 0.0;
        for (int i = 0; i < n; i++) eigenvalue_new += v_old[i] * v_new[i];
        error = fabs((eigenvalue_new - eigenvalue_old) / eigenvalue_new);
    }while (error > finalError);
    
    double eigenvalue_final = (1/eigenvalue_new);

    // ← único acréscimo
    for (int i = 0; i < n; i++) eigenvector[i] = v_old[i];

    cout << "Final Eigenvalue: " << eigenvalue_final << endl;
    cout << "Final Error: " << error << endl;
    return eigenvalue_final;
}

double convergeShifted(int n, double* M, double mu, double* eigenvector)
{
    double A_shifted[n*n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_shifted[i*n + j] = M[i*n + j] - (i == j ? mu : 0.0);

    // chama convergeInverse e pega o autovetor
    double lambda = convergeInverse(n, A_shifted, eigenvector);
    return lambda + mu;
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

    double A_final[n*n];
    double H[n*n];
    houseHolderMethod(n, (double*)A, A_final, H);

    // autovalores dominante de A#
    double v[n];
    converge(n, A_final, v);

    // demais autovalores via shifted — inspecionar diagonal de A_final pra escolher mus
    double mus[n];
    for (int i = 0; i < n; i++) mus[i] = A_final[i*n + i] * 0.95;

    for (int k = 0; k < 4; k++){
        double v_shifted[n];
        convergeShifted(n, A_final, mus[k], v_shifted);

        // autovetor de A original
        double x[n];
        matVec(n, H, v_shifted, x);
    }

    // Incial A eigenvalues
    converge(n, (double*)A, v);

    return 0;
}