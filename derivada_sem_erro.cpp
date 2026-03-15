/******************************************************************************

Welcome to GDB Online.
  GDB online is an online compiler and debugger tool for C, C++, Python, PHP, Ruby, 
  C#, OCaml, VB, Perl, Swift, Prolog, Javascript, Pascal, COBOL, HTML, CSS, JS
  Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip> 

using namespace std;

double fx(double x) 
{
    return pow(x, 4) + 2*pow(x, 3) - pow(x, 2) + 5*x + sin(x);
}

double forward(double x, double deltax, int ordem)
{
    if (ordem == 1)
    {
        return (fx(x+deltax)-fx(x))/deltax;
    }
    return (forward(x+deltax, deltax, ordem-1) - forward(x, deltax, ordem-1))/deltax;
}

double backward(double x, double deltax, int ordem)
{
    if (ordem == 1)
    {
        return (fx(x)-fx(x-deltax))/deltax;
    }
    return (backward(x, deltax, ordem-1) - backward(x-deltax, deltax, ordem-1))/deltax;
}

double central(double x, double deltax, int ordem)
{
    if (ordem == 1)
    {
        return (fx(x+deltax)-fx(x-deltax))/(2*deltax);
    }
    return (central(x+deltax, deltax, ordem-1) - central(x-deltax, deltax, ordem-1))/(2*deltax);
}

int main()
{
    cout << fixed << setprecision(5);
    int philosophy;
    cout << "Bem-vindo! Antes de começarmos, por favor escolha a sua filosofia de orientação: " << endl;
    cout << "(0) Forward; (1) Backward; (2) Central; (-1) Sair" << endl;
    cin >> philosophy;
    if (philosophy == -1) 
    {
        return 0;
    }
    double x;
    double deltaX;
    int ordem;
    cout << "Boa escolha! Qual X você gostaria de usar?" << endl;
    cin >> x;
    cout << "E qual será o delta x?" << endl;
    cin >> deltaX;
    cout << "Qual será a ordem da derivada?" << endl;
    cin >> ordem;
    cout << "Calculando..." << endl;
    double result;
    switch(philosophy)
    {
        case 0: 
            result = forward(x, deltaX, ordem);
            break;
        case 1: 
            result = backward(x, deltaX, ordem);
            break;
        case 2:
            result = central(x, deltaX, ordem);
            break;
        default:
            result = -1;
            break; 
    }
    cout << "Seu resultado é: " << result << endl;
    return 0;
}