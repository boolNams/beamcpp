//ОСНОВНОЙ МОДУЛЬ ПРОГРАММЫ ДЛЯ РАСЧЕТА БАЛКИ
//ДЛИНА БАЛКИ | 25

#include <iostream>
#include <fstream>
#include "header.h"
using namespace std;

namespace CONST
{
    extern const double M = 0.5;
    extern const double P = 0.0;
    extern const double k = 0.0;
    extern const double qA = 1.0;
    extern const double qB = 0.0;
    extern const double E = 300.0;
    extern const double J = 10.0;
}

extern const double L = 25.0;
extern const double xM = 13.0;
extern const double xP = 20.0;
extern const double xk = 23.0;
extern const double xq1 = 10.0, xq2 = 22.0;
extern const double xR1 = 1.0, xR2 = 8.0, xR3 = 16.0, xR4 = 25.0;

int main(void)
{
    //КОЛИЧЕСТВО УЗЛОВ
    int N = 10;

    Solution sol(N);

    cout << "STOP C++" << endl;
    return 0;
}
