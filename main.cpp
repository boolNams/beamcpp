//ОСНОВНОЙ МОДУЛЬ ПРОГРАММЫ ДЛЯ РАСЧЕТА БАЛКИ
//ДЛИНА БАЛКИ | 25

#include <iostream>
#include <fstream>
#include <string>
#include "header.h"
using namespace std;

namespace CONST
{
    double M;
    double P;
    double k;
    double qA;
    double qB;
    double E;
    double J;
}

extern const double L = 25.0;
extern const double xM = 13.0;
extern const double xP = 20.0;
extern const double xk = 23.0;
extern const double xq1 = 10.0, xq2 = 22.0;
extern const double xR1 = 1.0, xR2 = 8.0, xR3 = 16.0, xR4 = 25.0;

int main(void)
{
    //ЗАГРУЗКА ЗНАЧЕНИЙ
    load_params();

    //КОЛИЧЕСТВО УЗЛОВ
    int N = 20;

    Solution sol(N);

    cout << "STOP C++" << endl;
    return 0;
}

void load_params()
{
    //ЗАГРУЖАЕТ ЗНАЧЕНИЯ С ФАЙЛА "params.txt"

    //МАССИВ ДЛЯ ЗНАЧЕНИЙ
    float data[7];

    //ОТКРЫТИЕ ФАЙЛА
    ifstream file("params.txt");

    //СЧИТЫВАНИЕ ДАННЫХ
    string line;
    int i = 0;
    while(getline(file, line))
    {
        int pos = line.find('|');
        string str_value = line.substr(0, pos);
        data[i] = stof(str_value);
        i += 1;
    }

    //ПРИСВАИВАНИЕ ЗНАЧЕНИЙ
    CONST::P = data[0];
    CONST::M = data[1];
    CONST::k = data[2];
    CONST::qA = data[3];
    CONST::qB = data[4];
    CONST::E = data[5];
    CONST::J = data[6];
}
