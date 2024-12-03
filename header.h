#pragma once
using namespace std;

extern ofstream logfile;

namespace CONST
{
    extern const double M;
    extern const double P;
    extern const double k;
    extern const double qA;
    extern const double qB;
    extern const double E;
    extern const double J;
}

extern const double L;
extern const double xM;
extern const double xP;
extern const double xk;
extern const double xq1, xq2;
extern const double xR1, xR2, xR3, xR4;

class Solution
{
    //КЛАСС ДЛЯ РАСЧЕТА БАЛКИ 
    //N КОЛИЧЕСТВО УЗЛОВ
    //M КОЛИЧЕТСТВО КОНЕЧНЫХ ЭЛЕМЕНТОВ
    //x МАССИВ УЗЛОВ
    //A МАТРИЦА ЖЕСТКОСТИ
    //B ВЕКТОР НАГРУЗКИ
    //N2 ЧИСЛО СТЕПЕНЕЙ СВОБОДЫ
    //a МАТРИЦА ЖЕТСКОСТИ КОНЕЧНОГО ЭЛЕМЕНТА
    //b ВЕКТОР НАГРУЗКИ КОНЕЧНОГО ЭЛЕМЕНТА
    //c МАТРИЦА ПРОИЗВЕДЕНИЯ ВТОРЫХ ЧАСТНЫХ ПРОИЗВОДНЫХ
    //kk ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ПРУЖИНОЙ
    //km ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С МОМЕНТОМ
    //kp ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С СИЛОЙ
    //kq1 ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА В КОТОРОМ НАЧИНАЕТСЯ РАСПРЕДЕЛННАЯ НАГРУЗКА
    //kq2 ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА В КОТОРОМ ЗАКАНЧИВАЕТСЯ РАСПРЕДЕЛННАЯ НАГРУЗКА
    //Asol МАТРИЦА ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    //Bsol ВЕКТОР ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    //N3 РАЗМЕРНОСТЬ СИСТЕМЫ ПОДЛЕЖАЩЕЙ РЕШЕНИЮ
    //Ridx МАССИВ ИНДЕКСОВ ОПОР В СЕТКЕ
    //Csol ВЕКТОР РЕШЕНИЯ СИСТЕМЫ
    //C ВЕКТОР РАЗЛОЖЕНИЯ ПО БАЗИСУ

    int N;
    int M;
    int N2;
    int kk;
    int km;
    int kp;
    int kq1;
    int kq2;
    int N3;

    int Ridx[4];

    double *x;
    double *A;
    double *B;
    double *Asol;
    double *Bsol;
    double *Csol;
    double *C;

    double a[16];
    double c[16];
    double b[4];

public:

    Solution();
    Solution(int N);
    ~Solution();

    //info_base ВНОСИТ БАЗОВУЮ ИНФОРМАЦИЮ В LOG-ФАЙЛ
    //info_mesh ВНОСИТ ИНФОРМАЦИЮ О СЕТКЕ В LOG-ФАЙЛ
    //info_A ВНОСИТ ИНФОРМАЦИЮ О МАТРИЦЕ ЖЕСТКОСТИ В LOG-ФАЙЛ
    //info_B ВНОСИТ ИНФОРМАЦИЮ О ВЕКТОРЕ B
    //info_C ВНОСИТ ИНФОРМАЦИЮ О ВЕКТОРЕ C
    //info_a ВНОСИТ ИНФОРМАЦИЮ О МАТРИЦЕ ЖЕСТКОСТИ k-ГО КОНЕЧНОГО ЭЛЕМЕНТА
    //info_b ВНОСИТ ИНФОРМАЦИЮ О ВЕКТОРЕ НАГРУЗКИ k-ГО КОНЕЧНОГО ЭЛЕМЕНТА
    //info_Asol ВНОСИТ ИНФОРМАЦИЮ О МАТРИЦЕ Asol
    //info_Bsol ВНОСИТ ИНФОРМАЦИЮ О ВЕКТОРЕ Bsol
    //info_Csol ВНОСИТ ИНФОРМАЦИЮ О ВЕКТОРЕ Csol
    //nv_info ВНОСИТ ИНФОРМАЦИЮ О НОРМЕ ВЕКТОРА НЕВЯЗКИ

    void info_base();
    void info_mesh();
    void info_A();
    void info_B();
    void info_a(int k);
    void info_b(int k);
    void info_Asol();
    void info_Bsol();
    void info_Csol();
    void info_C();
    void info_nv();

    //create_mesh ПОСТРОЕНИЕ СЕТКИ
    //fill_A ЗАПОЛНЕНИЕ МАТРИЦЫ A
    //fill_B ЗАПОЛНЕНИЕ МАТРИЦЫ B
    //fill_c ЗАПОЛНЕНИЕ МАТРИЦЫ с
    //set_zeros_A КАЖДЫЙ ЭЛЕМЕНТ МАТРИЦЫ A ЭТО НОЛЬ
    //set_zeros_a КАЖДЫЙ ЭЛЕМЕНТ МАТРИЦЫ a ЭТО НОЛЬ
    //fill_Asol ЗАПОЛНЕНИЕ МАТРИЦЫ Asol
    //fill_Bsol ЗАПОЛНЕНИЕ ВЕКТОРА Bsol
    //solve РЕШАЕТ СИСТЕМУ УРАВНЕНИЙ
    //fill_txt ЗАПОЛНЯЕТ ФАЙЛ graph_val ЗНАЧЕНИЯМИ

    void create_mesh();
    void fill_A();
    void fill_B();
    void fill_c();
    void fill_Asol();
    void fill_Bsol();
    void set_zeros_A();
    void set_zeros_B();
    void set_zeros_a();
    void set_zeros_b();
    void solve();
    void fill_txt();

    //min_h МИНИМАЛЬНАЯ ДЛИНА В СЕТКЕ
    //opa БИЛИНЕЙНЫЙ ОПЕРАТОР НА ФУНКЦИЯХ С ИНДЕКСАМИ i И j НА ЭЛЕМЕНТЕ С ИНДЕКСОМ k
    //opb ИНТЕГРАЛ РАСПРЕДЕЛЕННОЙ НАГРУЗКИ НА ФУНКЦИИ С ИНДЕКСОМ i НА ЭЛЕМЕНТЕ С ИНДЕКСОМ k
    //w ВЫЧИСЛЯЕТ ЗНАЧЕНИЕ ПРОГИБА В ТОЧКЕ z ИЗ [0,L]
    //nv ВЫЧИСЛЯЕТ НОРМУ ВЕКТОРА НЕВЯЗКИ Asol*Csol - Bsol

    double min_h();
    double opa(int k, int i, int j);
    double opb(int k, int i);
    double w(double z);
    double nv();
};