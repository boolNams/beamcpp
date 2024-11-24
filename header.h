#pragma once
using namespace std;

extern ofstream logfile;

extern const double L;
extern const double M;
extern const double P;
extern const double k;
extern const double qA;
extern const double qB;
extern const double E;
extern const double J;

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
    //A ПЕРВОНАЧАЛЬНАЯ МАТРИЦА ЖЕСТКОСТИ
    //N2 ЧИСЛО СТЕПЕНЕЙ СВОБОДЫ
    //a МАТРИЦА ЖЕТСКОСТИ КОНЕЧНОГО ЭЛЕМЕНТА
    //c МАТРИЦА ПРОИЗВЕДЕНИЯ ВТОРЫХ ЧАСТНЫХ ПРОИЗВОДНЫХ
    //kk ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ПРУЖИНОЙ

    int N;
    int M;
    int N2;
    int kk;

    double *x;
    double *A;

    double a[16];
    double c[16];

public:

    Solution();
    Solution(int N);
    ~Solution();

    //info_base ВНОСИТ БАЗОВУЮ ИНФОРМАЦИЮ В LOG-ФАЙЛ
    //info_mesh ВНОСИТ ИНФОРМАЦИЮ О СЕТКЕ В LOG-ФАЙЛ
    //info_A ВНОСИТ ИНФОРМАЦИЮ О ПЕРВОНАЧАЛЬНОЙ МАТРИЦЕ ЖЕСТКОСТИ В LOG-ФАЙЛ
    //info_a ВНОСИТ ИНФОРМАЦИЮ О МАТРИЦЕ ЖЕСТКОСТИ k-ГО КОНЕЧНОГО ЭЛЕМЕНТА

    void info_base();
    void info_mesh();
    void info_A();
    void info_a(int k);

    //create_mesh ПОСТРОЕНИЕ СЕТКИ
    //fill_A ЗАПОЛНЕНИЕ МАТРИЦЫ A
    //fill_c ЗАПОЛНЕНИЕ МАТРИЦЫ с
    //set_zeros_A КАЖДЫЙ ЭЛЕМЕНТ МАТРИЦЫ A ЭТО НОЛЬ
    //set_zeros_a КАЖДЫЙ ЭЛЕМЕНТ МАТРИЦЫ a ЭТО НОЛЬ

    void create_mesh();
    void fill_A();
    void fill_c();
    void set_zeros_A();
    void set_zeros_a();

    //min_h МИНИМАЛЬНАЯ ДЛИНА В СЕТКЕ
    //op БИЛИНЕЙНЫЙ ОПЕРАТОР НА ФУНКЦИЯХ С ИНДЕКСАМИ i И j НА ЭЛЕМЕНТЕ С ИНДЕКСОМ k

    double min_h();
    double op(int k, int i, int j);
};