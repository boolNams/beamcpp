//МОДУЛЬ ФУНКЦИОНАЛА РЕШЕНИЯ
//ВСЕ ВЫЧИСЛЕНИЯ В LOG-ФАЙЛЕ beamlog.txt

#include <iostream>
#include <fstream>
#include "math.h"
#include "header.h"
using namespace std;

//ОПРЕДЕЛЕНИЕ ПОТОКА LOG-ФАЙЛА
ofstream logfile;

//БАЗИСНЫЕ ФУНКЦИИ НА [0,1]
double N1(double y){return 2.0*pow(y,3) - 3.0*pow(y,2) + 1.0;}
double L1(double y){return pow(y,3) - 2.0*pow(y,2) + y;}
double N2(double y){return -2.0*pow(y,3) + 3.0*pow(y,2);}
double L2(double y){return pow(y,3) - pow(y,2);}

//МАССИВ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*BASIS[4])(double) = {N1, L1, N2, L2};

Solution::Solution(){}

Solution::Solution(int N)
{
	//ОТКРЫТИЕ LOG ФАЙЛА
    logfile.open("beamlog.txt");

	//МИНИМАЛЬНОЕ ЧИСЛО УЗЛОВ 7 (2 УЗЛА РАСПРЕДЕЛЕННОЙ НАГРУЗКИ, 4 ОПОРЫ, 1 УЗЕЛ НА ЛЕВОМ ТОРЦЕ)
    if(N < 7)
    {
        logfile << "ЧИСЛО УЗЛОВ НЕ МОЖЕТ БЫТЬ МЕНЬШЕ 7!!!" << endl;
        logfile.close();
        throw -1;
    }

    //ПРИСВАИВАНИЕ ЗНАЧЕНИЙ И ВЫДЕЛЕНИЕ ПАМЯТИ
    this->N = N;
    this->M = N - 1;
    this->N2 = 2*N;
    this->x = new double[N];
    this->A = new double[N2*N2];

    //ЗАПОЛНЕНИЕ МАТРИЦЫ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ
    fill_c();

    //ВЫВОД ПЕРВОНАЧАЛЬНЫХ ДАННЫХ
    info_base();

    //ПОСТРОЕНИЕ СЕТКИ
    create_mesh();

    //ВЫВОД ДАННЫХ О СЕТКЕ
    info_mesh();

    //ЗАПОЛНЕНИЕ ПЕРВОНАЧАЛЬНОЙ МАТРИЦЫ ЖЕСТКОСТИ	
    fill_A();

    //ВЫВОД ДАННЫХ О ПЕРВОНАЧАЛЬНОЙ МАТРИЦЕ ЖЕСТКОСТИ
    info_A();

    //ЗАКРЫТИЕ LOG ФАЙЛА
    logfile.close();
}

Solution::~Solution()
{
	//ЧИСТКА ПАМЯТИ
	if(x != nullptr) delete x;
	if(A != nullptr) delete A;
}

void Solution::info_base()
{
	logfile << "### ДАННЫЕ О РЕШЕНИИ ЗАДАЧИ ###" << endl;
	logfile << "МИНИМАЛЬНО ДОПУСТИМОЕ ЧИСЛО УЗЛОВ | N_min = 7" << endl;
	logfile << "ЧИСЛО УЗЛОВ | N = " << N << endl;
	logfile << "ЧИСЛО КОНЕЧНЫХ ЭЛЕМЕНТОВ | M = " << M << endl;
}

void Solution::info_mesh()
{
	logfile << "МИНИМАЛЬНЫЙ РАЗМЕР ЭЛЕМЕНТА В СЕТКЕ | min_h = " << min_h() << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С ПРУЖИНОЙ | kk = " << kk + 1 << endl;

	logfile << "КООРДИНАТЫ ТОЧЕК В СЕТКЕ::" << endl;
	for(int i = 0; i < N; ++i)
	{
		logfile << x[i] << endl;
	}
}

void Solution::info_A()
{
	logfile << "ПЕРВОНАЧАЛЬНАЯ МАТРИЦА ЖЕСТКОСТИ A::" << endl;
	for(int i = 0; i < N2; ++i)
	{
		for(int j = 0; j < N2; ++j) logfile << A[i*N2 + j] << " ";
		logfile << endl; 
	}
}

void Solution::info_a(int k)
{
	logfile << "МАТРИЦА ЖЕСТКОСТИ " << k << "-ГО КОНЕЧНОГО ЭЛЕМЕНТА::" << endl;
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j) logfile << a[i*4+j] << " ";
		logfile << endl;
	}

}

double Solution::min_h()
{
	double min_h = L;
	for(int i = 0; i < N - 1; ++i)
	{
		if((x[i+1] - x[i]) < min_h) min_h = x[i+1] - x[i];
	}
	return min_h;
}

void Solution::create_mesh()
{	
	//МИНИМАЛЬНОЕ ЧИСЛО УЗЛОВ 7 (2 УЗЛА РАСПРЕДЕЛЕННОЙ НАГРУЗКИ, 4 ОПОРЫ, 1 УЗЕЛ НА ЛЕВОМ ТОРЦЕ)
	//СТРОИТСЯ РАВНОМЕРНАЯ СЕТКА, ВКЛЮЧАЮЩАЯ В СЕБЯ N-6 УЗЛОВ (БЕЗ 2 УЗЛОВ РАСПР. НАГР. И 4 УЗЛОВ ОПОРЫ)
	//УЗЕЛ НА ПРАВОМ ТОРЦЕ ЭТО УЗЕЛ ОДНОЙ ИЗ ОПОР

	double h = L / (N - 6);
    for(int i = 0; i < N - 6; ++i) x[i] = i*h;
    x[N-6] = xR1;
   	x[N-5] = xR2;
    x[N-4] = xR3;
	x[N-3] = xR4;
   	x[N-2] = xq1;
   	x[N-1] = xq2;

   	//СОРТИРОВКА ПУЗЫРЬКОМ
   	double tmp;
   	for(int i = 0; i < N - 1; ++i)
   	{
   		for(int j = 0; j < N - 1 - i; ++j)
   		{
   			if(x[j] > x[j+1])
   			{
   				tmp = x[j];
   				x[j] = x[j+1];
   				x[j+1] = tmp;
   			}
   		}
   	}

   	//УСТРАНЕНИЕ "СЛИПАНИЯ" ТОЧЕК 
   	for(int i = 1; i < N - 1; ++i)
   	{
   		if(x[i+1] - x[i] < 1.0e-16) x[i] = 0.5*(x[i-1] + x[i+1]);
   	}

   	//ПРОВЕРКА КАЧЕСТВА СЕТКИ
   	if(min_h() < 1.0e-16)
   	{
   		logfile << "НЕККОРЕТКТНО ПОСТРОЕНА СЕТКА!!!" << endl;
   		logfile.close();
        throw -1;
   	}

   	//ПРИСВАИВАНИЕ ЗНАЧЕНИЯ kk ИНДЕКСУ КОНЕЧНОГО ЭЛЕМЕНТА С ПРУЖИНОЙ
   	for(int k = 0; k < M; ++k)
   	{
   		if(xk < x[k+1] and xk > x[k]) kk = k;
   	}
}

void Solution::fill_c()
{
	//МАТРИЦА с ИСПОЛЬЗУЕТСЯ В МЕТОДЕ op ДЛЯ ВЫЧИСЛЕНИЯ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ
	//ДАННЫЙ МЕТОД ЗАПОЛНЯЕТ МАТРИЦУ с ПОЛУЧЕННЫМИ ЗНАЧЕНИЯМИ ИЗ SAGE ПРОГРАММЫ

	c[0]  =  12.0; c[1]  =  6.0; c[2]  = -12.0; c[3]  =  6.0;
	c[4]  =   6.0; c[5]  =  4.0; c[6]  =  -6.0; c[7]  =  2.0;
	c[8]  = -12.0; c[9]  = -6.0; c[10] =  12.0; c[11] = -6.0;
	c[12] =   6.0; c[13] =  2.0; c[14] =  -6.0; c[15] =  4.0;
}

void Solution::fill_A()
{
	//ВЫЧИСЛЕНИЕ ПЕРВОНАЧАЛЬНОЙ МАТРИЦЫ ЖЕСТКОСТИ
	set_zeros_A();

	//ЦИКЛ ПО КОНЕЧНЫМ ЭЛЕМЕНТАМ
	for(int k = 0; k < M; ++k)
	{
		//ЗАНУЛЕНИЕ МАТРИЦЫ ЖЕСТКОСТИ КОНЕЧНОГО ЭЛЕМЕНТА
		set_zeros_a();

		//ВЫЧИСЛЕНИЕ МАТРИЦЫ ЖЕТСКОСТИ КОНЕЧНОГО ЭЛЕМЕНТА
		for(int i = 0; i < 4; ++i)
		{
			//ФУНКЦИЯ op ПО ИНДЕКСАМ ФУНКЦИЙ ВЫЧИСЛИТ ИХ БИЛИНЕЙНЫЙ ОПЕРАТОР НА (k+1)-ОМ ЭЛЕМЕНТЕ
			for(int j = 0; j < i + 1; ++j) a[i*4+j] = op(k, i, j);

			//ЗАПОЛНЕНИЕ ВЕРХНЕЙ ЧАСТИ МАТРИЦЫ СИММЕТРИЧНО
			for(int j = 0; j < i; ++j) a[j*4+i] = a[i*4+j];
		}
		//ВЫЧИСЛЕНА МАТРИЦА ЖЕСТКОСТИ КОНЕЧНОГО ЭЛЕМЕНТА

		//ВЫВОД МАТРИЦЫ ЖЕСТКОСТИ ДЛЯ (k+1)-ГО КОНЕЧНОГО ЭЛЕМЕНТА
		info_a(k+1);

		//ЗАПОЛНЕНИЕ ЧАСТИ МАТРИЦЫ ЖЕСТКОСТИ СООТВЕСТВУЮЩЕЙ (k+1)-МУ КОНЕЧНОМУ ЭЛЕМЕНТУ
		for(int i = 0; i < 4; ++i)
		{
			//ГЛОБАЛЬНЫЕ НОМЕРА ФУНКЦИЙ ДЛЯ БИЛИНЕЙНОГО ОПЕРАТОРА (i,j) -> (2k+i, 2k+j)
			for(int j = 0; j < 4; ++j) A[(2*k+i)*N2 + 2*k+j] += a[i*4+j];
		}
		//ЗАПОЛНЕНА ЧАСТЬ МАТРИЦЫ ЖЕСТКОСТИ СООТВЕТСТВУЮЩАЯ (k+1)-МУ КОНЕЧНОМУ ЭЛЕМЕНТУ
	}

}

void Solution::set_zeros_A()
{
	for(int i = 0; i < N2*N2; ++i) A[i] = 0.0;
}

void Solution::set_zeros_a()
{
	for(int i = 0; i < 16; ++i) a[i] = 0.0;
}

double Solution::op(int k, int i, int j)
{
	// ВЫЧИСЛЕНИЕ БИЛИНЕЙНОГО ФУНКЦИОНАЛА
	// k ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА
	// (i,j) ИНДЕКСЫ БАЗИСНЫХ ФУНКЦИЙ

	// ДОМНОЖАЕМ ПРОИЗВЕДЕНИЕ ФУНКЦИЙ НА dot
	double dot = E * J / pow(x[k+1] - x[k], 3);

	// ПРИБАВЛЯЕМ СЛАГАЕМОЕ ОТ ПРУЖИНЫ ЕСЛИ РАБОТАЕМ НА КОНЕЧНОМ ЭЛЕМЕНТЕ С ПРУЖИНОЙ
	if(k == kk)
	{
		//ТОЧКА В КОТОРОЙ НЕОБХОДИМО ВЫЧИСЛИТЬ ЗНАЧЕНИЕ ФУНКЦИЙ
		double y = (xk - x[k])/(x[k+1] - x[k]);

		//ДОБАВЛЯЕМ К ПРОИЗВЕДЕНИЮ ЧЛЕН ОТВЕЧАЮЩИЙ ЗА ПРУЖИНУ
		return dot * c[i*4 + j] + k*BASIS[i](y)*BASIS[j](y);
	}

	// ИНАЧЕ ТОЛЬКО ПРОИЗВЕДЕНИЕ
	return dot * c[i*4 + j];

}