//МОДУЛЬ ФУНКЦИОНАЛА РЕШЕНИЯ
//ВСЕ ВЫЧИСЛЕНИЯ В LOG-ФАЙЛЕ beamlog.txt

#include <iostream>
#include <iomanip>
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

//ПЕРВЫЕ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double dN1(double y){return 6.0*pow(y,2) - 6.0*y;}
double dL1(double y){return 3.0*pow(y,2) - 4.0*y + 1.0;}
double dN2(double y){return -6.0*pow(y,2) + 6.0*y;}
double dL2(double y){return 3.0*pow(y,2) - 2.0*y;}

//МАССИВ ПЕРВЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*dBASIS[4])(double) = {dN1, dL1, dN2, dL2};

Solution::Solution(){}

Solution::Solution(int N)
{
	//ОТКРЫТИЕ LOG ФАЙЛА
    logfile.open("beamlog.txt");
    logfile.precision(3);

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
    this->N3 = N2 - 4;

    this->x = new double[N];
    this->A = new double[N2*N2];
    this->B = new double[N2];
    this->C = new double[N2];
    this->Asol = new double[N3*N3];
    this->Bsol = new double[N3];
    this->Csol = new double[N3];

    //ЗАПОЛНЕНИЕ МАТРИЦЫ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ
    fill_c();

    //ВЫВОД ПЕРВОНАЧАЛЬНЫХ ДАННЫХ
    info_base();

    //ПОСТРОЕНИЕ СЕТКИ
    create_mesh();

    //ВЫВОД ДАННЫХ О СЕТКЕ
    info_mesh();

    //ЗАПОЛНЕНИЕ МАТРИЦЫ ЖЕСТКОСТИ	
    fill_A();

    //ВЫВОД ДАННЫХ О МАТРИЦЕ ЖЕСТКОСТИ
    info_A();

    //ЗАПОЛНЕНИЕ ВЕКТОРА НАГРУЗКИ B
    fill_B();

    //ВЫВОД ДАННЫХ О ВЕКТОРЕ НАГРУЗКИ B
    info_B();

    //ЗАПОЛНЕНИЕ МАТРИЦЫ ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    fill_Asol();

    //ВЫВОД ДАННЫХ О МАТРИЦЕ ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    info_Asol();

    //ЗАПОЛНЕНИЕ ВЕКТОРА ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    fill_Bsol();

    //ВЫВОД ДАННЫХ О ВЕКТОРЕ ДЛЯ РЕШЕНИЯ СИСТЕМЫ
    info_Bsol();

    //РЕШЕНИЕ СИСТЕМЫ
    solve();

    //ВЫВОД ВЕКТОРА РЕШЕНИЯ СИСТЕМЫ Csol
    info_Csol();

    //ВЫВОД ВЕКТОРА РАЗЛОЖЕНИЯ ПО БАЗИСУ C
    info_C();

    //ЗАКРЫТИЕ LOG ФАЙЛА
    logfile.close();
}

Solution::~Solution()
{
	//ЧИСТКА ПАМЯТИ
	if(x != nullptr) delete x;
	if(A != nullptr) delete A;
	if(B != nullptr) delete B;
	if(C != nullptr) delete C;
	if(Asol != nullptr) delete Asol;
	if(Bsol != nullptr) delete Bsol;
	if(Csol != nullptr) delete Csol;
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
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С МОМЕНТОМ | km = " << km + 1 << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С СИЛОЙ | kp = " << kp + 1 << endl;
	logfile << "НОМЕР УЗЛА ПЕРВОЙ ОПОРЫ | Ridx = " << Ridx[0] + 1 << endl;
	logfile << "НОМЕР УЗЛА ВТОРОЙ ОПОРЫ | Ridx = " << Ridx[1] + 1 << endl;
	logfile << "НОМЕР УЗЛА ТРЕТЬЕЙ ОПОРЫ | Ridx = " << Ridx[2] + 1 << endl;
	logfile << "НОМЕР УЗЛА ЧЕТВЕРТОЙ ОПОРЫ | Ridx = " << Ridx[3] + 1 << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА В КОТОРОМ НАЧИНАЕТСЯ РАСПРЕДЕЛЕННАЯ НАГРУЗКА | kq1 = " << kq1 + 1 << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА В КОТОРОМ ЗАКАНЧИВАЕТСЯ РАСПРЕДЕЛЕННАЯ НАГРУЗКА | kq2 = " << kq2 + 1 << endl;

	logfile << "КООРДИНАТЫ ТОЧЕК В СЕТКЕ::" << endl;
	for(int i = 0; i < N; ++i)
	{
		logfile << x[i] << endl;
	}
}

void Solution::info_A()
{
	logfile << "ПОЛНАЯ МАТРИЦА ЖЕСТКОСТИ A::" << endl;
	for(int i = 0; i < N2; ++i)
	{
		for(int j = 0; j < N2; ++j)
		{
			logfile << setw(10) << A[i*N2 + j] << " ";
		}

		logfile << endl; 
	}
}

void Solution::info_B()
{
	logfile << "ПОЛНЫЙ ВЕКТОР НАГРУЗКИ B::" << endl;
	for(int i = 0; i < N2; ++i) logfile << B[i] << endl;
}

void Solution::info_C()
{
	logfile << "ВЕКТОР РАЗЛОЖЕНИЯ ПО БАЗИСУ C::" << endl;
	for(int i = 0; i < N2; ++i) logfile << C[i] << endl;
}

void Solution::info_a(int k)
{
	logfile << "МАТРИЦА ЖЕСТКОСТИ " << k << "-ГО КОНЕЧНОГО ЭЛЕМЕНТА::" << endl;
	for(int i = 0; i < 4; ++i)
	{
		for(int j = 0; j < 4; ++j)
		{
			logfile << setw(10) << a[i*4 + j] << " ";
		}

		logfile << endl;
	}

}

void Solution::info_b(int k)
{
	logfile << "ВЕКТОР НАГРУЗКИ " << k << "-ГО КОНЕЧНОГО ЭЛЕМЕНТА::" << endl;
	for(int i = 0; i < 4; ++i) logfile << b[i] << endl;
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

void Solution::info_Asol()
{
	logfile << "МАТРИЦА ДЛЯ РЕШЕНИЯ СИСТЕМЫ Asol::" << endl;
	for(int i = 0; i < N3; ++i) 
	{
		for(int j = 0; j < N3; ++j)
		{
			logfile << setw(10) << Asol[i*N3 + j] << " ";
		}

		logfile << endl;
	}

}

void Solution::info_Bsol()
{
	logfile << "ВЕКТОР ДЛЯ РЕШЕНИЯ СИСТЕМЫ Bsol::" << endl;
	for(int i = 0; i < N3; ++i) logfile << Bsol[i] << endl;
}

void Solution::info_Csol()
{
	logfile << "ВЕКТОР РЕШЕНИЯ СИСТЕМЫ Csol::" << endl;
	for(int i = 0; i < N3; ++i) logfile << Csol[i] << endl;
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

   	//ПРИСВАИВАНИЕ ЗНАЧЕНИЯ kk km kp ИНДЕКСАМ КОНЕЧНЫХ ЭЛЕМЕНТОВ С ПРУЖИНОЙ МОМЕНТОМ И СИЛОЙ
   	for(int k = 0; k < M; ++k)
   	{
   		//ЕСЛИ ОДНО ИЗ УСЛОВИЙ ЛЕЖИТ ПРЯМО НА УЗЛЕ ТО СЧИТАЕМ ЧТО УСЛОВИЕ ПРЕНАДЛЕЖИТ КОНЕЧНОМУ ЭЛЕМЕНТУ ЛЕВЫЙ УЗЕЛ КОТОРОГО ЯВЛЯЕТСЯ ДАННЫМ
   		if((xk < x[k+1] and xk > x[k]) or abs(xk - x[k]) < 1.0e-16) kk = k;
   		if((xM < x[k+1] and xM > x[k]) or abs(xM - x[k]) < 1.0e-16) km = k;
   		if((xP < x[k+1] and xP > x[k]) or abs(xP - x[k]) < 1.0e-16) kp = k;
   	}

   	//ПРИСВАИВАНИЕ ЗНАЧЕНИЯ kq1 kq2 ИНДЕКСАМ КОНЕЧНЫХ ЭЛЕМЕНТОВ ОПРЕДЕЛЯЮЩИХ НАЧАЛО И КОНЕЦ РАСПРЕДЕЛЕННОЙ НАГРУЗКИ
   	for(int k = 0; k < M; ++k)
   	{
   		if(abs(xq1 - x[k]) < 1.0e-16) kq1 = k;
   		if(abs(xq2 - x[k]) < 1.0e-16) kq2 = k;
   	}

   	//ПРИСВАИВАНИЕ ЗНАЧЕНИЙ МАССИВУ ИНДЕКСОВ ОПОР В СЕТКЕ
   	for(int i = 0; i < N; ++i)
   	{
   		if(abs(x[i] - xR1) < 1.0e-16) Ridx[0] = i;
   		if(abs(x[i] - xR2) < 1.0e-16) Ridx[1] = i;
   		if(abs(x[i] - xR3) < 1.0e-16) Ridx[2] = i;
   		if(abs(x[i] - xR4) < 1.0e-16) Ridx[3] = i;
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
	//ВЫЧИСЛЕНИЕ МАТРИЦЫ ЖЕСТКОСТИ
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
			for(int j = 0; j < i + 1; ++j) a[i*4+j] = opa(k, i, j);

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

void Solution::fill_B()
{
	//ВЫЧИСЛЕНИЕ ВЕКТОРА НАГРУЗКИ
	set_zeros_B();

	//ЛОГИЧЕСКАЯ ПЕРЕМЕННАЯ ДЛЯ ОТСЛЕЖИВАНИЯ КОНЕЧНЫХ ЭЛЕМЕНТОВ С РАСПРЕДЕЛЕННОЙ НАГРУЗКОЙ
	//ИЗНАЧАЛЬНО СЧИТАЕМ ЧТО НАГРУЗКИ НЕТ
	bool is_q = false;

	//ЦИКЛ ПО КОНЕЧНЫМ ЭЛЕМЕНТАМ
	for(int k = 0; k < M; ++k)
	{
		//ЗАНУЛЕНИЕ ВЕКТОРА КОНЕЧНОГО ЭЛЕМЕНТА
		set_zeros_b();

		//НАЧАЛАСЬ РАСПРЕДЕЛЕННАЯ НАГРУЗКА
		if(k == kq1) is_q = true;

		//ЗАКОНЧИЛАСЬ РАСПРЕДЕЛЕННАЯ НАГРУЗКА
		if(k == kq2) is_q = false;

		//ВЫЧИСЛЕНИЕ ВЕКТОРА КОНЕЧНОГО ЭЛЕМЕНТА

		//СЧИТАЕМ ИНТЕГРАЛ ЕСЛИ ЕСТЬ РАСПРЕДЕЛЕННАЯ НАГРУЗКА
		if(is_q)
		{
			for(int i = 0; i < 4; ++i) b[i] += opb(k, i);
		}

		//ДОБАВЛЯЕМ МОМЕНТ ЕСЛИ ОН ЕСТЬ НА КОНЕЧНОМ ЭЛЕМЕНТЕ
		if(k == km)
		{
			for(int i = 0; i < 4; ++i) b[i] += (1.0/(x[k+1] - x[k]))*M*dBASIS[i]((xM - x[k])/(x[k+1] - x[k]));
		}

		//ДОБАВЛЯЕМ СИЛЫ ЕСЛИ ОНА ЕСТЬ НА КОНЕЧНОМ ЭЛЕМЕНТЕ
		if(k == kp)
		{
			for(int i = 0; i < 4; ++i) b[i] -= P*BASIS[i]((xP - x[k])/(x[k+1] - x[k]));
		}

		//ВЫЧИСЛЕН ВЕКТОР КОНЕЧНОГО ЭЛЕМЕНТА

		//ВЫВОД ВЕКТОРА (k+1)-ГО КОНЕЧНОГО ЭЛЕМЕНТА
		info_b(k+1);

		//ЗАПОЛНЕНИЕ ЧАСТИ ВЕКТОРА НАГРУЗКИ СООТВЕТСТВУЮЩЕЙ (k+1)-МУ КОНЕЧНОМУ ЭЛЕМЕНТУ
		for(int i = 0; i < 4; ++i)
		{
			B[2*k+i] += b[i];
		}
		//ЧАСТЬ ВЕКТОРА СООТВЕТСТВУЮЩАЯ (k+1)-МУ КОНЕЧНОМУ ЭЛЕМЕНТУ ЗАПОЛНЕНА
	}

}

void Solution::fill_Asol()
{
	//ИМЕЕТСЯ 4 ОПОРЫ СЛЕДОВАТЕЛЬНО НЕОБХОДИМО УДАЛИТЬ ИЗ МАТРИЦЫ ЖЕСТКОСТИ A 4 СТРОКИ И 4 СТОЛБЦА
	//РАЗМЕРНОСТЬ МАТРИЦЫ Asol (N3xN3) ГДЕ N3 = N2 - 4

	//МАССИВ ИНДЕКСОВ СТОЛБЦОВ И СТРОК ДЛЯ УДАЛЕНИЯ
	int to_del[4] = {2*Ridx[0], 2*Ridx[1], 2*Ridx[2], 2*Ridx[3]};

	//СЧЕТЧИКИ УДАЛЕННЫХ СТРОК И СТОЛБЦОВ СООТВЕТСТВЕННО
	int counti = 0, countj = 0;

	//ЦИКЛ ПО СТРОКАМ МАТРИЦЫ A
	for(int i = 0; i < N2; ++i)
	{
		//ЕСЛИ СТРОКА МАТРИЦЫ A ПОДЛЕЖИТ УДАЛЕНИЮ ТО ПРОПУСКАЕМ ЕЕ 
		if(i == to_del[counti])
		{
			counti += 1;
			continue;
		}

		//ОБНУЛЕНИЕ СЧЕТЧИКА
		countj = 0;

		//ЕСЛИ СТРОКА НЕ ДЛЯ УДАЛЕНИЯ ЗАПУСКАЕТСЯ ЦИКЛ ПО СТОЛБЦАМ МАТРИЦЫ A
		for(int j = 0; j < N2; ++j)
		{

			//ПРОВЕРКА НА ТО ПОДЛЕДЖИТ ЛИ СТОЛБЕЦ МАТРИЦЫ A УДАЛЕНИЮ
			if(j == to_del[countj])
			{
				countj += 1;
				continue;
			}

			//ЗАПОЛНЕНИЕ МАТРИЦЫ Asol
			//НОМЕРА СТРОК И СТОЛБЦОВ ДЛЯ МАТРИЦЫ Asol ЕСТЬ (i - counti) И (j - countj)
			Asol[(i - counti)*N3 + (j - countj)] = A[i*N2 + j];
		}
	}	
}

void Solution::fill_Bsol()
{
	//РАЗМЕРНОСТЬ ВЕКТОРА Bsol N3
	//НЕОБХОДИМО УДАЛИТЬ ИЗ Bsol 4 КОМПОНЕНТЫ

	//МАССИВ ИНДЕКСОВ СТОЛБЦОВ И СТРОК ДЛЯ УДАЛЕНИЯ
	int to_del[4] = {2*Ridx[0], 2*Ridx[1], 2*Ridx[2], 2*Ridx[3]};

	//СЧЕТЧИК УДАЛЕННЫХ КОМПОНЕНТ
	int count = 0;

	//ЦИКЛ ПО КОМПОНЕНТАМ ВЕКТОРА B
	for(int i = 0; i < N2; ++i)
	{
		//ЕСЛИ КОМПОНЕНТА ПОДЛЕЖИТ УДАЛЕНИЮ ПРОПУСКАЕМ ЕЕ
		if(i == to_del[count])
		{
			count += 1;
			continue;
		}

		//НОМЕР КОМПОНЕНТЫ Bsol (i - count)
		Bsol[i-count] = B[i];
	}

}

void Solution::solve()
{
	//РЕШАЕТСЯ СИСТЕМА УРАВНЕНИЙ МЕТОДОМ СОПРЯЖЕННЫХ ГРАДИЕНТОВ (ЗАПОЛНЕНИЕ ВЕКТОРА Csol)
	//РАЗМЕРНОСТЬ СИСТЕМЫ N3
	//ЗАПОЛНЯЕТСЯ ВЕКТОР C

	//РЕШЕНИЕ СИСТЕМЫ МЕТОДОМ СОПРЯЖЕННЫХ ГРАДИЕНТОВ
	//ПЕРЕМЕННЫЕ НЕОБХОДИМЫЕ ДЛЯ РЕАЛИЗАЦИИ МЕТОДА
	double alpha;
	double beta;
	double *r = new double[N3];
	double *p = new double[N3];
	double *q = new double[N3];

	//ПОДГОТОВКА ЗНАЧЕНИЙ 
	for(int i = 0; i < N3; ++i)
	{
		Csol[i] = 0.0;
		r[i] = Bsol[i];
		p[i] = Bsol[i];
	}

	double rtmp = 0.0;
	for (size_t i = 0; i < N3; ++i) rtmp += r[i] * r[i];

	//ОСНОВНОЙ ЦИКЛ РЕШЕНИЯ
	double ptmp;
	for(size_t  k = 0; k < N3; ++k)
	{
		ptmp = 0.0;
		for(size_t i = 0; i < N3; ++i)
		{
			q[i] = 0.0;
			for(size_t j = 0; j < N3; ++j) {q[i] += Asol[i * N + j] * p[j];}
			ptmp += p[i] * q[i];
		}
		alpha = rtmp / ptmp;

		ptmp = 0.0;
		for(size_t i = 0; i < N3; ++i)
		{
			Csol[i] += alpha * p[i];
			r[i] -= alpha * q[i];
			ptmp += r[i] * r[i];
		}
		beta = ptmp / rtmp;
		rtmp = ptmp;
		for(size_t i = 0; i < N3; ++i) { p[i] = r[i] + beta * p[i];}
	}
	//СИСТЕМА РЕШЕНА
	//ЧИСТКА ПАМЯТИ ПОСЛЕ РЕШЕНИЯ СИСТЕМЫ
	delete r,p,q;

	//НАЧИНАЕТСЯ ЗАПОЛНЕНИЕ ВЕКТОРА С

	//ПРОПУЩЕННЫЕ ИНДЕКСЫ
	int miss[4] = {2*Ridx[0], 2*Ridx[1], 2*Ridx[2], 2*Ridx[3]};

	//СЧЕТЧИК ВЫСТАВЛЕННЫХ ЗНАЧЕНИЙ
	int count = 0;
	for(int i = 0; i < N2; ++i)
	{
		//В СЛУЧАЕ ПРОПУЩЕННОГО ЗНАЧЕНИЯ
		if(i == miss[count])
		{
			count += 1;
			C[i] = 0.0;
			continue;
		}

		//В СЛУЧАЕ НУЖНОГО ЗНАЧЕНИЯ
		C[i] = Csol[i - count];
	}
	//ВЕКТОР C ЗАПОЛНЕН
}

void Solution::set_zeros_A()
{
	for(int i = 0; i < N2*N2; ++i) A[i] = 0.0;
}

void Solution::set_zeros_a()
{
	for(int i = 0; i < 16; ++i) a[i] = 0.0;
}

void Solution::set_zeros_B()
{
	for(int i = 0; i < N2; ++i) B[i] = 0.0;
}

void Solution::set_zeros_b()
{
	for(int i = 0; i < 4; ++i) b[i] = 0.0;
}

double Solution::opa(int k, int i, int j)
{
	//ВЫЧИСЛЕНИЕ БИЛИНЕЙНОГО ФУНКЦИОНАЛА
	//k ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА
	//(i,j) ИНДЕКСЫ БАЗИСНЫХ ФУНКЦИЙ

	//ДОМНОЖАЕМ ПРОИЗВЕДЕНИЕ ФУНКЦИЙ НА dot
	double dot = E * J / pow(x[k+1] - x[k], 3);

	//ПРИБАВЛЯЕМ СЛАГАЕМОЕ ОТ ПРУЖИНЫ ЕСЛИ РАБОТАЕМ НА КОНЕЧНОМ ЭЛЕМЕНТЕ С ПРУЖИНОЙ
	if(k == kk)
	{
		//ТОЧКА В КОТОРОЙ НЕОБХОДИМО ВЫЧИСЛИТЬ ЗНАЧЕНИЕ ФУНКЦИЙ
		double y = (xk - x[k])/(x[k+1] - x[k]);

		//ДОБАВЛЯЕМ К ПРОИЗВЕДЕНИЮ ЧЛЕН ОТВЕЧАЮЩИЙ ЗА ПРУЖИНУ
		return dot * c[i*4 + j] + k*BASIS[i](y)*BASIS[j](y);
	}

	//ИНАЧЕ ТОЛЬКО ПРОИЗВЕДЕНИЕ
	return dot * c[i*4 + j];

}

double Solution::opb(int k, int i)
{
	//СЧИТАЕТ ИНТЕГРАЛ ОТ РАСПРЕДЛЕННОЙ НАГРУЗКИ

	//МАССИВЫ ИНТЕГРАЛОВ ОТ БАЗИСНЫХ ФУНКЦИЙ И Y*БАЗИСНЫЕ ФУНКЦИИ НА [0,1]
	double integral[4] = {0.5, 0.08333333, 0.5, -0.08333333};
	double integral_y[4] = {0.15, 0.03333333, 0.35, -0.05};

	//КОНСТАНТЫ РАСПРЕДЕЛЕННОЙ НАГРУЗКИ ДЛЯ УПРОЩЕНИЯ ЗАПИСИ ВЫРАЖЕНИЯ q(x) = k1 + k2*x
	double k2 = (qB - qA)/(xq2 - xq1);
	double k1 = qA - xq1*k2;

	//ШАГ
	double h = x[k+1] - x[k];

	//ИНТЕГРАЛ ВЫРАЖЕННЫЙ ЧЕРЕЗ ЛОКАЛЬНЫЕ КООРДИНАТЫ
	return h*(k1+k2*x[k])*integral[i] + pow(h,2)*k2*integral_y[i];
}