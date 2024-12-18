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
double N1(double y, double h){return 2.0*pow(y,3) - 3.0*pow(y,2) + 1.0;}
double L1(double y, double h){return h*(pow(y,3) - 2.0*pow(y,2) + y);}
double N2(double y, double h){return -2.0*pow(y,3) + 3.0*pow(y,2);}
double L2(double y, double h){return h*(pow(y,3) - pow(y,2));}

//ПЕРВЫЕ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double dN1(double y, double h){return 6.0*pow(y,2) - 6.0*y;}
double dL1(double y, double h){return h*(3.0*pow(y,2) - 4.0*y + 1.0);}
double dN2(double y, double h){return -6.0*pow(y,2) + 6.0*y;}
double dL2(double y, double h){return h*(3.0*pow(y,2) - 2.0*y);}

//ВТОРЫЕ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double d2N1(double y, double h){return 12.0*y - 6.0;}
double d2L1(double y, double h){return h*(6.0*y - 4.0);}
double d2N2(double y, double h){return -12.0*y + 6.0;}
double d2L2(double y, double h){return h*(6.0*y - 2.0);}

//ТРЕТЬИ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double d3N1(double y, double h){return 12.0;}
double d3L1(double y, double h){return h*6.0;}
double d3N2(double y, double h){return -12.0;}
double d3L2(double y, double h){return h*6.0;}

//МАССИВ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*BASIS[4])(double, double) = {N1, L1, N2, L2};

//МАССИВ ПЕРВЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*dBASIS[4])(double, double) = {dN1, dL1, dN2, dL2};

//МАССИВ ВТОРЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*d2BASIS[4])(double, double) = {d2N1, d2L1, d2N2, d2L2};

//МАССИВ ТРЕТЬИХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
double (*d3BASIS[4])(double, double) = {d3N1, d3L1, d3N2, d3L2};



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
    this->Rnv = new double[N2];

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

    //ВЫВОД НОРМЫ НЕВЯЗКИ
    info_nv();

    //ВЫЧИСЛЕНИЕ ВЕКТОРА НЕВЯЗКИ МАТРИЦЫ ЖЕСТКОСТИ
    fill_Rnv();

    //ВЫВОД ДАННЫХ О ВЕКТОРЕ НЕВЯЗКИ МАТРИЦЫ ЖЕСТКОСТИ
    info_Rnv();

    //ЗАПОЛНЕНИЕ TXT ФАЙЛА ДЛЯ ПОСТРОЕНИЯ ГРАФИКИ
    fill_txt();

    //ЗАПОЛНЕНИЕ TXT ФАЙЛА ДЛЯ ПОСТРОЕНИЯ ТОЧЕК ПРОВЕРКИ ПРОГИБА И УГЛОВ
    check_txt();

    //ЗАПОЛНЕНИЕ TXT ФАЙЛА УСИЛИЯМИ R
    Rtxt();

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
	logfile << "ЗАДАННАЯ СИЛА | P = " << CONST::P << endl;
	logfile << "ЗАДАННЫЙ МОМЕНТ | M = " << CONST::M << endl;
	logfile << "ЖЕСТКОСТЬ ПРУЖИНЫ | k = " << CONST::k << endl;
	logfile << "РАСПРЕДЕЛЕННАЯ В ТОЧКЕ A | qA = " << CONST::qA << endl;
	logfile << "РАСПРЕДЕЛЕННАЯ В ТОЧКЕ B | qB = " << CONST::qA << endl;
	logfile << "МОДУЛЬ ЮНГА | E = " << CONST::E << endl;
	logfile << "МОМЕНТ ИНЕРЦИИ | J = " << CONST::J << endl;
	logfile << "МИНИМАЛЬНО ДОПУСТИМОЕ ЧИСЛО УЗЛОВ | N_min = 10" << endl;
	logfile << "ЧИСЛО УЗЛОВ | N = " << N << endl;
	logfile << "ЧИСЛО КОНЕЧНЫХ ЭЛЕМЕНТОВ | M = " << M << endl;
}

void Solution::info_mesh()
{
	logfile << "МИНИМАЛЬНЫЙ РАЗМЕР ЭЛЕМЕНТА В СЕТКЕ | min_h = " << min_h() << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С ПРУЖИНОЙ | kk = " << kk + 1 << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С МОМЕНТОМ | km = " << km + 1 << endl;
	logfile << "НОМЕР КОНЕЧНОГО ЭЛЕМЕНТА С СИЛОЙ | kp = " << kp + 1 << endl;
	logfile << "НОМЕР УЗЛА ПЕРВОЙ ОПОРЫ | No1 = " << Ridx[0] + 1 << endl;
	logfile << "НОМЕР УЗЛА ВТОРОЙ ОПОРЫ | No2 = " << Ridx[1] + 1 << endl;
	logfile << "НОМЕР УЗЛА ТРЕТЬЕЙ ОПОРЫ | No3 = " << Ridx[2] + 1 << endl;
	logfile << "НОМЕР УЗЛА ЧЕТВЕРТОЙ ОПОРЫ | No4 = " << Ridx[3] + 1 << endl;
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

void Solution::info_nv()
{
	logfile << "НОРМА L1 НЕВЯЗКИ УРЕЗАННОЙ МАТРИЦЫ ЖЕСТКОСТИ Asol::" << endl;
	logfile << nv() << endl;
}

void Solution::info_Rnv()
{
	logfile << "ВЕКТОР НЕВЯЗКИ МАТРИЦЫ ЖЕТСКОСТИ Rnv::" << endl;
	for(int i = 0; i < N2; ++i) logfile << Rnv[i] << endl;
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

double Solution::nv()
{
	//ВЫЧИСЛЯЕТ НОРМУ L1 НЕВЯЗКИ Asol*Csol - Bsol
	double nv = 0.0;

	//i-Я КОМПОНЕНТА
	double nvi;
	for(int i = 0; i < N3; ++i)
	{
		nvi = 0.0;
		for(int j = 0; j < N3; ++j) nvi += Asol[i*N3 + j]*Csol[j];
		nv += abs(nvi - Bsol[i]);
	}
	return nv;
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

void Solution::sort(int idx)
{
	//ИДЕТ СОРТИРОВКА ДО ЭЛЕМЕНТА С НОМЕРОМ idx

	//ПУЗЫРЬКОВАЯ СОРТИРОВКА
	double tmp;
   	for(int i = 0; i < idx - 1; ++i)
   	{
   		for(int j = 0; j < idx - 1 - i; ++j)
   		{
   			if(x[j] > x[j+1])
   			{
   				tmp = x[j];
   				x[j] = x[j+1];
   				x[j+1] = tmp;
   			}
   		}
   	}
}

void Solution::create_mesh()
{	
	//МИНИМАЛЬНОЕ ЧИСЛО УЗЛОВ 10 (2 УЗЛА РАСПРЕДЕЛЕННОЙ НАГРУЗКИ, 4 ОПОРЫ, 1 УЗЕЛ НА ЛЕВОМ ТОРЦЕ, 1 ПРУЖИНА, 1 СИЛА, 1 МОМЕНТ)
	//СТРОИТСЯ РАВНОМЕРНАЯ СЕТКА, ВКЛЮЧАЮЩАЯ В СЕБЯ N-10 УЗЛОВ БЕЗ КРАЕВЫХ
	//УЗЕЛ НА ПРАВОМ ТОРЦЕ ЭТО УЗЕЛ ОДНОЙ ИЗ ОПОР

    x[0] = 0.0;
    x[1] = xk;
   	x[2] = xP;
    x[3] = xM;
    x[4] = xR1;
   	x[5] = xR2;
    x[6] = xR3;
	x[7] = xR4;
   	x[8] = xq1;
   	x[9] = xq2;

   	//КАЖДУЮ СЛЕДУЮЩУЮ ТОЧКУ КЛАДЕМ В СЕРЕДИНУ САМОГО БОЛЬШОГО КОНЕЧНОГО ЭЛЕМЕНТА
   	for(int i = 0; i < N - 10; ++i)
   	{
   		//ИНДЕКС КУДА НЕОБХОДИМО ПОЛОЖИТЬ ТОЧКУ
   		int idx = 10 + i;

   		//СОРТРОВКА УЖЕ ПОСТАВЛЕННЫХ УЗЛОВ
   		sort(idx);

   		//ДЛИНА САМОГО БОЛЬШОГО ОТРЕЗКА
   		double h_max = 0;
   		//ИНДЕКС ЛЕВОГО УЗЛА САМОГО БОЛЬШОГО ОТРЕЗКА
   		int idx_max = 0; 
   		//ПОИСК САМОГО БОЛЬШОГО ОТРЕЗКА
   		for(int j = 0; j < idx - 1; ++j)
   		{
   			if(x[j+1] - x[j] > h_max)
   			{
   				h_max = x[j+1] - x[j];
   				idx_max = j;
   			}
   		}

   		//КЛАДЕМ НОВУЮ ТОЧКУ
   		x[idx] = x[idx_max] + 0.5*h_max;
   	}

   	//СОРТИРУЕМ ВСЮ ПОСТРОЕННУЮ СЕТКУ
   	sort(N);

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

void Solution::fill_c(double h)
{
	//МАТРИЦА с ИСПОЛЬЗУЕТСЯ В МЕТОДЕ op ДЛЯ ВЫЧИСЛЕНИЯ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ
	//ДАННЫЙ МЕТОД ЗАПОЛНЯЕТ МАТРИЦУ с ПОЛУЧЕННЫМИ ЗНАЧЕНИЯМИ ИЗ SAGE ПРОГРАММЫ

	c[0]  =  12.0;   c[1]  =  h*6.0;        c[2]  =  -12.0;   c[3]  = h*6.0;
	c[4]  =   h*6.0; c[5]  =  pow(h,2)*4.0; c[6]  =  -h*6.0;  c[7]  = pow(h,2)*2.0;
	c[8]  = -12.0;   c[9]  = -h*6.0;        c[10] =  12.0;    c[11] = -h*6.0;
	c[12] =   h*6.0; c[13] =  pow(h,2)*2.0; c[14] =  -h*6.0;  c[15] = pow(h,2)*4.0;
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
			for(int i = 0; i < 4; ++i) b[i] -= (1.0/(x[k+1] - x[k]))*CONST::M*dBASIS[i]((xM - x[k])/(x[k+1] - x[k]), (x[k+1] - x[k]));
		}

		//ДОБАВЛЯЕМ СИЛУ ЕСЛИ ОНА ЕСТЬ НА КОНЕЧНОМ ЭЛЕМЕНТЕ
		if(k == kp)
		{
			for(int i = 0; i < 4; ++i) b[i] += CONST::P*BASIS[i]((xP - x[k])/(x[k+1] - x[k]), (x[k+1] - x[k]));
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
	//НЕОБХОДИМО УДАЛИТЬ ИЗ B 4 КОМПОНЕНТЫ

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
	for(size_t  k = 0; k < N3 + 30; ++k)
	{
		ptmp = 0.0;
		for(size_t i = 0; i < N3; ++i)
		{
			q[i] = 0.0;
			for(size_t j = 0; j < N3; ++j) q[i] += Asol[i * N3 + j] * p[j];
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
	double dot = CONST::E * CONST::J / pow(x[k+1] - x[k], 3);

	//ЗАПОЛНЕНИЕ МАТРИЦЫ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ
	fill_c(x[k+1] - x[k]);

	//ПРИБАВЛЯЕМ СЛАГАЕМОЕ ОТ ПРУЖИНЫ ЕСЛИ РАБОТАЕМ НА КОНЕЧНОМ ЭЛЕМЕНТЕ С ПРУЖИНОЙ
	if(k == kk)
	{
		//ТОЧКА В КОТОРОЙ НЕОБХОДИМО ВЫЧИСЛИТЬ ЗНАЧЕНИЕ ФУНКЦИЙ
		double y = (xk - x[k])/(x[k+1] - x[k]);

		//ДОБАВЛЯЕМ К ПРОИЗВЕДЕНИЮ ЧЛЕН ОТВЕЧАЮЩИЙ ЗА ПРУЖИНУ
		return dot * c[i*4 + j] + CONST::k*BASIS[i](y,x[k+1] - x[k])*BASIS[j](y,x[k+1] - x[k]);
	}

	//ИНАЧЕ ТОЛЬКО ПРОИЗВЕДЕНИЕ
	return dot * c[i*4 + j];

}

double Solution::opb(int k, int i)
{
	//СЧИТАЕТ ИНТЕГРАЛ ОТ РАСПРЕДЛЕННОЙ НАГРУЗКИ

	//ШАГ
	double h = x[k+1] - x[k];

	//МАССИВЫ ИНТЕГРАЛОВ ОТ БАЗИСНЫХ ФУНКЦИЙ И Y*БАЗИСНЫЕ ФУНКЦИИ НА [0,1]
	double integral[4] = {0.5, h*0.08333333, 0.5, -h*0.08333333};
	double integral_y[4] = {0.15, h*0.03333333, 0.35, -h*0.05};

	//КОНСТАНТЫ РАСПРЕДЕЛЕННОЙ НАГРУЗКИ ДЛЯ УПРОЩЕНИЯ ЗАПИСИ ВЫРАЖЕНИЯ q(x) = k1 + k2*x
	double k2 = (CONST::qB - CONST::qA)/(xq2 - xq1);
	double k1 = CONST::qA - xq1*k2;

	//ИНТЕГРАЛ ВЫРАЖЕННЫЙ ЧЕРЕЗ ЛОКАЛЬНЫЕ КООРДИНАТЫ
	return h*(k1+k2*x[k])*integral[i] + pow(h,2)*k2*integral_y[i];
}

int Solution::find_k(double z)
{
	//ВЫЧИСЛЕНИЕ ИНДЕКСА КОНЕЧНОГО ЭЛЕМЕНТА КОТОРОМУ ПРЕНАДЛЕЖИТ ТОЧКА z
	//ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА
	int kz;

	//ЕСЛИ z ПОПАЛ НА УЗЕЛ ТО БУДЕМ ПОЛАГАТЬ ЧТО z ПРЕНАДЛЕЖИТ  ЭЛЕМЕНТУ ЛЕВЫЙ УЗЕЛ КОТОРОГО ДАННЫЙ
	//ЦИКЛ ПО КОНЕЧНЫМ ЭЛЕМЕНТАМ
	for(int k = 0; k < M; ++k)
	{
		if((z > x[k] and z < x[k+1]) or abs(z - x[k]) < 1.0e-16) kz = k;
	}

	//ЕСЛИ ТОЧКА y = L - ДЛИНЕ БАЛКИ ТО ИНДЕКС ПОСЛЕДНЕГО КОНЕЧНОГО ЭЛЕМЕНТА
	//ТУТ УЗЕЛ КОНЕЧНОГО ЭЛЕМЕНТА ПРАВЫЙ
	if(abs(z - x[N - 1]) < 1.0e-16) kz = M - 1;

	return kz;
}

double Solution::w(double z)
{
	//НЕОБХОДИМО ВЫЧИСЛИТЬ ПРОГИБ В ТОЧКЕ z

	//ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ТОЧКОЙ z
	int kz = find_k(z);

	//ТОЧКА В ЛОКАЛЬНЫХ КООРДИНАТАХ СООТВЕТСТВУЮЩАЯ ТОЧКЕ z В ГЛОБАЛЬНЫХ
	double y = (z - x[kz]) / (x[kz + 1] - x[kz]);

	//ЗНАЧЕНИЕ ПРОГИБА
	double wz = 0.0;
	for(int i = 0; i < 4; ++i) wz += C[2*kz + i]*BASIS[i](y, x[kz+1] - x[kz]);
	//ВЫЧИСЛЕНО ЗНАЧЕНИЕ ПРОГИБА

	return wz;
}

double Solution::dw(double z)
{
	//НЕОБХОДИМО ВЫЧИСЛИТЬ ПРОИЗВОДНУЮ ПРОГИБА В ТОЧКЕ z

	//ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ТОЧКОЙ z
	int kz = find_k(z);

	//ТОЧКА В ЛОКАЛЬНЫХ КООРДИНАТАХ СООТВЕТСТВУЮЩАЯ ТОЧКЕ z В ГЛОБАЛЬНЫХ
	double y = (z - x[kz]) / (x[kz + 1] - x[kz]);

	//ЗНАЧЕНИЕ ПРОИЗВОДНОЙ ПРОГИБА
	double dwz = 0.0;
	for(int i = 0; i < 4; ++i) dwz += C[2*kz + i]*(1.0/(x[kz+1] - x[kz]))*dBASIS[i](y, x[kz+1] - x[kz]);
	//ВЫЧИСЛЕНО ЗНАЧЕНИЕ ПРОГИБА

	return dwz;
}

double Solution::d2w(double z)
{
	//НЕОБХОДИМО ВЫЧИСЛИТЬ ВТОРУЮ ПРОИЗВОДНУЮ ПРОГИБА В ТОЧКЕ z

	//ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ТОЧКОЙ z
	int kz = find_k(z);

	//ТОЧКА В ЛОКАЛЬНЫХ КООРДИНАТАХ СООТВЕТСТВУЮЩАЯ ТОЧКЕ z В ГЛОБАЛЬНЫХ
	double y = (z - x[kz]) / (x[kz + 1] - x[kz]);

	//ЗНАЧЕНИЕ ВТОРОЙ ПРОИЗВОДНОЙ ПРОГИБА
	double d2wz = 0.0;
	for(int i = 0; i < 4; ++i) d2wz += C[2*kz + i]*pow((1.0/(x[kz+1] - x[kz])),2)*d2BASIS[i](y, x[kz+1] - x[kz]);
	//ВЫЧИСЛЕНО ЗНАЧЕНИЕ ПРОГИБА

	return d2wz;
}

double Solution::d3w(double z)
{
	//НЕОБХОДИМО ВЫЧИСЛИТЬ ТРЕТЬЮ ПРОИЗВОДНУЮ ПРОГИБА В ТОЧКЕ z

	//ИНДЕКС КОНЕЧНОГО ЭЛЕМЕНТА С ТОЧКОЙ z
	int kz = find_k(z);

	//ТОЧКА В ЛОКАЛЬНЫХ КООРДИНАТАХ СООТВЕТСТВУЮЩАЯ ТОЧКЕ z В ГЛОБАЛЬНЫХ
	double y = (z - x[kz]) / (x[kz + 1] - x[kz]);

	//ЗНАЧЕНИЕ ВТОРОЙ ПРОИЗВОДНОЙ ПРОГИБА
	double d3wz = 0.0;
	for(int i = 0; i < 4; ++i) d3wz += C[2*kz + i]*pow((1.0/(x[kz+1] - x[kz])),3)*d3BASIS[i](y, x[kz+1] - x[kz]);
	//ВЫЧИСЛЕНО ЗНАЧЕНИЕ ПРОГИБА

	return d3wz;
}


void Solution::fill_txt()
{
	//НЕОБХОДИМО ЗАПОЛНИТЬ graph_val.txt ФАЙЛ ЗНАЧЕНИЯМИ
	ofstream val;

	//ФАЙЛ ДЛЯ СЧИТЫВАНИЯ ТОЧЕК В КОТОРЫХ НЕОБХОДИМО ПРОИЗВЕСТИ ВЫЧИСЛЕНИЕ graph_arg.txt
	ifstream argw;

	//ОТКРЫТИЕ graph_arg.txt ФАЙЛА ДЛЯ ЧТЕНИЯ
    argw.open("graph_arg.txt");

    //КОЛИЧЕСТВО ТОЧЕК
    int Nz;
    argw >> Nz;

    //МАССИВ ТОЧЕК
    double *z = new double[Nz];
    for(int i = 0; i < Nz; ++i) argw >> z[i];

    //ЗАКРЫТИЕ graph_arg.txt ФАЙЛА
    argw.close();

	//ОТКРЫТИЕ ФАЙЛА graph_val.txt ДЛЯ ЗАПИСИ
	val.open("graph_val.txt");
	
	for(int i = 0; i < Nz; ++i)
	{
		val << w(z[i]) << " " << dw(z[i]) << " " 
		<< d2w(z[i]) << " " << d3w(z[i]) << endl;
	}

	//ЗАКРЫТИЕ graph_val.txt ФАЙЛА
	val.close(); 

   	//ЧИСТКА ПАМЯТИ 
   	delete z;
}

void Solution::check_txt()
{
	//ЗАПОЛНЯЕТ ФАЙЛ check.txt ДЛЯ ПОСТРОЕНИЯ ТОЧЕК НА ГРАФИКЕ
	//ЗНАЧЕНИЯ ПРОГИБА И УГЛОВ В УЗЛАХ

	//ФАЙЛ ДЛЯ ЗАПОЛНЕНИЯ
	ofstream check;

	//ОТКРЫТИЕ ФАЙЛА
	check.open("check.txt");

	for(int i = 0; i < N; ++i)
	{
		check << x[i] << " " << C[2*i] << " " << C[2*i + 1] << endl;
	}

	//ЗАКРЫТИЕ ФАЙЛА
	check.close();
}

void Solution::fill_Rnv()
{
	//НЕОБХОДИМО ЗАПОЛНИТЬ ВЕКТОР Rnv НЕВЯЗКИ ПОЛНОЙ СИСТЕМЫ

	//i-АЯ КОМПОНЕНТА A*C
	double aci;

	//СЧИТАЕМ ПРОИЗВЕДЕНИЕ МАТРИЦЫ НА ВЕКТОР С И ВЫЧИТАЕМ B
	for(int i = 0; i < N2; ++i)
	{
		aci = 0.0;
		for(int j = 0; j < N2; ++j) aci += A[i*N2 + j] * C[j];

		//i-АЯ КОМПОНЕНТА НЕВЯЗКИ
		Rnv[i] = aci - B[i];
	}

}

void Solution::Rtxt()
{
	//ЗАПОЛНЕНИЕ ФАЙЛА "R.txt" ЗНАЧЕНИЯМИ УСИЛИЙ

	//ФАЙЛ ДЛЯ ЗАПОЛНЕНИЯ
	ofstream file;

	//ОТКРЫТИЕ ФАЙЛА
	file.open("R.txt", ios::app);

	file << "C++ R:" << endl;

	for(int i = 0; i < N2; ++i)
	{
		if(abs(Rnv[i]) > 1.0e-4) file << Rnv[i] << endl;
	}

	//ЗАКРЫТИЕ ФАЙЛА
	file.close();

}