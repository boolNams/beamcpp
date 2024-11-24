#МОДУЛЬ ДЛЯ ВЫЧИСЛЕНИЯ ЯВНЫХ ВЫРАЖЕНИЙ ОПЕРАТОРА НА БАЗИСНЫХ ФУНКЦИЯХ В ЛОКАЛЬНЫХ КООРДИНАТАХ
#ВЫЧИСЛЕНИЯ СОХРАНЯЮТСЯ В ФАЙЛ sageinfo.txt

import numpy as np

info = open(r"sageinfo.txt", mode='wt') 
print("### SAGE ВЫЧИСЛЕНИЯ ###", file = info)

var('y')

#БАЗИСНЫЕ ФУНКЦИИ НА ОТРЕЗКЕ [0,1]
N1(y) = 1 - 3*y^2 + 2*y^3 
L1(y) = y - 2*y^2 + y^3
N2(y) = 3*y^2 - 2*y^3
L2(y) = -y^2 + y^3 

#СПИСОК ФУНКЦИЙ
F = [N1, L1, N2, L2]

print("БАЗИСНЫЕ ФУНКЦИИ::", file = info)
for item in F:
	print(item, file = info)

#СПИСОК ИХ ПЕРВЫХ ПРОИЗВОДНЫХ
F1 = []
for item in F:
	F1.append(diff(item, y))

print("ПЕРВЫЕ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ::", file = info)
for item in F1:
	print(item, file = info)

#СПИСОК ИХ ВТОРЫХ ПРОИЗВОДНЫХ
F2 = []
for item in F1:
	F2.append(diff(item,y))

print("ВТОРЫЕ ПРОИЗВОДНЫЕ БАЗИСНЫХ ФУНКЦИЙ::", file = info)
for item in F2:
	print(item, file = info)

#МАТРИЦА ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ НА [0,1]
#ЭТО ЕСТЬ ИНТЕГРАЛ ОТ 0 ДО 1 ОТ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ БАЗИСНЫХ ФУНКЦИЙ
A = np.zeros([4,4])

print("МАТРИЦА ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ НА [0,1]::", file = info)

#ВЫЧИСЛЕНИЕ МАТРИЦЫ
for i in range(0,4):
	for j in range(0,4):
		#СЧИТАЕТСЯ ИНТЕГРАЛ ПРОИЗВЕДЕНИЯ ВТОРЫХ ПРОИЗВОДНЫХ
		A[i,j] = integral(F2[i]*F2[j],(y,0,1))

#ВЫВОД МАТРИЦЫ
print(A, file = info)