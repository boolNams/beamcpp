#
#ПРОГРАММА РАСЧЕТА БАЛКИ С ПОМОЩЬЮ СИСТЕМЫ КОМПЬЮТЕРНОЙ АЛГЕБРЫ | SAGE
#ДЛИНА БАЛКИ | 25
#
#ВЫЗЫВАЕТСЯ ПРОГРАММА C++
#СТРОЯТСЯ ГРАФИКИ МКЭ

import os

#ПОДГРУЗКА ЗАДАННЫХ ЗНАЧЕНИЙ ИЗ ФАЙЛА
data = []
params = open("params.txt", "r")
for line in params:
    str_value = line.split('|')[0].strip()
    data.append(float(str_value))

#================================================== SAGE 
#================================================== ПРОГИБ
M, P, k, qA, qB, E, J = var('M P k qA qB E J')							# переменные заданные по условию
P = data[0]; M = data[1]; k = data[2]
qA = data[3]; qB = data[4]
E = data[5]; J = data[6]
w0, wk, teta0, RA, RB, RC, RD = var('w0 wk teta0 R1 R2 R3 R4')			# переменные подлежащие определению

x = var('x')															# координата (0,25)     
prd = E*J
div = 1/(E*J)

a = var('a')															# параметр многочленов 
y2(x, a) = (1/2) * (x - a)^2
y3(x, a) = (1/6) * (x - a)^3
y4(x, a) = (1/24) * (x - a)^4
y5(x, a) = (1/120) * (x - a)^5                                                   

AR, BR, CR, DR = 1, 8, 16, 25											# координаты упругих опор
hAR(x,RA) = RA*heaviside(x - AR)										# функции хевисайда для реакций
hBR(x,RB) = RB*heaviside(x - BR)
hCR(x,RC) = RC*heaviside(x - CR)
hDR(x,RD) = RD*heaviside(x - DR)
wR = y3(x,AR)*hAR + y3(x,BR)*hBR + y3(x,CR)*hCR + y3(x,DR)*hDR			# прогиб (опоры)

AM = 13																	# координата приложенного момента M
hAM(x) = M*heaviside(x - AM)											# функция хевисайда для момента M
wM = y2(x,AM)*hAM														# прогиб (момент)

AP = 20																	# координата приложенной силы P
hAP(x) = P*heaviside(x - AP)											# функция хевисайда для силы P
wP = y3(x,AP)*hAP             											# прогиб (сила)

Aq, Bq = 10, 22 														# координаты распределенной нагрузки
tg = (qB - qA) / (Bq - Aq)												# тангенс угла наклона
hAq(x) = heaviside(x - Aq)												# функция хевисайда основной нагрузки (+)
hBq(x) = heaviside(x - Bq)												# функция хевисайда компенсирующей нагрузки (-)
wAq = (qA*y4(x,Aq) + tg*y5(x,Aq))*hAq									# прогиб (основная нагрузка +) 
wBq = -(qB*y4(x,Bq) + tg*y5(x,Bq))*hBq									# прогиб (компенсирующая нагрузка -)

Ak = 23																	# координата пружины жетскости k
Rk = -k*wk                                                              # сила от пружины
hAk(x, wk) = Rk*heaviside(x - Ak)									    # функция хевисайда пружины k
wAk = y3(x,Ak)*hAk  													# прогиб (пружина)

wConst(x, w0, teta0) = w0 + teta0*x										# прогиб (константы интегрирвания)

#==================================================
w = wConst + div * (wR + wM + wP + wAq + wBq + wAk)						# прогиб
#================================================== СИСТЕМА

eq1 = RA+RB+RC+RD + P + (1/2)*(Bq-Aq)*(qA+qB) + Rk == 0   			    # уравнение равновесия (1)

MR = RA*AR+RB*BR+RC*CR+RD*DR                                      		# момент создаваемый опорами
MQ = (1/2)*qA*(Bq^2-Aq^2) + tg*((1/3)*(Bq^3 - Aq^3) - (Aq/2)*(Bq^2 - Aq^2)) # момент создаваемый распределенной нагрузкой
eq2 = MR + P*AP - M + Rk*Ak + MQ == 0									# уравнение моментов (2)

eq3 = w(x = AR) == 0              										# уравнения на опорах (3-6)
eq4 = w(x = BR) == 0
eq5 = w(x = CR) == 0
eq6 = w(x = DR) == 0

eq7 = w(x = Ak) - wk == 0      											# уравнение на пружину (7)

eq = [eq1,eq2,eq3,eq4,eq5,eq6,eq7] 										# система уравнений
sol = solve(eq, w0, wk, teta0, RA, RB, RC, RD, solution_dict = True)[0]
w0, wk, teta0 = sol[w0], sol[wk], sol[teta0]
RA, RB, RC, RD = sol[RA], sol[RB], sol[RC], sol[RD]

dots = [AR, BR, CR, DR, AM, AP, Aq, Bq, Ak]

w = w(RA=RA, RB=RB, RC=RC, RD=RD, teta0=teta0, w0=w0, wk=wk)			# функция прогиба от x

teta = w.diff(x)                                                        # функция угла от x
for dot in dots:
    teta = teta.subs(dirac_delta(x - dot) == 0)

MM = teta.diff(x)                                                       # функция момента от x
for dot in dots:
    MM = MM.subs(dirac_delta(x - dot) == 0)

Q = MM.diff(x)
for dot in dots:
    Q = Q.subs(dirac_delta(x - dot) == 0)

#================================================== END SAGE
#================================================== C++
import numpy as np

#КОЛИЧЕСТВО ТОЧЕК НА ГРАФИКЕ
N = 500

#ФАЙЛ ДЛЯ ЗАПИСИ УСИЛИЙ ВЫЧИСЛЕННЫХ SAGE И С++
with open("R.txt",'w') as txt:
    txt.write("SAGE R:\n")
    txt.write(f"{RA.numerical_approx()}\n")
    txt.write(f"{RB.numerical_approx()}\n")
    txt.write(f"{RC.numerical_approx()}\n") 
    txt.write(f"{RD.numerical_approx()}\n")

#НЕОБХОДИМО ВЫЧИСЛИТЬ ЗНАЧЕНИЯ В ДАННЫХ ТОЧКАХ
arg = np.linspace(0,25.0,N)

#graph_arg.txt ФАЙЛ ДЛЯ ПЕРЕДАЧИ ЗНАЧЕНИЙ В ПРОГРАММУ C++
with open("graph_arg.txt",'w') as txt:
    txt.write(f"{arg.size}\n")
    np.savetxt(txt, arg)

#ВЫЗОВ ПРОГРАММЫ C++
os.system("bash run.sh")

#graph_w ФАЙЛ ДЛЯ ПЕРЕДАЧИ ЗНАЧЕНИЙ В ПРОГРАММУ SAGE
valCPP = np.genfromtxt("graph_val.txt")

#check_txt ФАЙЛ ДЛЯ ПЕРЕДАЧИ ЗНАЧЕНИЙ ТОЧЕК ПРОВЕРКИ
checkCPP = np.genfromtxt("check.txt")

#================================================== END C++
#================================================== ГРАФИКА

import matplotlib.pyplot as plt

fig, axs = plt.subplots(4,2)

axs[0,0].set_title("Сила")
axs[1,0].set_title("Момент")
axs[2,0].set_title("Угол")
axs[3,0].set_title("Прогиб")

width = 2                                                               # толщина графика

valSAGE = np.zeros(N)

for i in range(0,N):                                                    # построение силы
    valSAGE[i] = Q(x = arg[i]).subs(heaviside(0.0) == 0)
axs[0,0].plot(arg, valSAGE, linewidth = width)

axs[0,1].plot(arg, valCPP[:,3], linewidth = width)

for i in range(0,N):                                                    # построение моментов
	valSAGE[i] = MM(x = arg[i]).subs(heaviside(0.0) == 0)
axs[1,0].plot(arg, valSAGE, linewidth = width)

axs[1,1].plot(arg, valCPP[:,2], linewidth = width)

for i in range(0,N):                                                    # построение углов
	valSAGE[i] = teta(x = arg[i]).subs(heaviside(0.0) == 0)
axs[2,0].plot(arg, valSAGE, linewidth = width)

axs[2,1].plot(arg, valCPP[:,1], linewidth = width)
axs[2,1].plot(checkCPP[:,0], checkCPP[:,2], 'ko', markersize = 3)

for i in range(0,N):                                                    # построение прогиба
	valSAGE[i] = w(x = arg[i]).subs(heaviside(0.0) == 0)
axs[3,0].plot(arg, valSAGE, linewidth = width)

axs[3,1].plot(arg, valCPP[:,0], linewidth = width)
axs[3,1].plot(checkCPP[:,0], checkCPP[:,1], 'ko', markersize = 3)

for i in range(0,len(axs)):
    arg = np.array([AR, BR, CR, DR], dtype = float)							# построение опор
    valSAGE = np.zeros(4)
    axs[i,0].plot(arg, valSAGE, 'ro', label = "R")
    axs[i,1].plot(arg, valSAGE, 'ro', label = "R")

    arg = np.array([AP], dtype = float)										# построение силы P
    valSAGE = np.zeros(1)
    axs[i,0].plot(arg,valSAGE, 'bo', label = "P")
    axs[i,1].plot(arg, valSAGE, 'bo', label = "P")

    arg = np.array([AM], dtype = float)										# построение момента M
    valSAGE = np.zeros(1)
    axs[i,0].plot(arg,valSAGE, 'mo', label = "M")
    axs[i,1].plot(arg, valSAGE, 'mo', label = "M")

    arg = np.array([Ak], dtype = float)										# построение пружины k
    valSAGE = np.zeros(1)
    axs[i,0].plot(arg,valSAGE, 'yo', label = "k")
    axs[i,1].plot(arg, valSAGE, 'yo', label = "k")

    
    argq = np.array([Aq, Bq], dtype = float)									        # построение распределенной нагрузки (только на SAGE)
    valq = np.array([0.0, 0.0], dtype = float)
    axs[i,0].plot(argq,valq, color = 'g', alpha = 0.3, label = "q")
    
    argq = np.array([Aq, Bq], dtype = float)                                         # построение распределенной нагрузки (только на SAGE)
    valq = np.array([0.0, 0.0], dtype = float)
    axs[i,1].plot(argq,valq, color = 'g', alpha = 0.3,label = "q")

    axs[i,0].grid()
    axs[i,1].grid()
    axs[i,0].legend()
    axs[i,1].legend()
    axs[i,0].set_xlim(0,25)
    axs[i,1].set_xlim(0,25)

plt.show()

#==================================================
print("STOP SAGE")