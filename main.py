from numpy import zeros, linspace, sqrt
from matplotlib.pyplot import style, axes
import matplotlib.pyplot as plt

# Функция f подготавливает массив, содержащий элементы вектор-функции,
# определяющей правую часть решаемой системы ОДУ
def f(u) :
    f = zeros(3)
    f[0] = u[1]
    f[1] = (u[2]**4*u[1]**2 - 4*u[2]**3*u[0]*u[1] - 1)/(2*u[2]**4*u[0])
    f[2] = 1
    return f

def g(u):
    g = zeros(3)
    g[0] = u[1]
    g[1] = (u[2] ** 4 * u[1] ** 2 - 4 * u[2] ** 3 * u[0] * u[1] - 1) / (2 * u[2] ** 4 * u[0])
    g[2] = 1
    g = g/sqrt(1 + g[0]**2 + g[1]**2)
    return g

dl = 0.5
# Определение входных данных задачи
x_0 = 0.5; X = 1.2

# Определение числа интервалов сетки,
# на которой будет искаться приближённое решение
M = 30
J = 1000


# Определение сетки
tau = (X - x_0)/M
x = linspace(x_0, X, M + 1)
# Выделение памяти под массив сеточных значений решения системы ОДУ
# В строке с номером m этого массива хранятся сеточные значения решения,
# соответствующие моменту времени t_m
u1 = zeros((M + 1, 3))
u2 = zeros((M + 1, 3))
u3 = zeros((M + 1, 3))
u4 = zeros((J, 3))

# Задание начальных условий
#(записываются в строку с номером 0 массива V)
u1[0] = [-1.5, 8, 0.5]
u2[0] = [-1.5, 8, 0.5]
u3[0] = [-1.5, 8, 0.5]
u4[0] = [-1.5, 8, 0.5]


for m in range(M) :
    #Реализация схемы Эйлера
    u1[m + 1] = u1[m] + tau*f(u1[m])

    #Реализация схемы ERK2
    w_1_erk2 = f(u2[m])
    w_2_erk2 = f(u2[m] + tau * 2/3 * w_1_erk2)
    u2[m + 1] = u2[m] + tau * (1/4 * w_1_erk2 + 3/4 * w_2_erk2)

    #Реализация схемы ERK4
    w_1_erk4 = f(u3[m])
    w_2_erk4 = f(u3[m] + tau * 1/2 * w_1_erk4)
    w_3_erk4 = f(u3[m] + tau * 1/2 * w_2_erk4)
    w_4_erk4 = f(u3[m] + tau * 1 * w_3_erk4)
    u3[m + 1] = u3[m] + tau * (1/6 * w_1_erk4 + 1/3 * w_2_erk4 + 1/3 * w_3_erk4 + 1/6 * w_4_erk4)

#Реализация перехода на длину дуги кривой
j = 0
while u4[j,2] < X:
    w_1 = g(u4[j])
    w_2 = g(u4[j] + dl*2/3*w_1)
    u4[j + 1] = u4[j] + dl*(1/4*w_1 + 3/4*w_2)
    j = j + 1





#Отрисовка решения
style.use('dark_background')
fig = plt.figure()
ax = axes()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.scatter(u1[:,2], u1[:,0], color = 'blue')
ax.scatter(u2[:,2], u2[:,0], color = 'yellow')
ax.scatter(u4[0:j+1,2], u4[0:j+1,0], color = 'cyan')
ax.scatter(u3[:,2], u3[:,0], color = 'green')
ax.scatter(x,  (x**2 - 1)/(2*x**2), color = 'red')
plt.show()



