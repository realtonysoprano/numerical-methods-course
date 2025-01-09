import numpy as np
import matplotlib.pyplot as plt

def F(iA, iB, iC, x):
    return iC + iB * np.exp(iA * x)

def e_interp(x, y):

    A = np.zeros(len(x) - 2)
    B = np.zeros(len(x) - 2)
    C = np.zeros(len(x) - 2)

    for i in range(len(A)):
        dif = (y[i + 2] - y[i + 1]) * (x[i + 1] - x[i]) / ((y[i + 1] - y[i]) * (x[i + 2] - x[i + 1]))

        if abs(dif) == 1:
            A[i] = (y[i + 2] - y[i]) / (x[i + 2] - x[i])
            B[i] = y[i] - A[i] * x[i]
            print('Допустима только линейная интерполяция')
            C[i] = 0

        elif x[i + 2] + x[i] == x[i + 1] + x[i + 1]:
            z = (2 * y[i + 1] - y[i] - y[i + 2])
            r = (y[i + 2] - y[i + 1]) / (y[i + 1] - y[i])
            A[i] = np.log(r) / (x[i + 2] - x[i + 1])
            C[i] = (y[i + 1] * y[i + 1] - y[i] * y[i + 2]) / z
            B[i] = (y[i] - C[i]) * np.exp(-A[i] * x[i])

        else:

            Amin = np.log(dif) / (x[i + 2] - x[i])
            A0 = 2 * Amin
            while (1):
                u = np.exp(A0 * (x[i + 2] - x[i + 1]))
                v = np.exp(-A0 * (x[i + 1] - x[i]))
                F = (y[i + 1] - y[i]) * (u - 1) + (y[i + 2] - y[i + 1]) * (v - 1)
                FF = (y[i + 1] - y[i]) * (x[i + 2] - x[i + 1]) * u - (y[i + 2] - y[i + 1]) * (x[i + 1] - x[i]) * v
                dA = -F / FF
                A0 += dA
                if (abs(dA / Amin) < 0.000000001):
                    break

            A[i] = A0
            B[i] = (y[i] - y[i + 1]) / (np.exp(A0 * x[i]) - np.exp(A0 * x[i + 1]))
            C[i] = y[i] - B[i] * np.exp(A0 * x[i])

    return A, B, C


def get_graphics(A, B, C, x, y):

    n = len(A)
    x1 = np.arange(x[0], x[1] + 0.01, 0.01)
    y1 = C[0] + B[0] * np.exp(A[0] * x1)
    plt.title("Экспоненциальная интерполяция:")
    plt.plot(x1, y1, color="red")
    plt.scatter(x, y, c='b')

    for i in range(1, n):
        x1 = np.arange(x[i], x[i + 1] + 0.01, 0.01)
        y1 = ((x[i + 1] - x1) * (C[i - 1] + B[i - 1] * np.exp(A[i - 1] * x1)) + (x1 - x[i]) * (C[i] + B[i] * np.exp(A[i] * x1))) / (x[i + 1] - x[i])
        plt.plot(x1, y1, color="red")

    x1 = np.arange(x[-2], x[-1] + 0.01, 0.01)
    y1 = C[n - 1] + B[n - 1] * np.exp(A[n - 1] * x1)
    plt.plot(x1, y1, color="red")
    plt.grid()
    plt.show()


x = [0, 1, 2, 3, 4, 5]
y = [1, 4, 6, 10, 12, 15]

x1 = [0, 4, 7, 9, 10]
y1 = [2, 4, 6, 15, 20]

x3 = [0, 1, 2, 3, 3.5, 5, 7, 10]
y3 = [500, 1000, 2230, 3040, 3974, 4587, 5000, 5500]


A, B, C  = e_interp(x, y)
for i in range(len(x)-2):
    print(f'Экспоненциальная интерп. функция {i+1}-го сегмента: y = {C[i]} + {B[i]} * exp({A[i]} * x)')

get_graphics(A, B, C, x, y)
