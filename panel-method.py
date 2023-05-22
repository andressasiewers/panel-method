import numpy as np
import matplotlib.pyplot as plt

n = 20  # panels
phi1 = np.pi / 2  # 90°

X = np.zeros(n+1)
Y = np.zeros(n+1)
x = np.zeros(n)
y = np.zeros(n)
phi = np.zeros(n)
bet = np.zeros(n)
U = np.zeros(n)

alpha = (2 * np.pi) / n  # internal angle

for i in range(n):
    X[i] = np.cos(np.pi + (alpha / 2) - (i - 1) * alpha)
    Y[i] = np.sin(np.pi + (alpha / 2) - (i - 1) * alpha)

X[n] = X[0]
Y[n] = Y[0]

for i in range(n):
    x[i] = (X[i + 1] + X[i]) / 2
    y[i] = (Y[i + 1] + Y[i]) / 2
    phi[i] = phi1 - (i - 1) * alpha
    bet[i] = phi[i] + (np.pi / 2)
    U[i] = -2 * np.pi * np.cos(bet[i])

S = np.sqrt((X[1] - X[0]) ** 2 + (Y[1] - Y[0]) ** 2)  # panel length

A = np.zeros((n, n))
B = np.zeros((n, n))
C = np.zeros((n, n))
D = np.zeros((n, n))
E = np.zeros((n, n))
I1 = np.zeros((n, n))
I2 = np.zeros((n, n))
I3 = np.zeros((n, n))
Iijn = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        if i == j:
            Iijn[i, j] = np.pi
        else:
            A[i, j] = -((x[i] - X[j]) * np.cos(phi[j])) - ((y[i] - Y[j]) * np.sin(phi[j]))
            B[i, j] = (x[i] - X[j]) ** 2 + (y[i] - Y[j]) ** 2
            C[i, j] = np.sin(phi[i] - phi[j])
            D[i, j] = (y[i] - Y[j]) * np.cos(phi[i]) - (x[i] - X[j]) * np.sin(phi[i])
            E[i, j] = np.sqrt(B[i, j] - A[i, j] ** 2)
            I1[i, j] = (C[i, j] / 2) * np.log((S ** 2 + (2 * A[i, j] * S) + B[i, j]) / B[i, j])
            I2[i, j] = np.arctan((S + A[i, j]) / E[i, j])
            I3[i, j] = np.arctan(A[i, j] / E[i, j])
            Iijn[i, j] = I1[i, j] + ((D[i, j] - (A[i, j] * C[i, j])) / E[i, j]) * (I2[i, j] - I3[i, j])

lambida = np.linalg.inv(Iijn) @ U

Vs = np.zeros(n)
IS = np.zeros((n, n))
IS1 = np.zeros((n, n))
IS2 = np.zeros((n, n))
IS3 = np.zeros((n, n))
cp = np.zeros(n)

for i in range(n):
    for j in range(n):
        if i == j:
            IS[i, j] = 0
        else:
            IS1[i, j] = (D[i, j] - A[i, j] * C[i, j]) / (2 * E[i, j])
            IS2[i, j] = (S ** 2 + 2 * A[i, j] * S + B[i, j]) / B[i, j]
            IS3[i, j] = np.arctan((S + A[i, j]) / E[i, j]) - np.arctan(A[i, j] / E[i, j])
            IS[i, j] = IS1[i, j] * np.log(IS2[i, j]) - (C[i, j] * IS3[i, j])
        Vs[i] += lambida[j] * IS[i, j] / (2 * np.pi)

    Vs[i] += np.sin(bet[i])
    cp[i] = 1 - (Vs[i] ** 2)

teta = np.linspace(-np.pi, np.pi, 100)
cpt = 1 - 4 * (np.sin(teta)) ** 2  # theoretical values of Cp

# Plotting
plt.plot(teta, cpt, ':b', linewidth=3, label='Cp theoretical')
plt.plot(bet, cp, 'ro', markersize=3, label='Cp each panel')
plt.xlabel('teta', fontsize=14)
plt.ylabel('Cp', fontsize=14)
plt.title('Panel method Cp x Theoretical Cp', fontsize=16)
plt.legend()
plt.show()