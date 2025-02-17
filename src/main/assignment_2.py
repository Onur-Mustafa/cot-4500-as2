def neville(x_points, y_points, x):
    n = len(x_points)
    p = [[0] * n for _ in range(n)]

    for i in range(n):
        p[i][i] = y_points[i]

    for j in range(1, n):
        for i in range(n - j):
            p[i][i + j] = ((x - x_points[i]) * p[i + 1][i + j] - (x - x_points[i + j]) * p[i][i + j - 1]) / (x_points[i + j] - x_points[i])

    return p[0][n - 1]

x_points = [3.6, 3.8, 3.9]
y_points = [1.675, 1.436, 1.318]
x = 3.7

result = neville(x_points, y_points, x)
print(f"{result}\n")

def newton_forward(x_points, y_points):
    n = len(x_points)
    diff_table = [y_points[:]]

    for i in range(1, n):
        diff_table.append([
            (diff_table[i - 1][j + 1] - diff_table[i - 1][j]) / (x_points[j + i] - x_points[j])
            for j in range(n - i)
        ])

    def polynomial(x, degree):
        result = y_points[0]
        product = 1
        for i in range(1, degree + 1):
            product *= (x - x_points[i - 1])
            result += diff_table[i][0] * product
        return result

    return [polynomial(x_points[0], d) for d in range(1, n)], polynomial

from math import factorial

x_points = [7.2, 7.4, 7.5, 7.6]
y_points = [23.5492, 25.3913, 26.8224, 27.4589]

results, polynomial = newton_forward(x_points, y_points)
for degree, result in enumerate(results, start=1):
    print(f"{degree}: {result}")
print() 

# Approximate f(7.3)
approximation = polynomial(7.3, 3)
print(f"{approximation}\n")

def hermite_interpolation(x_points, y_points, y_derivatives):
    n = len(x_points)
    z = [0] * (2 * n)
    q = [[0] * (2 * n) for _ in range(2 * n)]

    for i in range(n):
        z[2 * i] = z[2 * i + 1] = x_points[i]
        q[2 * i][0] = q[2 * i + 1][0] = y_points[i]
        q[2 * i + 1][1] = y_derivatives[i]
        if i != 0:
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1])

    for i in range(2, 2 * n):
        for j in range(2, i + 1):
            q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / (z[i] - z[i - j])

    return q

x_points = [3.6, 3.8, 3.9]
y_points = [1.675, 1.436, 1.318]
y_derivatives = [-1.195, -1.188, -1.182]

hermite_matrix = hermite_interpolation(x_points, y_points, y_derivatives)
print()
for row in hermite_matrix:
    print(row)
print() 

import numpy as np

def cubic_spline(x_points, y_points):
    n = len(x_points) - 1
    h = [x_points[i+1] - x_points[i] for i in range(n)]
    
    A = np.zeros((n+1, n+1))
    A[0, 0] = 1
    A[n, n] = 1
    for i in range(1, n):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
    
    b = np.zeros(n+1)
    for i in range(1, n):
        b[i] = 3 * ((y_points[i+1] - y_points[i]) / h[i] - (y_points[i] - y_points[i-1]) / h[i-1])
    
    x = np.linalg.solve(A, b)
    
    return A, b, x

x_points = [2, 5, 8, 10]
y_points = [3, 5, 7, 9]

A, b, x = cubic_spline(x_points, y_points)

print("")
print(A)
print("\n")
print(b)
print("\n")
print(x)
