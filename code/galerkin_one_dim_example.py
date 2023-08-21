import numpy as np
import matplotlib.pyplot as plt

"""
Problem:
    -div(grad(u)) = 2 in (0,1)
    u(0) = 0, u(1) = 1

Exact solution:
    u(x) = x(1 - x)

Galerkin approach:
    - Divide [0,1] into 2^n sub-intervals that are equally spaced,
      i.e. create a partition of points that break [0,1] into 2^n sub-intervals.
        0 = x_0 < x_1 < ... < x_2^n = 1
        x_0 = 0, x_1 = 1/2^n, x_2 = 2/2^n, ..., x_j = j/2^n
    - To construct the nodal basis functions, now consider only the j = 1,..., 2^n - 1 points in the partition (excludes the boundary points).
    - Now construct the stiffness matrix A and the vector f to calculate the nodal basis coefficients U in the system AU = f.
    - After the coefficients are calculated, reconstruct the solution at some x in [0,1].
"""

number_of_figures = 4
fig, ax = plt.subplots(1, number_of_figures)
for n in range(1, number_of_figures + 1):
    partition = [j/(2**n) for j in range(0, 2**n + 1)]
    PNB = [i for i in range(1, 2**n)]

    A = [[0]*(2**n - 1) for j in range(0, 2**n - 1)]

    for j in range(1, 2**n):
        for k in range(1, 2**n):
            if(j == k):
                A[j-1][k-1] = 2**(n+1)
            if(k == j - 1 or k == j + 1):
                A[j-1][k-1] = -2**n

    f = [2**(1-n)]*(2**n - 1)

    A_np = np.array(A)
    f_np = np.array(f)
    U_np = np.linalg.solve(A_np, f_np)
    U = U_np.tolist()

    def phi_j(x, j):
        total = 0
        if(x > partition[j] and x <= partition[j+1]):
            total += (partition[j+1] - x)/(partition[j+1] - partition[j])
        if(x >= partition[j-1] and x < partition[j]):
            total += (x - partition[j-1])/(partition[j] - partition[j-1])
        if(x == partition[j]):
            total += 1
        return total

    x = 0
    x_coords = []
    y_coords = []
    while(x <= 1):
        galerkin_solution = 0
        for j in range(0, 2**n - 1):
            galerkin_solution += U[j]*phi_j(x, PNB[j])
        #print("U({}) = {}".format(x, galerkin_solution))
        x_coords.append(x)
        y_coords.append(galerkin_solution)
        x += 0.001

    x = np.linspace(0, 1)
    y = x*(1 - x)
    ax[n-1].plot(x, y, color='r', label="Exact Solution")
    ax[n-1].plot(x_coords, y_coords, color='b', label="Galerkin Approximation")
    ax[n-1].set_title("{} Elements".format(2**n))

plt.legend()
plt.show()
