import numpy as np
import matplotlib.pyplot as plt

"""
Problem:
    -div(grad(u)) = 2 in [0,1]
    u(0) = 0, u(1) = 1.

Exact solution:
    u(x) = x(1 - x)

Galerkin approach:
    - Divide [0,1] into 2^n sub-intervals that are equally spaced,
      i.e. create a partition of points that break [0,1] into 2^n sub-intervals
        0 = x_0 < x_1 < ... < x_2^n = 1
        x_0 = 0, x_1 = 1/2^n, x_2 = 2/2^n, ..., x_j = j/2^n
    - To construct the nodal basis functions, now consider only the j = 1,..., 2^n - 1 points in the partition (excludes the boundary points)
    - Now construct the stiffness matrix A and the vector f to calculate the nodal basis coefficients U in the system AU = f
        - To construct the stiffness matrix A, we will use lemma 3.9 from Bartels to calculate the gradients of each nodal basis function
        - In this case it is easy to calculate that each element in f will always be equal to 2^(1-n)
    - After the coefficients are calculated, reconstruct the solution at some x in [0,1] again using lemma 3.9 from Bartels, and compare against the exact solution
"""

n = int(input("Enter n value: "))

# Partition the domain
partition = [j/(2**n) for j in range(0, 2**n + 1)]
#print(partition)

# Keep track of all the sub-intervals
intervals = []
for i in range(0, 2**n):
    intervals.append([partition[i], partition[i+1]])
#print(intervals)

# Initialize the stiffness matrix
A = [[0]*(2**n - 1) for j in range(0, 2**n - 1)]
#print(A)

# Compute all the entries in the stiffness matrix
for j in range(1, 2**n):
    for k in range(1, 2**n):
        if(j == k):
            A[j-1][k-1] = 2**(n+1)
        if(k == j - 1 or k == j + 1):
            A[j-1][k-1] = -2**n
        """
        sum_Ajk = 0
        for interval in intervals:
            if(partition[j] in interval and partition[k] in interval):
                dot_grads = (-interval[0]/(interval[1] - interval[0])) * 1/(interval[1] - interval[0])
                length_interval = interval[1] - interval[0]
                sum_Ajk += dot_grads * length_interval
        A[j-1][k-1] = sum_Ajk
        #A[j-1][k-1] = 1
        """

f = [2**(1-n)]*(2**n - 1)

#for j in range(0, 2**n - 1):
#    print(A[j])
#print(f)

A_np = np.array(A)
f_np = np.array(f)
U_np = np.linalg.solve(A_np, f_np)
#print(np.allclose(np.dot(A_np,U_np), f_np))
U = U_np.tolist()
#print(U)

def phi_of_x(x1, x2, k, x):
    if(k == 0):
        if((x - x1)/(x2 - x1) >= 0 and (x - x1)/(x2 - x1) <= 1):
            return 1 - ((x - x1)/(x2 - x1))
        return 0
    if((x - x1)/(x2 - x1) >= 0 and (x - x1)/(x2 - x1) <= 1):
        return ((x - x1)/(x2 - x1))
    return 0

x = 0
x_coords = []
y_coords = []
while(x < 1):
    galerkin_solution = 0
    for j in range(1, 2**n):
        for k in range(0,1):
            galerkin_solution += U[j-1+k]*phi_of_x(intervals[j][0], intervals[j][1], k, x)
    #print("U({}) = {}".format(x, galerkin_solution))
    x_coords.append(x)
    y_coords.append(galerkin_solution)
    x += 0.001


fig, ax = plt.subplots()
ax.plot(x_coords, y_coords, color='b')
x = np.linspace(0, 1)
y = x*(1 - x)
ax.plot(x, y, color='r')

plt.show()
