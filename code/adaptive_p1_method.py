import numpy as np
import matplotlib.pyplot as plt
from shapely import geometry

"""
Consider the model problem
    -div(grad(u)) = f  in [0,1]x[0,1]
                u = 0  on the boundary of [0,1]x[0,1].

For certain f we have exact solutions for u, however, we wish to use AFEM to approximate the solution to u.
To do this we implement the algorithm Solve-Estimate-Mark-Refine.

The algorithm is described as follows:
Given a triangulation T_0 of [0,1]x[0,1] and parameters 0 < p < 1  and end_k > 0, set k = 0.
Then
    1. Compute the Galerkin approximation U satisfying the weak formulation of the model problem.
    2. For every element t in T_k compute the refinement indicator E_t(U).
    3. Choose a set of marked elements M_k in T_k such that E_t(U) > p*max(E_t'(U)) where t' also exists in T_k.
    4. Apply red-refinement to every t in M_k to produce a new triangulation.
    5. Set k = k + 1, and if k == end_k then stop, otherwise go back to 1.
"""
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def equals(self, q):
        if(self.x == q.x and self.y == q.y):
            return True
        return False
    def midpoint(self, q):
        return Point((q.x + self.x)/2, (q.y + self.y)/2)
    def distance(self, q):
        return np.sqrt((q.x - self.x)**2 + (q.y - self.y)**2)

class Triangle:
    def __init__(self, p, q, r):
        self.p = p
        self.q = q
        self.r = r
    def det_X_t(self):
        return np.linalg.det(np.array([[1, self.p.x, self.p.y], [1, self.q.x, t.q.y], [1, self.r.x, self.r.y]]))
    def vol_t(self):
        return self.det_X_t()/2
    def grad_phi(self, j):
        points = [self.p, self.q, self.r]
        return [(points[(j+1)%3].y - points[(j+2)%3].y)/(self.det_X_t()),
                (points[(j+2)%3].x - points[(j+1)%3].x)/(self.det_X_t())]
    def vertex_of_t(self, s):
        points = [self.p, self.q, self.r]
        for i in range(0, 3):
            if(points[i].equals(s)):
                return i
        return -1
    def x_in_closure_t(self, x):
        a = ( ((self.q.y - self.r.y)*(x.x - self.r.x) + (self.r.x - self.q.x)*(x.y - self.r.y)) / 
             ((self.q.y - self.r.y)*(self.p.x - self.r.x) + (self.r.x - self.q.x)*(self.p.y - self.r.y)) )
        b = ( ((self.r.y - self.p.y)*(x.x - self.r.x) + (self.p.x - self.r.x)*(x.y - self.r.y)) / 
             ((self.q.y - self.r.y)*(self.p.x - self.r.x) + (self.r.x - self.q.x)*(self.p.y - self.r.y)))
        c = 1 - a - b
        if((a >= 0 and a <= 1) and (b >= 0 and b <= 1) and (c >= 0 and c <= 1)):
            return True
        return False
    def red_refine(self):
        if(np.maximum(self.p.distance(self.q), np.maximum(self.q.distance(self.r), self.r.distance(self.p))) ==  self.p.distance(self.q)):
            points = [self.p, self.q, self.r]
        elif(np.maximum(self.p.distance(self.q), np.maximum(self.q.distance(self.r), self.r.distance(self.p))) ==  self.q.distance(self.r)):
            points = [self.q, self.r, self.p]
        else:
            points = [self.r, self.p, self.q]
        mp0 = points[0].midpoint(points[1]) 
        mp1 = points[1].midpoint(points[2]) 
        mp2 = points[2].midpoint(points[0]) 
        return [Triangle(points[0],mp0,mp2), Triangle(mp0,points[1],mp1), Triangle(mp0,mp1,mp2), Triangle(mp1, points[2], mp2)], [mp0, mp1, mp2]

def galerkin_basis_coefficients(T, vEb, scalar_f):
    local_f_vectors = []
    for t in T:
        local_f = [0]*len(vEb)
        for i in range(0, len(vEb)):
            #print("({},{})({},{})({},{}) ({},{})".format(t.p.x, t.p.y, t.q.x, t.q.y, t.r.x, t.r.y, vEb[i].x, vEb[i].y))
            #print(t.vertex_of_t(vEb[i]))
            if(t.vertex_of_t(vEb[i]) >= 0):
                local_f[i] = scalar_f*(t.vol_t()/3)
        local_f_vectors.append(local_f)
    f = [0]*len(vEb)
    for local_f in local_f_vectors:
        f = np.add(np.array(f), np.array(local_f)).tolist()
    local_A_matrices = []
    for t in T:
        local_A = [[0]*len(vEb) for i in range(0, len(vEb))]
        for i in range(0, len(vEb)):
            k = t.vertex_of_t(vEb[i])
            for j in range(0, len(vEb)):
                l = t.vertex_of_t(vEb[j])
                if(k != -1 and l != -1):
                    grad_phi_vi = t.grad_phi(k)
                    grad_phi_vj = t.grad_phi(l)
                    dot_grads = (grad_phi_vi[0] * grad_phi_vj[0]) + (grad_phi_vi[1] * grad_phi_vj[1])
                    local_A[i][j] = dot_grads * t.vol_t()
        local_A_matrices.append(local_A)
    A = [[0]*len(vEb) for i in range(0, len(vEb))]
    for local_A in local_A_matrices:
        A = np.add(np.array(A),np.array(local_A)).tolist()
    #print(A)
    #print(f)
    U = np.linalg.solve(np.array(A), np.array(f)).tolist()
    return U

def phi_of_x(T, v, x):
    for t in T:
        j = t.vertex_of_t(v)
        if(j != False):
            if(t.x_in_closure_t(x)):
                points = [t.p, t.q, t.r]
                det1 = np.linalg.det(np.array([
                    [1, x.x, x.y],
                    [1, points[(j+1)%3].x, points[(j+1)%3].y],
                    [1, points[(j+2)%3].x, points[(j+2)%3].y]]))
                det2 = np.linalg.det(np.array([
                    [1, points[j].x, points[j].y],
                    [1, points[(j+1)%3].x, points[(j+1)%3].y],
                    [1, points[(j+2)%3].x, points[(j+2)%3].y]]))
                return det1/det2
    return 0


def recreate_galerkin_solution_at_x(T, U, vEb, x):
    galerkin_solution = 0
    for i in range(0, len(vEb)):
        galerkin_solution += U[i] * phi_of_x(T, vEb[i], x)
    return galerkin_solution


p0 = Point(0,0)
p1 = Point(1,0)
p2 = Point(1,1)
p3 = Point(0,1)
vertices = [p0, p1, p2, p3]
domain = geometry.Polygon([(0,0),(1,0),(1,1),(0,1)])
boundary = domain.boundary

t0 = Triangle(p0,p1,p2)
t1 = Triangle(p0,p2,p3)

T = [t0, t1]

for k in range(0,6):
    T_copy = T.copy()
    for t in T_copy:
        new_t, new_vertices = t.red_refine()
        for t1 in new_t:
            T.append(t1)
        for v in  new_vertices:
            vertices.append(v)
        T.remove(t)

vEb = []
for v in vertices:
    if(boundary.contains(geometry.Point(v.x, v.y)) == False):
        flag = False
        for v1 in vEb:
            if(v.equals(v1)):
                flag = True
                break
        if(flag == False):
            vEb.append(v)

#print(len(T))
#print(len(vertices))
#print(len(vEb))
U = galerkin_basis_coefficients(T, vEb, 4)
u0 = recreate_galerkin_solution_at_x(T, U, vEb, Point(0.5,0.5))
u1 = recreate_galerkin_solution_at_x(T, U, vEb, Point(0.25,0.25))
u2 = recreate_galerkin_solution_at_x(T, U, vEb, Point(0.75,0.75))
u3 = recreate_galerkin_solution_at_x(T, U, vEb, Point(0.0,0.0))
u4 = recreate_galerkin_solution_at_x(T, U, vEb, Point(1,1))
#print(U)
print(u0, u1, u2, u3, u4)

fig, ax = plt.subplots()
for t in T:
    ax.plot([t.p.x,t.q.x],[t.p.y,t.q.y],color='black',linestyle='-')
    ax.plot([t.q.x,t.r.x],[t.q.y,t.r.y],color='black',linestyle='-')
    ax.plot([t.r.x,t.p.x],[t.r.y,t.p.y],color='black',linestyle='-')

plt.show()
