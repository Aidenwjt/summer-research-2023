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
        return np.linalg.det(np.array([[1, self.p.x, self.p.y], [1, self.q.x, self.q.y], [1, self.r.x, self.r.y]]))
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
    def A(self, i, j):
        grad_phi_i = self.grad_phi(i)
        grad_phi_j = self.grad_phi(j)
        dot_grads = (grad_phi_i[0] * grad_phi_j[0]) + (grad_phi_i[1] * grad_phi_j[1])
        return (dot_grads * self.vol_t())/2
        """
        b1 = self.r.y - self.p.y
        b2 = self.p.y - self.q.y
        c1 = self.p.x - self.r.x
        c2 = self.q.x - self.p.x
        if(i == 0 and j == 0):
            return ((-b1 - b2)**2 + (-c1 - c2)**2)/(4*self.vol_t())
        if((i == 1 and j == 0) or (i == 0 and j == 1)):
            return (b1*(-b1 - b2) + c1*(-c1 - c2))/(4*self.vol_t())
        if(i == 1 and j == 1):
            return(b1**2 + c1**2)/(4*self.vol_t())
        if((i == 2 and j == 1) or (i == 1 and j == 2)):
            return (b1*b2 + c1*c2)/(4*self.vol_t())
        if(i == 2 and j == 2):
            return(b2**2 + c2**2)/(4*self.vol_t())
        if((i == 2 and j == 0) or (i == 0 and j == 2)):
            return (b2*(-b1 - b2) + c2*(-c1 - c2))/(4*self.vol_t())
        """
    def red_refine(self, vertices, P_copy, P, k, bv, boundary):
        if(np.maximum(self.p.distance(self.q), np.maximum(self.q.distance(self.r), self.r.distance(self.p))) ==  self.p.distance(self.q)):
            v_indices = [P_copy[k][0], P_copy[k][1], P_copy[k][2]]
            v = [self.p, self.q, self.r]
        elif(np.maximum(self.p.distance(self.q), np.maximum(self.q.distance(self.r), self.r.distance(self.p))) ==  self.q.distance(self.r)):
            v_indices = [P_copy[k][1], P_copy[k][2], P_copy[k][0]]
            v = [self.q, self.r, self.p]
        else:
            v_indices = [P_copy[k][2], P_copy[k][0], P_copy[k][1]]
            v = [self.r, self.p, self.q]
        new_v = [v[0].midpoint(v[1]), v[1].midpoint(v[2]), v[2].midpoint(v[0])]
        index0 = -1
        index1 = -1
        index2 = -1
        for i in range(0,len(vertices)):
            if(new_v[0].equals(vertices[i])):
                index0 = i
            if(new_v[1].equals(vertices[i])):
                index1 = i
            if(new_v[2].equals(vertices[i])):
                index2 = i
        if(index0 == -1):
            vertices.append(new_v[0])
            index0 = len(vertices) - 1
            if(boundary.contains(geometry.Point(new_v[0].x, new_v[0].y))):
                bv.append(index0)
        if(index1 == -1):
            vertices.append(new_v[1])
            index1 = len(vertices) - 1
            if(boundary.contains(geometry.Point(new_v[1].x, new_v[1].y))):
                bv.append(index1)
        if(index2 == -1):
            vertices.append(new_v[2])
            index2 = len(vertices) - 1
            if(boundary.contains(geometry.Point(new_v[2].x, new_v[2].y))):
                bv.append(index2)
        P.append([v_indices[0], index0, index2])
        P.append([index0, v_indices[1], index1])
        P.append([index0, index1, index2])
        P.append([index1, v_indices[2], index2])
        return [Triangle(v[0], new_v[0], new_v[2]), Triangle(new_v[0], v[1], new_v[1]), Triangle(new_v[0], new_v[1], new_v[2]), Triangle(new_v[1], v[2], new_v[2])]

def phi_of_x(T, v, x):
    for t in T:
        j = t.vertex_of_t(v)
        #print("({},{}),({},{}),({},{})".format(t.p.x, t.p.y, t.q.x, t.q.y, t.r.x, t.r.y))
        #print("({},{})".format(v.x,v.y))
        #print(j)
        if(j > -1):
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


def recreate_galerkin_solution_at_x(T, U, vertices, x):
    galerkin_solution = 0
    for i in range(0, len(vertices)):
        if(U[i] != 0):
            #print(U[i])
            galerkin_solution += U[i] * phi_of_x(T, vertices[i], x)
    return galerkin_solution

scalar_f = 4
p0 = Point(0,0)
p1 = Point(1,0)
p2 = Point(1,1)
p3 = Point(0,1)
domain = geometry.Polygon([(0,0),(1,0),(1,1),(0,1)])
boundary = domain.boundary
t0 = Triangle(p0,p1,p2)
t1 = Triangle(p0,p2,p3)
T = [t0, t1]
vertices = [p0, p1, p2, p3]
P = [[0, 1, 2], [0, 2, 3]]
bv = [0, 1, 2, 3]

for l in range(0,6):
    T_copy = T.copy()
    P_copy = P.copy()
    for k in range(0, len(T_copy)):
    #for t in T_copy:
        new_t = T_copy[k].red_refine(vertices, P_copy, P, k, bv, boundary)
        for t1 in new_t:
            T.append(t1)
        T.remove(T_copy[k])
        P = P[1:]
    """
    for t in T:
        print("({},{}),({},{}),({},{})".format(t.p.x,t.p.y,t.q.x,t.q.y,t.r.x,t.r.y))
    print("---")
    for v in vertices:
        print("({},{})".format(v.x,v.y))
    print("---")
    for i in range(0, len(T)):
        print(P[i])
    print("---")
    """
    A = [[0]*len(vertices) for i in range(0, len(vertices))]
    f = [0]*len(vertices)
    for k in range(0, len(T)):
        for j in range(0, 3):
            for i in range(0, 3):
                A[P[k][i]][P[k][j]] += T[k].A(i, j)
            f[P[k][j]] += ((scalar_f*T[k].vol_t())/3)
    #for i in range(0, len(vertices)):
    #    print(A[i])
    #print("---")
    #print(f)
    #print("---")
    for i in range(0, len(bv)):
        for j in range(0, len(vertices)):
            if(j == bv[i]):
                A[j][bv[i]] = 1
            else:
                A[j][bv[i]] = 0
                A[bv[i]][j] = 0
            #A[vEb[i]-1][vEb[i]] = 0
            #A[vEb[i]][vEb[i]-1] = 0
            #A[vEb[i]+1][vEb[i]] = 0
            #A[vEb[i]][vEb[i]+1] = 0
        f[bv[i]] = 0
    """
    for i in range(0, len(vertices)):
        print(A[i])
    print("---")
    print(f)
    print("---")
    """
    U = np.linalg.solve(np.array(A), np.array(f)).tolist()
    #print(U)
    #print(len(T))
    #print(len(vertices))
    #print(len(P))
    #print(len(bv))
    u0 = recreate_galerkin_solution_at_x(T, U, vertices, Point(0.5,0.5))
    #u1 = recreate_galerkin_solution_at_x(T, U, vertices, Point(0,0))
    #u2 = recreate_galerkin_solution_at_x(T, U, vertices, Point(0.25,0.25))
    #u3 = recreate_galerkin_solution_at_x(T, U, vertices, Point(0.75,0.75))
    print(u0)
    #print(u1)
    #print(u2)
    #print(u3)

fig, ax = plt.subplots()
for t in T:
    ax.plot([t.p.x,t.q.x],[t.p.y,t.q.y],color='black',linestyle='-')
    ax.plot([t.q.x,t.r.x],[t.q.y,t.r.y],color='black',linestyle='-')
    ax.plot([t.r.x,t.p.x],[t.r.y,t.p.y],color='black',linestyle='-')

plt.show()
