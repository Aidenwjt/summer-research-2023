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
        return dot_grads * self.vol_t()
        # NOTE: below seems to be an equivalent implementation
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
    def shared_edges(self, T):
        v0 = [self.p, self.q, self.r]
        v1 = [T.p, T.q, T.r]
        for i in range(0,3):
            for j in range(0,3):
                if(( (v0[(i+1)%3].x == v1[(j+1)%3].x and v0[(i+1)%3].y == v1[(j+1)%3].y) and (v0[(i+2)%3].x == v1[(j+2)%3].x and v0[(i+2)%3].y == v1[(j+2)%3].y) )
                   or ( (v0[(i+1)%3].x == v1[(j+2)%3].x and v0[(i+1)%3].y == v1[(j+2)%3].y) and (v0[(i+2)%3].x == v1[(j+1)%3].x and v0[(i+2)%3].y == v1[(j+1)%3].y) )):
                    return (i, j)
        return None

class Node:
    def __init__(self, T, i, j, k, GT):
        self.parent = None
        self.left = None
        self.right = None
        self.T = T
        self.GT = GT
        self.P = [i, j, k] # NOTE: this keeps track of the indices of the elements vertices in the vertices list
        self.ET = None
        self.neighbor0 = None
        self.neighbor1 = None
        self.neighbor2 = None
        self.marked = False
    def find_refinement_edge(self):
        # NOTE: ET = 0 corresponds to the edge opposite p, ET = 1 edge opposite q, ET = 2 edge opposite r.
        if(np.maximum(self.T.q.distance(self.T.r), np.maximum(self.T.r.distance(self.T.p), self.T.p.distance(self.T.q))) == self.T.q.distance(self.T.r)):
            self.ET = 0
        if(np.maximum(self.T.q.distance(self.T.r), np.maximum(self.T.r.distance(self.T.p), self.T.p.distance(self.T.q))) == self.T.r.distance(self.T.p)):
            self.ET = 1
        if(np.maximum(self.T.q.distance(self.T.r), np.maximum(self.T.r.distance(self.T.p), self.T.p.distance(self.T.q))) == self.T.p.distance(self.T.q)):
            self.ET = 2
    def update_neighbor(self, elem):
        edges = self.T.shared_edges(elem.T)
        if(edges != None):
            i,j = edges
            if(i == 0):
                self.neighbor0 = elem
            if(i == 1):
                self.neighbor1 = elem
            if(i == 2):
                self.neighbor2 = elem
            if(j == 0):
                elem.neighbor0 = self
            if(j == 1):
                elem.neighbor1 = self
            if(j == 2):
                elem.neighbor2 = self
    def update_neighbors(self):
        if self.parent != None:
            parents_neighbors = [self.parent.neighbor0, self.parent.neighbor1, self.parent.neighbor2]
            for neighbor in parents_neighbors:
                if neighbor != None:
                    if neighbor.left != None:
                        children = [neighbor.left, neighbor.right]
                        for child in children:
                            self.update_neighbor(child)
                    else:
                        self.update_neighbor(neighbor)
    def green_bisect(self, T, vertices, bv, boundary):
        # Define a list of vertices of the element
        v_elem = [self.T.p, self.T.q, self.T.r]
        v_indices = [self.P[0], self.P[1], self.P[2]]
        # Define the new point, which is the midpoint on the refinement edge
        mp = v_elem[(self.ET+1)%3].midpoint(v_elem[(self.ET+2)%3])
        # Check if this point is already in vertices, and if it is on the boundary
        index = -1
        for i in range(0, len(vertices)):
            if(mp.equals(vertices[i])):
                index = i
                break
        if(index == -1):
            vertices.append(mp)
            index = len(vertices) - 1
            if(boundary.contains(geometry.Point(mp.x, mp.y))):
                bv.append(index)
        # Construct the new triangles
        child_triangle1 = Triangle(v_elem[self.ET], v_elem[(self.ET+1)%3], mp)
        child_triangle2 = Triangle(v_elem[self.ET], mp, v_elem[(self.ET+2)%3])
        # Construct the new nodes
        child1 = Node(child_triangle1, v_indices[self.ET], v_indices[(self.ET+1)%3], index, self.GT + 1)
        child2 = Node(child_triangle2, v_indices[self.ET], index, v_indices[(self.ET+2)%3], self.GT + 1)
        # Update their parents, neighbors, and refinement edges
        child1.parent = self
        child2.parent = self
        self.left = child1
        self.right = child2
        child1.update_neighbor(child2)
        child1.update_neighbors()
        child2.update_neighbors()
        child1.find_refinement_edge()
        child2.find_refinement_edge()
        # Remove the element from the triangulation and add the two new elements
        T.remove(self)
        T.append(child1)
        T.append(child2)
        return

def refine_recursive(T, elem, vertices, bv, boundary):
    # Create a list of the elements neighbors
    neighbors = [elem.neighbor0, elem.neighbor1, elem.neighbor2]
    # Create a reference to the neighbors sharing its refinement edge
    FT = neighbors[elem.ET]
    # If FT is not None and has an older generation, then it needs to be refined first
    if(FT != None and FT.GT < elem.GT):
        refine_recursive(T, FT, vertices, bv, boundary)
    # Else, we have a compatible refinement patch, so bisect the current element
    elem.green_bisect(T, vertices, bv, boundary)
    # If we just returned from a recursive call, and FT is still None, then just return from this function call
    if FT == None:
        return
    # Else, we also need to to bisect FT
    FT.green_bisect(T, vertices, bv, boundary)
    return

def refine(T, vertices, bv, boundary):
    # Create a shallow copy of T as we will be removing elements from it
    T_copy = T.copy()
    # Loop through all the elements
    #for k in range(0, len(T_copy)):
    for elem in T_copy:
        # If an element is marked for refinement, and it has not already been refined, then refine it
        if elem.marked == True and elem.left == None and elem.right == None:
            refine_recursive(T, elem, vertices, bv, boundary)

def phi_of_x(T, v, x):
    for elem in T:
        j = elem.T.vertex_of_t(v)
        if(j > -1):
            if(elem.T.x_in_closure_t(x)):
                points = [elem.T.p, elem.T.q, elem.T.r]
                det1 = np.linalg.det(np.array([
                    [1, x.x, x.y],
                    [1, points[(j+1)%3].x, points[(j+1)%3].y],
                    [1, points[(j+2)%3].x, points[(j+2)%3].y]]))
                det2 = np.linalg.det(np.array([
                    [1, points[j].x, points[j].y],
                    [1, points[(j+1)%3].x, points[(j+1)%3].y],
                    [1, points[(j+2)%3].x, points[(j+2)%3].y]]))
                print(det1/det2)
                return det1/det2
    return 0

# Define the scalar function f, which will just be constant in our case
scalar_f = 4

# Define the points/vertices of the polygonal domain
p0 = Point(0,0)
p1 = Point(1,0)
p2 = Point(1,1)
p3 = Point(0,1)

# Construct a boundary object to check if points are on the boundary in the future
domain = geometry.Polygon([(0,0),(1,0),(1,1),(0,1)])
boundary = domain.boundary

# Define the triangles in the initial triangulation
t0 = Triangle(p0,p1,p2)
t1 = Triangle(p0,p2,p3)

# Create a list that keeps track of all the individual vertices in the triangulation
vertices = [p0, p1, p2, p3]

# Create a conectivity matrix which will be necessary for Galerkin approximation
#P = [[0, 1, 2], [0, 2, 3]]

# Create a list to keep track of the indices in the vertices list that are on the boundary
bv = [0, 1, 2, 3]

# Use the initial triangles to construct roots of binary tree to keep track of neighbors, generation, etc...
root0 = Node(t0, 0, 1, 2, 0)
root1 = Node(t1, 0, 2, 3, 0)
T = [root0, root1]
for i in range(0, len(T)):
    T[i].marked = True
    T[i].find_refinement_edge()
    for j in range(0, len(T)):
        if(i != j):
            T[i].update_neighbor(T[j])

# Main loop that implements the Solve->Estimate->Mark->Refine algorithm
for l in range(0,3):
    # Our current triangulation requires us to do an initial refinement, which we will do here
    refine(T, vertices, bv, boundary) # NOTE: python will update the lists passed to a function, so no need for a return
    A = [[0]*len(vertices) for i in range(0, len(vertices))]
    f = [0]*len(vertices)
    #for k in range(0, len(T)):
    for elem in T:
        for j in range(0, 3):
            for i in range(0, 3):
                A[elem.P[i]][elem.P[j]] += elem.T.A(i, j)
            f[elem.P[j]] += ((scalar_f * elem.T.vol_t())/3)
    print("---")
    for i in range(0, len(A)):
        print(A[i])
    print("---")
    for i in range(0, len(bv)):
        for j in range(0, len(vertices)):
            if(j == bv[i]):
                A[j][bv[i]] = 1
            else:
                A[j][bv[i]] = 0
                A[bv[i]][j] = 0
        f[bv[i]] = 0
    U = np.linalg.solve(np.array(A), np.array(f)).tolist()
    for v in vertices:
        print("({},{})".format(v.x, v.y))
    print("---")
    for i in range(0, len(A)):
        print(A[i])
    print("---")
    print(U)
    print("---")
    print(f)
    print("---")
    galerkin_solution = 0
    for i in range(0, len(vertices)):
        if(U[i] != 0):
            galerkin_solution += U[i] * phi_of_x(T, vertices[i], Point(0.5, 0.5))
    print(galerkin_solution)
    print("---")
    for i in range(0, len(T)):
        T[i].marked = True

fig, ax = plt.subplots()
for elem in T:
    ax.plot([elem.T.p.x,elem.T.q.x],[elem.T.p.y,elem.T.q.y],color='black',linestyle='-')
    ax.plot([elem.T.q.x,elem.T.r.x],[elem.T.q.y,elem.T.r.y],color='black',linestyle='-')
    ax.plot([elem.T.r.x,elem.T.p.x],[elem.T.r.y,elem.T.p.y],color='black',linestyle='-')

plt.show()
