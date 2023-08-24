import numpy as np
import matplotlib.pyplot as plt
from shapely import geometry
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm

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
    def diam_of_triangle(self):
        return np.maximum(np.maximum(self.p.distance(self.q),self.q.distance(self.r)),self.r.distance(self.p))
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
        # NOTE: equivalent implementations below
        """
        G = np.matmul(np.linalg.inv(np.array([[1,1,1],[self.p.x, self.q.x, self.r.x],[self.p.y, self.q.y, self.r.y]])), np.array([[0,0],[1,0],[0,1]]))
        GG_tran = np.matmul(G, np.transpose(G))
        for row in GG_tran:
            row[0] *= self.det_X_t()/2
            row[1] *= self.det_X_t()/2
            row[2] *= self.det_X_t()/2
        print("---")
        for row in GG_tran:
            print(row)
        print("---")
        return GG_tran[i][j]
        """
        """
        points = [self.p, self.q, self.r]
        return ((points[(i+1)%3].y - points[(i+2)%3].y)*(points[(j+1)%3].y - points[(j+2)%3].y) + (points[(i+2)%3].x - points[(i+1)%3].x)*(points[(j+2)%3].x - points[(j+1)%3].x))/(4*self.vol_t())
        """
        """
        grad_phi_i = self.grad_phi(i)
        grad_phi_j = self.grad_phi(j)
        dot_grads = (grad_phi_i[0] * grad_phi_j[0]) + (grad_phi_i[1] * grad_phi_j[1])
        return dot_grads*self.vol_t()
        """
        b0 = self.q.y - self.r.y
        b1 = self.r.y - self.p.y
        b2 = self.p.y - self.q.y
        c0 = self.r.x - self.q.x
        c1 = self.p.x - self.r.x
        c2 = self.q.x - self.p.x
        b = [b0, b1, b2]
        c = [c0, c1, c2]
        return (b[i]*b[j] + c[i]*c[j])/(4*self.vol_t())
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
        self.eta = None
        self.neighbor0 = None
        self.neighbor1 = None
        self.neighbor2 = None
        self.marked = False
    def jump(self, boundary, vertices, U):
        # Variable to keep track of the overall sum
        jump_sum = 0
        # Vertices of current triangle
        v_elem = [self.T.p, self.T.q, self.T.r]
        # Loop through all three edges/sides
        for i in range(0, 3):
            # Check if the edge is on the boundary of the domain
            if(boundary.contains(geometry.Point(v_elem[(i+1)%3].x,v_elem[(i+1)%3].y)) == False 
               or boundary.contains(geometry.Point(v_elem[(i+2)%3].x,v_elem[(i+2)%3].y)) == False):
                # Compute the squared volume of the edge
                vol_S_squared = v_elem[(i+1)%3].distance(v_elem[(i+2)%3])**2
                # Find the neighbor that shares this edge, should exists if the edge is not on the boundary
                neighbors = [self.neighbor0, self.neighbor1, self.neighbor2]
                neighbor = neighbors[i]
                neighbor_v_elem = [neighbor.T.p, neighbor.T.q, neighbor.T.r]
                # Find the index of the shared edge for the neighbor
                _,j = self.T.shared_edges(neighbor.T)
                # Compute the outer unit normals on the side to the neighbor, then to the current triangle
                dx1 = v_elem[(i+2)%3].x - v_elem[(i+1)%3].x
                dy1 = v_elem[(i+2)%3].y - v_elem[(i+1)%3].y
                length1 = np.sqrt((-dy1)**2 + dx1**2)
                n1 = [-dy1/length1, dx1/length1]
                dx2 = neighbor_v_elem[(j+2)%3].x - neighbor_v_elem[(j+1)%3].x
                dy2 = neighbor_v_elem[(j+2)%3].y - neighbor_v_elem[(j+1)%3].y
                length2 = np.sqrt((-dy2)**2 + dx2**2)
                n2 = [-dy2/length2, dx2/length2]
                # Compute the gradients of the Galerkin solutions evaluted on each triangle, which are simplified in the nodal basis
                grad_u_self = [0, 0]
                grad_u_neighbor = [0, 0]
                for i in range(0, len(vertices)):
                    k = self.T.vertex_of_t(vertices[i])
                    l = neighbor.T.vertex_of_t(vertices[i])
                    if(k > -1):
                        grad_phi_k = self.T.grad_phi(k)
                        grad_u_self[0] += U[i]*grad_phi_k[0]
                        grad_u_self[1] += U[i]*grad_phi_k[1]
                    if(l > -1):
                        grad_phi_l = self.T.grad_phi(l)
                        grad_u_neighbor[0] += U[i]*grad_phi_l[0]
                        grad_u_neighbor[1] += U[i]*grad_phi_l[1]
                # Compute the current summand
                jump_sum += vol_S_squared*( np.absolute( (grad_u_self[0] * n1[0] + grad_u_self[1] * n1[1]) + (grad_u_neighbor[0] * n2[0] + grad_u_neighbor[1] * n2[1]) )**2 )
        return jump_sum
    # TODO: This function probabaly is no longer necessary
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
                    if neighbor.left != None and neighbor.right != None:
                        children = [neighbor.left, neighbor.right]
                        for child in children:
                            self.update_neighbor(child)
                    else:
                        self.update_neighbor(neighbor)
    def bisect(self, T, vertices, bv, boundary):
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
        child1.ET = 2
        child2.ET = 1
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
    neighbors = [elem.neighbor0, elem.neighbor1, elem.neighbor2]
    FT = neighbors[elem.ET]
    # Else, we have a compatible refinement patch, so bisect the current element
    elem.bisect(T, vertices, bv, boundary)
    # If we just returned from a recursive call, and FT is still None, then just return from this function call
    #if FT == None or FT.left != None:
    if FT == None:
        elem.left.update_neighbors()
        elem.right.update_neighbors()
        return
    # Else, we also need to to bisect FT
    FT.bisect(T, vertices, bv, boundary)
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
                return det1/det2
    return 0

# Define the scalar function f, which will just be constant in our case
scalar_f = 1

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

# Our current triangulation requires us to do an initial refinement, which we will do here
refine(T, vertices, bv, boundary) # NOTE: python will update the lists passed to a function, so no need for a return

# Define parameters for Solve-Estimate->Mark->Refine algorithm
theta = 0.5 # NOTE: theta in (0,1]
eps_stop = 0.05 # NOTE: eps_stop > 0
error_bound = 1 # NOTE: variable for computing the error estimate

test = 0
# Main loop that implements the Solve->Estimate->Mark->Refine algorithm
while(error_bound > eps_stop):
    # Construct the stiffness matrix and f vector
    A = [[0]*len(vertices) for i in range(0, len(vertices))]
    f = [0]*len(vertices)
    # Construct the preliminary stiffness matrix my adding up all of the local contributions
    #for v in vertices:
    #    print("({},{})".format(v.x,v.y))
    for elem in T:
        #print("---")
        #print("(({},{}),({},{}),({},{})".format(elem.T.p.x, elem.T.p.y, elem.T.q.x, elem.T.q.y, elem.T.r.x, elem.T.r.y))
        #G = np.matmul(np.linalg.inv(np.array([[1,1,1],[elem.T.p.x, elem.T.q.x, elem.T.r.x],[elem.T.p.y, elem.T.q.y, elem.T.r.y]])), np.array([[0,0],[1,0],[0,1]]))
        #GG_tran = np.matmul(G, np.transpose(G))
        #for row in GG_tran:
        #    row[0] *= elem.T.det_X_t()/2
        #    row[1] *= elem.T.det_X_t()/2
        #    row[2] *= elem.T.det_X_t()/2
        #for row in GG_tran:
        #    print(row)
        #print("---")
        for j in range(0, 3):
            for i in range(0, 3):
                A[elem.P[i]][elem.P[j]] = A[elem.P[i]][elem.P[j]] + elem.T.A(i, j)
            f[elem.P[j]] = f[elem.P[j]] + scalar_f*(elem.T.vol_t()/3)
    #print("---")
    ##for i in range(0, len(A)):
    #    print(A[i])
    #print("---")
    # After the preliminary matrix has been formed, enforce the boundary conditions
    for i in range(0, len(bv)):
        for j in range(0, len(vertices)):
            if(j == bv[i]):
                A[j][bv[i]] = 1
            else:
                A[j][bv[i]] = 0
                A[bv[i]][j] = 0
        f[bv[i]] = 0
    eigenvalues,_ = np.linalg.eig(A)
    lambda_min = np.amin(eigenvalues)
    lambda_max = np.amax(eigenvalues)
    print(lambda_max/lambda_min)
    #print("---")
    #for i in range(0, len(A)):
    #    print(A[i])
    #print("---")
    #print(f)
    # Now with the stiffness matrix and f vector, we can solve for the Galerkin coefficients
    U = np.linalg.solve(np.array(A), np.array(f)).tolist()
    #for v in vertices:
    #    print("({},{})".format(v.x, v.y))
    #print("---")
    #for i in range(0, len(A)):
    #    print(A[i])
    #print("---")
    #print(U)
    #print("---")
    #print(f)
    #print("---")
    # With the Galerkin coefficients, we can reconstruct the Galerkin solution at some point in our domain
    galerkin_solution = 0
    for i in range(0, len(vertices)):
        if(U[i] != 0):
            galerkin_solution += U[i] * phi_of_x(T, vertices[i], Point(0.5, 0.5))
    #print(galerkin_solution)
    #print("---")
    # With the Galerkin coefficients, we now compute the a posteriori error estimator, keeping track of the maximum
    max_eta = 0
    etas = []
    for elem in T:
        #elem.eta = np.sqrt((scalar_f**2)*(elem.T.vol_t()**2)) + elem.jump(boundary, vertices, U))
        #elem.eta = np.sqrt((scalar_f**2)*(elem.T.vol_t()*(elem.T.diam_of_triangle()**2)) + elem.jump(boundary, vertices, U))
        elem.eta = np.sqrt((scalar_f**2)*((elem.T.vol_t()**2)*(elem.T.diam_of_triangle()**2)) + elem.jump(boundary, vertices, U))
        etas.append(elem.eta)
        if(elem.eta > max_eta):
            max_eta = elem.eta
    error_bound = 0
    for eta in etas:
        error_bound += eta**2
    error_bound = np.sqrt(error_bound)
    #print(error_bound)
    for elem in T:
        if(elem.eta >= theta*max_eta):
            elem.marked = True
    test +=1
    #if test == 1:
    #    break
    if(error_bound > eps_stop):
        refine(T, vertices, bv, boundary)

fig, ax = plt.subplots(1, 3, subplot_kw={"projection": "3d"})
triangles = []
for elem in T:
    triangles.append((
        (elem.T.p.x, elem.T.p.y, 0), 
        (elem.T.q.x, elem.T.q.y, 0), 
        (elem.T.r.x, elem.T.r.y, 0) 
    ))
ax[0].add_collection(Poly3DCollection(triangles, edgecolor='black', facecolor='white'))
ax[0].set_zlim(0)
ax[0].set_title("2D Mesh")

triangles = []
for elem in T:
    triangles.append((
        (elem.T.p.x, elem.T.p.y, U[elem.P[0]] * phi_of_x(T, vertices[elem.P[0]], vertices[elem.P[0]])), 
        (elem.T.q.x, elem.T.q.y, U[elem.P[1]] * phi_of_x(T, vertices[elem.P[1]], vertices[elem.P[1]])), 
        (elem.T.r.x, elem.T.r.y, U[elem.P[2]] * phi_of_x(T, vertices[elem.P[2]], vertices[elem.P[2]])) 
    ))

ax[1].add_collection(Poly3DCollection(triangles, edgecolor='black', facecolor='white'))
ax[1].set_title("3D Mesh")

X = np.linspace(0,1)
Y = np.linspace(0,1)
X, Y = np.meshgrid(X, Y)
Z = X*(1-X) + Y*(1-Y)
surf = ax[2].plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax[2].set_title("Exact Solution")
ax[2].set_zlim(0,1)

#fig.colorbar(surf, shrink=0.2, aspect=5)


plt.show()
