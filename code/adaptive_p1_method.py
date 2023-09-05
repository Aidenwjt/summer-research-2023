#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from shapely import geometry
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# TODO:
#   - Put functions where they should be
#   - Better variable and function naming
#   - Maybe have a function to accept points for a simple polygon
#   - Comment everything

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
    def diam_of_T(self):
        return np.maximum(np.maximum(self.p.distance(self.q),self.q.distance(self.r)),self.r.distance(self.p))
    def det_X_T(self):
        return np.linalg.det(np.array([[1, self.p.x, self.p.y], [1, self.q.x, self.q.y], [1, self.r.x, self.r.y]]))
    def vol_T(self):
        return self.det_X_T()/2
    def grad_phi(self, j):
        points = [self.p, self.q, self.r]
        return [(points[(j+1)%3].y - points[(j+2)%3].y)/(self.det_X_T()),
                (points[(j+2)%3].x - points[(j+1)%3].x)/(self.det_X_T())]
    def contained_in_T(self, x):
        a = ( ((self.q.y - self.r.y)*(x.x - self.r.x) + (self.r.x - self.q.x)*(x.y - self.r.y)) / 
             ((self.q.y - self.r.y)*(self.p.x - self.r.x) + (self.r.x - self.q.x)*(self.p.y - self.r.y)) )
        b = ( ((self.r.y - self.p.y)*(x.x - self.r.x) + (self.p.x - self.r.x)*(x.y - self.r.y)) / 
             ((self.q.y - self.r.y)*(self.p.x - self.r.x) + (self.r.x - self.q.x)*(self.p.y - self.r.y)))
        c = 1 - a - b
        if((a >= 0 and a <= 1) and (b >= 0 and b <= 1) and (c >= 0 and c <= 1)):
            return True
        return False
    def vertex_of_T(self, s):
        points = [self.p, self.q, self.r]
        for i in range(0, 3):
            if(points[i].equals(s)):
                return i
        return False
    def A(self, i, j):
        grad_phi_i = self.grad_phi(i)
        grad_phi_j = self.grad_phi(j)
        dot_grads = (grad_phi_i[0] * grad_phi_j[0]) + (grad_phi_i[1] * grad_phi_j[1])
        return dot_grads * self.vol_T()
    def shared_edges(self, other_T):
        v0 = [self.p, self.q, self.r]
        v1 = [other_T.p, other_T.q, other_T.r]
        for i in range(0,3):
            for j in range(0,3):
                if(( (v0[(i+1)%3].x == v1[(j+1)%3].x and v0[(i+1)%3].y == v1[(j+1)%3].y) and (v0[(i+2)%3].x == v1[(j+2)%3].x and v0[(i+2)%3].y == v1[(j+2)%3].y) )
                   or ( (v0[(i+1)%3].x == v1[(j+2)%3].x and v0[(i+1)%3].y == v1[(j+2)%3].y) and (v0[(i+2)%3].x == v1[(j+1)%3].x and v0[(i+2)%3].y == v1[(j+1)%3].y) )):
                    return (i, j)
        return None

class Node:
    def __init__(self, T, P, GT):
        self.parent = None
        self.left = None
        self.right = None
        self.T = T
        self.GT = GT
        self.P = P # NOTE: local connectivity vector
        self.ET = None
        self.eta = None
        self.neighbor0 = None
        self.neighbor1 = None
        self.neighbor2 = None
        self.marked = False
    def find_refinement_edge(self):
        # TODO: This function will not be necessary later
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
    def bisect(self, mesh, vertices, boundary_vertices, boundary):
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
                boundary_vertices.append(index)
        # Construct the new triangles
        child_triangle1 = Triangle(v_elem[self.ET], v_elem[(self.ET+1)%3], mp)
        child_triangle2 = Triangle(v_elem[self.ET], mp, v_elem[(self.ET+2)%3])
        # Construct the new nodes
        child1 = Node(child_triangle1, [v_indices[self.ET], v_indices[(self.ET+1)%3], index], self.GT + 1)
        child2 = Node(child_triangle2, [v_indices[self.ET], index, v_indices[(self.ET+2)%3]], self.GT + 1)
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
        mesh.remove(self)
        mesh.append(child1)
        mesh.append(child2)
        return

def refine_recursive(mesh, elem, vertices, boundary_vertices, boundary):
    # Create a list of the elements neighbors
    neighbors = [elem.neighbor0, elem.neighbor1, elem.neighbor2]
    # Create a reference to the neighbors sharing its refinement edge
    FT = neighbors[elem.ET]
    # If FT is not None and has an older generation, then it needs to be refined first
    if(FT != None and FT.GT < elem.GT):
        refine_recursive(mesh, FT, vertices, boundary_vertices, boundary)
    neighbors = [elem.neighbor0, elem.neighbor1, elem.neighbor2]
    FT = neighbors[elem.ET]
    # Else, we have a compatible refinement patch, so bisect the current element
    elem.bisect(mesh, vertices, boundary_vertices, boundary)
    # If we just returned from a recursive call, and FT is still None, then just return from this function call
    #if FT == None or FT.left != None:
    if FT == None:
        elem.left.update_neighbors()
        elem.right.update_neighbors()
        return
    # Else, we also need to to bisect FT
    FT.bisect(mesh, vertices, boundary_vertices, boundary)
    return

def refine(mesh, vertices, boundary_vertices, boundary):
    # Create a shallow copy of T as we will be removing elements from it
    mesh_copy = mesh.copy()
    # Loop through all the elements
    for elem in mesh_copy:
        # If an element is marked for refinement, and it has not already been refined, then refine it
        if elem.marked == True and elem.left == None and elem.right == None:
            refine_recursive(mesh, elem, vertices, boundary_vertices, boundary)

def phi_of_x(mesh, v, x):
    for elem in mesh:
        j = elem.T.vertex_of_T(v)
        if(j is not False):
            if(elem.T.contained_in_T(x)):
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

def main():
    # Define the scalar function f, which will just be constant in our case
    f = 4
    # TODO: I am going to make a functionality for the initialization of some polygonal domain, given a list of vertices
    # ==================================================
    # Define the vertices of the simple polygon that will make up the domain
    p0 = Point(0,0)
    p1 = Point(1,0)
    p2 = Point(1,1)
    p3 = Point(0,1)
    p4 = Point(-1,1)
    p5 = Point(-1,0)
    p6 = Point(-1,-1)
    p7 = Point(0,-1)

    # Construct a boundary object that we can use to check if a point is on the boundary or not
    domain = geometry.Polygon([(0,0),(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1)])
    boundary = domain.boundary

    # Define a list of vertices and a list of indices to keep track of which vertices are on the boundary
    vertices = [p0, p1, p2, p3, p4, p5, p6, p7]
    boundary_vertices = [0,1,2,3,4,5,6,7] # NOTE: all initial vertices are on the boundary by default

    # Define the triangles to keep track of which elements contain which vertices
    t0 = Triangle(p0, p1, p2)
    t1 = Triangle(p0, p2, p3)
    t2 = Triangle(p0, p3, p5)
    t3 = Triangle(p3, p4, p5)
    t4 = Triangle(p5, p6, p0)
    t5 = Triangle(p6, p7, p0)

    # Define the elements themselves to keep track of connectivity to the vertices list, neighbors, and element generation
    root0 = Node(t0, [0, 1, 2], 0)
    root1 = Node(t1, [0, 2, 3], 0)
    root2 = Node(t2, [0, 3, 5], 0)
    root3 = Node(t3, [3, 4, 5], 0)
    root4 = Node(t4, [5, 6, 0], 0)
    root5 = Node(t5, [6, 7, 0], 0)

    # Define a list of all the elements so they can be iterates through
    mesh = [root0, root1, root2, root3, root4, root5]

    # Update the elements neighbors
    for i in range(0, len(mesh)):
        mesh[i].find_refinement_edge()
        for j in range(0, len(mesh)):
            if(i != j):
                mesh[i].update_neighbor(mesh[j])
    # ==================================================

    # Define parameters for Solve-Estimate->Mark->Refine algorithm
    theta = 0.5 # NOTE: theta is the parameter that dictates how many elements will be refined, where theta is in (0,1] 
    eps_stop = 0.8 # NOTE: eps_stop is the epsilon stopping/thresholding value used to halt the program at some desired error, where eps_stop > 0
    error_bound = 1 # NOTE: error_bound is the variable for computing the error estimate, we start at 1, but we could really start at any value greater than eps_stop

    # Main loop that implements the Solve->Estimate->Mark->Refine algorithm
    while(error_bound > eps_stop):
        # Construct the stiffness matrix and the nodal load vector F in the system AU = F, where U is the vector of Galerkin basis coefficients
        A = [[0]*len(vertices) for i in range(0, len(vertices))]
        F = [0]*len(vertices)
        # Construct the preliminary stiffness matrix by adding up all of the local contributions
        for elem in mesh:
            for j in range(0, 3):
                for i in range(0, 3):
                    A[elem.P[i]][elem.P[j]] = A[elem.P[i]][elem.P[j]] + elem.T.A(i, j)
                F[elem.P[j]] = F[elem.P[j]] + f*(elem.T.vol_T()/3)
        # Now enforce the purely Dirichlet boundary condition at the vertices on the boundary
        for i in range(0, len(boundary_vertices)):
            for j in range(0, len(vertices)):
                if(j == boundary_vertices[i]):
                    A[j][boundary_vertices[i]] = 1
                else:
                    A[j][boundary_vertices[i]] = 0
                    A[boundary_vertices[i]][j] = 0
            F[boundary_vertices[i]] = 0
        # Now with the stiffness matrix and nodal load vector F, we can solve for the Galerkin coefficients
        U = np.linalg.solve(np.array(A), np.array(F)).tolist()
        # With the Galerkin coefficients, we now compute the a posteriori error estimator, keeping track of the maximum so we can mark elements later
        max_eta = 0
        etas = []
        for elem in mesh:
            # 
            left_summand = (f**2)*(elem.T.vol_T()**2)
            #left_summand = (f)*(elem.T.vol_T()**2)
            grads_T = np.linalg.solve(np.array([[1,1,1],[elem.T.p.x,elem.T.q.x,elem.T.r.x],[elem.T.p.y,elem.T.q.y,elem.T.r.y]]), np.array([[0,0],[1,0],[0,1]]))
            nabla_u_T = np.matmul(np.transpose(grads_T), np.array([U[elem.P[0]], U[elem.P[1]], U[elem.P[2]]]))
            normal_times_area = np.multiply(-2*elem.T.vol_T(), grads_T)
            jump_vector = np.matmul(normal_times_area, nabla_u_T)
            right_summand = 0
            v_elem = [elem.T.p, elem.T.q, elem.T.r]
            for i in range(0, len(jump_vector)):
                if(boundary.contains(geometry.Point(v_elem[(i+1)%3].x,v_elem[(i+1)%3].y)) == False 
                   or boundary.contains(geometry.Point(v_elem[(i+2)%3].x,v_elem[(i+2)%3].y)) == False):
                    right_summand += np.absolute(v_elem[(i+1)%3].distance(v_elem[(i+2)%3])) * np.absolute(jump_vector[i])**2
            elem.eta = np.sqrt(left_summand + right_summand)
            etas.append(elem.eta)
            if(elem.eta > max_eta):
                max_eta = elem.eta
        error_bound = 0
        for eta in etas:
            error_bound += eta**2
        error_bound = np.sqrt(error_bound)
        print(error_bound) # NOTE: keep this to check if the error bound is actually going down
        for elem in mesh:
            if(elem.eta >= theta*max_eta):
                elem.marked = True
        if(error_bound > eps_stop):
            refine(mesh, vertices, boundary_vertices, boundary)

    fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    triangles1 = []
    triangles2 = []
    for elem in mesh:
        triangles1.append((
            (elem.T.p.x, elem.T.p.y, 0), 
            (elem.T.q.x, elem.T.q.y, 0), 
            (elem.T.r.x, elem.T.r.y, 0) 
        ))
        triangles2.append((
            (elem.T.p.x, elem.T.p.y, U[elem.P[0]] * phi_of_x(mesh, vertices[elem.P[0]], vertices[elem.P[0]])), 
            (elem.T.q.x, elem.T.q.y, U[elem.P[1]] * phi_of_x(mesh, vertices[elem.P[1]], vertices[elem.P[1]])), 
            (elem.T.r.x, elem.T.r.y, U[elem.P[2]] * phi_of_x(mesh, vertices[elem.P[2]], vertices[elem.P[2]])) 
        ))

    ax[0].add_collection(Poly3DCollection(triangles1, edgecolor='black', facecolor='white'))
    ax[0].set_title("2D Mesh")
    ax[0].set_xlim(-2, 2)
    ax[0].set_ylim(-2, 2)

    ax[1].add_collection(Poly3DCollection(triangles2, edgecolor='black', facecolor='white'))
    ax[1].set_title("3D Mesh")
    ax[1].set_xlim(-2, 2)
    ax[1].set_ylim(-2, 2)

    plt.show()

main()
