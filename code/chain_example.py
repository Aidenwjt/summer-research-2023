#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import helpers

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Edge:
    def __init__(self, p, q):
        self.p = p
        self.q = q
    def equals(self, e):
        if((self.p == e.p and self.q == e.q) or (self.p == e.q and self.q == e.p)):
            return True
        return False

class Triangle:
    def __init__(self, v0, v1, v2, e0, e1, e2):
        self.v = [v0, v1, v2]
        self.e = [e0, e1, e2]
    def shares_edge(self, t):
        for i in range(0,2):
            for j in range(0,2):
                if(self.e[i].equals(t.e[j]) == True):
                    return True
        return False

# https://www.tutorialspoint.com/python_data_structure/python_binary_tree.htm
class Node:
    def __init__(self, T, GT):
        self.parent = None
        self.left = None
        self.right = None
        self.T = T
        self.GT = GT
        self.neighbor0 = None # NOTE: This is always F(T)
        self.neighbor1 = None
        self.neighbor2 = None
        self.root = None
    def bisect(self):
        p = helpers.midpoint(self.T.v[1], self.T.v[2])
        child_triangle1 = Triangle(p, self.T.v[0], self.T.v[1], Edge(self.T.v[0],self.T.v[1]), Edge(p,self.T.v[1]), Edge(p,self.T.v[0]))
        child_triangle2 = Triangle(p, self.T.v[2], self.T.v[0], Edge(self.T.v[0],self.T.v[2]), Edge(p,self.T.v[0]), Edge(p,self.T.v[2]))
        child_node1 = Node(child_triangle1, self.GT + 1)
        child_node2 = Node(child_triangle2, self.GT + 1)
        self.left = child_node1
        self.right = child_node2
        child_node1.parent = self
        child_node2.parent = self
        return (child_node1, child_node2)

def update_neighbors(child1, child2, parents_neighbor):
    if(parents_neighbor.left != None):
        # Iterate through the neighbors left children
        for e in parents_neighbor.left.T.e:
            # Check if left child shares an edge with any of the child triangles
            if(child1.T.e[0].equals(e) == True):
                # Update neighbor relations
                child1.neighbor0 = parents_neighbor.left
                if(parents_neighbor.left.neighbor0 == child1.parent):
                    parents_neighbor.left.neighbor0 = child1
                elif(parents_neighbor.left.neighbor1 == child1.parent):
                    parents_neighbor.left.neighbor1 = child1
                else:
                    parents_neighbor.left.neighbor2 = child1
            if(child2.T.e[0].equals(e) == True):
                child2.neighbor0 = parents_neighbor.left
                if(parents_neighbor.left.neighbor0 == child1.parent):
                    parents_neighbor.left.neighbor0 = child2
                elif(parents_neighbor.left.neighbor1 == child1.parent):
                    parents_neighbor.left.neighbor1 = child2
                else:
                    parents_neighbor.left.neighbor2 = child2
        for e in parents_neighbor.right.T.e:
            # Check if right child shares an edge with any of the child triangles
            if(child1.T.e[0].equals(e) == True):
                # Update neighbor relations
                child1.neighbor0 = parents_neighbor.right
                if(parents_neighbor.right.neighbor0 == child1.parent):
                    parents_neighbor.right.neighbor0 = child1
                elif(parents_neighbor.right.neighbor1 == child1.parent):
                    parents_neighbor.right.neighbor1 = child1
                else:
                    parents_neighbor.right.neighbor2 = child1
            if(child2.T.e[0].equals(e) == True):
                child2.neighbor0 = parents_neighbor.right
                if(parents_neighbor.right.neighbor0 == child1.parent):
                    parents_neighbor.right.neighbor0 = child2
                elif(parents_neighbor.right.neighbor1 == child1.parent):
                    parents_neighbor.right.neighbor1 = child2
                else:
                    parents_neighbor.right.neighbor2 = child2
    else:
        for e in parents_neighbor.T.e:
            if(child1.T.e[0].equals(e) == True):
                # Update neighbor relations
                child1.neighbor0 = parents_neighbor
                if(parents_neighbor.neighbor0 == child1.parent):
                    parents_neighbor.neighbor0 = child1
                elif(parents_neighbor.neighbor1 == child1.parent):
                    parents_neighbor.neighbor1 = child1
                else:
                    parents_neighbor.neighbor2 = child1
            if(child2.T.e[0].equals(e) == True):
                child2.neighbor0 = parents_neighbor
                if(parents_neighbor.neighbor0 == child1.parent):
                    parents_neighbor.neighbor0 = child2
                elif(parents_neighbor.neighbor1 == child1.parent):
                    parents_neighbor.neighbor1 = child2
                else:
                    parents_neighbor.neighbor2 = child2

def refine_recursive(mesh, elem):
    if(elem.neighbor0 != None and elem.neighbor0.GT < elem.GT):
        mesh = refine_recursive(mesh, elem.neighbor0)
    child1, child2 = elem.bisect()
    mesh.remove(elem)
    mesh.append(child1)
    mesh.append(child2)
    child1.neighbor2 = child2 # NOTE: still dont know if edge 0 and edge 1 have triangles that share them
    child2.neighbor1 = child1 # NOTE: still dont konw if edge 0 and edge 2 have triangles that share them
    if(child1.parent.neighbor1 != None):
        update_neighbors(child1, child2, child1.parent.neighbor1)
    if(child1.parent.neighbor2 != None):
        update_neighbors(child1, child2, child1.parent.neighbor2)
    if elem.neighbor0 == None:
        return mesh
    child3, child4 = elem.neighbor0.bisect()
    child3.neighbor2 = child4
    child4.neighbor1 = child3
    child1.neighbor1 = child4
    child2.neighbor2 = child3
    child3.neighbor1 = child2
    child4.neighbor2 = child1
    if(child3.parent.neighbor1 != None):
        update_neighbors(child3, child4, child3.parent.neighbor1)
    if(child3.parent.neighbor2 != None):
        update_neighbors(child3, child4, child3.parent.neighbor2)
    mesh.remove(elem.neighbor0)
    mesh.append(child3)
    mesh.append(child4)
    return mesh

def refine(mesh, marked):
    for m in marked:
        mesh = refine_recursive(mesh, m)
    return mesh

# Initial mesh in example
p0 = Point(0,0)
p1 = Point(1,1)
p2 = Point(-1,1)
p3 = Point(-1,-1)
p4 = Point(1,-1)
T1 = Triangle(p0, p4, p1, Edge(p1 , p4), Edge(p0, p1), Edge(p0, p4))
T2 = Triangle(p0, p3, p4, Edge(p3 , p4), Edge(p0, p4), Edge(p3, p0))
T3 = Triangle(p0, p2, p3, Edge(p2 , p3), Edge(p0, p3), Edge(p0, p2))
T4 = Triangle(p0, p1, p2, Edge(p1 , p2), Edge(p2, p0), Edge(p0, p1))

# Creating forest of binary trees
root1 = Node(T1, 0)
root2 = Node(T2, 0)
root3 = Node(T3, 0)
root4 = Node(T4, 0)

root1.root = root1
root1.neighbor0 = None
root1.neighbor1 = root4
root1.neighbor2 = root2

root2.root = root2
root2.neighbor0 = None
root2.neighbor1 = root1
root2.neighbor2 = root3

root3.root = root3
root3.neighbor0 = None
root3.neighbor1 = root2
root3.neighbor2 = root4

root4.root = root4
root4.neighbor0 = None
root4.neighbor1 = root3
root4.neighbor2 = root1

# Initial mesh
mesh = [root1, root2, root3, root4]

# First iteration
marked = [root1, root4]
mesh = refine(mesh, marked)

# Second iteration
marked = [root1.right]
mesh = refine(mesh, marked)

# Third iteration
marked = [root1.right.right]
mesh = refine(mesh, marked)


fig, ax = plt.subplots()
for elem in mesh:
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
