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
        if(( (self.p.x == e.p.x and self.p.y == e.p.y) and (self.q.x == e.q.x and self.q.y == e.q.y) )
           or ( (self.p.x == e.q.x and self.p.y == e.q.y) and (self.q.x == e.p.x and self.q.y == e.p.y) )):
            return True
        return False

class Triangle:
    def __init__(self, v0, v1, v2, e0, e1, e2):
        self.v = [v0, v1, v2]
        self.e = [e0, e1, e2]
    def shared_edges(self, T):
        for i in range(0,3):
            for j in range(0,3):
                if(self.e[i].equals(T.e[j]) == True):
                    return (i,j)
        return None

class Node:
    def __init__(self, T, GT):
        self.parent = None
        self.left = None
        self.right = None
        self.T = T
        self.GT = GT
        self.neighbor0 = None # NOTE: F(T) for now
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
        child_node1.neighbor2 = child_node2
        child_node2.neighbor1 = child_node1
        return (child_node1, child_node2)
    def update_neighbor(self, elem):
        edges = self.T.shared_edges(elem.T)
        if(edges != None):
            i,j = edges
            if(i == 0):
                self.neighbor0 = elem
            elif(i == 1):
                self.neighbor1 = elem
            else:
                self.neighbor2 = elem
            if(j == 0):
                elem.neighbor0 = self
            elif(j == 1):
                elem.neighbor1 = self
            else:
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

def refine_recursive(mesh, elem):
    if(elem.neighbor0 != None and elem.neighbor0.GT < elem.GT):
        mesh = refine_recursive(mesh, elem.neighbor0)
    child1, child2 = elem.bisect()
    mesh.remove(elem)
    mesh.append(child1)
    mesh.append(child2)
    child1.update_neighbors()
    child2.update_neighbors()
    if elem.neighbor0 == None:
        return mesh
    child3, child4 = elem.neighbor0.bisect()
    child3.update_neighbors()
    child4.update_neighbors()
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

# Fourth iteration
marked = [root1.right.right.left]
mesh = refine(mesh, marked)

#marked = [root1.right.right.left.right]
#mesh = refine(mesh, marked)

#marked = [root1.right.right.left.right.left]
#mesh = refine(mesh, marked)

#marked = [root1.right.right.left.right.left.left]
#mesh = refine(mesh, marked)

fig, ax = plt.subplots()
for elem in mesh:
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
