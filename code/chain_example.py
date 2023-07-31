#!/usr/bin/env python3

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
        child_node1 = Node(child_triangle1, self.T.GT + 1)
        child_node2 = Node(child_triangle2, self.T.GT + 1)
        self.T.left = child_node1
        self.T.right = child_node2
        child_node1.parent = self.T
        child_node2.parent = self.T
        return (child_node1, child_node2)

def refine_recursive(mesh, elem):
    if (elem.neighbor0 != None):
        if(elem.neighbor0.GT < elem.GT):
            mesh = refine_recursive(mesh, elem.neighbor0)
    else:
        child1, child2 = elem.bisect()
        mesh.remove(elem)
        mesh.append(child1,child2)
        # TODO: Set new neighbors of child1 and child2
        #       - go through all three neighbors of parent node
        #           - if the neighbors have children, then go through the children instead
        #       - NEED TO ALSO UPDATE RELATION OF THE NEIGHBORING TRIANGLES TO THE NEW ONES AS WELL?
        child1.neighbor2 = child2 # NOTE: still dont know if edge 0 and edge 1 have triangles that share them
        child2.neighbor1 = child1 # NOTE: still dont konw if edge 0 and edge 2 have triangles that share them
        if(child1.parent.neighbor0 != None):
            if(child1.parent.neighbor0.left != None):
        if(child1.parent.neighbor1 != None):
            if(child1.parent.neighbor1.left != None):
        if(child1.parent.neighbor2 != None):
            if(child1.parent.neighbor2.left != None):
        if elem.neighbor0 == None:
            return mesh
        child3, child4 = neighbor2.bisect()
        child3.neighbor2 = child4
        child4.neighbor1 = child3
        # TODO: Set new neighbors of child3 child4
        #       - go through all three neighbors of parent node
        #           - if the neighbors have children, then go through the children instead
        mesh.remove(elem)
        mesh.remove(elem.neighbor0)
        mesh.append(child1,child2,child3,child4)
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
root1 = Node(T1, 0, None)
root2 = Node(T2, 0, None)
root3 = Node(T3, 0, None)
root4 = Node(T4, 0, None)

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

# TODO: Refine to get to same place as example
mesh = [root1, root2, root3, root4]
marked = [root1, root4]
