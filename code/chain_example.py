#!/usr/bin/env python3

import numpy as np
import helpers

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Edge:
    def __init__(self, p, q, label):
        self.p = p
        self.q = q
        self.label = label
    def equals(self, e):
        if((self.p == e.p and self.q == e.q) or (self.p == e.q and self.q == e.p)):
            return True
        return False

class Triangle:
    def __init__(self, v0, v1, v2, e0, e1, e2):
        self.v = [v0, v1, v2]
        self.e = [e0, e1, e2]
        self.refinement_edge = None
    def find_refinement_edge(self):
        refinement_edge = None
        for edge in self.e:
            if(refinement_edge == None)
                
    def shares_refinement_edge(self, t):
        refinement_edge = self.e[0]
        for i in range(1,2):
            if(refinement_edge.label > self.e[i].label):
                refinement_edge = self.e[i]
        for e in t.e:
            if(refinement_edge.equals(e) == True):
                return True
        return False

# https://www.tutorialspoint.com/python_data_structure/python_binary_tree.htm
class Node:
    # TODO: Need a bisection function?
    # TODO: Need funtion to return list of all leaf nodes
    # TODO: Need function to find leaf node with a given triangle
    def __init__(self, T, GT, FT, root):
        self.left = None
        self.right = None
        self.T = T
        self.GT = GT
        self.FT = FT
        self.ET = None
        self.root = root
    def bisect(self, mesh):
        if(np.minimum(np.minimum(self.T.e[0], self.T.e[1]), self.T.e[2]) == self.T.e[0]):
            P = helpers.midpoint(self.T.v[1],self.T.v[2])
            child_triangle1 = Triangle(P, self.T.v[2], self.T.v[0], self.T.e[1], self.T.e[2] + 1, self.T.e[0] + 2)
            child_triangle2 = Triangle(P, self.T.v[0], self.T.v[1], self.T.e[2], self.T.e[0] + 2, self.T.e[1] + 1)
            child_node1 = Node(child_triangle1, self.GT + 1, None, T.root)
            child_node2 = Node(child_triangle2, self.GT + 1, None, T.root)
            if((((self.T.v[2].x - self.T.v[0].x) * (self.T.v[1].x - self.T.v[0].x)) - ((self.T.v[2].y - self.T.v[0].y) * (self.T.v[1].y - self.T.v[0].y))) < 0):
                return (child_node1, child_node2)
            return (child_node2, child_node1)
        elif(np.minimum(np.minimum(self.T.e[0], self.T.e[1]), self.T.e[2]) == self.T.e[1]):
            P = helpers.midpoint(self.T.v[0],self.T.v[2])
            child1 = Triangle(P, self.T.v[0], self.T.v[1], self.T.e[2], self.T.e[0] + 1, self.T.e[1] + 2)
            child2 = Triangle(P, self.T.v[1], self.T.v[2], self.T.e[0], self.T.e[1] + 2, self.T.e[2] + 1)
            if((((self.T.v[0].x - self.T.v[1].x) * (self.T.v[2].x - self.T.v[1].x)) - ((self.T.v[0].y - self.T.v[1].y) * (self.T.v[2].y - self.T.v[1].y))) < 0):
                return (child1, child2)
            return (child2, child1)
        P = helpers.midpoint(self.T.v[0],self.T.v[1])
        child1 = Triangle(P, self.T.v[1], self.T.v[2], self.T.e[0], self.T.e[1] + 1, self.T.v[2] + 2)
        child2 = Triangle(P, self.T.v[2], self.T.v[0], self.T.e[1], self.T.e[2] + 2, self.T.v[0] + 1)
        if((((self.T.v[1].x - self.T.v[2].x) * (self.T.v[0].x - self.T.v[2].x)) - ((self.T.v[1].x - self.T.v[2].x) * (self.T.v[0].x - self.T.v[2].x))) < 0):
            return (child1, child2)
        return (child2, child1)

def refine(mesh, marked):
    pass

# Initial mesh in example
T1 = Triangle(Point(0,0), Point(1,-1), Point(1,1), 0, 1, 1)
T2 = Triangle(Point(0,0), Point(-1,-1), Point(1,-1), 0, 1, 1)
T3 = Triangle(Point(0,0), Point(-1,1), Point(-1,-1), 0, 1, 1)
T4 = Triangle(Point(0,0), Point(1,1), Point(-1,1), 0, 1, 1)

# Creating forest of binary trees
root1 = Node(T1, 0, None, None)
root2 = Node(T2, 0, None, None)
root3 = Node(T3, 0, None, None)
root4 = Node(T4, 0, None, None)
root1.root = root1
root2.root = root2
root3.root = root3
root4.root = root4

# TODO: Refine to get to same place as example
mesh = [root1, root2, root3, root4]
marked = [root1, root4]
