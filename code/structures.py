#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import helpers as h
from shapely import geometry
import copy

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.boundary = False
    def equals(self, v):
        if(self.x == v.x and self.y == v.y):
            return True
        return False

class Vertex:
    def __init__(self, coords):
        self.coords = coords
        self.interior_angle = None
        self.prev = None
        self.next = None
    def calculate_interior_angle(self):
        if(self.prev != None and self.next != None):
            a = Point(self.prev.coords.x - self.coords.x, self.prev.coords.y - self.coords.y)
            b = Point(self.next.coords.x - self.coords.x, self.next.coords.y - self.coords.y)
            a_magnitude = np.sqrt(a.x**2 + a.y**2)
            b_magnitude = np.sqrt(b.x**2 + b.y**2)
            if(h.convex_check(self.prev.coords, self.coords, self.next.coords) > 0):
                self.interior_angle = np.arccos((a.x*b.x + a.y*b.y)/(a_magnitude * b_magnitude))
            else:
                self.interior_angle = 2*np.pi - np.arccos((a.x*b.x + a.y*b.y)/(a_magnitude * b_magnitude))

class Polygon:
    """
    Cyclical doubly linked list that represents the structure of a polygon.
    
    ...
    
    Attributes:
    -----------
    head : Vertex
        the first vertex in the polygon (user defined).

    Functions:
    ---------
    is_empty()
        Check if the list is empty.
    
    append(coords)
        Add a new vertex with the given coords at the end of the list.
    
    prepend(coords)
        Add a new vertex with the given coords at the beginning of the list.
    
    insert_after(coords, after_coords)
        Insert a new vertex with the given coords after a specified vertex in the list.
    
    remove(coords)
        Remove the vertex with the specified coords from the list.
    
    display()
        Display the elements of the circular doubly linked list.
    """
    def __init__(self):
        self.head = None

    def is_empty(self):
        """
        Check if the list is empty.
        """
        return self.head is None

    def append(self, coords):
        """
        Add a new vertex with the given coords at the end of the list.
        """
        new_vertex = Vertex(coords)
        if self.is_empty():
            # If the list is empty, set the new vertex as the head
            self.head = new_vertex
        else:
            # If the list is not empty, update the references to maintain the circular structure
            last_vertex = self.head.prev
            last_vertex.next = new_vertex
            new_vertex.prev = last_vertex
        new_vertex.next = self.head
        self.head.prev = new_vertex

    def prepend(self, coords):
        """
        Add a new vertex with the given coords at the beginning of the list.
        """
        new_vertex = Vertex(coords)
        if self.is_empty():
            # If the list is empty, set the new vertex as the head
            self.head = new_vertex
        else:
            # If the list is not empty, update the references to maintain the circular structure
            last_vertex = self.head.prev
            new_vertex.next = self.head
            new_vertex.prev = last_vertex
            last_vertex.next = new_vertex
        self.head.prev = new_vertex
        self.head = new_vertex

    def insert_after(self, coords, after_coords):
        """
        Insert a new vertex with the given coords after a specified vertex in the list.
        """
        if self.is_empty():
            raise Exception("List is empty")
        current = self.head
        while current.coords != after_coords:
            current = current.next
            if current == self.head:
                raise Exception(f"{after_coords} not found in the list")
        new_vertex = Vertex(coords)
        new_vertex.next = current.next
        new_vertex.prev = current
        current.next.prev = new_vertex
        current.next = new_vertex

    def remove(self, coords):
        """
        Remove the vertex with the specified coords from the list.
        """
        if self.is_empty():
            raise Exception("List is empty")
        current = self.head
        while current.coords != coords:
            current = current.next
            if current == self.head:
                raise Exception(f"{coords} not found in the list")
        current.prev.next = current.next
        current.next.prev = current.prev
        if current == self.head:
            self.head = current.next

    def display(self):
        """
        Display the elements of the circular doubly linked list.
        """
        if self.is_empty():
            print("List is empty")
            return
        current = self.head
        while True:
            print(current.coords, end=" ")
            current = current.next
            if current == self.head:
                break
        print()

class Edge:
    def __init__(self, p, q):
        self.p = p
        self.q = q
        self.dist = h.distance(p, q)
        self.label = None
    def equals(self, e):
        if(( (self.p.x == e.p.x and self.p.y == e.p.y) and (self.q.x == e.q.x and self.q.y == e.q.y) )
           or ( (self.p.x == e.q.x and self.p.y == e.q.y) and (self.q.x == e.p.x and self.q.y == e.p.y) )):
            return True
        return False

class Triangle:
    def __init__(self, v0, v1, v2, e0, e1, e2):
        self.v = [v0, v1, v2]
        self.e = [e0, e1, e2]
        """
        v = [v0, v1, v2]
        e = [e0, e1, e2]
        # TODO: this is most likely not actually doing anything
        for i in range(0,3):
            tau1 = Point(v[1].x - v[0].x, v[1].y - v[0].y)
            tau1_perp = Point(-tau1.y, tau1.x)
            tau2 = Point(v[2].x - v[0].x, v[2].y - v[0].y)
            if( ((tau2.x * tau1_perp.x) + (tau2.y * tau1_perp.y)) > 0 ):
                break
            temp_v0 = v[0]
            temp_v1 = v[1]
            temp_v2 = v[2]
            v[0] = temp_v2
            v[1] = temp_v0
            v[2] = temp_v1
            temp_e0 = e[0]
            temp_e1 = e[1]
            temp_e2 = e[2]
            e[0] = temp_e2
            e[1] = temp_e0
            e[2] = temp_e1
        self.v = v.copy()
        self.e = e.copy()
        """
    def label_longest_edge(self):
        for i in range(0,3):
            if(np.maximum(self.e[i].dist, np.maximum(self.e[(i+1)%3].dist, self.e[(i+2)%3].dist)) == self.e[i].dist):
                self.e[i].label = 0
                self.e[(i+1)%3].label = 1
                self.e[(i+2)%3].label = 1
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
        self.ET = None # NOTE: F(T) corresponds to value of E(T)
        self.neighbor0 = None
        self.neighbor1 = None
        self.neighbor2 = None
        self.root = None
        self.marked = False
    def find_refinement_edge(self):
        for i in range(0,3):
            if(np.minimum(self.T.e[i].label, np.minimum(self.T.e[(i+1)%3].label, self.T.e[(i+2)%3].label)) == self.T.e[i].label):
                self.ET = i
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
    def bisect(self, p, vertices):
        #if(self.T.v[(self.ET+1)%3].boundary == True and self.T.v[(self.ET+2)%3].boundary == True):
        poly = geometry.Polygon(vertices)
        boundary = poly.boundary
        point = geometry.Point(p.x, p.y)
        if(boundary.contains(point) == True):
            p.boundary = True
        ordered_v = [self.T.v[self.ET], self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3]]
        ordered_e = [self.T.e[self.ET], self.T.e[(self.ET+1)%3], self.T.e[(self.ET+2)%3]]
        e0 = Edge(ordered_v[1], p)
        e1 = Edge(p, ordered_v[0])
        e2 = Edge(ordered_v[0], ordered_v[1])
        e0.label = ordered_e[0].label + 2
        e1.label = ordered_e[1].label + 1
        e2.label = ordered_e[2].label
        child_triangle1 = Triangle(ordered_v[0], ordered_v[1], p, e0, e1, e2)
        e3 = Edge(p, ordered_v[2])
        e4 = Edge(ordered_v[2], ordered_v[0])
        e5 = Edge(ordered_v[0], p)
        e3.label = ordered_e[0].label + 2
        e4.label = ordered_e[1].label
        e5.label = ordered_e[2].label + 1
        child_triangle2 = Triangle(ordered_v[0], p, ordered_v[2], e3, e4, e5)
        child_node1 = Node(child_triangle1, self.GT + 1)
        child_node2 = Node(child_triangle2, self.GT + 1)
        child_node1.root = self.root
        child_node2.root = self.root
        child_node1.find_refinement_edge()
        child_node2.find_refinement_edge()
        self.left = child_node1
        self.right = child_node2
        child_node1.parent = self
        child_node2.parent = self
        child_node1.update_neighbor(child_node2)
        return (child_node1, child_node2)
    def iterate_along_line(self, i, return_dict, sigma, p, t, step):
        if(i == 0):
            diam_triangle = h.diam_of_triangle(self.T.v[self.ET], self.T.v[(self.ET+1)%3], p)
            area_root2 = h.triangle_area_root2(self.T.v[self.ET], self.T.v[(self.ET+1)%3], p)
            diam_circle = h.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            if(diam_circle <= area_root2 and area_root2 <= diam_triangle and diam_triangle <= sigma * diam_circle):
                return_dict[i] = p
                return
            t = t - step
            p = Point(
                self.T.v[(self.ET+2)%3].x+(t/(h.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].x-self.T.v[(self.ET+2)%3].x),
                self.T.v[(self.ET+2)%3].y+(t/(h.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].y-self.T.v[(self.ET+2)%3].y)
            )
            if(h.barycentric_point_check(self.T.v[self.ET], self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3], p) == False):
                return_dict[i] = None
                return
        if(i == 1):
            diam_triangle = h.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            area_root2 = h.triangle_area_root2(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            diam_circle = h.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            if(diam_circle <= area_root2 and area_root2 <= diam_triangle and diam_triangle <= sigma * diam_circle):
                return_dict[i] = p
                return
            t = t + step
            p = Point(
                self.T.v[(self.ET+2)%3].x+(t/(h.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].x-self.T.v[(self.ET+2)%3].x),
                self.T.v[(self.ET+2)%3].y+(t/(h.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].y-self.T.v[(self.ET+2)%3].y)
            )
            if(h.barycentric_point_check(self.T.v[self.ET], self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3], p) == False):
                return_dict[i] = None
                return
        self.iterate_along_line(i, return_dict, sigma, p, t, step)
    def shape_regularity_bisect(self, sigma):
        p = h.midpoint(self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3])
        t = h.distance(self.T.v[(self.ET+1)%3], p)
        step = t/100 # TODO: smaller step?
        manager = mp.Manager()
        return_dict = manager.dict()
        processes = []
        for i in range(0,2):
            proc = mp.Process(target=self.iterate_along_line, args=(i, return_dict, sigma, p, t, step))
            processes.append(proc)
            proc.start()
        for proc in processes:
            proc.join()
        points = []
        for _,value in return_dict.items():
            points.append(value)
        if(points[0] == None and points[1] == None):
            print("Error") # NOTE: should never get here
        elif(points[0] != None and points[1] == None):
            child1, child2 = self.bisect(points[0])
            return (child1, child2)
        elif(points[0] == None and points[1] != None):
            child1, child2 = self.bisect(points[1])
            return (child1, child2)
        else:
            upper_bound1 = sigma*h.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            diam1 = h.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            upper_bound2 = sigma*h.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            diam2 = h.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            if(upper_bound1 - diam1 <= upper_bound2 - diam2):
                child1, child2 = self.bisect(points[1])
                return (child1, child2)
            else:
                child1, child2 = self.bisect(points[0])
                return (child1, child2)
