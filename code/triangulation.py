#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import helpers

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

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
            if(helpers.convex_check(self.prev.coords, self.coords, self.next.coords) > 0):
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
        self.dist = helpers.distance(p, q)
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
        for i in range(0,3):
            tao1 = Point(self.v[1].x - self.v[0].x, self.v[1].y - self.v[0].y)
            tao1_perp = Point(-tao1.y, tao1.x)
            tao2 = Point(self.v[2].x - self.v[0].x, self.v[2].y - self.v[0].y)
            if( ((tao2.x * tao1_perp.x) + (tao2.y * tao1_perp.y)) > 0 ):
                break
            temp_v0 = self.v[0]
            temp_v1 = self.v[1]
            temp_v2 = self.v[2]
            self.v[0] = temp_v2
            self.v[1] = temp_v0
            self.v[2] = temp_v1
            temp_e0 = self.e[0]
            temp_e1 = self.e[1]
            temp_e2 = self.e[2]
            self.e[0] = temp_e2
            self.e[1] = temp_e0
            self.e[2] = temp_e1
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
    def bisect(self, p):
        new_edge1 = Edge(p, self.T.v[(self.ET+1)%3])
        new_edge1.label = self.T.e[self.ET].label + 2
        new_edge2 = Edge(p, self.T.v[self.ET])
        new_edge2.label = self.T.e[(self.ET+1)%3].label + 1
        child_triangle1 = Triangle(self.T.v[self.ET], self.T.v[(self.ET+1)%3], p, new_edge1, new_edge2, self.T.e[(self.ET+2)%3])
        new_edge3 = Edge(p, self.T.v[(self.ET+2)%3])
        new_edge3.label = self.T.e[self.ET].label + 2
        new_edge4 = Edge(p, self.T.v[self.ET])
        new_edge4.label = self.T.e[(self.ET+2)%3].label + 1
        child_triangle2 = Triangle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3], new_edge3, self.T.e[(self.ET+1)%3], new_edge4)
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
            diam_triangle = helpers.diam_of_triangle(self.T.v[self.ET], self.T.v[(self.ET+1)%3], p)
            area_root2 = helpers.triangle_area_root2(self.T.v[self.ET], self.T.v[(self.ET+1)%3], p)
            diam_circle = helpers.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            print("Process 0: diam_triangle = {} area_root2 = {} diam_circle = {}".format(diam_triangle, area_root2, diam_circle))
            if(diam_circle <= area_root2 and area_root2 <= diam_triangle and diam_triangle <= sigma * diam_circle):
                return_dict[i] = p
                return
            t = t - step
            p = Point(
                self.T.v[(self.ET+2)%3].x+(t/(helpers.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].x-self.T.v[(self.ET+2)%3].x),
                self.T.v[(self.ET+2)%3].y+(t/(helpers.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].y-self.T.v[(self.ET+2)%3].y)
            )
            print("({},{})".format(p.x,p.y))
            #if(p.x == self.T.v[(self.ET+2)%3].x and p.y == self.T.v[(self.ET+2)%3].y):
            if(helpers.barycentric_point_check(self.T.v[self.ET], self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3], p) == False):
                print("Here1")
                return_dict[i] = None
                return
        if(i == 1):
            diam_triangle = helpers.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            area_root2 = helpers.triangle_area_root2(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            diam_circle = helpers.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            print("Process 1: diam_triangle = {} area_root2 = {} diam_circle = {}".format(diam_triangle, area_root2, diam_circle))
            if(diam_circle <= area_root2 and area_root2 <= diam_triangle and diam_triangle <= sigma * diam_circle):
                return_dict[i] = p
                return
            t = t + step
            p = Point(
                self.T.v[(self.ET+2)%3].x+(t/(helpers.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].x-self.T.v[(self.ET+2)%3].x),
                self.T.v[(self.ET+2)%3].y+(t/(helpers.distance(self.T.v[(self.ET+2)%3],self.T.v[(self.ET+1)%3])))*(self.T.v[(self.ET+1)%3].y-self.T.v[(self.ET+2)%3].y)
            )
            print("({},{})".format(p.x,p.y))
            #if(p.x == self.T.v[(self.ET+1)%3].x and p.y == self.T.v[(self.ET+1)%3].y):
            if(helpers.barycentric_point_check(self.T.v[self.ET], self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3], p) == False):
                print("Here2")
                return_dict[i] = None
                return
        self.iterate_along_line(i, return_dict, sigma, p, t, step)
    def shape_regularity_bisect(self, sigma):
        p = helpers.midpoint(self.T.v[(self.ET+1)%3], self.T.v[(self.ET+2)%3])
        print("p: ({},{})".format(p.x,p.y))
        t = helpers.distance(self.T.v[(self.ET+1)%3], p)
        step = t/100
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
            print("Error")
        elif(points[0] != None and points[1] == None):
            child1, child2 = self.bisect(points[0])
            return (child1, child2)
        elif(points[0] == None and points[1] != None):
            child1, child2 = self.bisect(points[1])
            return (child1, child2)
        else:
            upper_bound1 = sigma*helpers.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            diam1 = helpers.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+2)%3])
            upper_bound2 = sigma*helpers.diam_of_inscribed_circle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            diam2 = helpers.diam_of_triangle(self.T.v[self.ET], p, self.T.v[(self.ET+1)%3])
            if(upper_bound1 - diam1 <= upper_bound2 - diam2):
            #if(helpers.distance(self.T.v[(self.ET+2)%3], points[0]) >= helpers.distance(self.T.v[(self.ET+1)%3], points[1])):
                child1, child2 = self.bisect(points[1])
                return (child1, child2)
            else:
                child1, child2 = self.bisect(points[0])
                return (child1, child2)

def ear_tip_status(v,reflex_v):
    """ 
    Identify whether or not a vertex is an ear tip or not based on the conditions that ear tips
        - must be convex, and 
        - the closure of the triangle made up of the vertex and its two neighbor vertices does not contain any reflex vertices of the polygon,
            except possibly the vertex's neighbors.
    """
    if(len(reflex_v) > 0):
        p1 = v.coords
        p2 = v.prev.coords
        p3 = v.next.coords
        for i in reflex_v:
            p = i.coords
            if((p != p2) and (p != p3)):
                a = ( ((p2.y - p3.y)*(p.x - p3.x) + (p3.x - p2.x)*(p.y - p3.y)) / 
                     ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y)) )
                b = ( ((p3.y - p1.y)*(p.x - p3.x) + (p1.x - p3.x)*(p.y - p3.y)) / 
                     ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y)))
                c = 1 - a - b
                if((a >= 0 and a <= 1) and (b >= 0 and b <= 1) and (c >= 0 and c <= 1)):
                    return False
    return True

def sort_ear_tips(ear_tips):
    """ Sorts an array of ear tips from smallest to largest interior angle. """
    for i in range(0,len(ear_tips)):
        for j in range(i+1,len(ear_tips)):
            if(ear_tips[i].interior_angle > ear_tips[j].interior_angle):
                temp = ear_tips[i]
                ear_tips[i] = ear_tips[j]
                ear_tips[j] = temp

def update_ear_tip_status(v,convex_v,reflex_v,ear_tips):
    """ Updates lists of convex, reflex, and ear tip vertices based on the given vertex v. """
    if(v in reflex_v):
        if(v.interior_angle < np.pi):
            reflex_v.remove(v)
            convex_v.append(v)
    if(v in convex_v):
        if((ear_tip_status(v,reflex_v) == True) and (v not in ear_tips)):
            ear_tips.append(v)
            ear_tips = sort_ear_tips(ear_tips)
    if(((v in reflex_v) or (ear_tip_status(v,reflex_v) == False)) and (v in ear_tips)):
            ear_tips.remove(v)

def earclipping(vertices):
    """ Returns a triangulation of the given vertices of a simple polygon using the ear clipping method. """
    # Construct polygon
    polygon = Polygon()
    for v in vertices:
        polygon.append(Point(v[0],v[1]))

    # Calculate interior angles of each vertex in the simple polygon
    v = polygon.head
    convex_v = []
    reflex_v = []
    for i in range(0,len(vertices)):
        v.calculate_interior_angle()
        if(v.interior_angle < np.pi):
            convex_v.append(v)
        else:
            reflex_v.append(v)
        v = v.next

    # Find the ear tips
    ear_tips = []
    for v in convex_v:
        if(ear_tip_status(v,reflex_v) == True):
            ear_tips.append(v)
    sort_ear_tips(ear_tips)

    #  Start clipping the ears
    triangulation = []
    while(True):
        # Construct triangle and add to triangulation
        v = ear_tips[0]
        v_prev = v.prev
        v_next = v.next
        T = Triangle(v_prev.coords,v.coords,v_next.coords, Edge(v.coords,v_next.coords), Edge(v_prev.coords,v_next.coords), Edge(v_prev.coords,v.coords))
        T.label_longest_edge()
        elem = Node(T, 0)
        elem.find_refinement_edge()
        elem.root = elem
        triangulation.append(elem)
        # Update relationships in polygon
        polygon.remove(v.coords)
        # Let this vertex no longer be an ear tip
        ear_tips.remove(v)
        # Breaking out in the middle avoids runtime warning
        if(len(triangulation) == len(vertices)-2):
            break
        # Compute new interior angles of the neighbor vertices
        v_prev.calculate_interior_angle()
        v_next.calculate_interior_angle()
        # Update ear tip status of the neighbor vertices
        update_ear_tip_status(v_prev,convex_v,reflex_v,ear_tips)
        update_ear_tip_status(v_next,convex_v,reflex_v,ear_tips)
    
    sigma = 0
    for elem in triangulation:
        temp = helpers.diam_of_triangle(elem.T.v[0],elem.T.v[1],elem.T.v[2])/helpers.diam_of_inscribed_circle(elem.T.v[0],elem.T.v[1],elem.T.v[2]) 
        if(temp > sigma):
            sigma = temp

    initial_mesh = []
    for elem in triangulation:
        #child1, child2 = elem.shape_regularity_bisect(sigma)
        #grandchild1, grandchild2 = child1.shape_regularity_bisect(sigma)
        #grandchild3, grandchild4 = child2.shape_regularity_bisect(sigma)
        child1, child2 = elem.bisect(helpers.midpoint(elem.T.v[(elem.ET+1)%3], elem.T.v[(elem.ET+2)%3]))
        grandchild1, grandchild2 = child1.bisect(helpers.midpoint(child1.T.v[(child1.ET+1)%3], child1.T.v[(child1.ET+2)%3]))
        grandchild3, grandchild4 = child2.bisect(helpers.midpoint(child2.T.v[(child2.ET+1)%3], child2.T.v[(child2.ET+2)%3]))
        grandchildren = [grandchild1, grandchild2, grandchild3, grandchild4]
        for grandchild in grandchildren:
            neighbors = [grandchild.neighbor0, grandchild.neighbor1, grandchild.neighbor2]
            for i in range(0,3):
                if(neighbors[i] == None or neighbors[i].parent != grandchild.parent):
                    grandchild.T.e[i].label = 1
                else:
                    grandchild.T.e[i].label = 0
            initial_mesh.append(grandchild)

    for i in range(0, len(initial_mesh)):
        for j in range(0, len(initial_mesh)):
            if(i != j):
                initial_mesh[i].update_neighbor(initial_mesh[j])
    return initial_mesh
