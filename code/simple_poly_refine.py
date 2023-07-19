#!/usr/bin/env python3
"""
Program that demonstrates the recursive REFINE procedure described in Nochetto and Vesser's Primer of AFEM.

The user is prompted to enter vertices of a simple polygon with more than 4 vertices, where this polygon is triangulated using the
Ear Clipping method to generate a conforming mesh that is then passed into the REFINE procedure.
"""

import matplotlib.pyplot as plt
import numpy as np

# https://statisticsglobe.com/circular-doubly-linked-list-python
# TODO: Nodes might not need the convex identifier
class Node:
    def __init__(self, coords, inter_angle, convex):
        self.coords = coords
        self.inter_angle = inter_angle
        self.convex = convex
        self.next = None
        self.prev = None
class CircularDoublyLinkedList:
    def __init__(self):
        self.head = None

    def is_empty(self):
        """
        Check if the list is empty.
        """
        return self.head is None

    def append(self, coords, inter_angle, convex):
        """
        Add a new node with the given coords at the end of the list.
        """
        new_node = Node(coords, inter_angle, convex)
        if self.is_empty():
            # If the list is empty, set the new node as the head
            self.head = new_node
        else:
            # If the list is not empty, update the references to maintain the circular structure
            last_node = self.head.prev
            last_node.next = new_node
            new_node.prev = last_node
        new_node.next = self.head
        self.head.prev = new_node

    def prepend(self, coords, inter_angle, convex):
        """
        Add a new node with the given coords at the beginning of the list.
        """
        new_node = Node(coords, inter_angle, convex)
        if self.is_empty():
            # If the list is empty, set the new node as the head
            self.head = new_node
        else:
            # If the list is not empty, update the references to maintain the circular structure
            last_node = self.head.prev
            new_node.next = self.head
            new_node.prev = last_node
            last_node.next = new_node
        self.head.prev = new_node
        self.head = new_node

    def insert_after(self, coords, after_coords, inter_angle, convex):
        """
        Insert a new node with the given coords after a specified node in the list.
        """
        if self.is_empty():
            raise Exception("List is empty")
        current = self.head
        while current.coords != after_coords:
            current = current.next
            if current == self.head:
                raise Exception(f"{after_coords} not found in the list")
        new_node = Node(coords, inter_angle, convex)
        new_node.next = current.next
        new_node.prev = current
        current.next.prev = new_node
        current.next = new_node

    def remove(self, coords):
        """
        Remove the node with the specified coords from the list.
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
            print(current.coords, current.inter_angle, current.convex, end=" ")
            current = current.next
            if current == self.head:
                break
        print()

# Midpoint function
def midpoint(v1,v2):
    return [(v2[0] + v1[0])/2, (v2[1] + v1[1])/2]

# Distance function
def distance(v1,v2):
    return np.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2)

# REFINE procedure
# TODO: Needs to be updated
def refine(t,n,max_n):
    # Base case
    if(n == max_n):
        return
    # Calculate coordinate of new vertex using midpoint formula
    # v = [(t[0][1][0] + t[0][0][0])/2, (t[0][1][1] + t[0][0][1])/2]
    v = midpoint(t[0][0],t[0][1])
    # Bisect t based on coordinate of new vertex
    ax.plot([v[0],t[0][2][0]],[v[1],t[0][2][1]],color='black',linestyle='-')
    # Define child triangles based off the bisection
    child1 = [
            [ [t[0][0][0],t[0][0][1]], [t[0][2][0],t[0][2][1]], [v[0],v[1]] ],
            [ t[1][2]+2, t[1][2]+2, t[1][0] ]
    ]
    child2 = [
            [ [t[0][1][0],t[0][1][1]], [t[0][2][0],t[0][2][1]], [v[0],v[1]] ],
            [ t[1][2]+2, t[1][2]+2, t[1][0] ]
    ]
    # Recurse
    refine(child1,n+1,max_n)
    refine(child2,n+1,max_n)

def main():
    # Get input from user
    xcoords = input("Please enter the x-coordinates defining your simple polygon with no holes in counter-clockwise order (e.g. -1,-1,1,1): ")
    ycoords = input("Please enter the corresponding y-coordinates defining your simple polygon with no holes in counter-clockwise order (e.g. 1,-1,-1,1): ")
    iterations = int(input("How many refinement iterations would you like?: "))
    regularity = float(input("What value for the shape regularity constant would you like?: "))

    # Parse coordinates
    # TODO: Decide whether vertices should just be arrays, tuples, or objects
    xs = []
    ys = []
    xtokens = xcoords.split(',')
    ytokens = ycoords.split(',')
    n = len(xtokens)
    # TODO: make sure n >= 4 and that len(ytokens) is the same value
    for i in range(0,n):
        xs.append(float(xtokens[i]))
        ys.append(float(ytokens[i]))

    # Plot points
    # TODO: see if ax can just be passed as an argument to avoid global variables
    global ax
    fig, ax = plt.subplots()
    ax.scatter(xs,ys,color='white')

    # Draw lines connecting points and fill doubly linked list with vertices
    polygon = CircularDoublyLinkedList()
    for i in range(0,n):
        # Draw lines between initial vertices
        ax.plot([xs[i],xs[(i+1)%n]],[ys[i],ys[(i+1)%n]],color='black',linestyle='-')
        # Add x,y coordinates of each vertex to polygon structure
        polygon.append([xs[i],ys[i]],None,None)

    # Compute interior angles of each vertex in the simple polygon
    convex_vi = []
    reflex_vi = []
    vi = polygon.head
    for i in range(0,n):
        # Compute interior angle of vi
        l1 = distance(vi.prev.coords,vi.coords)
        l2 = distance(vi.coords,vi.next.coords)
        tri_area_x2 = (vi.next.coords[0]-vi.prev.coords[0])*(vi.prev.coords[1]-vi.coords[1]) - (vi.prev.coords[0]-vi.coords[0])*(vi.next.coords[1]-vi.prev.coords[1])
        l3 = np.absolute(tri_area_x2 / distance(vi.prev.coords,vi.next.coords))
        vi.inter_angle = np.arccos(l3/l1) + np.arccos(l3/l2) # in radians
        # Check convexity
        if(tri_area_x2 > 0):
            vi.convex = True
            convex_vi.append(vi)
        else:
            vi.convex = False
            reflex_vi.append(vi)
        temp = vi.next # TODO: is this necessary?
        vi = temp
    
    """
    Identify ear tip vertices by the conditions that
        - they must be convex, and 
        - the closure of the triangle made up of a vertex and its two neighbor vertices does not contain any reflex vertices of the polygon.
    The second point can be tested very easily using a Barycentric Coordinate System.
    """
    # TODO: This should probably be a function
    ear_tips = []
    for i in convex_vi:
        # Check if there are any reflex vertices
        if(len(reflex_vi) > 0):
            p1 = i.coords
            p2 = i.prev.coords
            p3 = i.next.coords
            for j in reflex_vi:
                p = j.coords
                a = ( ((p2[1] - p3[1])*(p[0] - p3[0]) + (p3[0] - p2[0])*(p[1] - p3[1])) / 
                     ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1])))
                b = ( ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) / 
                     ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1])))
                c = 1 - a - b
                if((a < 0 or a > 1) and (b < 0 or b > 1) and (c < 0 or c > 1)):
                    ear_tips.append(i)
        # If not, then every vertex in the polygon is convex, futher implying that every vertex is an ear tip
        else:
            ear_tips.append(i)

    # Sort ear tips from smallest to largest interior angle
    for i in range(0,len(ear_tips)):
        for j in range(i+1,len(ear_tips)):
            if(ear_tips[i].inter_angle > ear_tips[j].inter_angle):
                temp = ear_tips[i]
                ear_tips[i] = ear_tips[j]
                ear_tips[j] = temp

    # TODO: Start clipping the ears and adding them to the initial triangulation/mesh
    mesh_0 = []
    
    plt.show()

main()
