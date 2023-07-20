#!/usr/bin/env python3
"""
Program that demonstrates the recursive REFINE procedure described in Nochetto and Vesser's Primer of AFEM.

In this case, the user is prompted to enter vertices of a simple polygon with more than 4 vertices, where this polygon is triangulated using the
Ear Clipping Method to generate an initial mesh that is then passed into the REFINE procedure.
"""

import matplotlib.pyplot as plt
import numpy as np

class Triangle:
    def __init__(self,v1,v2,v3,e1,e2,e3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3

# https://statisticsglobe.com/circular-doubly-linked-list-python
class Vertex:
    """
    Object that represents a vertex in a polygon.
    
    ...
    
    Attributes:
    -----------
    coords : tuple
        a tuple containing the x,y coordinates of the vertex, i.e. coords = (x,y).
    
    inter_angle : float
        the interior angle the vertex shares with its two neighbor vertices.
    
    prev : Vertex
        neighbor vertex.
    
    next : Vertex
        neighbor vertex.
    """
    def __init__(self, coords, inter_angle):
        self.coords = coords
        self.inter_angle = inter_angle
        self.prev = None
        self.next = None

class CircularDoublyLinkedList:
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
    
    append(coords, inter_angle)
        Add a new vertex with the given coords at the end of the list.
    
    prepend(coords, inter_angle)
        Add a new vertex with the given coords at the beginning of the list.
    
    insert_after(coords, after_coords, inter_angle)
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

    def append(self, coords, inter_angle):
        """
        Add a new vertex with the given coords at the end of the list.
        """
        new_vertex = Vertex(coords, inter_angle)
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

    def prepend(self, coords, inter_angle):
        """
        Add a new vertex with the given coords at the beginning of the list.
        """
        new_vertex = Vertex(coords, inter_angle)
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

    def insert_after(self, coords, after_coords, inter_angle):
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
        new_vertex = Vertex(coords, inter_angle)
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
            print(current.coords, current.inter_angle, end=" ")
            current = current.next
            if current == self.head:
                break
        print()

def midpoint(v1,v2):
    """ Return the midpoint of two points in 2D. """
    return [(v2.coords[0] + v1.coords[0])/2, (v2.coords[1] + v1.coords[1])/2]

def distance(v1,v2):
    """ Return the distance between two vertices in 2D. """
    return np.sqrt((v2.coords[0]-v1.coords[0])**2 + (v2.coords[1]-v1.coords[1])**2)

def interior_angle(vi):
    """ Calculate the interior angle that a vertex shares with its two neighbor vertices in a polygon. """
    l1 = distance(vi.prev,vi)
    l2 = distance(vi,vi.next)
    tri_area_x2 = (vi.next.coords[0]-vi.prev.coords[0])*(vi.prev.coords[1]-vi.coords[1]) - (vi.prev.coords[0]-vi.coords[0])*(vi.next.coords[1]-vi.prev.coords[1])
    l3 = np.absolute(tri_area_x2 / distance(vi.prev,vi.next)) # TODO: RuntimeWarning: invalid value encountered in scalar divide
    if(tri_area_x2 > 0):
        vi.inter_angle = np.arccos(l3/l1) + np.arccos(l3/l2) # in radians
    else:
        vi.inter_angle = 2*(np.pi) - (np.arccos(l3/l1) + np.arccos(l3/l2)) # in radians

def ear_tip_status(vi,reflex_vi):
    """ 
    Identify whether or not a vertex is an ear tip or not based on the conditions that ear tips
        - must be convex, and 
        - the closure of the triangle made up of the vertex and its two neighbor vertices does not contain any reflex vertices of the polygon,
            except possibly the vertex's neighbors.
    The second point can be tested very easily using a Barycentric Coordinate System.
    """
    if(len(reflex_vi) > 0):
        p1 = vi.coords
        p2 = vi.prev.coords
        p3 = vi.next.coords
        for j in reflex_vi:
            if((j.coords != vi.prev.coords) and (j.coords != vi.next.coords)):
                p = j.coords
                a = ( ((p2[1] - p3[1])*(p[0] - p3[0]) + (p3[0] - p2[0])*(p[1] - p3[1])) / 
                     ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1])))
                b = ( ((p3[1] - p1[1])*(p[0] - p3[0]) + (p1[0] - p3[0])*(p[1] - p3[1])) / 
                     ((p2[1] - p3[1])*(p1[0] - p3[0]) + (p3[0] - p2[0])*(p1[1] - p3[1])))
                c = 1 - a - b
                if((a >= 0 and a <= 1) and (b >= 0 and b <= 1) and (c >= 0 and c <= 1)):
                    return False
    return True

def sort_ear_tips(ear_tips):
    """ Sort an array of ear tips from smallest to largest interior angle. """
    for i in range(0,len(ear_tips)):
        for j in range(i+1,len(ear_tips)):
            if(ear_tips[i].inter_angle > ear_tips[j].inter_angle):
                temp = ear_tips[i]
                ear_tips[i] = ear_tips[j]
                ear_tips[j] = temp
    return ear_tips

def update_ear_tip_status(vi,convex_vi,reflex_vi,ear_tips):
    if(vi in reflex_vi):
        if(vi.inter_angle < np.pi):
            reflex_vi.remove(vi)
            convex_vi.append(vi)
    if(vi in convex_vi):
        if((ear_tip_status(vi,reflex_vi) == True) and (vi not in ear_tips)):
            ear_tips.append(vi)
            ear_tips = sort_ear_tips(ear_tips)
    if(((vi in reflex_vi) or (ear_tip_status(vi,reflex_vi) == False)) and (vi in ear_tips)):
            ear_tips.remove(vi)
    return convex_vi,reflex_vi,ear_tips

# REFINE procedure
# TODO: Needs to be updated
def refine(t,iteration,max_iterations,regularity,ax):
    # Base case
    if(iteration == max_max_iterations):
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
    # Get input from user and parse coordinates
    xcoords = input("Please enter the x-coordinates defining your simple polygon of 4 or more vertices with no holes in counter-clockwise order (e.g. -1,-1,1,1): ")
    ycoords = input("Please enter the corresponding y-coordinates (e.g. 1,-1,-1,1): ")
    xtokens = xcoords.split(',')
    ytokens = ycoords.split(',')
    n = len(xtokens)
    while((len(ytokens) != n) or (n < 4)):
        xcoords = input("Please enter the x-coordinates defining your simple polygon of 4 or more vertices with no holes in counter-clockwise order (e.g. -1,-1,1,1): ")
        ycoords = input("Please enter the corresponding y-coordinates (e.g. 1,-1,-1,1): ")
        xtokens = xcoords.split(',')
        ytokens = ycoords.split(',')
        n = len(xtokens)
    iterations = int(input("How many refinement iterations would you like? (whole number greater than or equal to zero): "))
    regularity = float(input("What value for the shape regularity constant would you like? (any number greater than one): "))
    xs = []
    ys = []
    for i in range(0,n):
        xs.append(float(xtokens[i]))
        ys.append(float(ytokens[i]))

    # Plot points
    fig, ax = plt.subplots()
    ax.scatter(xs,ys,color='white')

    # Draw lines connecting points and fill doubly linked list with vertices
    polygon = CircularDoublyLinkedList()
    for i in range(0,n):
        # Draw lines between initial vertices
        ax.plot([xs[i],xs[(i+1)%n]],[ys[i],ys[(i+1)%n]],color='black',linestyle='-')
        # Add x,y coordinates of each vertex to polygon structure
        polygon.append((xs[i],ys[i]),None)

    # Compute interior angles of each vertex in the simple polygon
    vi = polygon.head
    convex_vi = []
    reflex_vi = []
    for i in range(0,n):
        interior_angle(vi)
        if(vi.inter_angle < np.pi):
            convex_vi.append(vi)
        else:
            reflex_vi.append(vi)
        vi = vi.next

    # Find ear tips
    ear_tips = []
    for vi in convex_vi:
        if(ear_tip_status(vi,reflex_vi) == True):
            ear_tips.append(vi)
    ear_tips = sort_ear_tips(ear_tips)

    #  Start clipping the ears and adding them to the initial mesh
    mesh_0 = []
    while(len(mesh_0) < n-2):
        # Construct triangle
        vi = ear_tips[0]
        vi_prev = vi.prev
        vi_next = vi.next
        # mesh_0.append([[vi_prev.coords, vi.coords, vi_next.coords], [0,0,0]]) # [[vi_prev,vi,vi_next],[edge opposite vi_prev,edge opposite vi,edge opposite vi_next]]
        mesh_0.append(Triangle(vi_prev,vi,vi_next,0,0,0))
        # Update relationships in polygon
        polygon.remove(vi.coords)
        # Let this vertex no longer be an ear tip
        ear_tips.remove(vi)
        # Compute new interior angles of the neighbor vertices
        interior_angle(vi_prev)
        interior_angle(vi_next)
        # Update ear tip status of the neighbor vertices
        convex_vi, reflex_vi, ear_tips = update_ear_tip_status(vi_prev,convex_vi,reflex_vi,ear_tips)
        convex_vi, reflex_vi, ear_tips = update_ear_tip_status(vi_next,convex_vi,reflex_vi,ear_tips)

    # TODO: Implement shape regularity for bisections
    # TODO: Create distance from point to line function (also returns point on the line?)
    # TODO: Generate the special labeling of the initial mesh and draw the initial mesh.
    for t in mesh_0:
        ax.plot([t.v1.coords[0],t.v2.coords[0]],[t.v1.coords[1],t.v2.coords[1]],color='black',linestyle='-')
        ax.plot([t.v2.coords[0],t.v3.coords[0]],[t.v2.coords[1],t.v3.coords[1]],color='black',linestyle='-')
        ax.plot([t.v3.coords[0],t.v1.coords[0]],[t.v3.coords[1],t.v1.coords[1]],color='black',linestyle='-')
   
    """
    for t in mesh_0:
        refine(t,0,iterations,regularity,ax)
    """
    plt.show()

main()
