#!/usr/bin/env python3

import numpy as np
import structures as s
import helpers as h

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

# TODO: separate earclipping functionality so I can implement other triangulation methods
def create(vertices):
    """ Returns a triangulation of the given vertices of a simple polygon using the ear clipping method. """
    # Construct polygon
    polygon = s.Polygon()
    for v in vertices:
        vertex = s.Point(v[0],v[1])
        vertex.boundary = True
        polygon.append(vertex)

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
        T = s.Triangle(v_prev.coords,v.coords,v_next.coords, s.Edge(v.coords,v_next.coords), s.Edge(v_prev.coords,v_next.coords), s.Edge(v_prev.coords,v.coords))
        T.label_longest_edge()
        elem = s.Node(T, 0)
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
   
    # NOTE: might not need this block
    sigma = 0
    for elem in triangulation:
        temp = h.diam_of_triangle(elem.T.v[0],elem.T.v[1],elem.T.v[2])/h.diam_of_inscribed_circle(elem.T.v[0],elem.T.v[1],elem.T.v[2]) 
        if(temp > sigma):
            sigma = temp

    initial_mesh = []
    for elem in triangulation:
        child1, child2 = elem.bisect(h.midpoint(elem.T.v[(elem.ET+1)%3], elem.T.v[(elem.ET+2)%3]), vertices)
        grandchild1, grandchild2 = child1.bisect(h.midpoint(child1.T.v[(child1.ET+1)%3], child1.T.v[(child1.ET+2)%3]), vertices)
        grandchild3, grandchild4 = child2.bisect(h.midpoint(child2.T.v[(child2.ET+1)%3], child2.T.v[(child2.ET+2)%3]), vertices)
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
