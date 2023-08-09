#!/usr/bin/env python3

import numpy as np

class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y

def midpoint(P,Q):
    """ Returns the midpoint of two points in 2D. """
    return Point((Q.x + P.x)/2, (Q.y + P.y)/2)

def distance(P,Q):
    """ Returns the distance between two points in 2D. """
    return np.sqrt((Q.x-P.x)**2 + (Q.y-P.y)**2)

def convex_check(A,B,C):
    return ((C.x-A.x)*(A.y-B.y) - (A.x-B.x)*(C.y-A.y))

def tri_area_x2(A,B,C):
    """ Returns twice the area of the triangle ABC. """
    return np.absolute(convex_check(A,B,C))

def nearest_distance_to_line(A,B,P):
    """ Returns the nearest distance from the point P to the line AB. """
    return tri_area_x2(A,B,P) / distance(A,B)

def diam_of_triangle(A,B,C):
    """ Returns the diameter of the triangle ABC. """
    return np.maximum(np.maximum(distance(A,B),distance(B,C)),distance(C,A))

def triangle_area_root2(A,B,C):
    """ Returns the area of the triangle ABC root 2. """
    return np.sqrt(0.5*(tri_area_x2(A,B,C)))

def section_formula(A,B,m,n):
    """ Returns the point on the line AB defined by the section formula, given ratios m and n. """
    return Point((m*B.x + n*A.x)/(m+n),(m*B.y + n*A.y)/(m+n))

def diam_of_inscribed_circle(A,B,C):
    """ Returns the diameter of the circle inscribed in the triangle ABC. """
    bisector1 = section_formula(A, C, distance(A, B), distance(B, C))
    m1 = (bisector1.y - B.y) / (bisector1.x - B.x)
    b1 =  B.y - m1 * B.x
    bisector2 = section_formula(B, A,distance(B, C),distance(C, A))
    m2 = (bisector2.y - C.y)/(bisector2.x - C.x)
    b2 =  C.y - m2 * C.x
    x = (b2 - b1)/(m1 - m2)
    y = m2*x + b2
    radius = nearest_distance_to_line(A,B,Point(x,y))
    return 2*radius
