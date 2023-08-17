#!/usr/bin/env python3

import numpy as np
import structures as s

def midpoint(P,Q):
    """ Returns the midpoint of two points in 2D. """
    return s.Point((Q.x + P.x)/2, (Q.y + P.y)/2)

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
    return s.Point((m*B.x + n*A.x)/(m+n),(m*B.y + n*A.y)/(m+n))

def diam_of_inscribed_circle(A,B,C):
    """ Returns the diameter of the circle inscribed in the triangle ABC. """
    D = section_formula(A, C, distance(A, B), distance(B, C))
    F = section_formula(B, D, distance(B, C), distance(D, C))
    return 2*nearest_distance_to_line(A,B,F)

def barycentric_point_check(A,B,C,P):
    a = ( ((B.y - C.y)*(P.x - C.x) + (C.x - B.x)*(P.y - C.y)) / 
         ((B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y)) )
    b = ( ((C.y - A.y)*(P.x - C.x) + (A.x - C.x)*(P.y - C.y)) / 
         ((B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y)))
    c = 1 - a - b
    if((a >= 0 and a <= 1) and (b >= 0 and b <= 1) and (c >= 0 and c <= 1)):
        return True
    return False
