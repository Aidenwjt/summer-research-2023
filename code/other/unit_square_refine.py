#!/usr/bin/env python3
"""
Program that demonstrates a really simple version of the recursive REFINE procedure described in Nochetto and Vesser's Primer of AFEM.

The domain used here is the unit square, and the initial conforming mesh is the set of 4 triangles that equally split up the unit square.
"""

import matplotlib.pyplot as plt
import numpy as np

def refine(t,n,max_n):
    # Base case
    if(n == max_n):
        return
    # Calculate coordinate of new vertex using midpoint formula
    v = [(t[0][1][0] + t[0][0][0])/2, (t[0][1][1] + t[0][0][1])/2]
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

iterations = int(input("How many refinement iterations would you like?: "))

# x,y coordinates for vertices of triangles
xs = [-1,-1,0,1,1]
ys = [-1,1,0,1,-1]

# Plot vertices
global ax
fig, ax = plt.subplots()
ax.scatter(xs,ys,color='white')

# Draw lines connecting vertices of initial conforming mesh
ax.plot([xs[0],xs[1]],[ys[0],ys[1]], color='black',linestyle='-')
ax.plot([xs[0],xs[2]],[ys[0],ys[2]], color='black',linestyle='-')
ax.plot([xs[1],xs[2]],[ys[1],ys[2]], color='black',linestyle='-')
ax.plot([xs[1],xs[3]],[ys[1],ys[3]], color='black',linestyle='-')
ax.plot([xs[3],xs[2]],[ys[3],ys[2]], color='black',linestyle='-')
ax.plot([xs[3],xs[4]],[ys[3],ys[4]], color='black',linestyle='-')
ax.plot([xs[4],xs[2]],[ys[4],ys[2]], color='black',linestyle='-')
ax.plot([xs[4],xs[0]],[ys[4],ys[0]], color='black',linestyle='-')

# Define triangles in the initial conforming mesh
t1 = [
        [ [xs[0],ys[0]] ,[xs[1],ys[1]] , [xs[2],ys[2]] ],
        [1,1,0]
]
t2 = [
        [ [xs[1],ys[1]] ,[xs[3],ys[3]] , [xs[2],ys[2]] ],
        [1,1,0]
]
t3 = [
        [ [xs[3],ys[3]] ,[xs[4],ys[4]] , [xs[2],ys[2]] ],
        [1,1,0]
]
t4 = [
        [ [xs[4],ys[4]] ,[xs[0],ys[0]] , [xs[2],ys[2]] ],
        [1,1,0]
]

# Recursively refine the initial mesh
mesh_0 = [t1,t2,t3,t4]
for t in mesh_0:
    refine(t,0,iterations)
plt.show()
