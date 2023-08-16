#!/usr/bin/env python3

import numpy as np
import structures as s
import helpers as h

# TODO:
# - unique_vertices_excluding_boundary(mesh)
#       -> returns a list of unique vertices in the mesh not including the vertices on the boundary
# - grad_phi(elem)
#       -> returns a 2D list of 3 gradient vectors corresponding to the 3 nodal basis functions of the 3 corresponding vertices in the element/triangle
# - galerkin_basis_coeffcients(mesh, vertices)
#       -> returns a list of Galerkin basis coefficients found in the TU = f system
# - phi_of_x(z, x)
#       -> returns the value of the nodal basis function corresponding to some vertex z evaluated at some value x in the domain
# - recreate_galerkin_solution(U, vertices, x)
#       -> returns the value of the Galerkin solution, recreated with the nodal basis functions, evaluated at some value x in the domain

def unique_vertices_excluding_boundary(mesh):
    vertices = []
    for elem in mesh:
        for v in elem.T.v:
            if(v.boundary == False):
                flag = False
                for v1 in vertices:
                    if(v.equals(v1) == True):
                        flag = True
                if(flag == False):
                    vertices.append(v)
    return vertices

def grad_phi(elem):
    det_XT = (elem.T.v[1].x * elem.T.v[2].y) + (elem.T.v[2].x * elem.T.v[0].y) + (elem.T.v[0].x * elem.T.v[1].y)
            - (elem.T.v[1].x * elem.T.v[0].y) - (elem.T.v[0].x * elem.T.v[2].y) - (elem.T.v[2].x * elem.T.v[1].y)
    grad_phi_v0 = [((elem.T.v[2].x * elem.T.v[0].y)-(elem.T.v[0].x * elem.T.v[2].y))/det_XT, ((elem.T.v[0].x * elem.T.v[1].y)-(elem.T.v[1].x * elem.T.v[0].y))/det_XT]
    grad_phi_v1 = [(elem.T.v[2].y - elem.T.v[0].y)/det_XT, (elem.T.v[0].y  - elem.T.v[1].y)/det_XT]
    grad_phi_v2 = [(elem.T.v[0].x - elem.T.v[2].x)/det_XT, (elem.T.v[1].x  - elem.T.v[0].x)/det_XT]
    return [grad_phi_v0, grad_phi_v1, grad_phi_v2]

def galerkin_basis_coeffcients(mesh, vertices, f):
    # Compute T matrix and f vector
    T = [[0]*len(vertices) for i in range(0, len(vertices))]
    f = [0]*len(vertices)
    for i in range(0, len(vertices)):
        sum_fi = 0
        for elem in mesh:
            if(elem.T.v[i] in vertices):
                sum_fi += f*((h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)/3)
        f[i] = sum_fi
        for j in range(0, len(vertices)):
            sum_Tij = 0
            for elem in mesh:
                if((elem.T.v[i] in vertices) and (elem.T.v[j] in vertices)):
                    grads = grad_phi(elem)
                    grad_phi_vi = grads[i]
                    grad_phi_vj = grads[j]
                    dot_grads = (grad_phi_vi[0] * grad_phi_vj[0]) + (grad_phi_vi[1] * grad_phi_vj[1])
                    area_of_T = h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2
                    sum_Tij += dot_grads * area_of_T
            T[i][j] = sum_Tij
    return np.linalg.solve(np.array(T), np.array(f)).tolist()

def phi_of_x(mesh, v, x):
    for elem in mesh:
        if(h.barycentric_point_check(elem.T.v[0], elem.T.v[1], elem.T.v[2], x) == True):
            for j in range(0,3):
                if(v == elem.T.v[j]):
                    det = (elem.T.v[(j+1)%3].x * elem.T.v[(j+2)%3].y) + (elem.T.v[(j+2)%3].x * x.y) + (x.x * elem.T.v[(j+1)%3].x)
                            - (elem.T.v[(j+1)%3].x * x.y) - (x.x * elem.T.v[(j+2)%3].x) - (elem.T.v[(j+2)%3].x * elem.T.v[(j+1)%3].y)
                    return (1/(2 * (h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)))*det
    return 0 # NOTE: I guess return nothing if the x value is not in the domain at all?

def recreate_galerkin_solution_at_x(mesh, U, vertices, x):
    galerkin_solution = 0
    for i in range(0, len(vertices)):
        galerkin_solution += U[i] * phi_of_x(mesh, vertices[i], x)
    return galerkin_solution
