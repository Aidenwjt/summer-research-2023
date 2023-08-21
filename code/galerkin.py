#!/usr/bin/env python3

import numpy as np
import structures as s
import helpers as h

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

"""
def grad_phi(elem):
    #print("({},{})({},{})({},{})".format(elem.T.v[0].x, elem.T.v[0].y, elem.T.v[1].x, elem.T.v[1].y, elem.T.v[2].x, elem.T.v[2].y))
    XT = [[1, 1, 1],[elem.T.v[0].x, elem.T.v[1].x, elem.T.v[2].x],[elem.T.v[0].y, elem.T.v[1].y, elem.T.v[2].y]]
    XT_inv = np.linalg.inv(np.array(XT))
    grads = np.matmul(XT_inv, np.array([[0,0],[1,0],[0,1]])).tolist()
    #print(grads[0], grads[1], grads[2])
    return grads
    det_XT = (elem.T.v[1].x * elem.T.v[2].y) + (elem.T.v[2].x * elem.T.v[0].y) + (elem.T.v[0].x * elem.T.v[1].y) - (elem.T.v[1].x * elem.T.v[0].y) - (elem.T.v[0].x * elem.T.v[2].y) - (elem.T.v[2].x * elem.T.v[1].y)
    grad_phi_v0 = [((elem.T.v[2].x * elem.T.v[0].y)-(elem.T.v[0].x * elem.T.v[2].y))/det_XT, ((elem.T.v[0].x * elem.T.v[1].y)-(elem.T.v[1].x * elem.T.v[0].y))/det_XT]
    grad_phi_v1 = [(elem.T.v[2].y - elem.T.v[0].y)/det_XT, (elem.T.v[0].y  - elem.T.v[1].y)/det_XT]
    grad_phi_v2 = [(elem.T.v[0].x - elem.T.v[2].x)/det_XT, (elem.T.v[1].x  - elem.T.v[0].x)/det_XT]
    print(grad_phi_v0, grad_phi_v1, grad_phi_v2)
    return [grad_phi_v0, grad_phi_v1, grad_phi_v2]
"""
def grad_phi(elem, j):
    return [(elem.T.v[(j+1)%3].y - elem.T.v[(j+2)%3].y)/(2*(h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)),
            (elem.T.v[(j+2)%3].x - elem.T.v[(j+1)%3].x)/(2*(h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2))]

def galerkin_basis_coefficients(mesh, vertices, scalar_f):
    # Compute T matrix and f vector
    T = [[0]*len(vertices) for i in range(0, len(vertices))]
    f = [0]*len(vertices)
    for i in range(0, len(vertices)):
        sum_fi = 0
        for elem in mesh:
            for k in range(0, 3):
                if(vertices[i].equals(elem.T.v[k]) == True):
                    sum_fi += scalar_f*((h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)/3)
                    #print(np.linalg.det(np.array([ [elem.T.v[1].x - elem.T.v[0].x, elem.T.v[2].x - elem.T.v[0].x], [elem.T.v[1].y - elem.T.v[0].y, elem.T.v[2].y - elem.T.v[0].y] ]))/6)
        f[i] = sum_fi
        for j in range(0, len(vertices)):
            sum_Tij = 0
            for elem in mesh:
                for k in range(0,3):
                    if((vertices[i].equals(elem.T.v[k]) == True and vertices[j].equals(elem.T.v[k]) == True) or
                       (vertices[i].equals(elem.T.v[k]) == True and vertices[j].equals(elem.T.v[(k+1)%3]) == True) or
                       (vertices[i].equals(elem.T.v[(k+1)%3]) == True and vertices[j].equals(elem.T.v[k]) == True)):
                        for l in range(0, 3):
                            if(vertices[i].equals(elem.T.v[l])):
                                grad_phi_vi = grad_phi(elem, l)
                        for m in range(0, 3):
                            if(vertices[j].equals(elem.T.v[m])):
                                grad_phi_vj = grad_phi(elem, m)
                        dot_grads = (grad_phi_vi[0] * grad_phi_vj[0]) + (grad_phi_vi[1] * grad_phi_vj[1])
                        area_of_T = h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2
                        sum_Tij += dot_grads * area_of_T
                        break
            T[i][j] = sum_Tij
    for i in range(0, len(vertices)):
        print(T[i])
    U = np.linalg.solve(np.array(T), np.array(f)).tolist()
    print(U)
    print(f)
    return U

def phi_of_x(mesh, v, x):
    for elem in mesh:
        for j in range(0, 3):
            if(v.equals(elem.T.v[j]) == True):
                if(h.barycentric_point_check(elem.T.v[0], elem.T.v[1], elem.T.v[2], x) == True):
                    det1 = np.linalg.det(np.array([[1, x.x, x.y], [1, elem.T.v[(j+1)%3].x, elem.T.v[(j+1)%3].y], [1, elem.T.v[(j+2)%3].x, elem.T.v[(j+2)%3].y]]))
                    det2 = np.linalg.det(np.array([[1, elem.T.v[j].x, elem.T.v[j].y], [1, elem.T.v[(j+1)%3].x, elem.T.v[(j+1)%3].y], [1, elem.T.v[(j+2)%3].x, elem.T.v[(j+2)%3].y]]))
                    #det = np.linalg.det(np.array([[1,1,1], [x.x, elem.T.v[(j+1)%3].x, elem.T.v[(j+2)%3].x], [x.y, elem.T.v[(j+1)%3].y, elem.T.v[(j+2)%3].y]]))
                    #print(det)
                    #det2 = (elem.T.v[(j+1)%3].x * elem.T.v[(j+2)%3].y) + (elem.T.v[(j+2)%3].x * x.y) + (x.x * elem.T.v[(j+1)%3].y) - (elem.T.v[(j+1)%3].x * x.y) - (x.x * elem.T.v[(j+2)%3].y) - (elem.T.v[(j+2)%3].x * elem.T.v[(j+1)%3].y)
                    #print(det2)
                    #total += (1/(2 * (h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)))*det
                    #print(total)
                    #return (1/(2 * (h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)))*det
                    #print(det1/det2)
                    #print((1/(2 * (h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)))*det)
                    #print(det1/det2 == (1/(2 * (h.triangle_area_root2(elem.T.v[0], elem.T.v[1], elem.T.v[2])**2)))*det)
                    return det1/det2
    return 0

def recreate_galerkin_solution_at_x(mesh, U, vertices, x):
    galerkin_solution = 0
    for i in range(0, len(vertices)):
        galerkin_solution += U[i] * phi_of_x(mesh, vertices[i], x)
    return galerkin_solution
