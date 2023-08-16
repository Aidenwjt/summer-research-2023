#!/usr/bin/env python3

import numpy as np
import structures as s
import helpers as h

# TODO:
# - unique_vertices_excluding_boundary(initial_mesh)
#       -> returns a list of unique vertices in the mesh not including the vertices on the boundary
# - grad_phi(elem)
#       -> returns a 2D list of 3 gradient vectors corresponding to the 3 nodal basis functions of the 3 corresponding vertices in the element/triangle
# - galerkin_basis_coeffcients(initial_mesh, vertices)
#       -> returns a list of Galerkin basis coefficients found in the TU = f system
# - phi_of_x(z, x)
#       -> returns the value of the nodal basis function corresponding to some vertex z evaluated at some value x in the domain
# - recreate_galerkin_solution(U, vertices, x)
#       -> returns the value of the Galerkin solution, recreated with the nodal basis functions, evaluated at some value x in the domain
