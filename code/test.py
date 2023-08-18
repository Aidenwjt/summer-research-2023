import matplotlib.pyplot as plt
import structures as s
import initial_mesh as im
import galerkin as g
import helpers as h

# TODO: this is broken for some reason
def refine_recursive(mesh, elem, poly):
    if(elem.FT() != None and elem.FT().GT < elem.GT):
        mesh = refine_recursive(mesh, elem.FT(), poly)
    child1, child2 = elem.bisect(h.midpoint(elem.T.v[(elem.ET+1)%3], elem.T.v[(elem.ET+2)%3]), poly)
    mesh.remove(elem)
    mesh.append(child1)
    mesh.append(child2)
    child1.update_neighbors()
    child2.update_neighbors()
    if elem.FT() == None:
        return mesh
    child3, child4 = elem.FT().bisect(h.midpoint(elem.FT().T.v[(elem.FT().ET+1)%3], elem.FT().T.v[(elem.FT().ET+2)%3]), poly)
    child3.update_neighbors()
    child4.update_neighbors()
    mesh.remove(elem.FT())
    mesh.append(child3)
    mesh.append(child4)
    return mesh

def refine(mesh, poly):
    marked = mesh.copy()
    for elem in marked:
        if elem.left == None:
            mesh = refine_recursive(mesh, elem, poly)
    return mesh

poly = [(0,0),(1,0),(1,1),(0,1)]
#poly = [(0,0),(4,0),(4,3),(6,3),(6,4),(3,4),(3,2),(0,2)]
mesh = im.create(poly)

vertices = g.unique_vertices_excluding_boundary(mesh)
for v in vertices:
    print("({},{})".format(v.x,v.y))
U = g.galerkin_basis_coefficients(mesh, vertices, 4)
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.5, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.5,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(1, 1, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(1,1))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.25))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.75, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.75,0.75))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.75))))

mesh = refine(mesh, poly)

vertices = g.unique_vertices_excluding_boundary(mesh)
for v in vertices:
    print("({},{})".format(v.x,v.y))
U = g.galerkin_basis_coefficients(mesh, vertices, 4)
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.5, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.5,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(1, 1, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(1,1))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.25))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.75, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.75,0.75))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.75))))

for elem in mesh:
    print("({},{}),({},{}),({},{}) ET = {}".format(elem.T.v[0].x, elem.T.v[0].y, elem.T.v[1].x, elem.T.v[1].y, elem.T.v[2].x, elem.T.v[2].y, elem.ET))

mesh = refine(mesh, poly)

vertices = g.unique_vertices_excluding_boundary(mesh)
for v in vertices:
    print("({},{})".format(v.x,v.y))
U = g.galerkin_basis_coefficients(mesh, vertices, 4)
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.5, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.5,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(1, 1, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(1,1))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.25))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.75, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.75,0.75))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.75))))
"""
"""

fig, ax = plt.subplots()
for elem in mesh:
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
