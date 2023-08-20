import matplotlib.pyplot as plt
import structures as s
import initial_mesh as im
import galerkin as g
import helpers as h

def refine_recursive(mesh, elem, poly, sigma):
    neighbors = [elem.neighbor0, elem.neighbor1, elem.neighbor2]
    FT = neighbors[elem.ET]
    if(FT != None and FT.GT < elem.GT):
        mesh = refine_recursive(mesh, FT, poly, sigma)
    child1, child2 = elem.bisect(h.midpoint(elem.T.v[(elem.ET+1)%3], elem.T.v[(elem.ET+2)%3]), poly)
    #child1, child2 = elem.shape_regularity_bisect(sigma, poly)
    mesh.remove(elem)
    mesh.append(child1)
    mesh.append(child2)
    child1.update_neighbors()
    child2.update_neighbors()
    if FT == None:
        return mesh
    child3, child4 = FT.bisect(h.midpoint(FT.T.v[(FT.ET+1)%3], FT.T.v[(FT.ET+2)%3]), poly)
    #child3, child4 = FT.shape_regularity_bisect(sigma, poly)
    child3.update_neighbors()
    child4.update_neighbors()
    mesh.remove(FT)
    mesh.append(child3)
    mesh.append(child4)
    return mesh

def refine(mesh, poly, sigma):
    marked = mesh.copy()
    for elem in marked:
        if elem.left == None:
            mesh = refine_recursive(mesh, elem, poly, sigma)
    return mesh

poly = [(0,0),(1,0),(1,1),(0,1)]
#poly = [(0,0),(4,0),(4,3),(6,3),(6,4),(3,4),(3,2),(0,2)]
#poly = [(1,0),(2,0),(3,1),(3,2),(2,3),(1,3),(0,2),(0,1)]
mesh = im.create(poly)

sigma = 0
for elem in mesh:
    temp = h.diam_of_triangle(elem.T.v[0],elem.T.v[1],elem.T.v[2])/h.diam_of_inscribed_circle(elem.T.v[0],elem.T.v[1],elem.T.v[2]) 
    if(temp > sigma):
        sigma = temp

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

for k in range(0,0):
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
    mesh = refine(mesh, poly, sigma)

fig, ax = plt.subplots()
for elem in mesh:
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
