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
"""
p0 = s.Point(0.5,0.5)
p1 = s.Point(1,0)
p1.boundary = True
p2 = s.Point(1,1)
p2.boundary = True
p3 = s.Point(0,1)
p3.boundary = True
p4 = s.Point(0,0)
p4.boundary = True
points = [p0, p1, p2, p3, p4]
root0 = s.Node(s.Triangle(points[0], points[1], points[2], s.Edge(points[1], points[2]), s.Edge(points[0], points[2]), s.Edge(points[0], points[1])), 0)
root1 = s.Node(s.Triangle(points[0], points[4], points[1], s.Edge(points[1], points[4]), s.Edge(points[0], points[1]), s.Edge(points[0], points[4])), 0)
root2 = s.Node(s.Triangle(points[0], points[3], points[4], s.Edge(points[3], points[2]), s.Edge(points[0], points[4]), s.Edge(points[0], points[3])), 0)
root3 = s.Node(s.Triangle(points[0], points[2], points[3], s.Edge(points[3], points[2]), s.Edge(points[0], points[3]), s.Edge(points[0], points[2])), 0)
root0.root = root0
root1.root = root1
root2.root = root2
root3.root = root3
root0.T.label_longest_edge()
root1.T.label_longest_edge()
root2.T.label_longest_edge()
root3.T.label_longest_edge()
root0.find_refinement_edge()
root1.find_refinement_edge()
root2.find_refinement_edge()
root3.find_refinement_edge()

mesh = [root0, root1, root2, root3]
for i in range(0, len(mesh)):
    for j in range(0, len(mesh)):
        if(i != j):
            mesh[i].update_neighbor(mesh[j])
"""
sigma = 0
for elem in mesh:
    temp = h.diam_of_triangle(elem.T.v[0],elem.T.v[1],elem.T.v[2])/h.diam_of_inscribed_circle(elem.T.v[0],elem.T.v[1],elem.T.v[2]) 
    if(temp > sigma):
        sigma = temp

vertices = g.unique_vertices_excluding_boundary(mesh)
U = g.galerkin_basis_coefficients(mesh, vertices, 4)
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.5, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.5,0.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.25))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.2, 0.2, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.2,0.2))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.75, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.75,0.75))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.7, 0.7, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.7,0.7))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0,0, 0.0, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.0,0.0))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(1.5, 1.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(1.5,1.5))))
print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0.25))))

for k in range(0,0):
    mesh = refine(mesh, poly, sigma)
    vertices = g.unique_vertices_excluding_boundary(mesh)
    for v in vertices:
        print("({},{})".format(v.x,v.y))
    U = g.galerkin_basis_coefficients(mesh, vertices, 4)
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.5, 0.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.5,0.5))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.25, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.25,0.25))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.2, 0.2, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.2,0.2))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.75, 0.75, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.75,0.75))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0.7, 0.7, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.7,0.7))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0,0, 0.0, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0.0,0.0))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(1.5, 1.5, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(1.5,1.5))))
    print("Galerkin solution at x = ({},{}) with f = 4: {}".format(0, 0.25, g.recreate_galerkin_solution_at_x(mesh, U, vertices, s.Point(0,0.25))))

fig, ax = plt.subplots()
for elem in mesh:
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
