import matplotlib.pyplot as plt
import triangulation

#initial_mesh = triangulation.earclipping([(0,0),(4,0),(4,3),(6,3),(6,4),(3,4),(3,2),(0,2)])
initial_mesh = triangulation.earclipping([(0,0),(1,0),(1,1),(0,1)])

fig, ax = plt.subplots()
print("v0           v1          v2        l0 l1 l2")
print("-------------------------------------------")
for elem in initial_mesh:
    print("({},{})    ({},{})       ({},{})  {}  {}  {}".format(elem.T.v[0].x, elem.T.v[0].y, elem.T.v[1].x, elem.T.v[1].y, elem.T.v[2].x, elem.T.v[2].y, elem.T.e[0].label, elem.T.e[1].label, elem.T.e[2].label))
    ax.plot([elem.T.v[0].x,elem.T.v[1].x],[elem.T.v[0].y,elem.T.v[1].y],color='black',linestyle='-')
    ax.plot([elem.T.v[1].x,elem.T.v[2].x],[elem.T.v[1].y,elem.T.v[2].y],color='black',linestyle='-')
    ax.plot([elem.T.v[2].x,elem.T.v[0].x],[elem.T.v[2].y,elem.T.v[0].y],color='black',linestyle='-')

plt.show()
