import matplotlib.pyplot as plt
import earclipping

triangulation = earclipping.triangulate([(0,0),(4,0),(4,3),(6,3),(6,4),(3,4),(3,2),(0,2)])

fig, ax = plt.subplots()
for t in triangulation:
    ax.plot([t.A.x,t.B.x],[t.A.y,t.B.y],color='black',linestyle='-')
    ax.plot([t.B.x,t.C.x],[t.B.y,t.C.y],color='black',linestyle='-')
    ax.plot([t.C.x,t.A.x],[t.C.y,t.A.y],color='black',linestyle='-')

plt.show()
