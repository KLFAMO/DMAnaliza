import matplotlib.pyplot as plt
import numpy as np
from parameters import inf, labs


R = 6378137.0
u = np.linspace(0, 2*np.pi, 100)
v = np.linspace(0, np.pi, 100)

xs = R*np.outer(np.cos(u), np.sin(v))
ys = R*np.outer(np.sin(u), np.sin(v))
zs = R*np.outer(np.ones_like(u), np.cos(v))


fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs, color="lightblue", alpha=0.2, linewidth=0)

ax.scatter(inf['UMK1']['X'], inf['UMK1']['Y'], inf['UMK1']['Z'], label='UMK1', color='lightgreen', marker="*", s=40)
ax.scatter(inf['UMK2']['X'], inf['UMK2']['Y'], inf['UMK2']['Z'], label='UMK2', color='lightgreen', marker="*", s=40)

ax.scatter(inf['NPLSr']['X'], inf['NPLSr']['Y'], inf['NPLSr']['Z'], label='NPLSr', color='gray')
ax.scatter(inf['NPLYb']['X'], inf['NPLYb']['Y'], inf['NPLYb']['Z'], label='NPLYb', color='gray')

#exclude = []
#exclude = ['UMK1', 'UMK2']
exclude = ['UMK1', 'UMK2','NPLSr','NPLYb','NICT','NMIJ']

for lab in labs:
    if lab in exclude:
        pass
    else:
        ax.scatter(inf[lab]['X'], inf[lab]['Y'], inf[lab]['Z'], label=lab), #color=inf[lab]['col'])
    

ax.scatter(inf['NICT']['X'], inf['NICT']['Y'], inf['NICT']['Z'], label='NICT', color='pink')
ax.scatter(inf['NMIJ']['X'], inf['NMIJ']['Y'], inf['NMIJ']['Z'], label='NMIJ', color='magenta')

L = R*1.2
ax.plot([-L,L],[0,0],[0,0], color='black')
ax.plot([0,0],[-L,L],[0,0], color='black')
ax.plot([0,0],[0,0],[-L,L], color='black')

ax.text(L,0,0,"X")
ax.text(0,L,0,"Y")
ax.text(0,0,L,"Z")

ax.set_xlim(-L,L)
ax.set_ylim(-L,L)
ax.set_zlim(-L,L)

ax.xaxis.pane.fill = False
#ax.view_init(elev=20, azim=30, roll=10)
ax.view_init(elev=90, azim=0, roll=0)

ax.set_box_aspect([1,1,1])
ax.grid(False)
plt.legend()

plt.show()