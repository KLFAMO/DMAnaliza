import matplotlib.pyplot as plt
import numpy as np


labs = ['UMK1','UMK2', 'NIST', 'NPLSr', 'NPLYb', 'NICT', 'SYRTE']
#labs = ['UMK1', 'UMK2', 'NICT']
inf = { 'UMK1': {'col':'green', 'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736}, 
        'UMK2': {'col':'red',   'atom':'88Sr',
                 'X':3644273,  'Y':1226649,  'Z':5071736},
        'NIST': {'col':'blue',  'atom':'171Yb',
                 'X':-1288363, 'Y':-4721684, 'Z':4078659},
        'NPLSr':{'col':'cyan',  'atom':'87Sr',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NPLYb':{'col':'black', 'atom':'171Yb+',
                 'X':3985500,  'Y':-23625,   'Z':4962941},
        'NICT': {'col':'gray',  'atom':'87Sr',
                 'X':-3941931, 'Y':3368182,  'Z':3702068},
        'SYRTE':{'col':'brown', 'atom':'87Sr',
                 'X':4202777,  'Y':171368,   'Z':4778660}
}

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 5300000 * np.outer(np.cos(u), np.sin(v))
y = 5300000 * np.outer(np.sin(u), np.sin(v))
z = 5300000 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')


for lab in labs:
    ax.scatter(inf[lab]['X'], inf[lab]['Y'], inf[lab]['Z'])
    #ax.quiver(0,0,0,inf[lab]['X'], inf[lab]['Y'], inf[lab]['Z'])

plt.show()