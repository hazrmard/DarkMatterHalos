from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
r = 2.0
ev = [2.0,1.0,1.0]
x = (ev[0]*r/max(ev)) * np.outer(np.cos(u), np.sin(v))
y = (ev[1]*r/max(ev)) * np.outer(np.sin(u), np.sin(v))
z = (ev[2]*r/max(ev)) * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.5)
ax.set_xlim3d(-2, 2)
ax.set_ylim3d(-2, 2)
ax.set_zlim3d(-2, 2)

plt.show()
