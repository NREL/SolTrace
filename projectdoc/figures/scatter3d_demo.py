import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
n = 50
#for c, m, zl, zh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:

c,m,zl,zh = ('k','s',0,10)
xs = randrange(n, 0, 100)
ys = randrange(n, 0, 100)
zs = randrange(n, zl, zh)

ax.set_axis_on()
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.scatter(xs, ys, zs, marker=m, s=100)

#ax.set_zlim(-50,50)
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

ax = fig.add_subplot(122, projection='3d')
ax.set_axis_off()
ax.scatter(xs, ys, zs, marker=m, s=100)

fig.subplots_adjust(hspace=0.02, wspace=0.02, left=0.01, right=0.99, top=0.99, bottom=0.01)

plt.show()
