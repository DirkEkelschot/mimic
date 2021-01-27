import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


filename = 'boundary_nodes.dat'
with open(filename) as f:
    A = f.readlines()[0:]
    
L = len(A)
x = np.zeros((L,1))
y = np.zeros((L,1))
z = np.zeros((L,1))
for i in range(0,L):
    s=A[i][:].split()
    x[i] = float(s[0])
    y[i] = float(s[1])
    z[i] = float(s[2])
print(L)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x[0:10000], y[0:10000], z[0:10000], c='r', marker='.')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
