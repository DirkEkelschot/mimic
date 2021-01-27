import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
import scipy.linalg

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p=np.zeros((7*3,1))
M = np.zeros((7,1))
p[3*0+0] = 0;
p[3*0+1] = 0;
p[3*0+2] = 0;
M[0] = 2.0;
ax.scatter(p[3*0+0], p[3*0+1], p[3*0+2], marker="o")





p[3*1+0] = 1;
p[3*1+1] = 0.01;
p[3*1+2] = -0.1;
M[1] = 1.0;
ax.scatter(p[3*1+0], p[3*1+1], p[3*1+2], marker="o")
p[3*2+0] = 0.01;
p[3*2+1] = 1.01;
p[3*2+2] = 0.1;
M[2] = 1.0;
ax.scatter(p[3*2+0], p[3*2+1], p[3*2+2], marker="o")
p[3*3+0] = 0.01;
p[3*3+1] = 0.01;
p[3*3+2] = 1.1;
M[3] = 1.0;
ax.scatter(p[3*3+0], p[3*3+1], p[3*3+2], marker="o")
p[3*4+0] = -1.01;
p[3*4+1] = 0.01;
p[3*4+2] = 0.1;
M[4] = 1.0;
ax.scatter(p[3*4+0], p[3*4+1], p[3*4+2], marker="o")
p[3*5+0] = 0.01;
p[3*5+1] = -1.01;
p[3*5+2] = 0.1;
M[5] = 1.0;
ax.scatter(p[3*5+0], p[3*5+1], p[3*5+2], marker="o")
p[3*6+0] = 0.01;
p[3*6+1] = 0.01;
p[3*6+2] = -1.1;
M[6] = 1.0;
ax.scatter(p[3*6+0], p[3*6+1], p[3*6+2], marker="o")

x0 = 0
y0 = 0
z0 = 0

dX=np.zeros((6,3))

for i in range(1,7):
    x = p[3*i+0];
    y = p[3*i+1];
    z = p[3*i+2];
    dx = x-x0
    dy = y-y0
    dz = z-z0
    dX[i-1,0] = dx;
    dX[i-1,1] = dy;
    dX[i-1,2] = dz;

dXnew = np.zeros((4,3))
dXnew[0,0] = -1.0;dXnew[0,1] = 0.0;dXnew[0,2] = 1.0;
dXnew[1,0] = 0.0;dXnew[1,1] = -1.0;dXnew[1,2] = 1.0;
dXnew[2,0] = 1.0;dXnew[2,1] = 0.0;dXnew[2,2] = 1.0;
dXnew[3,0] = 0.0;dXnew[3,1] = 1.0;dXnew[3,2] = 1.0;

b = np.zeros((4,1))
b[0] = -100;b[1] = -100;b[2] = 100;b[3] = 100;

print("b = ", b)
Q, R = np.linalg.qr(dXnew)
print(R)
p = np.dot(Q.T, b)

Rinv = np.linalg.inv(R)
print("Rinv*p = ", np.dot(Rinv,p))
result = np.matmul(Rinv,p)
print("result = ", result)
fin = np.matmul(dXnew,result)
print(fin)
LQ = 0.0
for i in range(0,4):
    print("loc_error ", (fin[i]-b[i])*(fin[i]-b[i]))
    LQ = LQ + (fin[i]-b[i])*(fin[i]-b[i]);

print("LQ=",LQ, "sqrt(LQ)",np.sqrt(LQ))
