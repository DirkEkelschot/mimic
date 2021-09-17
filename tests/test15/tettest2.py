import numpy as np;
import matplotlib.pyplot as plt;


def PlotTriangle(f0,ax,c):
    ax.plot([f0[0][0],f0[1][0]],[f0[0][1],f0[1][1]],[f0[0][2],f0[1][2]],c);
    ax.plot([f0[1][0],f0[2][0]],[f0[1][1],f0[2][1]],[f0[1][2],f0[2][2]],c);
    ax.plot([f0[2][0],f0[0][0]],[f0[2][1],f0[0][1]],[f0[2][2],f0[0][2]],c);
    
def PlotQuad(f0,ax,c):
    ax.plot([f0[0][0],f0[1][0]],[f0[0][1],f0[1][1]],[f0[0][2],f0[1][2]],c);
    ax.plot([f0[1][0],f0[2][0]],[f0[1][1],f0[2][1]],[f0[1][2],f0[2][2]],c);
    ax.plot([f0[2][0],f0[3][0]],[f0[2][1],f0[3][1]],[f0[2][2],f0[3][2]],c);
    ax.plot([f0[3][0],f0[0][0]],[f0[3][1],f0[0][1]],[f0[3][2],f0[0][2]],c);


p0 = [-0.678234, 0.212849, 0.0471629]
p1 = [-0.675937, 0.224903, 0.0682225]
p2 = [-0.671863, 0.23039, 0.0323036];
p3 = [-1.323, -0.0751015, 0.0406843];

f0 = [p1,p3,p2]
f1 = [p0,p3,p2]
f2 = [p0,p1,p3]
f3 = [p0,p1,p2]

fig = plt.figure()
ax = plt.axes(projection='3d')
c = '-r';
PlotTriangle(f0,ax,'-r')
PlotTriangle(f1,ax,'-r')
PlotTriangle(f2,ax,'-r')
PlotTriangle(f3,ax,'-r')
plt.show()
