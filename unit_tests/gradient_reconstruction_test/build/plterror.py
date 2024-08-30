import numpy as np
import matplotlib.pyplot as plt
import sys


s = 0;
conv1 = np.loadtxt('errors_or_'+str(10)+'.dat')
conv2 = np.loadtxt('errors_extended_quad_'+str(30)+'.dat')
dh = [1/10,1/15,1/20,1/25,1/30,1/40]
dhp = [10,15,20,25,30,40]
# dh = [10,15,20]
dudx_error = []
dudy_error = []
dudz_error = []
du2dx2_error = []
du2dxy_error = []
du2dxz_error = []

dudx_error_or_plus_vrt = []
dudy_error_or_plus_vrt = []
dudz_error_or_plus_vrt = []
du2dx2_error_or_plus_vrt = []
du2dxy_error_or_plus_vrt = []
du2dxz_error_or_plus_vrt = []

dudx_error_extended = []
dudy_error_extended = []
dudz_error_extended = []
du2dx2_error_extended = []
du2dxy_error_extended = []
du2dxz_error_extended = []

dudx_error_extended_quad = []
dudy_error_extended_quad = []
dudz_error_extended_quad = []
du2dx2_error_extended_quad = []
du2dxy_error_extended_quad = []
du2dxz_error_extended_quad = []
dhn1=[];
dhq1=[];
dhn2=[];
dhq2=[];
dhn=[];

for i in dhp:
    print('errors'+str(i)+'.dat')
    errors = np.loadtxt('errors_or_'+str(i)+'.dat');
    errors_or_plus_vrt = np.loadtxt('errors_or_plus_vrt_'+str(i)+'.dat');
    errors_extended = np.loadtxt('errors_extended_'+str(i)+'.dat');
    errors_extended_quad = np.loadtxt('errors_extended_quad_'+str(i)+'.dat');

    dudx_error.append(errors[0]);
    dudy_error.append(errors[1]);
    dudz_error.append(errors[2]);
    du2dx2_error.append(errors[3]);
    du2dxy_error.append(errors[4]);
    du2dxz_error.append(errors[5]);

    dudx_error_or_plus_vrt.append(errors_or_plus_vrt[0]);
    dudy_error_or_plus_vrt.append(errors_or_plus_vrt[1]);
    dudz_error_or_plus_vrt.append(errors_or_plus_vrt[2]);
    du2dx2_error_or_plus_vrt.append(errors_or_plus_vrt[3]);
    du2dxy_error_or_plus_vrt.append(errors_or_plus_vrt[4]);
    du2dxz_error_or_plus_vrt.append(errors_or_plus_vrt[5]);


    dudx_error_extended.append(errors_extended[0]);
    dudy_error_extended.append(errors_extended[1]);
    dudz_error_extended.append(errors_extended[2]);
    du2dx2_error_extended.append(errors_extended[3]);
    du2dxy_error_extended.append(errors_extended[4]);
    du2dxz_error_extended.append(errors_extended[5]);

    dudx_error_extended_quad.append(errors_extended_quad[0]);
    dudy_error_extended_quad.append(errors_extended_quad[1]);
    dudz_error_extended_quad.append(errors_extended_quad[2]);
    du2dx2_error_extended_quad.append(errors_extended_quad[3]);
    du2dxy_error_extended_quad.append(errors_extended_quad[4]);
    du2dxz_error_extended_quad.append(errors_extended_quad[5]);

    dhn1.append(conv1[-1]*(dh[s]/dh[-1])**1)
    dhq1.append(conv2[-1]*(dh[s]/dh[-1])**2)

    dhn2.append(conv1[-1]*(dh[s]/dh[-1])**1)
    dhq2.append(conv2[-1]*(dh[s]/dh[-1])**2)

    
    

    s = s + 1

plt.figure(1)
plt.loglog(dh,dudx_error,'-o',label='dudx_error');
plt.loglog(dh,dudy_error,'-o',label='dudy_error');
plt.loglog(dh,dudz_error,'-o',label='dudz_error');
plt.loglog(dh,dudx_error_or_plus_vrt,'-o',label='dudx_error_or_plus_vrt');
plt.loglog(dh,dudy_error_or_plus_vrt,'-o',label='dudy_error_or_plus_vrt');
plt.loglog(dh,dudz_error_or_plus_vrt,'-o',label='dudz_error_or_plus_vrt');
plt.loglog(dh,dudx_error_extended,'-o',label='dudx_error_extended');
plt.loglog(dh,dudy_error_extended,'-o',label='dudy_error_extended');
plt.loglog(dh,dudz_error_extended,'-o',label='dudz_error_extended');
plt.loglog(dh,dudx_error_extended_quad,'--o',label='dudx_error_extended quad');
plt.loglog(dh,dudy_error_extended_quad,'--o',label='dudy_error_extended quad');
plt.loglog(dh,dudz_error_extended_quad,'--o',label='dudz_error_extended quad');
plt.loglog(dh,dhn1,'--k',linewidth=2,label="1nd order ref.")
plt.loglog(dh,dhq1,'--b',linewidth=2,label="2nd order ref.")
# plt.loglog(dh,dhq,'-b')
plt.legend()
plt.grid()



plt.figure(2)
plt.loglog(dh,du2dx2_error,'-o',label='du2dx2_error');
plt.loglog(dh,du2dxy_error,'-o',label='du2dxy_error');
plt.loglog(dh,du2dxz_error,'-o',label='du2dxz_error');
plt.loglog(dh,du2dx2_error_or_plus_vrt,'-o',label='du2dx2_error_or_plus_vrt');
plt.loglog(dh,du2dxy_error_or_plus_vrt,'-o',label='du2dxy_error_or_plus_vrt');
plt.loglog(dh,du2dxz_error_or_plus_vrt,'-o',label='du2dxz_error_or_plus_vrt');
plt.loglog(dh,du2dx2_error_extended,'-o',label='du2dx2_error_extended');
plt.loglog(dh,du2dxy_error_extended,'-o',label='du2dxy_error_extended');
plt.loglog(dh,du2dxz_error_extended,'-o',label='du2dxz_error_extended');
plt.loglog(dh,du2dx2_error_extended_quad,'--o',label='du2dx2_error_extended_quad');
plt.loglog(dh,du2dxy_error_extended_quad,'--o',label='du2dxy_error_extended_quad');
plt.loglog(dh,du2dxz_error_extended_quad,'--o',label='du2dxz_error_extended_quad');
# plt.loglog(dh,dhn,'-k')
plt.loglog(dh,dhn2,'--k',linewidth=2,label="1st order ref.")
plt.loglog(dh,dhq2,'--b',linewidth=2,label="2st order ref.")
plt.legend()
plt.grid()
plt.show()