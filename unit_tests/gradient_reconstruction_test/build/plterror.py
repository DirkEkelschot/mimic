import numpy as np
import matplotlib.pyplot as plt
import sys


s = 0;
conv1 = np.loadtxt('errors_or_'+str(10)+'.dat')
conv2 = np.loadtxt('errors_extended_quad_'+str(10)+'.dat')
# dh = [1/10,1/15,1/20,1/25,1/30,1/40,1/50]
# dhp = [10,15,20,25,30,40,50]

dh = [1/10,1/15,1/20,1/25,1/30,1/40,1/50]
dhp = [10,15,20,25,30,40,50]
# dh = [10,15,20]
dudx_error = []
dudy_error = []
dudz_error = []
du2dx2_error = []
du2dxy_error = []
du2dxz_error = []
du3dx3_error = []
du3dz3_error = []

dudx_error_or_plus_vrt = []
dudy_error_or_plus_vrt = []
dudz_error_or_plus_vrt = []
du2dx2_error_or_plus_vrt = []
du2dxy_error_or_plus_vrt = []
du2dxz_error_or_plus_vrt = []
du3dx3_error_or_plus_vrt = []
du3dz3_error_or_plus_vrt = []

dudx_error_extended = []
dudy_error_extended = []
dudz_error_extended = []
du2dx2_error_extended = []
du2dxy_error_extended = []
du2dxz_error_extended = []
du3dx3_error_extended = []
du3dz3_error_extended = []

dudx_error_extended_quad = []
dudy_error_extended_quad = []
dudz_error_extended_quad = []
du2dx2_error_extended_quad = []
du2dxy_error_extended_quad = []
du2dxz_error_extended_quad = []
dhz1=[];
dhn1=[];
dhq1=[];
dhc1=[];

dhz2=[];
dhn2=[];
dhq2=[];
dhc2=[];
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
    # du3dx3_error.append(errors[6]);
    # du3dz3_error.append(errors[7]);

    dudx_error_or_plus_vrt.append(errors_or_plus_vrt[0]);
    dudy_error_or_plus_vrt.append(errors_or_plus_vrt[1]);
    dudz_error_or_plus_vrt.append(errors_or_plus_vrt[2]);
    du2dx2_error_or_plus_vrt.append(errors_or_plus_vrt[3]);
    du2dxy_error_or_plus_vrt.append(errors_or_plus_vrt[4]);
    du2dxz_error_or_plus_vrt.append(errors_or_plus_vrt[5]);
    # du3dx3_error_or_plus_vrt.append(errors_or_plus_vrt[6]);
    # du3dz3_error_or_plus_vrt.append(errors_or_plus_vrt[7]);


    dudx_error_extended.append(errors_extended[0]);
    dudy_error_extended.append(errors_extended[1]);
    dudz_error_extended.append(errors_extended[2]);
    du2dx2_error_extended.append(errors_extended[3]);
    du2dxy_error_extended.append(errors_extended[4]);
    du2dxz_error_extended.append(errors_extended[5]);
    # du3dx3_error_extended.append(errors_extended[6]);
    # du3dz3_error_extended.append(errors_extended[7]);

    dudx_error_extended_quad.append(errors_extended_quad[0]);
    dudy_error_extended_quad.append(errors_extended_quad[1]);
    dudz_error_extended_quad.append(errors_extended_quad[2]);
    du2dx2_error_extended_quad.append(errors_extended_quad[3]);
    du2dxy_error_extended_quad.append(errors_extended_quad[4]);
    du2dxz_error_extended_quad.append(errors_extended_quad[5]);
    print(conv1)
    print(conv1[-1],conv1[0],conv2[0],i,dh[s],dh[-1],"error ",conv1[0]*(dh[s]/dh[0])**1,dh[s],dh[0],dudx_error_extended_quad[0])

    dhz1.append(20*conv2[0]*(dh[s]/dh[0])**0)
    dhn1.append(20*conv2[0]*(dh[s]/dh[0])**1)
    dhq1.append(20*conv2[0]*(dh[s]/dh[0])**2)
    dhc1.append(20*conv2[0]*(dh[s]/dh[0])**3)

    dhz2.append(conv2[4]*(dh[s]/dh[0])**0)
    dhn2.append(conv2[4]*(dh[s]/dh[0])**1)
    dhq2.append(conv2[4]*(dh[s]/dh[0])**2)
    dhc2.append(conv2[4]*(dh[s]/dh[0])**3)
    
    

    s = s + 1

plt.figure(1,figsize=(10,8))
plt.loglog(dh,dudx_error,'-o',label='dudx LWLSQ');
plt.loglog(dh,dudy_error,'-o',label='dudy LWLSQ');
plt.loglog(dh,dudz_error,'-o',label='dudz LWLSQ');
plt.loglog(dh,dudx_error_or_plus_vrt,'-o',label='dudx LWLSQ+V');
plt.loglog(dh,dudy_error_or_plus_vrt,'-o',label='dudy LWLSQ+V');
plt.loglog(dh,dudz_error_or_plus_vrt,'-o',label='dudz LWLSQ+V');
# plt.loglog(dh,dudx_error_extended,'-o',label='dudx QWLSQ');
# plt.loglog(dh,dudy_error_extended,'-o',label='dudy QWLSQ');
# plt.loglog(dh,dudz_error_extended,'-o',label='dudz QWLSQ');
plt.loglog(dh,dudx_error_extended_quad,'--o',label='dudx QWLSQ');
plt.loglog(dh,dudy_error_extended_quad,'--o',label='dudy QWLSQ');
plt.loglog(dh,dudz_error_extended_quad,'--o',label='dudz QWLSQ');
plt.loglog(dh,dhn1,'--k',label="1st order ref.",linewidth=3.0)
plt.loglog(dh,dhq1,'--b',label="2nd order ref.",linewidth=3.0)
plt.loglog(dh,dhc1,'--r',linewidth=2,label="3rd order ref.")
plt.xlabel(r"$\frac{1}{h}$ [-]",fontsize=18)
plt.ylabel(r"$\|U_{rec}-U_{exact}\|_{L_2}$ [-]",fontsize=18)
# plt.loglog(dh,dhq,'-b')
plt.legend()
filename = 'dudxi_error'
plt.savefig(filename+'.png')
# plt.grid()



plt.figure(2,figsize=(10,8))
plt.loglog(dh,du2dx2_error,'-o',label='du2dx2 recursive LWLSQ');
plt.loglog(dh,du2dxy_error,'-o',label='du2dy2 recursive LWLSQ');
plt.loglog(dh,du2dxz_error,'-o',label='du2dz2 recursive LWLSQ');
plt.loglog(dh,du2dx2_error_or_plus_vrt,'-o',label='du2dx2 recursive LWLSQ+V');
plt.loglog(dh,du2dxy_error_or_plus_vrt,'-o',label='du2dy2 recursive LWLSQ+V');
plt.loglog(dh,du2dxz_error_or_plus_vrt,'-o',label='du2dz2 recursive LWLSQ+V');
plt.loglog(dh,du2dx2_error_extended,'-o',label='du2dx2 recursive QWLSQ');
plt.loglog(dh,du2dxy_error_extended,'-o',label='du2dy2 recursive QWLSQ');
plt.loglog(dh,du2dxz_error_extended,'-o',label='du2dz2 recursive QWLSQ');
plt.loglog(dh,du2dx2_error_extended_quad,'--o',label='du2dx2 QWLSQ');
plt.loglog(dh,du2dxy_error_extended_quad,'--o',label='du2dy2 QWLSQ');
plt.loglog(dh,du2dxz_error_extended_quad,'--o',label='du2dz2 QWLSQ');
# plt.loglog(dh,dhn,'-k')
plt.loglog(dh,dhz2,'--k',label="0th order ref.",linewidth=3.0)
plt.loglog(dh,dhn2,'--b',label="1st order ref.",linewidth=3.0)
plt.loglog(dh,dhc2,'--r',linewidth=2,label="3rd order ref.")
plt.xlabel(r"$\frac{1}{h}$ [-]",fontsize=18)
plt.ylabel(r"$\|U_{rec}-U_{exact}\|_{L_2}$ [-]",fontsize=18)
plt.legend()
filename = 'du2dxi2_error'
plt.savefig(filename+'.png')


# plt.figure(3,figsize=(10,8))
# plt.loglog(dh,du3dx3_error,'-o',label='du3dx3 recursive LWLSQ');
# plt.loglog(dh,du3dz3_error,'-o',label='du3dz3 recursive LWLSQ');
# plt.loglog(dh,du3dx3_error_or_plus_vrt,'-o',label='du3dx3 recursive LWLSQ+V');
# plt.loglog(dh,du3dz3_error_or_plus_vrt,'-o',label='du3dz3 recursive LWLSQ+V');
# plt.loglog(dh,du3dx3_error_extended,'-o',label='du3dx3 recursive QWLSQ');
# plt.loglog(dh,du3dz3_error_extended,'-o',label='du3dz3 recursive QWLSQ');
# # plt.loglog(dh,dhn,'-k')
# plt.loglog(dh,dhz2,'--k',label="0th order ref.",linewidth=3.0)
# plt.loglog(dh,dhn2,'--b',label="1st order ref.",linewidth=3.0)
# # plt.loglog(dh,dhc2,'--r',linewidth=2,label="3rd order ref.")
# plt.xlabel(r"$\frac{1}{h}$ [-]",fontsize=18)
# plt.ylabel(r"$\|U_{rec}-U_{exact}\|_{L_2}$ [-]",fontsize=18)
# plt.legend()
# filename = 'du3dxi3_error'
# plt.savefig(filename+'.png')
# plt.grid()
plt.show()