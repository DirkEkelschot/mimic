import numpy as np
import matplotlib.pyplot as plt
import sys


s = 0;
conv = np.loadtxt('errors'+str(10)+'.dat')
dh = [10,15,20,25,30,40,50]
# dh = [10,15,20]
dudx_error = []
dudy_error = []
dudz_error = []
du2dx2_error = []
du2dxy_error = []
du2dxz_error = []

dudx_error_old = []
dudy_error_old = []
dudz_error_old = []
du2dx2_error_old = []
du2dxy_error_old = []
du2dxz_error_old = []

for i in dh:
    print('errors'+str(i)+'.dat')
    errors = np.loadtxt('errors'+str(i)+'.dat');
    errors_old = np.loadtxt('errors_old'+str(i)+'.dat');

    dudx_error.append(errors[0]);
    dudy_error.append(errors[1]);
    dudz_error.append(errors[2]);
    du2dx2_error.append(errors[3]);
    du2dxy_error.append(errors[4]);
    du2dxz_error.append(errors[5]);

    dudx_error_old.append(errors_old[0]);
    dudy_error_old.append(errors_old[1]);
    dudz_error_old.append(errors_old[2]);
    du2dx2_error_old.append(errors_old[3]);
    du2dxy_error_old.append(errors_old[4]);
    du2dxz_error_old.append(errors_old[5]);
    print(errors)
    s = s + 1

plt.figure(1)
plt.loglog(dh,dudx_error,'-o',label='dudx_error');
plt.loglog(dh,dudy_error,'-o',label='dudy_error');
plt.loglog(dh,dudz_error,'-o',label='dudz_error');
plt.loglog(dh,dudx_error_old,'-o',label='dudx_error_old');
plt.loglog(dh,dudy_error_old,'-o',label='dudy_error_old');
plt.loglog(dh,dudz_error_old,'-o',label='dudz_error_old');
plt.legend()
plt.grid()



plt.figure(2)
plt.loglog(dh,du2dx2_error,'-o',label='du2dx2_error');
plt.loglog(dh,du2dxy_error,'-o',label='du2dxy_error');
plt.loglog(dh,du2dxz_error,'-o',label='du2dxz_error');
plt.loglog(dh,du2dx2_error_old,'-o',label='du2dx2_error_old');
plt.loglog(dh,du2dxy_error_old,'-o',label='du2dxy_error_old');
plt.loglog(dh,du2dxz_error_old,'-o',label='du2dxz_error_old');
plt.legend()
plt.grid()
plt.show()