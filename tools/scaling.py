import numpy as np
import matplotlib.pyplot as plt

nproc = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2056, 4000]);
nproc = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2056]);
nproc = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]);
time  = np.array([161.773,81.6285,40.78,20.5611,10.7906,5.73447,3.38092,2.09747,1.56928,1.48158,7.64706,37.603]);
time  = np.array([161.773,81.6285,40.78,20.5611,10.7906,5.73447,3.38092,2.09747,1.56928,1.48158,7.64706]);
time  = np.array([161.773,81.6285,40.78,20.5611,10.7906,5.73447,3.38092,2.09747,1.56928,1.48158]);

plt.figure(1)
plt.loglog(nproc/2,time/time[0],'-ok')
plt.ylabel("$t/t_{np=2}$")
plt.xlabel("np")
plt.grid()
plt.figure(2)
plt.loglog(nproc,time,'-ok')
plt.ylabel("$t [s]$")
plt.xlabel("np [-]")
plt.grid()
#plt.figure(2)
#plt.semilogx(nproc,time/time[0]*100,'-ok')
#plt.ylabel("$t/t_{np=2}*100$")
#plt.xlabel("nproc")
#plt.figure(3)
#plt.plot(nproc,time[0]/(nproc*time),'-ok')
#plt.ylabel("$t_{np=2}/(np*t)$")
#plt.xlabel("nproc")
plt.show()
