import numpy as np
import matplotlib.pyplot as plt

nproc = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 840]);
time  = np.array([294.766, 193.118, 105.906, 55.8185, 28.9824, 14.6696, 7.92044, 4.24703, 2.60059, 1.97075, 1.97305]);

plt.figure(1)
plt.loglog(nproc,time/time[0],'-ob')

plt.figure(2)
plt.semilogx(nproc,time/time[0]*100,'-ob')

plt.figure(3)
plt.plot(nproc,time[0]/(nproc*time),'-ob')

plt.show()
