import numpy as np
import matplotlib.pyplot as plt
import sys

conv = np.loadtxt('convergence.dat')
#print(len(conv[:,0]),len(conv[0,:]))
#l1,=plt.loglog(conv[:,11-3],conv[:,1],'-r',label="res vs #steps")
#l1,=plt.loglog(conv[:,0],conv[:,1],'-b',label="res vs time")
#
#plt.legend(handles=[l1])
#
#plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
#plt.grid()




fig, ax1 = plt.subplots()
t = conv[:,1];
#for i in range(0,len(t)):
#    t[i] = t[i]- t[0]
data1=conv[:,11-3]
color = 'tab:red'
ax1.set_xlabel('time (s)', color=color)
ax1.set_ylabel('residual')
ax1.loglog(data1,t, color=color)
ax1.tick_params(axis='x', labelcolor=color)

ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
data2=conv[:,0]
color = 'tab:blue'
ax2.set_xlabel('# steps', color=color)  # we already handled the x-label with ax1
ax2.loglog(data2,t, color=color)
ax2.set_ylabel('residual')
ax2.tick_params(axis='x', labelcolor=color)


fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()

