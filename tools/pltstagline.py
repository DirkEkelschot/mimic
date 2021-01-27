import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
with open('wvars_adapted_stagline2.dat') as f:
    A = f.readlines()[0:]
  
var = []
lines = [];
of=16; #skipping header

for i in range(0,10):
    vn=[]
    for j in range(0,37):
        s=A[of+j][:].split()
        print(s)
        for k in range(0,len(s)):
            vn.append(float(s[k]))
    of=of+37;
    var.append(vn);
print(var[4])
for i in range(of,len(A)):
    line=[]
    s=A[i][:].split()
    for k in range(0,len(s)):
        line.append(int(s[k]))
    lines.append(line)
    print(s)
    
l1,=plt.plot([var[0][lines[0][0]-1],var[0][lines[0][1]-1]],[var[4][lines[0][0]-1],var[4][lines[0][1]-1]],'-b',linewidth=2,label='adapted')
plt.legend(handles=[l1])
for i in range(0,int(len(lines)/2-20)):
    id0 = lines[i][0]
    id1 = lines[i][1]
    plt.plot([var[0][id0-1],var[0][id1-1]],[var[4][id0-1],var[4][id1-1]],'-b')
#===========================================================================
#===========================================================================
#===========================================================================
#===========================================================================
with open('wvars_original_stagline_original.dat') as f:
    A = f.readlines()[0:]
var = []
lines = [];
of=16; #skipping header
nl = 30;
for i in range(0,10):
    vn=[]
    for j in range(0,nl):
        s=A[of+j][:].split()
        print(s)
        for k in range(0,len(s)):
            vn.append(float(s[k]))
    of=of+nl;
    var.append(vn);
print(var[4])
for i in range(of,len(A)):
    line=[]
    s=A[i][:].split()
    for k in range(0,len(s)):
        line.append(int(s[k]))
    lines.append(line)
    print(s)
    
print(var[0])

l2,=plt.plot([var[0][lines[0][0]-1],var[0][lines[0][1]-1]],[var[4][lines[0][0]-1],var[4][lines[0][1]-1]],'-r',linewidth=2,label='original')
plt.legend(handles=[l2])
for i in range(0,len(lines)):
    id0 = lines[i][0]
    id1 = lines[i][1]
    plt.plot([var[0][id0-1],var[0][id1-1]],[var[4][id0-1],var[4][id1-1]],'-r')
plt.grid()

plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
plt.show()
