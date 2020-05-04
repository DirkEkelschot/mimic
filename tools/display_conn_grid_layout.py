
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def ComputeNormal(face):
    print(face)
    P1x=face[0,0];P1y=face[0,1];P1z=face[0,2];
    P2x=face[1,0];P2y=face[1,1];P2z=face[1,2];
    P3x=face[2,0];P3y=face[2,1];P3z=face[2,2];
    P4x=face[3,0];P4y=face[3,1];P4z=face[3,2];
    
    QP = np.zeros((3,1))
    QS = np.zeros((3,1))
    QP[0] = P2x-P1x;
    QP[1] = P2y-P1y;
    QP[2] = P2z-P1z;
    
    QS[0] = P4x-P1x;
    QS[1] = P4y-P1y;
    QS[2] = P4z-P1z;
    print("QP=",QP)
    print("QS=",QS)
    N=np.zeros((3,1))
    N[0] =   QP[1]*QS[2]-QS[1]*QP[2];
    N[1] = -(QP[0]*QS[2]-QS[0]*QP[2]);
    N[2] =   QP[0]*QS[1]-QS[0]*QP[1];
    return N;
    
    
'''
L = np.zeros((20,1))

f = h5.File("grid_sm.h5", "r")

datasetNames = [n for n in f.keys()]
for n in datasetNames:
    print(n)
print(" ")
f2 = h5.File("conn_sm.h5", "r")

datasetNames = [n for n in f2.keys()]
for n in datasetNames:
    print(n)
    
    
print("layout of the GRID.H5")
iet = reflectance = f['iet']
print(iet)
ifn = reflectance = f['ifn']
print(ifn)
info = reflectance = f['info']
print(info)
xcn = reflectance = f['xcn']
print(xcn)
zones = reflectance = f['zones']
print(zones)


print("layout of the CONN.H5")
iee = f2['iee']
print(iee)
if isinstance(iee, h5.Group):
    datasetNames = [n for n in iee.keys()]
    for n in datasetNames:
        print(n)


ief = f2['ief']
print(ief)
if isinstance(ief, h5.Group):
    datasetNames = [n for n in ief.keys()]
    for n in datasetNames:
        print(n)
    
    
ien = f2['ien']
print(ien)
if isinstance(ien, h5.Group):
    datasetNames = [n for n in ien.keys()]
    for n in datasetNames:
        print(n)


ife  = f2['ife']
print(ife.value)
if isinstance(ife, h5.Group):
    datasetNames = [n for n in ife.keys()]
    for n in datasetNames:
        print(n)


info = f2['info']
print(info)
if isinstance(info, h5.Group):
    print("yes!")
    datasetNames = [n for n in info.keys()]
    for n in datasetNames:
        print("grid=",n)


print("layout of the GRID.H5")
f3 = h5.File("grid_sm.h5", "r")
SetNames = [n for n in f3.keys()]
for n in SetNames:
    print(f3[n])
    if isinstance(n, h5.Group):
        datasetNames = [q for q in f3[n].keys()]
        for q in datasetNames:
            print(q)
print("layout of the CONN.H5")
f3 = h5.File("conn_sm.h5", "r")
SetNames = [n for n in f3.keys()]
for n in SetNames:
    print(n,f3[n])
    if isinstance(f3[n], h5.Group):
        datasetNames = [q for q in f3[n].keys()]
        for q in datasetNames:
            ds = f3[n][q]
            print("Is it a dataset ->",isinstance(ds, h5.Dataset),ds)

print(ien)
'''

print("========================================================")
print("==================== Print Conn data ===================")
print("========================================================")
filename = '../grids/piston/conn.h5'
data = h5.File(filename, 'r')
for group in data.keys():
    print("DATATA === ", data[group])
    if isinstance(data[group], h5.Dataset):
        ds = data[group]
        
        print("ds - ", ds)
    
    if isinstance(data[group], h5.Group):
       for dset in data[group].keys():
           print("groepjes ",h5.Group,dset)
           if dset == 'zdefs':
              zdefs = data[group][dset]
           if dset == 'znames':
              znames = data[group][dset]
           
           
           #for dset2 in data[group][dset].keys():
           #    print('gaatie2 = ',data[group][dset][dset2])


    else:
        #print('dataset =',data[group],group,len(data[group][:]))
        ds = data[group][:]
        #print(ds)
        if group=='ien':
            ien = ds
        if group=='iee':
            iee = ds
        if group=='ief':
            ief = ds
        if group=='ife':
            ife = ds

#print(znames[:])
#print(zdefs[:])
print(ief)
print(ien)

print("========================================================")
print("==================== Print Grid data ===================")
print("========================================================")
filename = '../grids/piston/grid.h5'
data = h5.File(filename, 'r')
for group in data.keys():
    if isinstance(data[group], h5.Group):
       print("1.=",data[group])
       for dset in data[group].keys():
           print(data[group][dset])
           ds = data[group][dset]
           
           print(data[group],group,ds)
    else:
        #print("print data !!! = ", data[group],group,len(data[group][:]))
        
        ds = data[group][:]
        if group=='iet':
            iet = ds
        if group=='xcn':
            xcn = ds
        if group=='ifn':
            ifn = ds
        if group=='zdefs':
            zdefs = ds

print("========= ifn ==========")
print(ifn)

print("========================================================")
print("==================== Print Data data ===================")
print("========================================================")
filename = 'data.h5.old'
data = h5.File(filename, 'r')
for group in data.keys():
    if isinstance(data[group], h5.Group):
       print("data.h5 => ",data[group])
       print(data[group].keys())
       for dset in data[group].keys():
           if isinstance(data[group][dset], h5.Group):
              print("dset = ", dset)
              print(data[group][dset].keys())
              for dset2 in data[group][dset].keys():
                  print("dset2 = ", dset2)
                  ds = data[group][dset][dset2]
                  if dset2=='boundaries':
                     bnds=ds
                  if dset2=='interior':
                     interior=ds
                  if dset2=='bvnames':
                     bvnames=dset2
           #print("dataset = ",data[group][dset],data[group])
           #ds = data[group][dset]
                   
           #print("new = ",data[group],group,ds)
           
           #print("ds",ds)
    else:
        print("print data !!! = ", data[group],group,len(data[group][:]))
        ds = data[group][:]

print("bnds=",len(bnds[0,:]),len(bnds[:,0]))
print("bnds array",bnds[0,:])
print("bvnames = ",bvnames)

#print("interior=",len(interior[0,:]),len(interior[:,0]))
#Points = np.zeros((8,3));
#print("interior array",interior[:])

#ax.scatter(xcn[:,0],xcn[:,1],xcn[:,2],'.',color='black')
#print("ien ", ien)
'''
for i in range(0,len(ien)):
    if(ien[i,0]==4):
        fig = plt.figure(i)
        ax = Axes3D(fig)
        ax = fig.add_subplot(111, projection='3d')
    
        Points[0,0] = xcn[ien[i,1]-1,0];      Points[0,1] = xcn[ien[i,1]-1,1];      Points[0,2] = xcn[ien[i,1]-1,2];
        Points[1,0] = xcn[ien[i,2]-1,0];      Points[1,1] = xcn[ien[i,2]-1,1];      Points[1,2] = xcn[ien[i,2]-1,2];
        Points[2,0] = xcn[ien[i,3]-1,0];      Points[2,1] = xcn[ien[i,3]-1,1];      Points[2,2] = xcn[ien[i,3]-1,2];
        Points[3,0] = xcn[ien[i,4]-1,0];      Points[3,1] = xcn[ien[i,4]-1,1];      Points[3,2] = xcn[ien[i,4]-1,2];
        Points[4,0] = xcn[ien[i,5]-1,0];      Points[4,1] = xcn[ien[i,5]-1,1];      Points[4,2] = xcn[ien[i,5]-1,2];
        Points[5,0] = xcn[ien[i,6]-1,0];      Points[5,1] = xcn[ien[i,6]-1,1];      Points[5,2] = xcn[ien[i,6]-1,2];
        Points[6,0] = xcn[ien[i,7]-1,0];      Points[6,1] = xcn[ien[i,7]-1,1];      Points[6,2] = xcn[ien[i,7]-1,2];
        Points[7,0] = xcn[ien[i,8]-1,0];      Points[7,1] = xcn[ien[i,8]-1,1];      Points[7,2] = xcn[ien[i,8]-1,2];

        print("iefn ", ief[i,1],ief[i,2],ief[i,3],ief[i,4],ief[i,5],ief[i,6])
        
        idn0_0 = ifn[np.abs(ief[i,1])-1,1]
        idn0_1 = ifn[np.abs(ief[i,1])-1,2]
        idn0_2 = ifn[np.abs(ief[i,1])-1,3]
        idn0_3 = ifn[np.abs(ief[i,1])-1,4]
        
        idn1_0 = ifn[np.abs(ief[i,2])-1,1]
        idn1_1 = ifn[np.abs(ief[i,2])-1,2]
        idn1_2 = ifn[np.abs(ief[i,2])-1,3]
        idn1_3 = ifn[np.abs(ief[i,2])-1,4]
        
        idn2_0 = ifn[np.abs(ief[i,3])-1,1]
        idn2_1 = ifn[np.abs(ief[i,3])-1,2]
        idn2_2 = ifn[np.abs(ief[i,3])-1,3]
        idn2_3 = ifn[np.abs(ief[i,3])-1,4]
        
        idn3_0 = ifn[np.abs(ief[i,4])-1,1]
        idn3_1 = ifn[np.abs(ief[i,4])-1,2]
        idn3_2 = ifn[np.abs(ief[i,4])-1,3]
        idn3_3 = ifn[np.abs(ief[i,4])-1,4]
        
        idn4_0 = ifn[np.abs(ief[i,5])-1,1]
        idn4_1 = ifn[np.abs(ief[i,5])-1,2]
        idn4_2 = ifn[np.abs(ief[i,5])-1,3]
        idn4_3 = ifn[np.abs(ief[i,5])-1,4]
        
        idn5_0 = ifn[np.abs(ief[i,6])-1,1]
        idn5_1 = ifn[np.abs(ief[i,6])-1,2]
        idn5_2 = ifn[np.abs(ief[i,6])-1,3]
        idn5_3 = ifn[np.abs(ief[i,6])-1,4]
        
        face = np.zeros((4,3))
    
        face[0,0] = xcn[idn0_0-1,0];face[0,1] = xcn[idn0_0-1,1];face[0,2] = xcn[idn0_0-1,2];
        face[1,0] = xcn[idn0_1-1,0];face[1,1] = xcn[idn0_1-1,1];face[1,2] = xcn[idn0_1-1,2];
        face[2,0] = xcn[idn0_2-1,0];face[2,1] = xcn[idn0_2-1,1];face[2,2] = xcn[idn0_2-1,2];
        face[3,0] = xcn[idn0_3-1,0];face[3,1] = xcn[idn0_3-1,1];face[3,2] = xcn[idn0_3-1,2];
        N=ComputeNormal(face);
        print("Normal1=",N)
        face[0,0] = xcn[idn1_0-1,0];face[0,1] = xcn[idn1_0-1,1];face[0,2] = xcn[idn1_0-1,2];
        face[1,0] = xcn[idn1_1-1,0];face[1,1] = xcn[idn1_1-1,1];face[1,2] = xcn[idn1_1-1,2];
        face[2,0] = xcn[idn1_2-1,0];face[2,1] = xcn[idn1_2-1,1];face[2,2] = xcn[idn1_2-1,2];
        face[3,0] = xcn[idn1_3-1,0];face[3,1] = xcn[idn1_3-1,1];face[3,2] = xcn[idn1_3-1,2];
        N=ComputeNormal(face);
        print("Normal2=",N)
    
        
    
        ax.scatter(Points[0,0],Points[0,1],Points[0,2],'.',color='magenta')
        ax.scatter(Points[1,0],Points[1,1],Points[1,2],'.',color='green')
        ax.scatter(Points[2,0],Points[2,1],Points[2,2],'.',color='red')
        ax.scatter(Points[3,0],Points[3,1],Points[3,2],'.',color='black')
    
        ax.scatter(Points[4,0],Points[4,1],Points[4,2],'.',color='magenta')
        ax.scatter(Points[5,0],Points[5,1],Points[5,2],'.',color='green')
        ax.scatter(Points[6,0],Points[6,1],Points[6,2],'.',color='red')
        ax.scatter(Points[7,0],Points[7,1],Points[7,2],'.',color='black')
    
    
        plt.plot([xcn[idn0_0-1,0],xcn[idn0_1-1,0]],[xcn[idn0_0-1,1],xcn[idn0_1-1,1]],[xcn[idn0_0-1,2],xcn[idn0_1-1,2]],'--k')
        plt.plot([xcn[idn0_1-1,0],xcn[idn0_2-1,0]],[xcn[idn0_1-1,1],xcn[idn0_2-1,1]],[xcn[idn0_1-1,2],xcn[idn0_2-1,2]],'--k')
        plt.plot([xcn[idn0_2-1,0],xcn[idn0_3-1,0]],[xcn[idn0_2-1,1],xcn[idn0_3-1,1]],[xcn[idn0_2-1,2],xcn[idn0_3-1,2]],'--k')
        plt.plot([xcn[idn0_3-1,0],xcn[idn0_0-1,0]],[xcn[idn0_3-1,1],xcn[idn0_0-1,1]],[xcn[idn0_3-1,2],xcn[idn0_0-1,2]],'--k')
    
        plt.plot([xcn[idn1_0-1,0],xcn[idn1_1-1,0]],[xcn[idn1_0-1,1],xcn[idn1_1-1,1]],[xcn[idn1_0-1,2],xcn[idn1_1-1,2]],'--g')
        plt.plot([xcn[idn1_1-1,0],xcn[idn1_2-1,0]],[xcn[idn1_1-1,1],xcn[idn1_2-1,1]],[xcn[idn1_1-1,2],xcn[idn1_2-1,2]],'--g')
        plt.plot([xcn[idn1_2-1,0],xcn[idn1_3-1,0]],[xcn[idn1_2-1,1],xcn[idn1_3-1,1]],[xcn[idn1_2-1,2],xcn[idn1_3-1,2]],'--g')
        plt.plot([xcn[idn1_3-1,0],xcn[idn1_0-1,0]],[xcn[idn1_3-1,1],xcn[idn1_0-1,1]],[xcn[idn1_3-1,2],xcn[idn1_0-1,2]],'--g')
    
        plt.plot([xcn[idn2_0-1,0],xcn[idn2_1-1,0]],[xcn[idn2_0-1,1],xcn[idn2_1-1,1]],[xcn[idn2_0-1,2],xcn[idn2_1-1,2]],'--b')
        plt.plot([xcn[idn2_1-1,0],xcn[idn2_2-1,0]],[xcn[idn2_1-1,1],xcn[idn2_2-1,1]],[xcn[idn2_1-1,2],xcn[idn2_2-1,2]],'--b')
        plt.plot([xcn[idn2_2-1,0],xcn[idn2_3-1,0]],[xcn[idn2_2-1,1],xcn[idn2_3-1,1]],[xcn[idn2_2-1,2],xcn[idn2_3-1,2]],'--b')
        plt.plot([xcn[idn2_3-1,0],xcn[idn2_0-1,0]],[xcn[idn2_3-1,1],xcn[idn2_0-1,1]],[xcn[idn2_3-1,2],xcn[idn2_0-1,2]],'--b')
    
        plt.plot([xcn[idn3_0-1,0],xcn[idn3_1-1,0]],[xcn[idn3_0-1,1],xcn[idn3_1-1,1]],[xcn[idn3_0-1,2],xcn[idn3_1-1,2]],'--r')
        plt.plot([xcn[idn3_1-1,0],xcn[idn3_2-1,0]],[xcn[idn3_1-1,1],xcn[idn3_2-1,1]],[xcn[idn3_1-1,2],xcn[idn3_2-1,2]],'--r')
        plt.plot([xcn[idn3_2-1,0],xcn[idn3_3-1,0]],[xcn[idn3_2-1,1],xcn[idn3_3-1,1]],[xcn[idn3_2-1,2],xcn[idn3_3-1,2]],'--r')
        plt.plot([xcn[idn3_3-1,0],xcn[idn3_0-1,0]],[xcn[idn3_3-1,1],xcn[idn3_0-1,1]],[xcn[idn3_3-1,2],xcn[idn3_0-1,2]],'--r')
    
        plt.plot([xcn[idn4_0-1,0],xcn[idn4_1-1,0]],[xcn[idn4_0-1,1],xcn[idn4_1-1,1]],[xcn[idn4_0-1,2],xcn[idn4_1-1,2]],'--y')
        plt.plot([xcn[idn4_1-1,0],xcn[idn4_2-1,0]],[xcn[idn4_1-1,1],xcn[idn4_2-1,1]],[xcn[idn4_1-1,2],xcn[idn4_2-1,2]],'--y')
        plt.plot([xcn[idn4_2-1,0],xcn[idn4_3-1,0]],[xcn[idn4_2-1,1],xcn[idn4_3-1,1]],[xcn[idn4_2-1,2],xcn[idn4_3-1,2]],'--y')
        plt.plot([xcn[idn4_3-1,0],xcn[idn4_0-1,0]],[xcn[idn4_3-1,1],xcn[idn4_0-1,1]],[xcn[idn4_3-1,2],xcn[idn4_0-1,2]],'--y')
    
        plt.plot([xcn[idn5_0-1,0],xcn[idn5_1-1,0]],[xcn[idn5_0-1,1],xcn[idn5_1-1,1]],[xcn[idn5_0-1,2],xcn[idn5_1-1,2]],'--m')
        plt.plot([xcn[idn5_1-1,0],xcn[idn5_2-1,0]],[xcn[idn5_1-1,1],xcn[idn5_2-1,1]],[xcn[idn5_1-1,2],xcn[idn5_2-1,2]],'--m')
        plt.plot([xcn[idn5_2-1,0],xcn[idn5_3-1,0]],[xcn[idn5_2-1,1],xcn[idn5_3-1,1]],[xcn[idn5_2-1,2],xcn[idn5_3-1,2]],'--m')
        plt.plot([xcn[idn5_3-1,0],xcn[idn5_0-1,0]],[xcn[idn5_3-1,1],xcn[idn5_0-1,1]],[xcn[idn5_3-1,2],xcn[idn5_0-1,2]],'--m')

        plt.show()
        

    if(ien[i,0]==6):
        fig = plt.figure(i)
        ax = Axes3D(fig)
        ax = fig.add_subplot(111, projection='3d')
        Points[0,0] = xcn[ien[i,1]-1,0];      Points[0,1] = xcn[ien[i,1]-1,1];      Points[0,2] = xcn[ien[i,1]-1,2];
        Points[1,0] = xcn[ien[i,2]-1,0];      Points[1,1] = xcn[ien[i,2]-1,1];      Points[1,2] = xcn[ien[i,2]-1,2];
        Points[2,0] = xcn[ien[i,3]-1,0];      Points[2,1] = xcn[ien[i,3]-1,1];      Points[2,2] = xcn[ien[i,3]-1,2];
        Points[3,0] = xcn[ien[i,4]-1,0];      Points[3,1] = xcn[ien[i,4]-1,1];      Points[3,2] = xcn[ien[i,4]-1,2];
        Points[4,0] = xcn[ien[i,5]-1,0];      Points[4,1] = xcn[ien[i,5]-1,1];      Points[4,2] = xcn[ien[i,5]-1,2];
        Points[5,0] = xcn[ien[i,6]-1,0];      Points[5,1] = xcn[ien[i,6]-1,1];      Points[5,2] = xcn[ien[i,6]-1,2];
        Points[6,0] = xcn[ien[i,7]-1,0];      Points[6,1] = xcn[ien[i,7]-1,1];      Points[6,2] = xcn[ien[i,7]-1,2];
        Points[7,0] = xcn[ien[i,8]-1,0];      Points[7,1] = xcn[ien[i,8]-1,1];      Points[7,2] = xcn[ien[i,8]-1,2];
        print(Points)
        
        ax.scatter(Points[0,0],Points[0,1],Points[0,2],'.',color='magenta')
        ax.scatter(Points[1,0],Points[1,1],Points[1,2],'.',color='green')
        ax.scatter(Points[2,0],Points[2,1],Points[2,2],'.',color='red')
        
        ax.scatter(Points[3,0],Points[3,1],Points[3,2],'.',color='magenta')
        ax.scatter(Points[4,0],Points[4,1],Points[4,2],'.',color='green')
        ax.scatter(Points[5,0],Points[5,1],Points[5,2],'.',color='red')
        
        idn0_0 = ifn[np.abs(ief[i,1])-1,1]
        idn0_1 = ifn[np.abs(ief[i,1])-1,2]
        idn0_2 = ifn[np.abs(ief[i,1])-1,3]
        
        idn1_0 = ifn[np.abs(ief[i,2])-1,1]
        idn1_1 = ifn[np.abs(ief[i,2])-1,2]
        idn1_2 = ifn[np.abs(ief[i,2])-1,3]
        
        
        plt.plot([xcn[idn0_0-1,0],xcn[idn0_1-1,0]],[xcn[idn0_0-1,1],xcn[idn0_1-1,1]],[xcn[idn0_0-1,2],xcn[idn0_1-1,2]],'--k')
        plt.plot([xcn[idn0_1-1,0],xcn[idn0_2-1,0]],[xcn[idn0_1-1,1],xcn[idn0_2-1,1]],[xcn[idn0_1-1,2],xcn[idn0_2-1,2]],'--k')
        plt.plot([xcn[idn0_2-1,0],xcn[idn0_0-1,0]],[xcn[idn0_2-1,1],xcn[idn0_0-1,1]],[xcn[idn0_2-1,2],xcn[idn0_0-1,2]],'--k')
        
        plt.plot([xcn[idn1_0-1,0],xcn[idn1_1-1,0]],[xcn[idn1_0-1,1],xcn[idn1_1-1,1]],[xcn[idn1_0-1,2],xcn[idn1_1-1,2]],'--k')
        plt.plot([xcn[idn1_1-1,0],xcn[idn1_2-1,0]],[xcn[idn1_1-1,1],xcn[idn1_2-1,1]],[xcn[idn1_1-1,2],xcn[idn1_2-1,2]],'--k')
        plt.plot([xcn[idn1_2-1,0],xcn[idn1_0-1,0]],[xcn[idn1_2-1,1],xcn[idn1_0-1,1]],[xcn[idn1_2-1,2],xcn[idn1_0-1,2]],'--k')
        
        
        #plt.plot([xcn[idn1_0-1,0],xcn[idn1_1-1,0]],[xcn[idn1_0-1,1],xcn[idn1_1-1,1]],[xcn[idn1_0-1,2],xcn[idn1_1-1,2]],'--k')
        #plt.plot([xcn[idn1_1-1,0],xcn[idn1_2-1,0]],[xcn[idn1_1-1,1],xcn[idn1_2-1,1]],[xcn[idn1_1-1,2],xcn[idn1_2-1,2]],'--k')
        #plt.plot([xcn[idn1_2-1,0],xcn[idn1_0-1,0]],[xcn[idn1_2-1,1],xcn[idn1_0-1,1]],[xcn[idn1_2-1,2],xcn[idn1_0-1,2]],'--k')
        
        
        plt.show()
    
'''
    




'''
elnum = 132960

print(ien[elnum,:])
fig = plt.figure(7)
ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')
points = []
for i in range(1,5):
    id = ien[elnum,i]
    print(id)
    ax.scatter(xcn[id,0],xcn[id,1],xcn[id,2],'.',color='red')
    points.append([xcn[id,0],xcn[id,1],xcn[id,2]])
print(points)

for i in range(5,len(ien[elnum,:])):
    id = ien[elnum,i]
    print(id)
    ax.scatter(xcn[id,0],xcn[id,1],xcn[id,2],'.',color='blue')
    points.append([xcn[id,0],xcn[id,1],xcn[id,2]])
print(points)
plt.show()

'''
'''
fig = plt.figure(7)
ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')

for i in range(0,8):
    ax.scatter(Points[i][0],Points[i][1],Points[i][2],'.',color='red')
plt.show()
'''
'''


test = [13428254, 1, 28462350, 1740, 29,
30, 0, 1741, 28462351, 13428373, 2,
1742, 1, 31, 28462352, 13428492, 3,
32, 1743, 2, 28462353, 13428611, 4,
1744, 3, 33, 28462354, 13428730, 5,
1745, 4, 34, 28462355, 13428849, 6,
5, 35, 1746, 28462356, 13428968, 7,
36, 1747, 6, 28462357, 13429087, 8,
1748, 37, 7, 28462358, 13429206, 9,
1749, 38, 8, 28462359, 13429325, 10,
9, 1750, 39, 28462360, 13429444, 11,
1751, 40, 10, 28462361, 13429563, 12,
11, 1752, 41, 28462362, 13429682, 13,
1753, 12, 42, 28462363, 13429801, 14,
43, 13, 1754, 28462364, 13429920, 15,
44, 1755, 14, 28462365, 13430039, 16,
1756, 15, 45, 28462366, 13430158, 17,
1757, 46, 16, 28462367, 13430277, 18,
1758, 47, 17, 28462368, 13430396, 19,
48, 18, 1759, 28462369, 13430515, 20]

fig = plt.figure(10)
ax = plt.axes(projection='3d')
for i in range(0,len(test)):

    #print("element = ", i )
    el_id = test[i]
    #print("element _id ", el_id)
    print(ien[el_id][:])
    vid0 = ien[el_id][1]-1
    vid1 = ien[el_id][2]-1
    vid2 = ien[el_id][3]-1
    vid3 = ien[el_id][4]-1
    
    vid4 = ien[el_id][5]-1
    vid5 = ien[el_id][6]-1
    vid6 = ien[el_id][7]-1
    vid7 = ien[el_id][8]-1
    
    # print(vid0,vid1,vid2,vid3,vid4,vid5,vid6,vid7)
    
    ax.plot([xcn[vid0][0],xcn[vid1][0]],[xcn[vid0][1],xcn[vid1][1]],[xcn[vid0][2],xcn[vid1][2]],'-r')
    ax.plot([xcn[vid1][0],xcn[vid2][0]],[xcn[vid1][1],xcn[vid2][1]],[xcn[vid1][2],xcn[vid2][2]],'-r')
    ax.plot([xcn[vid2][0],xcn[vid3][0]],[xcn[vid2][1],xcn[vid3][1]],[xcn[vid2][2],xcn[vid3][2]],'-r')
    ax.plot([xcn[vid3][0],xcn[vid0][0]],[xcn[vid3][1],xcn[vid0][1]],[xcn[vid3][2],xcn[vid0][2]],'-r')

    ax.plot([xcn[vid4][0],xcn[vid5][0]],[xcn[vid4][1],xcn[vid5][1]],[xcn[vid4][2],xcn[vid5][2]],'-r')
    ax.plot([xcn[vid5][0],xcn[vid6][0]],[xcn[vid5][1],xcn[vid6][1]],[xcn[vid5][2],xcn[vid6][2]],'-r')
    ax.plot([xcn[vid6][0],xcn[vid7][0]],[xcn[vid6][1],xcn[vid7][1]],[xcn[vid6][2],xcn[vid7][2]],'-r')
    ax.plot([xcn[vid7][0],xcn[vid4][0]],[xcn[vid7][1],xcn[vid4][1]],[xcn[vid7][2],xcn[vid4][2]],'-r')
    
    ax.plot([xcn[vid0][0],xcn[vid4][0]],[xcn[vid0][1],xcn[vid4][1]],[xcn[vid0][2],xcn[vid4][2]],'-r')
    ax.plot([xcn[vid1][0],xcn[vid5][0]],[xcn[vid1][1],xcn[vid5][1]],[xcn[vid1][2],xcn[vid5][2]],'-r')
    ax.plot([xcn[vid2][0],xcn[vid6][0]],[xcn[vid2][1],xcn[vid6][1]],[xcn[vid2][2],xcn[vid6][2]],'-r')
    ax.plot([xcn[vid3][0],xcn[vid7][0]],[xcn[vid3][1],xcn[vid7][1]],[xcn[vid3][2],xcn[vid7][2]],'-r')

    el_id = i

    vid0 = ien[el_id][1]-1
    vid1 = ien[el_id][2]-1
    vid2 = ien[el_id][3]-1
    vid3 = ien[el_id][4]-1

    vid4 = ien[el_id][5]-1
    vid5 = ien[el_id][6]-1
    vid6 = ien[el_id][7]-1
    vid7 = ien[el_id][8]-1

    # print(vid0,vid1,vid2,vid3,vid4,vid5,vid6,vid7)

    ax.plot([xcn[vid0][0],xcn[vid1][0]],[xcn[vid0][1],xcn[vid1][1]],[xcn[vid0][2],xcn[vid1][2]],'-b')
    ax.plot([xcn[vid1][0],xcn[vid2][0]],[xcn[vid1][1],xcn[vid2][1]],[xcn[vid1][2],xcn[vid2][2]],'-b')
    ax.plot([xcn[vid2][0],xcn[vid3][0]],[xcn[vid2][1],xcn[vid3][1]],[xcn[vid2][2],xcn[vid3][2]],'-b')
    ax.plot([xcn[vid3][0],xcn[vid0][0]],[xcn[vid3][1],xcn[vid0][1]],[xcn[vid3][2],xcn[vid0][2]],'-b')

    ax.plot([xcn[vid4][0],xcn[vid5][0]],[xcn[vid4][1],xcn[vid5][1]],[xcn[vid4][2],xcn[vid5][2]],'-b')
    ax.plot([xcn[vid5][0],xcn[vid6][0]],[xcn[vid5][1],xcn[vid6][1]],[xcn[vid5][2],xcn[vid6][2]],'-b')
    ax.plot([xcn[vid6][0],xcn[vid7][0]],[xcn[vid6][1],xcn[vid7][1]],[xcn[vid6][2],xcn[vid7][2]],'-b')
    ax.plot([xcn[vid7][0],xcn[vid4][0]],[xcn[vid7][1],xcn[vid4][1]],[xcn[vid7][2],xcn[vid4][2]],'-b')

    ax.plot([xcn[vid0][0],xcn[vid4][0]],[xcn[vid0][1],xcn[vid4][1]],[xcn[vid0][2],xcn[vid4][2]],'-b')
    ax.plot([xcn[vid1][0],xcn[vid5][0]],[xcn[vid1][1],xcn[vid5][1]],[xcn[vid1][2],xcn[vid5][2]],'-b')
    ax.plot([xcn[vid2][0],xcn[vid6][0]],[xcn[vid2][1],xcn[vid6][1]],[xcn[vid2][2],xcn[vid6][2]],'-b')
    ax.plot([xcn[vid3][0],xcn[vid7][0]],[xcn[vid3][1],xcn[vid7][1]],[xcn[vid3][2],xcn[vid7][2]],'-b')





'''
plt.show()
