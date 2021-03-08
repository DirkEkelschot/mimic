#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe


#with open('wvars_original_wall.dat') as f:



def PlotSpanWall(fname,type,hoffset,nEl,Npo,n_z,fign,c,l,lab,pl):

    A=[]
    with open(fname) as f:
        A = f.readlines()[0:]
    print("sizeA",len(A))
    ncol    = 5
    nstep   = int(Npo/n_z)
    print("nstep",nstep)
    noff    = nstep
    Nv      = int(Npo/5)
    print("Nv",Nv)
    Nvr     = Npo % 5
    if Nvr > 0:
        Nv  = int(Npo/5)+1
        
    Nvel    = int(nEl/5.0)
    Nvelr   = nEl % 5
    if Nvelr>0:
        Nvel = int(nEl/5)+1
    
    of = hoffset;
    x=[]
    for i in range(0,1):
        vn=[]
        for j in range(0,Nv):
            s=A[of+j][:].split()
            for k in range(0,len(s)):
                vn.append(float(s[k]))
        of=of+Nv;
        x.append(vn);


    y = []
    of=hoffset+Nv; #skipping header
    for i in range(0,1):
        vn=[]
        for j in range(0,Nv):
            s=A[of+j][:].split()
            for k in range(0,len(s)):
                vn.append(float(s[k]))
        of=of+Nv;
        y.append(vn);
    z = []
    var2 = []
    lines = [];
    of=hoffset+2*Nv; #skipping header
    for i in range(0,1):
        vn=[]
        for j in range(0,Nv):
            s=A[of+j][:].split()
            for k in range(0,len(s)):
                vn.append(float(s[k]))
        of=of+Nv;
        z.append(vn);

    qw = []
    of=hoffset+3*Nv; #skipping header
    Nvn = Nvel;
    for i in range(0,1):
        vn=[]
        for j in range(0,Nvn):
            s=A[of+j][:].split()
            for k in range(0,len(s)):
                vn.append(float(s[k]))
        of=of+Nvn;
        qw.append(vn);
    el = []
    of= hoffset+3*Nv+Nvel; #skipping header
    Nvn = nEl;
    print(Nv,hoffset+3*Nv,of,Nv)
    for i in range(0,1):
        vn=[]
        for j in range(0,Nvn):
            s=A[of+j][:].split()
            for k in range(0,len(s)):
                vn.append(float(s[k]))
        of=of+Nvn;
        el.append(vn);
    print(len(el[0]))
    qwno = np.zeros((Npo,1))
    qwtel = np.zeros((Npo,1))
    qwr = np.zeros((Npo,1))
    print("len(el)",len(el[0])/4,Npo)
    for i in range(0,nEl):
        for j in range(0,type):
            
            qwno[int(el[0][i*4+j])-1] = qwno[int(el[0][i*4+j])-1]+qw[0][i]
            qwtel[int(el[0][i*4+j])-1] = qwtel[int(el[0][i*4+j])-1]+1
            
    for i in range(0,Npo):
        qwr[i] = qwno[i]/qwtel[i];


    xn0 = np.zeros((noff,1));
    yn0 = np.zeros((noff,1));
    zn0 = np.zeros((noff,1));
    qwn0 = np.zeros((noff,1));

    xn0_s = np.zeros((n_z,noff));
    yn0_s = np.zeros((n_z,noff));
    zn0_s = np.zeros((n_z,noff));
    qwn0_s = np.zeros((n_z,noff));

    plt.figure(fign,figsize=(9,8))
    #plt.figure(fign)

    
    plt.xticks(fontsize=20, rotation=0)
    plt.yticks(fontsize=20, rotation=0)
    #plt.title("Comparing $q_w$", fontsize=20)
    cn = []
    cn.append(c[0])
    cn.append(c[1])
    cn.append(c[2])
    for i in range(0,n_z):
        for s in range(0,noff):
            xn0[s]  = xn0[s]+x[0][i*noff+s];
            yn0[s]  = yn0[s]+y[0][i*noff+s];
            zn0[s]  = zn0[s]+z[0][i*noff+s];
            qwn0[s] = qwn0[s]+qwr[i*noff+s];
            
            xn0_s[i,s]  = x[0][i*noff+s];
            yn0_s[i,s]  = y[0][i*noff+s];
            zn0_s[i,s]  = z[0][i*noff+s];
            qwn0_s[i,s] = qwr[i*noff+s];
        cn[0] = (i+1)*c[0]
        cn[1] = (i+1)*c[1]
        cn[2] = (i+1)*c[2]
        if cn[0]>=1:
           cn[0] = 1
        if cn[1]>=1:
           cn[1] = 1
        if cn[2]>=1:
           cn[2] = 1
#        cn[0] = i*c[0];
#        cn[1] = i*c[1];
#        cn[2] = i*c[2];
        pl+=plt.plot(yn0_s[i,:],qwn0_s[i,:]/(100*100),'-',color=[cn[0],cn[1],cn[2]],linewidth=l,label= lab+' z = '+str(zn0_s[i,0]))
        #print(zn0_s[i,:])
        
    plt.legend(handles=pl,loc='upper right')
    plt.xlabel("y [$m$]", fontsize=20)
    plt.ylabel("q [$W$/$m^2$]", fontsize=20)
    plt.grid()
    for i in range(0,noff):
        xn0[i] = xn0[i]/n_z;
        yn0[i] = yn0[i]/n_z;
        zn0[i] = zn0[i]/n_z;
        qwn0[i] = qwn0[i]/n_z;
    








fnamer0 = 'qw_base.dat'
fname_a1 = 'qw_it1.dat'
#fname_a2 = 'qw_it2.dat'

pl13 = []
PlotSpanWall(fname_a2,3,11,1188,700,7,14,[0.5,0.0,0.1],2,'tess +',pl13);
PlotSpanWall(fnamer,4,11,2394,2800,7,14,[0.01,0.01,0.01],2,'ref hex -> ',pl13);

plt.grid()
plt.show()
