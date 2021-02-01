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
    






fname1 = 'qw_Cgrid_h_0p005.dat'
fname2 = 'qw_Cgrid_h_0p010.dat'
fname3 = 'qw_Cgrid_h_0p0125.dat'
fname33 = 'qw_turb1.dat'
fname4 = 'qw_Cgrid_h_0p0150.dat'
fname5 = 'qw_Cgrid_h_0p020.dat'
fname6 = 'qw_L_f1p0m5.dat'
fname7 = 'qw_tess.dat'
fname8 = 'qw_tess_plus.dat'

fname9 = 'qw_Cgrid_f_tess_plus.dat'
fname10 = 'qw_Cgrid_f.dat'

fname11 = 'qw_fail.dat'
fname12 = 'qqw_hz_0p05.dat'

fnamer0 = 'qw_base.dat'
fnamer1 = 'qw_Cgrid_f.dat'
fnamer = 'qw_reference.dat'

pl1 = []
PlotSpanWall(fname1,3,11,1188,700,7,1,[0.5,0.0,0.1],2,'$f$ = 0.005 -> ',pl1);
PlotSpanWall(fnamer,4,11,2394,2800,7,1,[0.01,0.01,0.01],2,'ref hex -> ',pl1);
plt.grid()
pl2 = []
PlotSpanWall(fname2,3,11,1188,700,7,2,[0.5,0.0,0.1],2,'$f$ = 0.01 -> ',pl2);
PlotSpanWall(fnamer,4,11,2394,2800,7,2,[0.01,0.01,0.01],2,'ref hex -> ',pl2);
plt.grid()
pl3 = []
PlotSpanWall(fname3,3,11,1188,700,7,3,[0.5,0.0,0.1],2,'$f$ = 0.0125 -> ',pl3);
PlotSpanWall(fnamer,4,11,2394,2800,7,3,[0.01,0.01,0.01],2,'ref hex -> ',pl3);
plt.grid();
pl4 = []
PlotSpanWall(fname33,3,11,1188,700,7,4,[0.5,0.0,0.1],2,'$f$ = 0.015 ->',pl4);
PlotSpanWall(fnamer,4,11,2394,2800,7,4,[0.01,0.01,0.01],2,'ref hex -> ',pl4);
pl5 = []
PlotSpanWall(fname5,3,11,1188,700,7,5,[0.5,0.0,0.1],2,'$f$ = 0.02 ->',pl5);
PlotSpanWall(fnamer,4,11,2394,2800,7,5,[0.01,0.01,0.01],2,'ref hex -> ',pl5);
pl6 = []
PlotSpanWall(fname6,3,11,1188,700,7,6,[0.5,0.0,0.1],2,'length var ->',pl6);
PlotSpanWall(fnamer,4,11,2394,2800,7,6,[0.01,0.01,0.01],2,'ref hex -> ',pl6);
pl7 = []
PlotSpanWall(fname7,3,11,1188,700,7,7,[0.5,0.0,0.1],2,'tess ->',pl7);
PlotSpanWall(fnamer,4,11,2394,2800,7,7,[0.01,0.01,0.01],2,'ref hex -> ',pl7);
pl8 = []
PlotSpanWall(fname8,3,11,1188,700,7,8,[0.5,0.0,0.1],2,'tess +',pl8);
PlotSpanWall(fnamer,4,11,2394,2800,7,8,[0.01,0.01,0.01],2,'ref hex -> ',pl8);

pl9 = []
PlotSpanWall(fname9,3,11,2388,1400,7,9,[0.5,0.0,0.1],2,'tess +',pl9);
PlotSpanWall(fname10,4,11,1194,1400,7,9,[0.01,0.01,0.01],2,'ref hex -> ',pl9);


pl10 = []
PlotSpanWall(fnamer0,4,11,594,700,7,11,[0.0,0.0,0.11],2,'ref 0',pl10);
PlotSpanWall(fnamer1,4,11,1194,1400,7,11,[0.00,0.11,0.00],2,'ref 1',pl10);
PlotSpanWall(fnamer, 4,11,2394,2800,7,11,[0.11,0.00,0.00],2,'ref 2',pl10);

pl11 = []
PlotSpanWall(fname11,3,11,2388,1400,7,12,[0.5,0.0,0.1],2,'tess +',pl11);
PlotSpanWall(fnamer,4,11,2394,2800,7,12,[0.01,0.01,0.01],2,'ref hex -> ',pl11);

pl12 = []
PlotSpanWall(fname12,3,11,1188,700,7,13,[0.5,0.0,0.1],2,'tess +',pl12);
PlotSpanWall(fnamer,4,11,2394,2800,7,13,[0.01,0.01,0.01],2,'ref hex -> ',pl12);

plt.grid()
plt.show()
