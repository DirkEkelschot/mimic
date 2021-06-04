#include "src/adapt_io.h"
#include "src/adapt_recongrad.h"
#include "src/adapt_output.h"
#include "src/adapt_geometry.h"
#include "src/adapt_bltopology.h"
#include "src/adapt_parops.h"
#include "src/hex2tet.h"
#include "src/adapt_boundary.h"
#include <iomanip>
#include "main_api.h"

int madam_api( const char *fn_grid, const char *fn_conn, const char *fn_data, MPI_Comm comm, int world_rank, int world_size) {
    
    std::cout << "Starting MADAM_API";
    std::cout << " comm="<<comm<<" world_rank="<<world_rank<<" world_size="<<world_size<< std::endl;
    MPI_Info info = MPI_INFO_NULL;
    
    int debug = 0;
    char buffer[20];
    //========================================================================
    //========================================================================
    //========================================================================
        
    const char* fn_metric = "metric.inp";
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    int varia = 4;
        
    int ReadFromStats = 0;
    if(metric_inputs.size()==6)
    {
        ReadFromStats=metric_inputs[5];
    }
    std::cout << "grid="<<fn_grid <<"  conn="<<fn_conn<<"  data="<<fn_data << std::endl;
    std::cout << "HACK JMB!! fixing grid, and data to known values... something is weird with read in from US3D"<< std::endl;
    // fn_grid = "grid.h5";
    // fn_conn = "conn.h5";
    // fn_data = "data.h5";
    US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,world_rank,world_size,info);
    int Nve = us3d->xcn->getNglob();
        
    int Nel_part = us3d->ien->getNrow();
    

    std::cout << "HACK JMB!! Done with ReadUS3DData";
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(), comm, world_rank, world_size);
    std::cout << "HACK JMB!! Done with ParallelState 1"<< std::endl;
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(), comm, world_rank, world_size);
    std::cout << "HACK JMB!! Done with ParallelState 2"<< std::endl;
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv, comm, world_rank, world_size);
    std::cout << "HACK JMB!! Done with ParallelState 2"<< std::endl;
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(), comm, world_rank, world_size);
    std::cout << "HACK JMB!! Done with ParallelState stuff"<< std::endl;
    Array<double>* Uivar = new Array<double>(Nel_part,1);
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    for(int i=0;i<Nel_part;i++)
    {
        rhoState = us3d->interior->getVal(i,0);
        uState   = us3d->interior->getVal(i,1);
        vState   = us3d->interior->getVal(i,2);
        wState   = us3d->interior->getVal(i,3);
        TState   = us3d->interior->getVal(i,4);
        VtotState = sqrt(uState*uState+vState*vState+wState*wState);
        aState   = sqrt(1.4*287.05*TState);
        MState = VtotState/aState;
        Uivar->setVal(i,0,MState);
    }
        
    delete us3d->interior;
    
    clock_t t,t1;
    double tmax = 0.0;
    double tn = 0.0;
    t = clock();

    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Uivar, comm, world_rank, world_size);
        
        
    double duration = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Ptime = 0.0;
    MPI_Allreduce(&duration, &Ptime, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
    std::cout << "Timing partitioning: " << duration << std::endl;
    }
        
    std::vector<int> LocElem = P->getLocElem();
    std::vector<double> Uvaria  = P->getLocElemVaria();
    std::map<int,Array<double>*> Uvaria_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid = LocElem[i];
        UvariaV   = Uvaria[i];
        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map[gid] = Uarr;
    }
        
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm, world_rank, world_size);
        
    std::map<int,double> Volumes = meshTopo->getVol();
        
    P->AddStateVecForAdjacentElements(Uvaria_map,1,comm, world_rank, world_size);
        
    std::map<int,Array<double>* > var_vmap = P->ReduceStateVecToAllVertices(Uvaria_map,1);
    std::map<int,Array<double>* >::iterator vm;
    for(vm=var_vmap.begin();vm!=var_vmap.end();vm++)
    {
        if(vm->second->getVal(0,0)<0.0)
        {
            std::cout << "Mach vert " << vm->first << " " << vm->second->getVal(0,0) << std::endl;
        }
        
    }
        
        
    if(world_rank == 0)
    {
        std::cout << "Started creating mesh topology object... " << std::endl;
    }

        
    if(world_rank == 0)
    {
        std::cout << "Finished creating mesh topology object... "  << std::endl;
    }
        
    if(world_rank == 0)
    {
        std::cout << "Setting the ghost element data... "  << std::endl;
    }
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }
        
    delete us3d->ghost;
        
    t = clock();
    if(world_rank == 0)
    {
        std::cout << "Started reconstructing the gradient... " << std::endl;
    }
        
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,Uvaria_map,gB,comm, world_rank, world_size);
        
    double Gtiming = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Gmax_time = 0.0;
    MPI_Allreduce(&Gtiming, &Gmax_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        std::cout << "Finished reconstructing the gradient... " << std::endl;
        std::cout << "Timing gradient reconstruction... " << Gmax_time << std::endl;
    }
        
    P->AddStateVecForAdjacentElements(dUdXi,3,comm, world_rank, world_size);
        
    std::map<int,Array<double>* >::iterator grit;
    std::map<int,Array<double>* > dUidxi_map;
    std::map<int,Array<double>* > dUidyi_map;
    std::map<int,Array<double>* > dUidzi_map;
    std::map<int,Array<double>* > dUidXi_map;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
            
        Array<double>* dUdx_E = new Array<double>(1,1);
        dUdx_E->setVal(0,0,grit->second->getVal(0,0));
        Array<double>* dUdy_E = new Array<double>(1,1);
        dUdy_E->setVal(0,0,grit->second->getVal(1,0));
        Array<double>* dUdz_E = new Array<double>(1,1);
        dUdz_E->setVal(0,0,grit->second->getVal(2,0));
            
        dUidxi_map[grit->first]=dUdx_E;
        dUidyi_map[grit->first]=dUdy_E;
        dUidzi_map[grit->first]=dUdz_E;
            
        delete grit->second;
            
    }

    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUidxi_map,gB,comm,world_rank,world_size);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUidyi_map,gB,comm,world_rank,world_size);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUidzi_map,gB,comm,world_rank,world_size);
   
    std::map<int,Array<double>*> Hess_map;

    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;
        
    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
            
        Array<double>* Hess = new Array<double>(6,1);
            
        Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
        Hess->setVal(1,0,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(2,0,dU2dXi2[gid]->getVal(2,0));

        Hess->setVal(3,0,dU2dYi2[gid]->getVal(1,0));
        Hess->setVal(4,0,dU2dYi2[gid]->getVal(2,0));
        Hess->setVal(5,0,dU2dZi2[gid]->getVal(2,0));
            
        Hess_map[gid] = Hess;
            
        delete dU2dXi2[gid];
        delete dU2dYi2[gid];
        delete dU2dZi2[gid];
            
        t++;
    }
        
    dU2dXi2.clear();
    dU2dYi2.clear();
    dU2dZi2.clear();
        
        
    double* Hessie = new double[9];
    double * WRn = new double[3];
    Array<double>* DR  = new Array<double>(3,3);
    Array<double>* UR  = new Array<double>(3,3);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++
    double hmin         = metric_inputs[1];
    double hmax         = metric_inputs[2];
    double f            = metric_inputs[3];
    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double cmplxty = 0.0;
    double cmplxty2 = 0.0;
    double Volu=0.0,cmplxty_red=0.0,cmplxty_tmp=0.0,cmplxty_tmp2=0.0;
    double po = 6.0;
    for(int j=0;j<3;j++)
    {
        for(int k=0;k<3;k++)
        {
            DR->setVal(j,k,0.0);
        }
    }
    for(itgg=Hess_map.begin();itgg!=Hess_map.end();itgg++)
    {
        Hessie[0] = itgg->second->getVal(0,0);
        Hessie[1] = itgg->second->getVal(1,0);
        Hessie[2] = itgg->second->getVal(2,0);
            
        Hessie[3] = itgg->second->getVal(1,0);
        Hessie[4] = itgg->second->getVal(3,0);
        Hessie[5] = itgg->second->getVal(4,0);
            
        Hessie[6] = itgg->second->getVal(2,0);
        Hessie[7] = itgg->second->getVal(4,0);
        Hessie[8] = itgg->second->getVal(5,0);
            
        Eig* eig = ComputeEigenDecomp(3, Hessie);
            
        WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[2] = std::min(std::max(f*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
            
        DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);
            
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        double detRf = sqrt( Rf->getVal(0,0)*(Rf->getVal(1,1)*Rf->getVal(2,2)-Rf->getVal(2,1)*Rf->getVal(1,2))
                            -Rf->getVal(0,1)*(Rf->getVal(1,0)*Rf->getVal(2,2)-Rf->getVal(2,0)*Rf->getVal(1,2))
                            +Rf->getVal(0,2)*(Rf->getVal(1,0)*Rf->getVal(2,1)-Rf->getVal(2,0)*Rf->getVal(1,1)));
            
        Volu = Volumes[itgg->first];
        cmplxty_tmp = detRf;
        cmplxty_tmp2 = Volu;
        //std::pow(cmplxty_tmp,(po+1.0)/(2.0*po+3.0));
        cmplxty = cmplxty + cmplxty_tmp;
        cmplxty2 = cmplxty2 + cmplxty_tmp2;
        delete iVR;
        delete Rs;
        delete Rf;
        

        MPI_Allreduce(&cmplxty, &cmplxty_red, 1, MPI_DOUBLE, MPI_SUM, comm);
        //cmplxty_red = std::pow(cmplxty_red,(po+1)/(2.0*po+3.0));
        
        P->AddStateVecForAdjacentElements(Hess_map,6,comm,world_rank,world_size);

        std::map<int,Array<double>* > hess_vmap = P->ReduceStateVecToAllVertices(Hess_map,6);
    if(world_rank == 0)
    {
        std::cout << "complexity = " << cmplxty_red << " " << Nve << std::endl;
    }
        cmplxty_red = Nve/cmplxty_red;
        
        if(world_rank == 0)
    {
        std::cout << "complexity = " << cmplxty_red << std::endl;
        }
        ComputeMetric(P, metric_inputs, comm, world_rank, world_size, var_vmap, hess_vmap, 1.0, po);
        
        if(world_rank==0)
        {
            std::cout << "Started gathering metric data on rank 0..." <<std::endl;
        }
        Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, hess_vmap, comm,world_rank,world_size);
        
        hess_vmap.clear();
        
        
        if(world_rank==0)
        {
            std::cout << "Finished gathering metric data on rank 0..." <<std::endl;
        }
        
        Array<double>*  xcn_g;
        Array<int>*     ief_g;
        Array<int>*     ien_g;
        Array<int>*     ifn_g;
        Array<int>*     if_ref_g;
        Array<int>*     ife_g;
        
        if(world_rank == 0)
        {
            xcn_g       = new Array<double>(us3d->xcn->getNglob(),3);
            ief_g       = new Array<int>(us3d->ief->getNglob(),6);
            ien_g       = new Array<int>(us3d->ien->getNglob(),8);
            if_ref_g    = new Array<int>(us3d->ifn->getNglob(),1);
            ifn_g       = new Array<int>(us3d->ifn->getNglob(),4);
            ife_g       = new Array<int>(us3d->ifn->getNglob(),2);
        }
        else
        {
            xcn_g    = new Array<double>(1,1);
            ief_g    = new Array<int>(1,1);
            ien_g    = new Array<int>(1,1);
            if_ref_g = new Array<int>(1,1);
            ifn_g    = new Array<int>(1,1);
            ife_g    = new Array<int>(1,1);
        }

        int* ien_nlocs      = new int[world_size];
        int* ien_offsets    = new int[world_size];
        int* ief_nlocs      = new int[world_size];
        int* ief_offsets    = new int[world_size];
        int* xcn_nlocs      = new int[world_size];
        int* xcn_offsets    = new int[world_size];
        int* ifn_nlocs      = new int[world_size];
        int* ifn_offsets    = new int[world_size];
        int* if_ref_nlocs   = new int[world_size];
        int* if_ref_offsets = new int[world_size];
        int* ife_nlocs      = new int[world_size];
        int* ife_offsets    = new int[world_size];
        
        for(int i=0;i<world_size;i++)
        {
            xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]  *3;
            xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;

            ief_nlocs[i]   = ien_pstate->getNlocs()[i]  *6;
            ief_offsets[i] = ien_pstate->getOffsets()[i]*6;

            ien_nlocs[i]   = ien_pstate->getNlocs()[i]  *8;
            ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
            
            ifn_nlocs[i]   = ife_pstate->getNlocs()[i]  *4;
            ifn_offsets[i] = ife_pstate->getOffsets()[i]*4;
            
            if_ref_nlocs[i]   = ife_pstate->getNlocs()[i]  *1;
            if_ref_offsets[i] = ife_pstate->getOffsets()[i]*1;
            
            ife_nlocs[i]   = ife_pstate->getNlocs()[i]  *2;
            ife_offsets[i] = ife_pstate->getOffsets()[i]*2;
        }

        MPI_Gatherv(&us3d->xcn->data[0],
                    us3d->xcn->getNrow()*3,
                    MPI_DOUBLE,
                    &xcn_g->data[0],
                    xcn_nlocs,
                    xcn_offsets,
                    MPI_DOUBLE, 0, comm);

        MPI_Gatherv(&us3d->ief->data[0],
                    us3d->ief->getNrow()*6,
                    MPI_INT,
                    &ief_g->data[0],
                    ief_nlocs,
                    ief_offsets,
                    MPI_INT, 0, comm);

        MPI_Gatherv(&us3d->ien->data[0],
                    us3d->ien->getNrow()*8,
                    MPI_INT,
                    &ien_g->data[0],
                    ien_nlocs,
                    ien_offsets,
                    MPI_INT, 0, comm);

        MPI_Gatherv(&us3d->ifn->data[0],
                    us3d->ifn->getNrow()*4,
                    MPI_INT,
                    &ifn_g->data[0],
                    ifn_nlocs,
                    ifn_offsets,
                    MPI_INT, 0, comm);
        
        MPI_Gatherv(&us3d->if_ref->data[0],
                    us3d->if_ref->getNrow()*1,
                    MPI_INT,
                    &if_ref_g->data[0],
                    if_ref_nlocs,
                    if_ref_offsets,
                    MPI_INT, 0, comm);

        MPI_Gatherv(&us3d->ife->data[0],
                    us3d->ife->getNrow()*2,
                    MPI_INT,
                    &ife_g->data[0],
                    ife_nlocs,
                    ife_offsets,
                    MPI_INT, 0, comm);
        
        delete P;
        
        dUdXi.clear();
        
        dUidxi_map.clear();
        dUidyi_map.clear();
        dUidzi_map.clear();
        
        delete us3d->ien;
        delete us3d->ief;
        delete us3d->iee;
        delete us3d->xcn;
        delete us3d->iet;
        delete us3d->ife;
        delete us3d->if_ref;
        delete us3d->ifn;
        //delete us3d->ghost;
        delete us3d->ie_Nv;
        delete us3d->ie_Nf;
        
         
        if(world_rank == 0)
        {
            BoundaryMap* bmap = new BoundaryMap(ifn_g, if_ref_g);
            int wall_id = 3;
            int nLayer  = metric_inputs[4];
            if(nLayer>0)
            {
                int counter = 0;
                
                std::map<int,std::vector<int> > bnd_face_map = bmap->getBfaceMap();
                std::map<std::set<int>,int> tria_ref_map = bmap->getTriaRefMap();
                std::map<std::set<int>,int> quad_ref_map = bmap->getQuadRefMap();
                std::map<int,int> vert_ref_map = bmap->getNodeRefMap();
                
                BLShellInfo* BLshell = FindOuterShellBoundaryLayerMesh(wall_id, nLayer,xcn_g,ien_g,ief_g,ife_g,ifn_g,xcn_pstate,ien_pstate,bnd_face_map,vert_ref_map,comm,world_rank,world_size);
   
                //========================================================================
                //========================================================================
                //            struct BLShellInfo
                //            {
                //                std::map<int,int> ShellFace2BFace;
                //                std::map<std::set<int>,int> ShellTri2FaceID;
                //                std::map<std::set<int>,std::vector<int> > ShellFaceID2TriID
                //                std::map<int,int> FaceID2TopoType;
                //                std::map<int,std::map<int,int> > ShellFace2ShellVert2OppositeBoundaryVerts;
                //                Array<int>* ShellRef;
                //            };
                //========================================================================
                //========================================================================
                
                int nbHex       =  ien_g->getNrow();
                int nbPrisms    =  bnd_face_map[wall_id].size()*(nLayer)*2;
                //int nbTets      = (nbHex-bnd_face_map[wall_id].size()*(nLayer))*6;
                int nbHexsNew   = (nbHex-bnd_face_map[wall_id].size()*(nLayer));
                
                
            
                int ith = 0;
                std::set<int> u_tet_vert;
                std::map<int,int> lv2gv_tet_mesh;
                std::map<int,int> gv2lv_tet_mesh;
                std::map<int,double*> metric_hex2tet;
                Array<int>* ien_hex2tet = new Array<int>(nbHexsNew,8);
                int sv = 0;
                std::vector<int> locTet_verts;
                for(int i=0;i<ien_g->getNrow();i++)
                {
                    if(BLshell->elements_set.find(i)==BLshell->elements_set.end())
                    {
                        for(int j=0;j<8;j++)
                        {
                            int val = ien_g->getVal(i,j);
                            
                            if(u_tet_vert.find(val)==u_tet_vert.end())
                            {
                                u_tet_vert.insert(val);
                                locTet_verts.push_back(val);
                                gv2lv_tet_mesh[val]=sv;
                                lv2gv_tet_mesh[sv]=val;
                                ien_hex2tet->setVal(ith,j,sv);
                                double* met = new double[6];
                                
                                met[0] = mv_g->getVal(val,0);
                                met[1] = mv_g->getVal(val,1);
                                met[2] = mv_g->getVal(val,2);
                                met[3] = mv_g->getVal(val,3);
                                met[4] = mv_g->getVal(val,4);
                                met[5] = mv_g->getVal(val,5);
                                
                                
                                metric_hex2tet[val]=met;
                                sv++;
                            }
                            else
                            {
                                int lv = gv2lv_tet_mesh[val];
                                ien_hex2tet->setVal(ith,j,lv);
                            }
                        }
                        ith++;
                    }
                }
        
                MMG5_pMesh mmgMesh_TET = NULL;
                MMG5_pSol mmgSol_TET   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppMet,&mmgSol_TET,
                MMG5_ARG_end);
                
                int nbVerts_TET=locTet_verts.size();

                if ( MMG3D_Set_meshSize(mmgMesh_TET,nbVerts_TET,nbHexsNew*6,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh_TET,mmgSol_TET,MMG5_Vertex,mmgMesh_TET->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                
                int cshell = 0;
                int nshell = 0;
                
                for(int i=0;i<nbVerts_TET;i++)
                {
                    int gv=lv2gv_tet_mesh[i];
                    
                    mmgMesh_TET->point[i+1].c[0] = xcn_g->getVal(gv,0);
                    mmgMesh_TET->point[i+1].c[1] = xcn_g->getVal(gv,1);
                    mmgMesh_TET->point[i+1].c[2] = xcn_g->getVal(gv,2);
                    
                    mmgMesh_TET->point[i+1].ref  = BLshell->ShellRef->getVal(gv,0);
                    
                    if(BLshell->ShellRef->getVal(gv,0)==0)
                    {
                        std::cout << i << " " << gv << " " << xcn_g->getNrow() << " zero here already" <<std::endl;
                    }
                    if(BLshell->ShellRef->getVal(gv,0)==-1)
                    {
                        cshell++;
                    }
                    if(BLshell->ShellRef->getVal(gv,0)==-3)
                    {
                        nshell++;
                    }
                    
                    double m11 = mv_g->getVal(gv,0);
                    double m12 = mv_g->getVal(gv,1);
                    double m13 = mv_g->getVal(gv,2);
                    double m22 = mv_g->getVal(gv,3);
                    double m23 = mv_g->getVal(gv,4);
                    double m33 = mv_g->getVal(gv,5);
                    
                    if ( MMG3D_Set_tensorSol(mmgSol_TET, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
                }
                int ref = 20;
                
                int* hexTabNew = new int[9*(nbHexsNew+1)];
                for(int i=0;i<nbHexsNew;i++)
                {
                    int hexTabPosition = 9*(i+1);
                    for(int j=0;j<8;j++)
                    {
                        int val = ien_hex2tet->getVal(i,j)+1;
                        hexTabNew[hexTabPosition+j] = val;
                    }
                    hexTabNew[hexTabPosition+8] = ref;
                }
                std::cout << "Check the orientation and modify when necessary so that we can cut each hex up into 6 tetrahedra..."<<std::endl;
                int num = H2T_chkorient(mmgMesh_TET,hexTabNew,nbHexsNew);
                int* adjahexNew = NULL;
                adjahexNew = (int*)calloc(6*nbHexsNew+7,sizeof(int));
                assert(adjahexNew);
                 
                if(!H2T_hashHexa(hexTabNew,adjahexNew,nbHexsNew))
                {
                    std::cout << "Error :: setting up the new adjacency for the hexes after reorientation." << std::endl;
                }
                
                Hedge       hed22;
                hed22.size  = 6*nbHexsNew;
                hed22.hnxt  = 6*nbHexsNew;
                hed22.nhmax = (int)(16*6*nbHexsNew);
                hed22.item  = NULL;
                hed22.item  = (hedge*)calloc(hed22.nhmax+1,sizeof(hedge));

                for (int k=6*nbHexsNew; k<hed22.nhmax; k++)
                {
                    hed22.item[k].nxt = k+1;
                }
                
                std::cout << "Cut each hexahedral up into 6 tetrahedra..."<<std::endl;
                //we need to do this in order to determine the orientation of the triangles at the shell interface and trace back how we need to tesselate the wall boundary into triangles so that they match eachothers orientation.
                int ret = H2T_cuthex(mmgMesh_TET, &hed22, hexTabNew, adjahexNew, nbHexsNew);
                int nel_tets = mmgMesh_TET->ne;
                std::map<int,std::vector<int> > unique_shell_tri_map;
                
                std::set<std::set<int> > unique_shell_tris;
                int shell_T_id=0;
                int shell_T_id2 = 0;
                int tellertOr = 0;
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref == 20)
                    {
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                            
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        if( mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        
                        tellertOr++;
                    }
                }
                
                std::map<std::set<int>,int > shelltri2fid=BLshell->ShellTri2FaceID;
                std::map<int,int> shellFace2bFace=BLshell->ShellFace2BFace;
                std::set<std::set<int> >::iterator itset;
                std::map<int,std::vector<int> > TriID2ShellFaceID;
                std::vector<std::vector<int> > u_tris;
                std::vector<std::vector<int> > u_tris_loc;
                int teller=0;
                std::set<int> unique_shell_vert;
                std::vector<int> U_shell_vert;
                std::map<int,int> glob2loc_shell_vert;
                std::map<int,int> loc2glob_shell_vert;
                int loc_s_v=0;
                for(itset=unique_shell_tris.begin();itset!=unique_shell_tris.end();itset++)
                {
                    std::set<int> shell_tri = *itset;
                    std::set<int>::iterator itsh;
                    std::vector<int> u_tri_vec(3);
                    std::vector<int> u_tri_vec_loc(3);
                    int bb = 0;
                    for(itsh=shell_tri.begin();itsh!=shell_tri.end();itsh++)
                    {
                        u_tri_vec[bb]=*itsh;
                        
                        if(unique_shell_vert.find(*itsh)==unique_shell_vert.end())
                        {
                            //u_tri_vec_loc.push_back(loc_s_v);
                            unique_shell_vert.insert(*itsh);
                            U_shell_vert.push_back(*itsh);
                            loc2glob_shell_vert[loc_s_v]=*itsh;
                            glob2loc_shell_vert[*itsh]=loc_s_v;
                            u_tri_vec_loc[bb]= loc_s_v;
                            loc_s_v++;
                        }
                        else
                        {
                            int local_s_vId = glob2loc_shell_vert[*itsh];
                            u_tri_vec_loc[bb]=local_s_vId;
                        }
                        bb++;
                    }
                    
                    u_tris.push_back(u_tri_vec);
                    u_tris_loc.push_back(u_tri_vec_loc);
                    
                    int shell_faceid        = shelltri2fid[shell_tri];
                    int bfaceID             = shellFace2bFace[shell_faceid];
                    std::map<int,int> v2v   = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid];
                    TriID2ShellFaceID[teller].push_back(shell_faceid);
                    BLshell->ShellFaceID2TriID[shell_faceid].push_back(teller);
                    teller++;
                }
                
                std::cout << "Extracting the prismatic boundary layer mesh..." <<std::endl;
                Mesh_Topology_BL* mesh_topo_bl2 =  ExtractBoundaryLayerMeshFromShell(u_tris, BLshell, wall_id, 
                                                                                     nLayer, 
                                                                                     xcn_g, ien_g, ief_g, ife_g, ifn_g, 
                                                                                     xcn_pstate, ien_pstate, 
                                                                                     bnd_face_map, tria_ref_map, quad_ref_map, 
                                                                                     comm, world_rank, world_size);
            
                // counting the number of boundary triangles and quads in the BL mesh.
                int nTriangles_BL  = 0;
                int nQuads_BL      = 0;
                std::map<int,std::vector<std::vector<int> > >::iterator itrr;
                for(itrr=mesh_topo_bl2->bcTria.begin();itrr!=mesh_topo_bl2->bcTria.end();itrr++)
                {
                    nTriangles_BL  = nTriangles_BL+itrr->second.size();

                }
                for(itrr=mesh_topo_bl2->bcQuad.begin();itrr!=mesh_topo_bl2->bcQuad.end();itrr++)
                {
                    nQuads_BL  = nQuads_BL+itrr->second.size();
                }
                
                std::map<int,int> loc2glob_final_verts;
                std::map<int,int> glob2loc_final_verts;
                std::map<int,std::vector<std::vector<int> > >::iterator iter;
                std::set<int> unique_prism_verts;
                int locp = 0;
                for(iter=mesh_topo_bl2->BLlayersPrisms.begin();
                   iter!=mesh_topo_bl2->BLlayersPrisms.end();iter++)
                {
                    int numit=iter->second.size();
                    
                    for(int p=0;p<numit;p++)
                    {
                        std::vector<int> prism = iter->second[p];
                        for(int q=0;q<prism.size();q++)
                        {
                            int pp = prism[q];
                            if(unique_prism_verts.find(pp)==unique_prism_verts.end())
                            {
                                unique_prism_verts.insert(pp);
                                loc2glob_final_verts[locp]=pp;
                                glob2loc_final_verts[pp]=locp;
                                locp++;
                            }
                        }
                        
                    }
                }
                
                //OutputMesh_MMG(mmgMesh_TET,offset_el,offset_el,"OuterVolume_TET.dat");
                
                int refer;
                int tellert = 0;
                int tellert2 = 0;
                std::map<int,std::vector< int* > > bound_tet;
                double m11,m12,m13,m22,m23,m33;
                int ut = 0;
                int wtel = 0;
                std::map<int,int> loc2glob_hyb;
                std::map<int,int> glob2loc_hyb;
                
                std::map<int,std::set<int> > newvert2elem;
                std::map<int,std::vector<int> > newvert2vert;
                std::map<int,int*> bndtrisVol;
                std::map<int,int> bndtrisVolRef;
                int tra = 0;
                int st = 0;
                int stt = 0;
                int gnt = 0;
                int yte = 0;
                int missing = 0;

                double vxc = 0;
                double vyc = 0;
                double vzc = 0;
                
                double vxf = 0;
                double vyf = 0;
                double vzf = 0;

                std::vector<double> oris;
                
                std::cout << "Defining the tetrahedra mesh..."<<std::endl;
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref != 0)
                    {
                        // Determine the vertices that are newly introduced by the tesselation of the hexes into tets.
                        //std::vector<int> which;
                        //std::vector<int> indexes;
                        vxc=0;vyc=0;vzc=0;
                        for(int s=0;s<4;s++)
                        {
                            vxc = vxc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[0];
                            vyc = vyc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[1];
                            vzc = vzc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[2];
                            
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == 0)
                            {
                                newvert2elem[mmgMesh_TET->tetra[i].v[s]].insert(i);
                                
                                for(int t=0;t<4;t++)
                                {
                                    if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[t]].ref != 0
                                    && mmgMesh_TET->tetra[i].v[t]!=mmgMesh_TET->tetra[i].v[s])
                                    {
                                        newvert2vert[mmgMesh_TET->tetra[i].v[s]].push_back(mmgMesh_TET->tetra[i].v[t]);
                                    }
                                }
                            }
                            
                        }
                        

                        vxc = vxc/4;
                        vyc = vyc/4;
                        vzc = vzc/4;

                        int* tetras= new int[4];
                        for(int s=0;s<4;s++)
                        {
                            if(lv2gv_tet_mesh.find(mmgMesh_TET->tetra[i].v[s]-1)!=lv2gv_tet_mesh.end())
                            {
                                tetras[s] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[s]-1];
                            }
                            else
                            {
                                tetras[s] = -(s+1);
                            }
                        }

                        
                        int v0=tetras[0];
                        int v1=tetras[1];
                        int v2=tetras[2];
                        int v3=tetras[3];
                        
                        std::set<int> bface;
                        
                        std::set<int> tria0;
                        
                        tria0.insert(v0);
                        tria0.insert(v2);
                        tria0.insert(v1);
                        //
                        
                        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria0];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v2+1;
                            tria[2] = v1+1;
                        
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        std::set<int> tria1;
                            
                        tria1.insert(v1);
                        tria1.insert(v2);
                        tria1.insert(v3);
                            
                        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria1];
                            int* tria = new int[3];
                            tria[0] = v1+1;
                            tria[1] = v2+1;
                            tria[2] = v3+1;
                            
                            
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        
                        std::set<int> tria2;
                        
                        tria2.insert(v0);
                        tria2.insert(v3);
                        tria2.insert(v2);
                        
                        if(tria_ref_map.find(tria2)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria2];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v3+1;
                            tria[2] = v2+1;
                            
                            
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra] = tria;
                            bndtrisVolRef[tra] = refer;
                            tra++;
                        }
                        
                        std::set<int> tria3;
                                                   
                        tria3.insert(v0);
                        tria3.insert(v1);
                        tria3.insert(v3);
                       
                        if(tria_ref_map.find(tria3)!=tria_ref_map.end())
                        {
                           refer = tria_ref_map[tria3];
                           int* tria = new int[3];
                           tria[0] = v0+1;
                           tria[1] = v1+1;
                           tria[2] = v3+1;
                            
                           bound_tet[refer].push_back(tria);
                           bndtrisVol[tra]      = tria;
                           bndtrisVolRef[tra]   = refer;
                           tra++;
                        }
                        
                        tria0.clear();
                        tria1.clear();
                        tria2.clear();
                        tria3.clear();
                        delete[] tetras;
                        tellert++;
                    }
                }

                int nTriangles_Vol = 0;
                std::map<int,std::vector< int* > >::iterator itrrb;
                for(itrrb=bound_tet.begin();itrrb!=bound_tet.end();itrrb++)
                {
                    nTriangles_Vol  = nTriangles_Vol+itrrb->second.size();
                }
                
                
                std::map<int,std::vector<int> >::iterator nve;
                double m11n=0.0;double m12n=0.0;double m13n=0.0;
                double m22n=0.0;double m23n=0.0;
                double m33n=0.0;
                int vgg=0;
                std::map<int,double*> newvert2metric;
                for(nve=newvert2vert.begin();nve!=newvert2vert.end();nve++)
                {
                    m11n=0.0;m12n=0.0;m13n=0.0;
                    m22n=0.0;m23n=0.0;
                    m33n=0.0;

                    for(int q=0;q<nve->second.size();q++)
                    {
                        vgg  = lv2gv_tet_mesh[nve->second[q]-1];

                        m11n = m11n + mv_g->getVal(vgg,0);
                        m12n = m12n + mv_g->getVal(vgg,1);
                        m13n = m13n + mv_g->getVal(vgg,2);
                        m22n = m22n + mv_g->getVal(vgg,3);
                        m23n = m23n + mv_g->getVal(vgg,4);
                        m33n = m33n + mv_g->getVal(vgg,5);
                    }
                    
                    m11n = m11n/nve->second.size();
                    m12n = m12n/nve->second.size();
                    m13n = m13n/nve->second.size();
                    m22n = m22n/nve->second.size();
                    m23n = m23n/nve->second.size();
                    m33n = m33n/nve->second.size();
                    
                    double* newM = new double[6];
                    
                    newM[0] = m11n;newM[1] = m12n;newM[2] = m13n;
                    newM[3] = m22n;newM[4] = m23n;newM[5] = m33n;
                    
                    newvert2metric[nve->first]=newM;
                }
                

                //====================================================================
                //====================================================================
                //====================================================================
                //====================================================================
                //              Begin Create hybrid mesh
                //====================================================================
                //====================================================================
                //====================================================================
                //====================================================================
                
                MMG5_pMesh mmgMesh_hyb = NULL;
                MMG5_pSol mmgSol_hyb   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
                MMG5_ARG_end);
                
                int nVertices_New  = mmgMesh_TET->np+unique_prism_verts.size()-cshell;
                int nbTets_New     = tellert;
                
                if ( MMG3D_Set_meshSize(mmgMesh_hyb,nVertices_New,nbTets_New,nbPrisms,nTriangles_BL+nTriangles_Vol,nQuads_BL,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh_hyb,mmgSol_hyb,MMG5_Vertex,mmgMesh_hyb->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);

                int i     = 0;
                int tt    = 1;
                int qt    = 1;
                
                std::cout << "Set the prisms in the mmgMesh..."<<std::endl;
                for(iter=mesh_topo_bl2->BLlayersPrisms.begin();
                   iter!=mesh_topo_bl2->BLlayersPrisms.end();iter++)
                {
                    int numit=iter->second.size();
                    
                    for(int p=0;p<numit;p++)
                    {
                        std::vector<int> prism = iter->second[p];

                        mmgMesh_hyb->prism[i+1].v[0] = prism[0]+1;
                        m11 = mv_g->getVal(prism[0],0);
                        m12 = mv_g->getVal(prism[0],1);
                        m13 = mv_g->getVal(prism[0],2);
                        m22 = mv_g->getVal(prism[0],3);
                        m23 = mv_g->getVal(prism[0],4);
                        m33 = mv_g->getVal(prism[0],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[0]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[1] = prism[1]+1;
                        m11 = mv_g->getVal(prism[1],0);
                        m12 = mv_g->getVal(prism[1],1);
                        m13 = mv_g->getVal(prism[1],2);
                        m22 = mv_g->getVal(prism[1],3);
                        m23 = mv_g->getVal(prism[1],4);
                        m33 = mv_g->getVal(prism[1],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[1]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[2] = prism[2]+1;
                        m11 = mv_g->getVal(prism[2],0);
                        m12 = mv_g->getVal(prism[2],1);
                        m13 = mv_g->getVal(prism[2],2);
                        m22 = mv_g->getVal(prism[2],3);
                        m23 = mv_g->getVal(prism[2],4);
                        m33 = mv_g->getVal(prism[2],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[2]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[3] = prism[3]+1;
                        m11 = mv_g->getVal(prism[3],0);
                        m12 = mv_g->getVal(prism[3],1);
                        m13 = mv_g->getVal(prism[3],2);
                        m22 = mv_g->getVal(prism[3],3);
                        m23 = mv_g->getVal(prism[3],4);
                        m33 = mv_g->getVal(prism[3],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[3]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[4] = prism[4]+1;
                        m11 = mv_g->getVal(prism[4],0);
                        m12 = mv_g->getVal(prism[4],1);
                        m13 = mv_g->getVal(prism[4],2);
                        m22 = mv_g->getVal(prism[4],3);
                        m23 = mv_g->getVal(prism[4],4);
                        m33 = mv_g->getVal(prism[4],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[4]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[5] = prism[5]+1;
                        m11 = mv_g->getVal(prism[5],0);
                        m12 = mv_g->getVal(prism[5],1);
                        m13 = mv_g->getVal(prism[5],2);
                        m22 = mv_g->getVal(prism[5],3);
                        m23 = mv_g->getVal(prism[5],4);
                        m33 = mv_g->getVal(prism[5],5);
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[5]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].ref  = 0;

                        std::set<int> tria0;
                        std::set<int> tria1;
                        std::set<int> quad0;
                        std::set<int> quad1;
                        std::set<int> quad2;
                        
                        tria0.insert(prism[0]);
                        tria0.insert(prism[1]);
                        tria0.insert(prism[2]);
                        
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria0];
                            mmgMesh_hyb->tria[tt].v[0] = prism[0]+1;
                            mmgMesh_hyb->tria[tt].v[1] = prism[1]+1;
                            mmgMesh_hyb->tria[tt].v[2] = prism[2]+1;
                            mmgMesh_hyb->tria[tt].ref  = refer;
                            tt++;
                        }
                        tria1.insert(prism[3]);
                        tria1.insert(prism[4]);
                        tria1.insert(prism[5]);

                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria1];
                            mmgMesh_hyb->tria[tt].v[0] = prism[3]+1;
                            mmgMesh_hyb->tria[tt].v[1] = prism[4]+1;
                            mmgMesh_hyb->tria[tt].v[2] = prism[5]+1;
                            mmgMesh_hyb->tria[tt].ref  = refer;
                            tt++;
                        }
                        
                        quad0.insert(prism[0]);
                        quad0.insert(prism[2]);
                        quad0.insert(prism[4]);
                        quad0.insert(prism[3]);
                       
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(quad_ref_map.find(quad0)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad0];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                        }
                        quad1.insert(prism[1]);
                        quad1.insert(prism[5]);
                        quad1.insert(prism[4]);
                        quad1.insert(prism[2]);
                        
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(quad_ref_map.find(quad1)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad1];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[1]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                        }
                        quad2.insert(prism[0]);
                        quad2.insert(prism[3]);
                        quad2.insert(prism[5]);
                        quad2.insert(prism[1]);
                        
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(quad_ref_map.find(quad2)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad2];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[1]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                        }
                       
                        tria0.clear();
                        tria1.clear();
                        quad0.clear();
                        quad1.clear();
                        quad2.clear();
                        i++;
                    }
                }
                
                for(int i=0;i<bndtrisVol.size();i++)
                {
                    mmgMesh_hyb->tria[tt].v[0] = bndtrisVol[i][0];
                    mmgMesh_hyb->tria[tt].v[1] = bndtrisVol[i][1];
                    mmgMesh_hyb->tria[tt].v[2] = bndtrisVol[i][2];
                    
                    mmgMesh_hyb->tria[tt].ref = bndtrisVolRef[i];
                    tt++;
                }
                
                delete mesh_topo_bl2;
                bndtrisVol.clear();
                
                // tets first then prisms.
                // verts building up the tets first then the verts building up the prisms.
                 
                int tet   = 0;
                int off_v = 0;
                int www   = 0;
                int nbVertices = xcn_g->getNrow();
                for(int i=0;i<nbVertices;i++)
                {
                    mmgMesh_hyb->point[i+1].c[0] = xcn_g->getVal(i,0);
                    mmgMesh_hyb->point[i+1].c[1] = xcn_g->getVal(i,1);
                    mmgMesh_hyb->point[i+1].c[2] = xcn_g->getVal(i,2);
                    
                    m11 = mv_g->getVal(i,0);
                    m12 = mv_g->getVal(i,1);
                    m13 = mv_g->getVal(i,2);
                    m22 = mv_g->getVal(i,3);
                    m23 = mv_g->getVal(i,4);
                    m33 = mv_g->getVal(i,5);
                    
                    if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
                }
                
                std::map<int,double*>::iterator v2m;
                //std::map<int,int> locNew2globNew;
                std::map<int,int> globNew2locNew;
                int p = 0;
                for(v2m=newvert2metric.begin();v2m!=newvert2metric.end();v2m++)
                {
                    mmgMesh_hyb->point[nbVertices+p+1].c[0] = mmgMesh_TET->point[v2m->first].c[0];
                    mmgMesh_hyb->point[nbVertices+p+1].c[1] = mmgMesh_TET->point[v2m->first].c[1];
                    mmgMesh_hyb->point[nbVertices+p+1].c[2] = mmgMesh_TET->point[v2m->first].c[2];
                    
                    m11 = v2m->second[0];
                    m12 = v2m->second[1];
                    m13 = v2m->second[2];
                    m22 = v2m->second[3];
                    m23 = v2m->second[4];
                    m33 = v2m->second[5];
                    
                    //locNew2globNew[nbVertices+p+1] = v2m->first;
                    globNew2locNew[v2m->first]     = nbVertices+p+1;
                    if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,nbVertices+p+1) != 1 ) exit(EXIT_FAILURE);

                    p++;
                }

                std::cout << "Set the tetrahedra in the mmgMesh..."<<std::endl;

                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref == 20)
                    {
                        tet++;
                        
                        for(int s=0;s<4;s++)
                        {
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == 0)
                            {
                                m11 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][0];
                                m12 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][1];
                                m13 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][2];
                                m22 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][3];
                                m23 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][4];
                                m33 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][5];
                                
                                int locNew = globNew2locNew[mmgMesh_TET->tetra[i].v[s]];

                                mmgMesh_hyb->tetra[tet].v[s] = locNew;

                                if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,locNew) != 1 ) exit(EXIT_FAILURE);
                                
                                //here++;
                            }
                            else
                            {
                                int vg = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[s]-1];
                                mmgMesh_hyb->tetra[tet].v[s] = vg+1;

                                m11 = mv_g->getVal(vg,0);
                                m12 = mv_g->getVal(vg,1);
                                m13 = mv_g->getVal(vg,2);
                                m22 = mv_g->getVal(vg,3);
                                m23 = mv_g->getVal(vg,4);
                                m33 = mv_g->getVal(vg,5);
                                
                                if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,vg+1) != 1 ) exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
                
                delete xcn_g;
                delete ien_g;
                delete ifn_g;
                delete ife_g;
                delete if_ref_g;
                delete ief_g;
                delete mv_g;
                  
                lv2gv_tet_mesh.clear(); 
                globNew2locNew.clear();
                newvert2elem.clear();
                newvert2vert.clear();
                newvert2metric.clear();
                
                MMG3D_Free_all(MMG5_ARG_start,
                               MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppSols,&mmgSol_TET,
                               MMG5_ARG_end);
             
                //MMG3D_Set_handGivenMesh(mmgMesh_hyb);
                if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, metric_inputs[0]) != 1 )    exit(EXIT_FAILURE);

                //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
                MMG3D_Set_dparameter( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_DPARAM_hgradreq , -1 );
                std::cout<<"Start the adaptation of the tetrahedra..."<<std::endl;
                int ier = MMG3D_mmg3dlib(mmgMesh_hyb,mmgSol_hyb);
                std::cout<<"Finished the adaptation of the tetrahedra..."<<std::endl;

                MMG5_pMesh mmgMesh_TETCOPY = NULL;
                MMG5_pSol mmgSol_TETCOPY   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppMet,&mmgSol_TETCOPY,
                MMG5_ARG_end);
      
                if ( MMG3D_Set_meshSize(mmgMesh_TETCOPY,mmgMesh_hyb->np,mmgMesh_hyb->ne,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
                
                //if ( MMG3D_Set_solSize(mmgMesh_TETCOPY,mmgSol_TETCOPY,MMG5_Vertex,mmgMesh_TETCOPY->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                for(int i=0;i<mmgMesh_hyb->np;i++)
                {
                    mmgMesh_TETCOPY->point[i+1].c[0] = mmgMesh_hyb->point[i+1].c[0];
                    mmgMesh_TETCOPY->point[i+1].c[1] = mmgMesh_hyb->point[i+1].c[1];
                    mmgMesh_TETCOPY->point[i+1].c[2] = mmgMesh_hyb->point[i+1].c[2];
                    
                    mmgMesh_TETCOPY->point[i+1].ref  = 0;
                }
                
                for(int i=0;i<mmgMesh_hyb->ne;i++)
                {
                    mmgMesh_TETCOPY->tetra[i+1].v[0] = mmgMesh_hyb->tetra[i+1].v[0];
                    mmgMesh_TETCOPY->tetra[i+1].v[1] = mmgMesh_hyb->tetra[i+1].v[1];
                    mmgMesh_TETCOPY->tetra[i+1].v[2] = mmgMesh_hyb->tetra[i+1].v[2];
                    mmgMesh_TETCOPY->tetra[i+1].v[3] = mmgMesh_hyb->tetra[i+1].v[3];
                    mmgMesh_TETCOPY->tetra[i+1].ref  = 0;
                }
                
                
                std::cout<<"Started writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
                OutputMesh_MMG_Slice(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
                std::cout<<"Finished writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
                
                delete bmap;

                MMG3D_Free_all(MMG5_ARG_start,
                               MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppSols,&mmgSol_TETCOPY,
                               MMG5_ARG_end);

                std::cout<<"Started writing the adapted hybrid mesh in US3D format..."<<std::endl;
                //WriteUS3DGridFromMMG_it0(mmgMesh_hyb, us3d, bnd_face_map);
                WriteUS3DGridFromMMG_it0(mmgMesh_hyb, mmgSol_hyb, us3d);
                std::cout<<"Finished writing the adapted hybrid mesh in US3D format..."<<std::endl;
                //
            }
        }
        
        /**/
        MPI_Finalize();
        
    }
     
    return 0;
}