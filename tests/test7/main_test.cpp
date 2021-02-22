#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include <iomanip>

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
//    const char* fn_grid="../test_mesh/cyl_tess/grid.h5";
//    const char* fn_conn="../test_mesh/cyl_tess/conn.h5";
//    const char* fn_data="../test_mesh/cyl_tess/data.h5";
    
    const char* fn_grid="../test_mesh/cyl_a/grid.h5";
    const char* fn_conn="../test_mesh/cyl_a/conn.h5";
    const char* fn_data="../test_mesh/cyl_a/data.h5";
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);
    const char* fn_metric = "metric.inp";
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
        
    int Nel_part = us3d->ien->getNrow();

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    for(int i=0;i<Nel_part;i++)
    {
        Ui->setVal(i,0,us3d->interior->getVal(i,varia));
    }
    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
    
//    integer(KIND=US3D_GINT), parameter :: ET_TRI= 1      ! Triangle
//    integer(KIND=US3D_GINT), parameter :: ET_TET= 2      ! Tetrahedron
//    integer(KIND=US3D_GINT), parameter :: ET_QAD= 3      ! Quadrilateral
//    integer(KIND=US3D_GINT), parameter :: ET_HEX= 4      ! Hexahedral
//    integer(KIND=US3D_GINT), parameter :: ET_PYR= 5      ! Pyramid
//    integer(KIND=US3D_GINT), parameter :: ET_PRS= 6      ! Prism
//    integer(KIND=US3D_GINT), parameter :: ET_SEG= 7      ! Line segment
                                                                         
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, comm);
    
    double duration = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Ptime = 0.0;
    MPI_Allreduce(&duration, &Ptime, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
    std::cout << "Timing partitioning: " << duration << std::endl;
    }
    //Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);

//
    
    std::vector<int> LocElem        = P->getLocElem();
    std::vector<int> LocElemNv      = P->getLocElemNv();
    std::map<int,int> LocElem2Nv      = P->getLocElem2Nv();
    std::vector<double> Uvaria      = P->getLocElemVaria();

    std::map<int,double> Ui_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid     = LocElem[i];
        UvariaV     = Uvaria[i];
        Ui_map[gid] = UvariaV;
    }

    Domain* pDom = P->getPartitionDomain();
    std::vector<Vert> Verts = P->getLocalVerts();
    
    std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Ui_map,comm);

    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }

    t = clock();
    
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,UauxNew,gB,comm);
    
    //std::map<int,Array<double>* > dUdXi = ComputedUdx_MGG(P,UauxNew,meshTopo,gB,comm);
    
    double Gtiming = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Gmax_time = 0.0;
    MPI_Allreduce(&Gtiming, &Gmax_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        std::cout << "Finished reconstructing the gradient... " << std::endl;
        std::cout << "Timing gradient reconstruction... " << Gmax_time << std::endl;
    }
    
    Array<int>* lE2gE     = new Array<int>(dUdXi.size(),1);
    Array<double>* dUidxi = new Array<double>(dUdXi.size(),1);
    Array<double>* dUidyi = new Array<double>(dUdXi.size(),1);
    Array<double>* dUidzi = new Array<double>(dUdXi.size(),1);

    std::map<int,double> dUidxi_map;
    std::map<int,double> dUidyi_map;
    std::map<int,double> dUidzi_map;

    std::map<int,Array<double>* >::iterator grit;
    int ti=0;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        lE2gE->setVal(ti,0,grit->first);
        dUidxi->setVal(ti,0,grit->second->getVal(0,0));
        dUidyi->setVal(ti,0,grit->second->getVal(1,0));
        dUidzi->setVal(ti,0,grit->second->getVal(2,0));
//        dUidxi_map[grit->first]=grit->second->getVal(0,0);
//        dUidyi_map[grit->first]=grit->second->getVal(1,0);
//        dUidzi_map[grit->first]=grit->second->getVal(2,0);

        ti++;
    }
    
    int nlElem = us3d->ien->getNrow();
    int nElem  = us3d->ien->getNglob();
    int nvg    = us3d->xcn->getNglob();
    
    std::vector<double> GuX_loc;
    std::vector<int> vids;
    std::vector<double> Uids;
    std::vector<double> Hids;
    int gid,lid;
    int nval = 6;
    std::set<int> vdone;
//
    Array<int>*  lE2gE_g;
    Array<double>*  GuX_g;
    Array<double>*  GuX_gr;
    Array<double>*  GuY_g;
    Array<double>*  GuZ_g;

    int nEl_glob = us3d->ien->getNglob();

    if(world_rank == 0)
    {
        lE2gE_g     = new Array<int>(nEl_glob,1);
        GuX_g       = new Array<double>(nEl_glob,1);
        GuX_gr      = new Array<double>(nEl_glob,1);
        GuY_g       = new Array<double>(nEl_glob,1);
        GuZ_g       = new Array<double>(nEl_glob,1);
    }
    else
    {
        lE2gE_g     = new Array<int>(1,1);
        GuX_g       = new Array<double>(1,1);
        GuX_gr      = new Array<double>(1,1);
        GuY_g       = new Array<double>(1,1);
        GuZ_g       = new Array<double>(1,1);
    }

    int* G_nlocs      = new int[world_size];
    int* red_G_nlocs  = new int[world_size];
    int* G_offsets    = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        G_nlocs[i] = 0;
        
        if(i==world_rank)
        {
            G_nlocs[i] = dUidxi->getNrow();
        }
        else
        {
            G_nlocs[i] = 0;
        }
    }

    MPI_Allreduce(G_nlocs, red_G_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    
    int offset = 0;
    for(i=0;i<world_size;i++)
    {
        G_offsets[i] = offset;
        offset = offset+red_G_nlocs[i];
    }

    MPI_Gatherv(&lE2gE->data[0],
                lE2gE->getNrow(),
                MPI_INT,
                &lE2gE_g->data[0],
                red_G_nlocs,
                G_offsets,
                MPI_INT, 0, comm);
    
    
    MPI_Gatherv(&dUidxi->data[0],
                dUidxi->getNrow(),
                MPI_DOUBLE,
                &GuX_g->data[0],
                red_G_nlocs,
                G_offsets,
                MPI_DOUBLE, 0, comm);
    
    
    if(world_rank == 0)
    {
        
        std::vector<double> GuX_ref = ReadReferenceData();
        
        for(int i=0;i<nEl_glob;i++)
        {
            int gid = lE2gE_g->getVal(i,0);
            GuX_gr->setVal(gid,0,GuX_g->getVal(i,0));
        }
        
        int flip = 0;
        double err = 1.0e-08;
        for(int i=0;i<nEl_glob;i++)
        {
            double diff = fabs(GuX_gr->getVal(i,0)-GuX_ref[i]);

            if(diff>err)
            {
                std::cout << std::setprecision(16) << i << " " << diff << " " << GuX_gr->getVal(i,0) << " " << GuX_ref[i] << std::endl;

                flip = 1;
            }
        }

        if(flip == 1)
        {
            std::cout << " --::-- Parallel gradient reconstruction test has FAILED. --::-- " << std::endl;
        }
        if(flip == 0)
        {
            std::cout << " --::-- Parallel gradient reconstruction test has PASSED. --::-- " << std::endl;
        }
    }
    
    
    
    
    std::map<int,double > dUdxauxNew  = P->CommunicateAdjacentDataUS3D(dUidxi_map,comm);
    std::map<int,double > dUdyauxNew  = P->CommunicateAdjacentDataUS3D(dUidyi_map,comm);
    std::map<int,double > dUdzauxNew  = P->CommunicateAdjacentDataUS3D(dUidzi_map,comm);

    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUdxauxNew,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUdyauxNew,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUdzauxNew,gB,comm);

    std::map<int,double> d2udx2_map,d2udxy_map,d2udxz_map,
                         d2udyx_map,d2udy2_map,d2udyz_map,
                         d2udzx_map,d2udzy_map,d2udz2_map;

    std::map<int,Array<double>*> Hess_map;

    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;

    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
        Array<double>* Hess = new Array<double>(3,3);

        Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
        Hess->setVal(0,1,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(0,2,dU2dXi2[gid]->getVal(2,0));

        Hess->setVal(1,0,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(1,1,dU2dYi2[gid]->getVal(1,0));
        Hess->setVal(1,2,dU2dYi2[gid]->getVal(2,0));

        Hess->setVal(2,0,dU2dXi2[gid]->getVal(2,0));
        Hess->setVal(2,1,dU2dYi2[gid]->getVal(2,0));
        Hess->setVal(2,2,dU2dZi2[gid]->getVal(2,0));
//        if(world_rank == 0)
//        {
//            std::cout << "============="<<std::endl;
//            std::cout << Hess->getVal(0,0) << " " << Hess->getVal(0,1) << " " << Hess->getVal(0,2) << std::endl;
//            std::cout << Hess->getVal(1,0) << " " << Hess->getVal(1,1) << " " << Hess->getVal(1,2) << std::endl;
//            std::cout << Hess->getVal(2,0) << " " << Hess->getVal(2,1) << " " << Hess->getVal(2,2) << std::endl;
//            std::cout << "============="<<std::endl;
//        }
        Hess_map[gid] = Hess;
        t++;
    }
    
    std::map<int,Array<double>* > Hess_vm = P->ReduceMetricToVertices(Hess_map);

    std::map<int,double> u_vmap = P->ReduceFieldToVertices(Ui_map);

    std::map<int,double> dudx_vmap = P->ReduceFieldToVertices(dUidxi_map);
    
    std::map<int,Array<double>* > metric = ComputeMetric(P,metric_inputs, Hess_vm, comm);

    //Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, metric, comm);
    
    std::vector<std::vector<int> > Elements = pDom->Elements;
    std::vector<std::vector<int> > Hexes    = pDom->Hexes;
    std::vector<std::vector<int> > Tetras   = pDom->Tetras;
    std::vector<std::vector<int> > Prisms   = pDom->Prisms;
    
    int nTetras_loc   = Tetras.size();
    int nHexes_loc    = Hexes.size();
    int nPrisms_loc   = Prisms.size();
    int nHexes_glob   = 0;
    int nTetras_glob  = 0;
    int nPrisms_glob  = 0;
    
    MPI_Allreduce(&nTetras_loc, &nTetras_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nHexes_loc,  &nHexes_glob,  1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nPrisms_loc, &nPrisms_glob, 1, MPI_INT, MPI_SUM, comm);

    
    std::cout << "tetras = " << nTetras_loc << " " << nTetras_glob << std::endl;
    std::cout << "prisms = " << nPrisms_loc << " " << nPrisms_glob << std::endl;
    std::cout << "hexes  = " << nHexes_loc  << " " << nHexes_glob  << std::endl;

    Array<int>* Ate = new Array<int>(nTetras_loc,4);
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    int gv = 0;
    for(int i=0;i<nTetras_loc;i++)
    {
        for(int j=0;j<4;j++)
        {
            gv = lpartv2gv[Tetras[i][j]];
            Ate->setVal(i,j,gv);
        }
    }

    Array<int>* Apr = new Array<int>(nPrisms_loc,6);

    for(int i=0;i<nPrisms_loc;i++)
    {
        for(int j=0;j<6;j++)
        {
            gv = lpartv2gv[Prisms[i][j]];
            Apr->setVal(i,j,Prisms[i][j]);
        }
    }
    
    Array<int>* gAteRoot   = GatherArrayOnRoot<int>(Ate,comm,info);
    Array<int>* gAprRoot   = GatherArrayOnRoot<int>(Apr,comm,info);
    
    Mesh* us3dRoot = ReduceMeshToRoot(us3d->ien,
                                      us3d->ief,
                                      us3d->xcn,
                                      us3d->ifn,
                                      us3d->ife,
                                      us3d->if_ref,
                                      comm,info);
    
//    for(int i=0;i<nrow_ifn;i++)
//    {
//        ref = if_ref_g->getVal(i,0);
//        faceid = i;
//        if(ref != 2)
//        {
//            bnd_face_map[ref].push_back(faceid);
//        }
//    }

    //std::cout << "# tets global " << gAteRoot->getNrow() << " #prisms global " << gAprRoot->getNrow() << std::endl;
    
    //Array<double>* xcn_g = ReadDataSetFromFile<double>(fn_grid,"xcn");
    
    
    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
    //std::map<int,int> lpartv2gv   = pDom->lpartv2gv;
    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    
    std::ofstream myfilet;
    myfilet.open("parttet_" + std::to_string(world_rank) + ".dat");
    myfilet << "TITLE=\"new_volume.tec\"" << std::endl;
    myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
    myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << Tetras.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
    
    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid  = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfilet << Verts[loc_vid].x << " " << Verts[loc_vid].y << " " << Verts[loc_vid].z
         << " " << metric[glob_vid]->getVal(0,0)<< " " << metric[glob_vid]->getVal(0,1)<<
            " " << metric[glob_vid]->getVal(0,2)<< " " << metric[glob_vid]->getVal(1,1)<<
            " " << metric[glob_vid]->getVal(1,2)<< " " << metric[glob_vid]->getVal(2,2) <<  std::endl;
//        myfilet << Verts[loc_vid].x << " " << Verts[loc_vid].y << " " << Verts[loc_vid].z
//         << " " << u_vmap[glob_vid] << " " << dudx_vmap[glob_vid]<< " " <<  std::endl;
    }
    
    for(int i=0;i<Tetras.size();i++)
    {
        myfilet << Tetras[i][0]+1 << " " << Tetras[i][1]+1 << " " << Tetras[i][2]+1 << " " << Tetras[i][3]+1 << std::endl;
    }

    
    myfilet.close();
    
    std::cout << world_rank << " # prisms = " << Prisms.size() << "  " << " # tetras = " << Tetras.size() << std::endl;
    /*
    
    Array<int>*  lE2gE_g;
    Array<double>*  GuX_g;
    Array<double>*  GuX_gr;
    Array<double>*  GuY_g;
    Array<double>*  GuZ_g;

    int nEl_glob = us3d->ien->getNglob();

    if(world_rank == 0)
    {
        lE2gE_g     = new Array<int>(nEl_glob,1);
        GuX_g       = new Array<double>(nEl_glob,1);
        GuX_gr      = new Array<double>(nEl_glob,1);
        GuY_g       = new Array<double>(nEl_glob,1);
        GuZ_g       = new Array<double>(nEl_glob,1);
    }
    else
    {
        lE2gE_g     = new Array<int>(1,1);
        GuX_g       = new Array<double>(1,1);
        GuX_gr      = new Array<double>(1,1);
        GuY_g       = new Array<double>(1,1);
        GuZ_g       = new Array<double>(1,1);
    }
    
    
//    MMG5_pMesh mmgMesh_hyb = NULL;
//    MMG5_pSol mmgSol_hyb   = NULL;
//
//    MMG3D_Init_mesh(MMG5_ARG_start,
//    MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
//    MMG5_ARG_end);
    
    
    
//    //MMG3D_Set_handGivenMesh(mmgMesh_hyb);
//    if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, 3.0) != 1 )    exit(EXIT_FAILURE);
//
//    //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
//    MMG3D_Set_dparameter( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_DPARAM_hgradreq , -1 );
//    std::cout<<"Start the adaptation of the tetrahedra..."<<std::endl;
//    int ier = MMG3D_mmg3dlib(mmgMesh_hyb,mmgSol_hyb);
//    std::cout<<"Finished the adaptation of the tetrahedra..."<<std::endl;
//
//    MMG5_pMesh mmgMesh_TETCOPY = NULL;
//    MMG5_pSol mmgSol_TETCOPY   = NULL;
//
//    MMG3D_Init_mesh(MMG5_ARG_start,
//    MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppMet,&mmgSol_TETCOPY,
//    MMG5_ARG_end);
//
//    if ( MMG3D_Set_meshSize(mmgMesh_TETCOPY,mmgMesh_hyb->np,mmgMesh_hyb->ne,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
//
//    //if ( MMG3D_Set_solSize(mmgMesh_TETCOPY,mmgSol_TETCOPY,MMG5_Vertex,mmgMesh_TETCOPY->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
//
//    for(int i=0;i<mmgMesh_hyb->np;i++)
//    {
//        mmgMesh_TETCOPY->point[i+1].c[0] = mmgMesh_hyb->point[i+1].c[0];
//        mmgMesh_TETCOPY->point[i+1].c[1] = mmgMesh_hyb->point[i+1].c[1];
//        mmgMesh_TETCOPY->point[i+1].c[2] = mmgMesh_hyb->point[i+1].c[2];
//
//        mmgMesh_TETCOPY->point[i+1].ref  = 0;
//    }
//
//    for(int i=0;i<mmgMesh_hyb->ne;i++)
//    {
//        mmgMesh_TETCOPY->tetra[i+1].v[0] = mmgMesh_hyb->tetra[i+1].v[0];
//        mmgMesh_TETCOPY->tetra[i+1].v[1] = mmgMesh_hyb->tetra[i+1].v[1];
//        mmgMesh_TETCOPY->tetra[i+1].v[2] = mmgMesh_hyb->tetra[i+1].v[2];
//        mmgMesh_TETCOPY->tetra[i+1].v[3] = mmgMesh_hyb->tetra[i+1].v[3];
//        mmgMesh_TETCOPY->tetra[i+1].ref  = 0;
//    }
//
//
//    std::cout<<"Started writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
//    OutputMesh_MMG(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
//    std::cout<<"Finished writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
//
//    MMG3D_Free_all(MMG5_ARG_start,
//                   MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppSols,&mmgSol_TETCOPY,
//                   MMG5_ARG_end);
//
//    MMG3D_Free_all(MMG5_ARG_start,
//                   MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppSols,&mmgSol_TET,
//                   MMG5_ARG_end);
//
//    std::cout<<"Started writing the adapted hybrid mesh in US3D format..."<<std::endl;
//    WriteUS3DGridFromMMG(mmgMesh_hyb, us3d, bnd_face_map);
//    std::cout<<"Finished writing the adapted hybrid mesh in US3D format..."<<std::endl;
//    //
//
    */
    MPI_Finalize();
}
