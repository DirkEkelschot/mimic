#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include <iomanip>

struct ReferenceMesh{
    std::vector<std::vector<double> > Nodes;
    std::vector<std::vector<double> > Elements;
};

ReferenceMesh* ReadReferenceMesh()
{
    
    ReferenceMesh* refmesh = new ReferenceMesh;
    std::ifstream fin_v;
    fin_v.open("Nodes.ref");
    // Read the file row by row
    std::vector<std::vector<double> > Vref;
    
    std::vector<double> row_v(3);
    while(fin_v >> row_v[0] >> row_v[1] >> row_v[2] )
    {
        refmesh->Nodes.push_back(row_v);
    }
    fin_v.close();
    
    
    
    
    
    std::ifstream fin_e;
    fin_e.open("Elements.ref");
    // Read the file row by row
    std::vector<std::vector<double> > Eref;
    
    std::vector<double> row_e(4);
    while(fin_e >> row_e[0] >> row_e[1] >> row_e[2] >> row_e[3] )
    {
        refmesh->Elements.push_back(row_e);
    }
    fin_e.close();
    
    return refmesh;
    
}


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
    
    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    
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
    
    std::vector<int> LocElem        = P->getLocElem();
    std::vector<int> LocElemNv      = P->getLocElemNv();
    std::map<int,int> LocElem2Nv    = P->getLocElem2Nv();
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
        dUidxi_map[grit->first]=grit->second->getVal(0,0);
        dUidyi_map[grit->first]=grit->second->getVal(1,0);
        dUidzi_map[grit->first]=grit->second->getVal(2,0);

        ti++;
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
        
        Hess_map[gid] = Hess;
        
        t++;
    }
    
    std::map<int,Array<double>* > Hess_vm = P->ReduceMetricToVertices(Hess_map);

    std::map<int,double> u_vmap = P->ReduceFieldToVertices(Ui_map);

    std::map<int,double> dudx_vmap = P->ReduceFieldToVertices(dUidxi_map);
    
    std::map<int,Array<double>* > metric = ComputeMetric(P,metric_inputs, Hess_vm, comm);

    Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, metric, comm);
    
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

    Array<int>* Ate = new Array<int>(nTetras_loc,4);
    std::map<int,int> lv2gpv     = pDom->lv2gpv;
    int gv = 0;
    int lv = 0;
    

    for(int i=0;i<nTetras_loc;i++)
    {
        for(int j=0;j<4;j++)
        {
            gv = lv2gpv[Tetras[i][j]];
            Ate->setVal(i,j,gv);
        }
    }

    Array<int>* Apr = new Array<int>(nPrisms_loc,6);

    for(int i=0;i<nPrisms_loc;i++)
    {
        for(int j=0;j<6;j++)
        {
            gv = lv2gpv[Prisms[i][j]];
            Apr->setVal(i,j,gv);
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
    
    if(world_rank == 0)
    {
        MMG5_pMesh mmgMesh_hyb = NULL;
        MMG5_pSol mmgSol_hyb   = NULL;
        
        MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
        MMG5_ARG_end);
        
        //int nVerts    = us3dRoot->xcn->getNrow();
        int nTets     = gAteRoot->getNrow();
        int nPrisms   = gAprRoot->getNrow();
        
        std::map<int,int> l2g_tetv;
        std::map<int,int> g2l_tetv;
        std::set<int> uv_tet_set;
        int lv_uv = 0;
        for(int i=0;i<gAteRoot->getNrow();i++)
        {
            for(int j=0;j<gAteRoot->getNcol();j++)
            {
                int gv = gAteRoot->getVal(i,j);
                if(uv_tet_set.find(gv)==uv_tet_set.end())
                {
                    uv_tet_set.insert(gv);
                    l2g_tetv[lv_uv]=gv;
                    g2l_tetv[gv]=lv_uv;
                    lv_uv++;
                }
            }
        }

        
        
        for(int i=0;i<gAprRoot->getNrow();i++)
        {
            for(int j=0;j<gAprRoot->getNcol();j++)
            {
                int gv = gAprRoot->getVal(i,j);
                if(uv_tet_set.find(gv)==uv_tet_set.end())
                {
                    uv_tet_set.insert(gv);
                    l2g_tetv[lv_uv]=gv;
                    g2l_tetv[gv]=lv_uv;
                    lv_uv++;
                }
            }
        }
        
        
        
        
        int nVerts  = l2g_tetv.size();
        if ( MMG3D_Set_meshSize(mmgMesh_hyb,nVerts,nTets,nPrisms,0,0,0) != 1 )  exit(EXIT_FAILURE);
        
        if ( MMG3D_Set_solSize(mmgMesh_hyb,mmgSol_hyb,MMG5_Vertex,mmgMesh_hyb->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
        
        
        for(int i=0;i<gAprRoot->getNrow();i++)
        {
            for(int j=0;j<gAprRoot->getNcol();j++)
            {
                int gv = gAprRoot->getVal(i,j);
                mmgMesh_hyb->prism[i+1].v[j] = g2l_tetv[gv]+1;
            }
            mmgMesh_hyb->tetra[i+1].ref = 3;
        }
        
        for(int i=0;i<gAteRoot->getNrow();i++)
        {
            for(int j=0;j<gAteRoot->getNcol();j++)
            {
                int gv = gAteRoot->getVal(i,j);
                mmgMesh_hyb->tetra[i+1].v[j] = g2l_tetv[gv]+1;
            }
            mmgMesh_hyb->tetra[i+1].ref = 3;
        }
        
        
        std::map<int,int>::iterator itm;
        
        int i = 0;
        for(itm=l2g_tetv.begin();itm!=l2g_tetv.end();itm++)
        {
            int lvv = itm->first;
            int gvv = itm->second;
                        
            mmgMesh_hyb->point[i+1].c[0] = us3dRoot->xcn->getVal(gvv,0);
            mmgMesh_hyb->point[i+1].c[1] = us3dRoot->xcn->getVal(gvv,1);
            mmgMesh_hyb->point[i+1].c[2] = us3dRoot->xcn->getVal(gvv,2);
            mmgMesh_hyb->point[i+1].ref  = 1;
            
            double m11 = mv_g->getVal(gvv,0);
            double m12 = mv_g->getVal(gvv,1);
            double m13 = mv_g->getVal(gvv,2);
            double m22 = mv_g->getVal(gvv,3);
            double m23 = mv_g->getVal(gvv,4);
            double m33 = mv_g->getVal(gvv,5);
            
            if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
            
            i++;
        }
        
        if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, metric_inputs[0]) != 1 )    exit(EXIT_FAILURE);
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
        
        int flip = 0;
        double tol = 1.0e-05;
        
        ReferenceMesh* refmesh = ReadReferenceMesh();
//        std::ofstream myfile_nv;
//        myfile_nv.open("Nodes.ref");
        for(int i=0;i<mmgMesh_hyb->np;i++)
        {
            mmgMesh_TETCOPY->point[i+1].c[0] = mmgMesh_hyb->point[i+1].c[0];
            mmgMesh_TETCOPY->point[i+1].c[1] = mmgMesh_hyb->point[i+1].c[1];
            mmgMesh_TETCOPY->point[i+1].c[2] = mmgMesh_hyb->point[i+1].c[2];
            mmgMesh_TETCOPY->point[i+1].ref  = 0;
            
//            myfile_nv << mmgMesh_hyb->point[i+1].c[0] << " " << mmgMesh_hyb->point[i+1].c[1] << " " << mmgMesh_hyb->point[i+1].c[2] << std::endl;
            
            
            if(fabs(mmgMesh_hyb->point[i+1].c[0]-refmesh->Nodes[i][0])>tol || fabs(mmgMesh_hyb->point[i+1].c[1]-refmesh->Nodes[i][1])>tol || fabs(mmgMesh_hyb->point[i+1].c[2]-refmesh->Nodes[i][2])>tol)
            {
                flip = 1;
            }
        }
        
//        myfile_nv.close();
//        std::ofstream myfile_ne;
//        myfile_ne.open("Elements.ref");
        for(int i=0;i<mmgMesh_hyb->ne;i++)
        {
            mmgMesh_TETCOPY->tetra[i+1].v[0] = mmgMesh_hyb->tetra[i+1].v[0];
            mmgMesh_TETCOPY->tetra[i+1].v[1] = mmgMesh_hyb->tetra[i+1].v[1];
            mmgMesh_TETCOPY->tetra[i+1].v[2] = mmgMesh_hyb->tetra[i+1].v[2];
            mmgMesh_TETCOPY->tetra[i+1].v[3] = mmgMesh_hyb->tetra[i+1].v[3];
            mmgMesh_TETCOPY->tetra[i+1].ref  = 0;
            
//            myfile_ne << mmgMesh_hyb->tetra[i+1].v[0] << " " << mmgMesh_hyb->tetra[i+1].v[1] << " " << mmgMesh_hyb->tetra[i+1].v[2] << " " << mmgMesh_hyb->tetra[i+1].v[3] << std::endl;
            
            if( fabs(mmgMesh_hyb->tetra[i+1].v[0]-refmesh->Elements[i][0])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[1]-refmesh->Elements[i][1])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[2]-refmesh->Elements[i][2])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[3]-refmesh->Elements[i][3])>tol)
            {
                std::cout << fabs(mmgMesh_hyb->tetra[i+1].v[0]-refmesh->Elements[i][0]) << " "
                          << fabs(mmgMesh_hyb->tetra[i+1].v[1]-refmesh->Elements[i][1]) << " "
                          << fabs(mmgMesh_hyb->tetra[i+1].v[2]-refmesh->Elements[i][2]) << " "
                          << fabs(mmgMesh_hyb->tetra[i+1].v[3]-refmesh->Elements[i][3]) << std::endl;
                flip = 1;
            }
        }
        //myfile_ne.close();
        
        if(flip == 1)
        {
            std::cout << " --::-- Adaptation based on gradient reconstruction for a hybrid mesh has FAILED. --::-- " << std::endl;
        }
        if(flip == 0)
        {
            std::cout << " --::-- Adaptation based on gradient reconstruction for a hybrid mesh has PASSED. --::-- " << std::endl;
        }
        
//      std::cout << "tets = " << mmgMesh_TETCOPY->ne << " " << nTets << " " << nPrisms << std::endl;
        
        std::cout<<"Started writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
        OutputMesh_MMG(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
        std::cout<<"Finished writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
         
    }

    MPI_Finalize();
}
