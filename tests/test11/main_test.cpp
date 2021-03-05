#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
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

    
    
//    const char* fn_grid="../test_mesh/cylinder_hex/grid.h5";
//    const char* fn_conn="../test_mesh/cylinder_hex/conn.h5";
//    const char* fn_data="../test_mesh/cylinder_hex/data.h5";
    
//    const char* fn_grid="itn1/grid.h5";
//    const char* fn_conn="itn1/conn.h5";
//    const char* fn_data="itn1/data.h5";
    
    const char* fn_grid="../test_mesh/it1n/grid.h5";
    const char* fn_conn="../test_mesh/it1n/conn.h5";
    const char* fn_data="../test_mesh/it1n/data.h5";
    
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
    std::map<int,Array<double>*> Ui_map_arr;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid     = LocElem[i];
        UvariaV     = Uvaria[i];
        Array<double>* Uarr = new Array<double>(1,1);
        
        Uarr->setVal(0,0,UvariaV);
        Ui_map[gid] = UvariaV;
        Ui_map_arr[gid] = Uarr;
    }
    LocElem.clear();
    LocElemNv.clear();
    LocElem2Nv.clear();
    Uvaria.clear();
    
    std::vector<Vert*> Verts = P->getLocalVerts();
    std::map<int,std::map<int,double> > n2n = P->getNode2NodeMap();
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }

    
    
    P->AddStateVecForAdjacentElements(Ui_map_arr,1,comm);
    std::map<int,Array<double>* > u_vmap = P->ReduceStateVecToAllVertices(Ui_map_arr, 1);
    P->AddStateVecForAdjacentVertices(u_vmap, 1, comm);
            
    Domain* pDom = P->getPartitionDomain();
    
    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    
    std::map<int,std::map<int,double> >::iterator n2nit;
    double sum      = 0.0;
    double sum_dist = 0.0;
    double di = 0.0;
    double uval = 0.0;
    
    std::map<int,double > u_vmap_new;
    int vid,gvid;
    for(n2nit=n2n.begin();n2nit!=n2n.end();n2nit++)
    {
        gvid     = n2nit->first;
        sum      = 0.0;
        sum_dist = 0.0;
        uval     = 0.0;
        std::map<int,double>::iterator n2dit;
        for(n2dit=n2nit->second.begin();n2dit!=n2nit->second.end();n2dit++)
        {
            vid       = n2dit->first;
            di        = n2dit->second;
            sum_dist  = sum_dist+di;
            uval      = uval+u_vmap[vid]->getVal(0,0)*di;
        }
        
        u_vmap_new[gvid] = uval/sum_dist;
    }
    
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_Vrt_US3D(P,Ui_map,u_vmap_new,meshTopo,gB,comm);
    
    double Gtiming = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Gmax_time = 0.0;
    MPI_Allreduce(&Gtiming, &Gmax_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        std::cout << "Finished reconstructing the gradient... " << std::endl;
        std::cout << "Timing gradient reconstruction... " << Gmax_time << std::endl;
    }
    
    //============================================================
    std::cout << "Rank " << world_rank << " element gradient " << dUdXi.size() << std::endl;
    P->AddStateVecForAdjacentElements(dUdXi,3,comm);
    std::cout << "Rank " << world_rank << " element+adj gradient " << dUdXi.size() << std::endl;

    P->AddStateVecForAdjacentElements(dUdXi,3,comm);
    std::map<int,Array<double>* > dUdXi_vmap = P->ReduceStateVecToAllVertices(dUdXi,3);
    P->AddStateVecForAdjacentVertices(dUdXi_vmap,3,comm);

    
    
    //P->AddAdjacentVertexDataUS3D(dUdXi_vmap, comm);

    //============================================================
    //OOOOOOOOOOOOOOOOOOooooooooooooooOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOooooooooo
    //O             Begin reducing and smoothing data over vertices                O
    //OOOOOOOOOOOOOOOOOOooooooooooooooOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOooooooooo


    std::map<int,Array<double>* >::iterator itarr;
    int size_b = dUdXi_vmap.size();

    sum      = 0.0;
    sum_dist = 0.0;
    di       = 0.0;
    
    double duvalx = 0.0;
    double duvaly = 0.0;
    double duvalz = 0.0;
    
    std::map<int,double> dudx_vmap_new;
    std::map<int,double> dudy_vmap_new;
    std::map<int,double> dudz_vmap_new;
    
    for(n2nit=n2n.begin();n2nit!=n2n.end();n2nit++)
    {
        gvid       = n2nit->first;
        
        sum        = 0.0;
        duvalx     = 0.0;
        duvaly     = 0.0;
        duvalz     = 0.0;
        sum_dist   = 0.0;
        
        std::map<int,double>::iterator n2dit;
        for(n2dit=n2nit->second.begin();n2dit!=n2nit->second.end();n2dit++)
        {
            vid      = n2dit->first;
            di       = n2dit->second;
            sum_dist = sum_dist+di;
            
            duvalx     = duvalx+dUdXi_vmap[vid]->getVal(0,0)*di;
            duvaly     = duvaly+dUdXi_vmap[vid]->getVal(1,0)*di;
            duvalz     = duvalz+dUdXi_vmap[vid]->getVal(2,0)*di;
        }
        
        dudx_vmap_new[gvid] = duvalx/sum_dist;
        dudy_vmap_new[gvid] = duvaly/sum_dist;
        dudz_vmap_new[gvid] = duvalz/sum_dist;
        
    }
   
    //OOOOOOOOOOOOOOOOOOooooooooooooooOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOooooooooo
    //O               End reducing and smoothing data over vertices                O
    //OOOOOOOOOOOOOOOOOOooooooooooooooOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOooooooooo
    
    std::map<int,Array<double>* >::iterator grit;
    std::map<int,double> dUidxi_map;
    std::map<int,double> dUidyi_map;
    std::map<int,double> dUidzi_map;
    std::map<int,Array<double>* > dUidXi_map;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        dUidxi_map[grit->first]=grit->second->getVal(0,0);
        dUidyi_map[grit->first]=grit->second->getVal(1,0);
        dUidzi_map[grit->first]=grit->second->getVal(2,0);
        
        delete grit->second;
        
        i++;
    }
    
    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidxi_map,dudx_vmap_new,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidyi_map,dudy_vmap_new,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidzi_map,dudz_vmap_new,meshTopo,gB,comm);

//    dudx_vmap_new.clear();
//    dudy_vmap_new.clear();
//    dudz_vmap_new.clear();
    
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
    
    P->AddStateVecForAdjacentElements(Hess_map,6,comm);
    std::map<int,Array<double>* > hess_vmap = P->ReduceStateVecToAllVertices(Hess_map,6);
    P->AddStateVecForAdjacentVertices(hess_vmap,6,comm);
    
   
    sum      = 0.0;
    sum_dist = 0.0;
    di       = 0.0;
    
    double d2udx2val = 0.0;
    double d2udxyval = 0.0;
    double d2udxzval = 0.0;
    double d2udy2val = 0.0;
    double d2udyzval = 0.0;
    double d2udz2val = 0.0;

    std::map<int,Array<double>* > hess_vmap_new;
    std::map<int,Array<double>* > hess2_vmap_new;

    for(n2nit=n2n.begin();n2nit!=n2n.end();n2nit++)
    {
        gvid          = n2nit->first;
        sum           = 0.0;
        d2udx2val     = 0.0;
        d2udxyval     = 0.0;
        d2udxzval     = 0.0;
        d2udy2val     = 0.0;
        d2udyzval     = 0.0;
        d2udz2val     = 0.0;
        
        sum_dist   = 0.0;
        std::map<int,double>::iterator n2dit;

        for(n2dit=n2nit->second.begin();n2dit!=n2nit->second.end();n2dit++)
        {
            vid           = n2dit->first;
            di            = n2dit->second;
            sum_dist      = sum_dist+di;
            d2udx2val     = d2udx2val+hess_vmap[vid]->getVal(0,0)*di;
            d2udxyval     = d2udxyval+hess_vmap[vid]->getVal(1,0)*di;
            d2udxzval     = d2udxzval+hess_vmap[vid]->getVal(2,0)*di;
            d2udy2val     = d2udy2val+hess_vmap[vid]->getVal(3,0)*di;
            d2udyzval     = d2udyzval+hess_vmap[vid]->getVal(4,0)*di;
            d2udz2val     = d2udz2val+hess_vmap[vid]->getVal(5,0)*di;
        }
        
        Array<double>* hess_new = new Array<double>(6,1);
        hess_new->setVal(0,0,d2udx2val/sum_dist);
        hess_new->setVal(1,0,d2udxyval/sum_dist);
        hess_new->setVal(2,0,d2udxzval/sum_dist);
        hess_new->setVal(3,0,d2udy2val/sum_dist);
        hess_new->setVal(4,0,d2udyzval/sum_dist);
        hess_new->setVal(5,0,d2udz2val/sum_dist);
        
        hess_vmap_new[gvid] = hess_new;
        hess2_vmap_new[gvid] = hess_new;

        
    }
    
    hess_vmap.clear();
    
    // updates hess_vmap_new such that is contains the metric.
    ComputeMetric(P,metric_inputs, comm, hess_vmap_new);

    Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, hess_vmap_new, comm);
    
//    Mesh* us3dRoot = ReduceMeshToRoot(us3d->ien,
//                                      us3d->ief,
//                                      us3d->xcn,
//                                      us3d->ifn,
//                                      us3d->ife,
//                                      us3d->if_ref,
//                                      comm,info);
    
    
    //====================================================================
//
    std::vector<std::vector<int> > Elements = pDom->Elements;

    std::vector<std::vector<int> > tetras;
    std::vector<std::vector<int> > prisms;
    std::vector<std::vector<int> > hexes;

    for(int i=0;i<Elements.size();i++)
    {
        if(Elements[i].size() == 4)
        {
            std::vector<int> Et(4);
            for(int j=0;j<4;j++)
            {
                int nidt = Elements[i][j];
                Et[j] = nidt;

            }
            tetras.push_back(Et);
            Et.clear();
        }

        if(Elements[i].size() == 6)
        {
            std::vector<int> Ep(6);
            for(int j=0;j<6;j++)
            {
                int nidp = Elements[i][j];
                Ep[j]   = nidp;
            }
            prisms.push_back(Ep);
            Ep.clear();
        }

        if(Elements[i].size() == 8)
        {
            std::vector<int> Eh(8);
            for(int j=0;j<8;j++)
            {
                int nidh = Elements[i][j];
                Eh[j]   = nidh;
            }
            hexes.push_back(Eh);
            Eh.clear();
        }

    }

    std::cout << "Rank = " << world_rank << " N_hexes = " << hexes.size() << " N_prisms = " << prisms.size() << " N_tetras = " << tetras.size() << std::endl;
    std::ofstream myfilet;
    myfilet.open("output_" + std::to_string(world_rank) + ".dat");
    myfilet << "TITLE=\"new_volume.tec\"" << std::endl;
    myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"Up\", \"dUdx\", \"dUdy\", \"dUdz\", \"dUdx_v2\", \"dUdy_v2\", \"dUdz_v2\"" << std::endl;
    myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetras.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfilet << Verts[loc_vid]->x << " " << Verts[loc_vid]->y << " " << Verts[loc_vid]->z << " " << u_vmap[glob_vid]->getVal(0,0) << " " << u_vmap_new[glob_vid] << " " << dudx_vmap_new[glob_vid] << " " << dudy_vmap_new[glob_vid] << " " << dudz_vmap_new[glob_vid] << " " <<
            hess2_vmap_new[glob_vid]->getVal(0,0) << " " <<
            hess2_vmap_new[glob_vid]->getVal(0,1) << " " <<
            hess2_vmap_new[glob_vid]->getVal(0,2) <<std::endl;
    }

    for(int i=0;i<tetras.size();i++)
    {
        myfilet << tetras[i][0]+1 << " " << tetras[i][1]+1 << " "
                << tetras[i][2]+1 << " " << tetras[i][3]+1 <<  std::endl;
    }

    myfilet.close();
  // ====================================================================
    
    /*
    //std::map<int,Array<double>* > Hess_vm = P->ReduceMetricToVertices(Hess_map);
    //std::map<int,double> u_vmap           = P->ReduceFieldToVertices(Ui_map);
    std::map<int,Array<double>* > metric  = ComputeMetric(P,metric_inputs, hess_vmap_new, comm);

    Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, metric, comm);
    
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
    
    delete P;
    delete us3d;
    metric.clear();
    delete Ate;
    delete Apr;
    lv2gpv.clear();
    */
    
    if(world_rank == 0)
    {
        /*
        MMG5_pMesh mmgMesh_hyb = NULL;
        MMG5_pSol mmgSol_hyb   = NULL;
        
        MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
        MMG5_ARG_end);
        
        BoundaryMap* bmap = new BoundaryMap(us3dRoot->ifn, us3dRoot->if_ref);
        std::map<int,std::vector<int> > bnd_face_map = bmap->getBfaceMap();
        std::map<std::set<int>,int> tria_ref_map     = bmap->getTriaRefMap();
        std::map<std::set<int>,int> quad_ref_map     = bmap->getQuadRefMap();
        std::map<int,int> vert_ref_map               = bmap->getNodeRefMap();
        
        std::map<int,std::vector<int> > bndTriVol;
        std::map<int,int> bndTriVolRef;
        std::map<int,std::vector<int> > bndQuadVol;
        std::map<int,int> bndQuadVolRef;
        
        std::map<int,std::vector<int> >::iterator itb;
        
        for(itb=bnd_face_map.begin();itb!=bnd_face_map.end();itb++)
        {
            std::cout << " bnd_map = " << itb->first << " " << itb->second.size() << std::endl;
        }
        
        std::cout << "quad_ref_map.size() " << quad_ref_map.size() << std::endl;
        
        int tra = 0;
        int refer;
        int tt = 0;
        int qt = 0;
        int t3 = 0;
        for(int i=0;i<gAprRoot->getNrow();i++)
        {
            int gv0 = gAprRoot->getVal(i,0);
//            mmgMesh_hyb->prism[i+1].v[0] = gv0+1;
            int gv1 = gAprRoot->getVal(i,1);
//            mmgMesh_hyb->prism[i+1].v[1] = gv1+1;
            int gv2 = gAprRoot->getVal(i,2);
//            mmgMesh_hyb->prism[i+1].v[2] = gv2+1;
            int gv3 = gAprRoot->getVal(i,3);
//            mmgMesh_hyb->prism[i+1].v[3] = gv3+1;
            int gv4 = gAprRoot->getVal(i,4);
//            mmgMesh_hyb->prism[i+1].v[4] = gv4+1;
            int gv5 = gAprRoot->getVal(i,5);
//            mmgMesh_hyb->prism[i+1].v[5] = gv5+1;
//            mmgMesh_hyb->prism[i+1].ref = 3;
            
            std::set<int> tria0;
            tria0.insert(gv0);
            tria0.insert(gv1);
            tria0.insert(gv2);
            if(tria_ref_map.find(tria0)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria0];
                std::vector<int> tria(3);
                tria[0] = gv0+1;
                tria[1] = gv1+1;
                tria[2] = gv2+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                if(refer == 3)
                {
                    t3++;
                }
                tt++;
            }

            std::set<int> tria1;
            tria1.insert(gv3);
            tria1.insert(gv4);
            tria1.insert(gv5);
            if(tria_ref_map.find(tria1)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria1];
                std::vector<int> tria(3);
                tria[0] = gv3+1;
                tria[1] = gv4+1;
                tria[2] = gv5+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                if(refer == 3)
                {
                    t3++;
                }
                tt++;
            }
            
            std::set<int> quad0;
            quad0.insert(gv0);
            quad0.insert(gv1);
            quad0.insert(gv4);
            quad0.insert(gv3);
            if(quad_ref_map.find(quad0)!=quad_ref_map.end())
            {
                refer = quad_ref_map[quad0];
                std::vector<int> quad(4);
                quad[0] = gv0+1;
                quad[1] = gv1+1;
                quad[2] = gv4+1;
                quad[3] = gv3+1;
                bndQuadVol[qt] = quad;
                bndQuadVolRef[qt] = refer;
                qt++;
            }
            
            
            std::set<int> quad1;
            quad1.insert(gv1);
            quad1.insert(gv2);
            quad1.insert(gv5);
            quad1.insert(gv4);
            if(quad_ref_map.find(quad1)!=quad_ref_map.end())
            {
                refer = quad_ref_map[quad1];
                std::vector<int> quad(4);
                quad[0] = gv1+1;
                quad[1] = gv2+1;
                quad[2] = gv5+1;
                quad[3] = gv4+1;
                bndQuadVol[qt] = quad;
                bndQuadVolRef[qt] = refer;
                qt++;
            }
            
            
            std::set<int> quad2;
            quad2.insert(gv0);
            quad2.insert(gv3);
            quad2.insert(gv5);
            quad2.insert(gv2);
            if(quad_ref_map.find(quad2)!=quad_ref_map.end())
            {
                refer = quad_ref_map[quad2];
                std::vector<int> quad(4);
                quad[0] = gv0+1;
                quad[1] = gv3+1;
                quad[2] = gv5+1;
                quad[3] = gv2+1;
                bndQuadVol[qt] = quad;
                bndQuadVolRef[qt] = refer;
                qt++;
            }
            
            tria0.clear();
            tria1.clear();
            quad0.clear();
            quad1.clear();
            quad2.clear();

        }
        std::cout << "t3 - quadref " << t3 << " " << qt << " " << bndQuadVolRef.size()  << std::endl;
        for(int i=0;i<gAteRoot->getNrow();i++)
        {
            
            int gv0 = gAteRoot->getVal(i,0);
//            mmgMesh_hyb->tetra[i+1].v[0] = gv0+1;
            int gv1 = gAteRoot->getVal(i,1);
//            mmgMesh_hyb->tetra[i+1].v[1] = gv1+1;
            int gv2 = gAteRoot->getVal(i,2);
//            mmgMesh_hyb->tetra[i+1].v[2] = gv2+1;
            int gv3 = gAteRoot->getVal(i,3);
//            mmgMesh_hyb->tetra[i+1].v[3] = gv3+1;
//            mmgMesh_hyb->tetra[i+1].ref = 3;
//
            std::set<int> tria0;
            
            tria0.insert(gv0);
            tria0.insert(gv2);
            tria0.insert(gv1);
            
            if(tria_ref_map.find(tria0)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria0];
                std::vector<int> tria(3);
                tria[0] = gv0+1;
                tria[1] = gv2+1;
                tria[2] = gv1+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                
                if(refer == 3)
                {
                    t3++;
                }
                
                
                tt++;
            }
            std::set<int> tria1;
                
            tria1.insert(gv1);
            tria1.insert(gv2);
            tria1.insert(gv3);
            
            if(tria_ref_map.find(tria1)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria1];
                std::vector<int> tria(3);
                tria[0] = gv1+1;
                tria[1] = gv2+1;
                tria[2] = gv3+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                
                if(refer == 3)
                {
                    t3++;
                }
                tt++;
            }
            
            std::set<int> tria2;
            
            tria2.insert(gv0);
            tria2.insert(gv3);
            tria2.insert(gv2);
            
            if(tria_ref_map.find(tria2)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria2];
                std::vector<int> tria(3);
                tria[0] = gv0+1;
                tria[1] = gv3+1;
                tria[2] = gv2+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                
                if(refer == 3)
                {
                    t3++;
                }
                tt++;
            }
            
            std::set<int> tria3;
                                       
            tria3.insert(gv0);
            tria3.insert(gv1);
            tria3.insert(gv3);
            
            if(tria_ref_map.find(tria3)!=tria_ref_map.end())
            {
                refer = tria_ref_map[tria3];
                std::vector<int> tria(3);
                tria[0] = gv0+1;
                tria[1] = gv1+1;
                tria[2] = gv3+1;
                bndTriVol[tt]    = tria;
                bndTriVolRef[tt] = refer;
                
                if(refer == 3)
                {
                    t3++;
                }
                tt++;
            }
            
            tria0.clear();
            tria1.clear();
            tria2.clear();
            tria3.clear();
        }
        
        
        int nQuad     = bndQuadVolRef.size();
        int nTri      = bndTriVolRef.size();
        int nTets     = gAteRoot->getNrow();
        int nPrisms   = gAprRoot->getNrow();
        int nVerts    = us3dRoot->xcn->getNrow();
        
        
        if ( MMG3D_Set_meshSize(mmgMesh_hyb,nVerts,nTets,nPrisms,nTri,nQuad,0) != 1 )  exit(EXIT_FAILURE);
        
        if ( MMG3D_Set_solSize(mmgMesh_hyb,mmgSol_hyb,MMG5_Vertex,mmgMesh_hyb->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
        
        
        std::ofstream myfile_met;
        myfile_met.open("parttet_" + std::to_string(world_rank) + ".dat");
        myfile_met << "TITLE=\"new_volume.tec\"" << std::endl;
        myfile_met <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
        myfile_met <<"ZONE N = " << nVerts << ", E = " << nTets << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
        
        
        for(int i = 0;i < us3dRoot->xcn->getNrow();i++)
        {
            
            mmgMesh_hyb->point[i+1].c[0] = us3dRoot->xcn->getVal(i,0);
            mmgMesh_hyb->point[i+1].c[1] = us3dRoot->xcn->getVal(i,1);
            mmgMesh_hyb->point[i+1].c[2] = us3dRoot->xcn->getVal(i,2);
            
            
            mmgMesh_hyb->point[i+1].ref  = 1;
            
            double m11 = mv_g->getVal(i,0);
            double m12 = mv_g->getVal(i,1);
            double m13 = mv_g->getVal(i,2);
            double m22 = mv_g->getVal(i,3);
            double m23 = mv_g->getVal(i,4);
            double m33 = mv_g->getVal(i,5);
            
            myfile_met << us3dRoot->xcn->getVal(i,0) << " "
                       << us3dRoot->xcn->getVal(i,1) << " "
                       << us3dRoot->xcn->getVal(i,2) << " "
                       << m11 << " " << m12 << " " << m13 << " "
                       << m22 << " " << m23 << " " << m33 << std::endl;
            
            if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
            
        }
        
        
        
        for(int i=0;i<gAprRoot->getNrow();i++)
        {
            for(int j=0;j<gAprRoot->getNcol();j++)
            {
                int gv = gAprRoot->getVal(i,j);
                mmgMesh_hyb->prism[i+1].v[j] = gv+1;

            }
            
            mmgMesh_hyb->tetra[i+1].ref = -3;
        }
        
        for(int i=0;i<gAteRoot->getNrow();i++)
        {
            for(int j=0;j<gAteRoot->getNcol();j++)
            {
                int gv = gAteRoot->getVal(i,j);
                mmgMesh_hyb->tetra[i+1].v[j] = gv+1;
            }
            
            myfile_met << mmgMesh_hyb->tetra[i+1].v[0] << " " << mmgMesh_hyb->tetra[i+1].v[1] << " " << mmgMesh_hyb->tetra[i+1].v[2] << " " << mmgMesh_hyb->tetra[i+1].v[3] << std::endl;

            mmgMesh_hyb->tetra[i+1].ref = -3;
        }
        
        int st = 1;
        
        int t33 = 0;
        std::map<int,std::vector<int> >::iterator itbnd;
        for(itbnd=bndTriVol.begin();itbnd!=bndTriVol.end();itbnd++)
        {
            int itb = itbnd->first;
            
            mmgMesh_hyb->tria[st].v[0] = bndTriVol[itb][0];
            mmgMesh_hyb->tria[st].v[1] = bndTriVol[itb][1];
            mmgMesh_hyb->tria[st].v[2] = bndTriVol[itb][2];
            mmgMesh_hyb->tria[st].ref  = bndTriVolRef[itb];
            
            if(bndTriVolRef[itb] == 3)
            {
                t33++;
            }
            st++;
        }
        
//        for(int i=0;i<bndTriVolRef.size();i++)
//        {
//            mmgMesh_hyb->tria[i+1].v[0] = bndTriVol[i][0];
//            mmgMesh_hyb->tria[i+1].v[1] = bndTriVol[i][1];
//            mmgMesh_hyb->tria[i+1].v[2] = bndTriVol[i][2];
//            mmgMesh_hyb->tria[i+1].ref  = bndTriVolRef[i];
//        }
        
        int sq = 1;
        
        for(itbnd=bndQuadVol.begin();itbnd!=bndQuadVol.end();itbnd++)
        {
            int itb = itbnd->first;

            mmgMesh_hyb->quadra[sq].v[0] = bndQuadVol[itb][0];
            mmgMesh_hyb->quadra[sq].v[1] = bndQuadVol[itb][1];
            mmgMesh_hyb->quadra[sq].v[2] = bndQuadVol[itb][2];
            mmgMesh_hyb->quadra[sq].v[3] = bndQuadVol[itb][3];
            mmgMesh_hyb->quadra[sq].ref  = bndQuadVolRef[itb];
            sq++;
        }
        
        
        std::map<int,int>::iterator itm;
        
        
        
        myfile_met.close();
        
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
        
        //ReferenceMesh* refmesh = ReadReferenceMesh();
        
//        std::ofstream myfile_nv;
//        myfile_nv.open("Nodes.ref");
        for(int i=0;i<mmgMesh_hyb->np;i++)
        {
            mmgMesh_TETCOPY->point[i+1].c[0] = mmgMesh_hyb->point[i+1].c[0];
            mmgMesh_TETCOPY->point[i+1].c[1] = mmgMesh_hyb->point[i+1].c[1];
            mmgMesh_TETCOPY->point[i+1].c[2] = mmgMesh_hyb->point[i+1].c[2];
            mmgMesh_TETCOPY->point[i+1].ref  = 0;
            
//            myfile_nv << mmgMesh_hyb->point[i+1].c[0] << " " << mmgMesh_hyb->point[i+1].c[1] << " " << mmgMesh_hyb->point[i+1].c[2] << std::endl;
            
//            if(fabs(mmgMesh_hyb->point[i+1].c[0]-refmesh->Nodes[i][0])>tol ||
//               fabs(mmgMesh_hyb->point[i+1].c[1]-refmesh->Nodes[i][1])>tol ||
//               fabs(mmgMesh_hyb->point[i+1].c[2]-refmesh->Nodes[i][2])>tol)
//            {
//                flip = 1;
//            }
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
            
//          myfile_ne << mmgMesh_hyb->tetra[i+1].v[0] << " "
//                        << mmgMesh_hyb->tetra[i+1].v[1] << " "
//                        << mmgMesh_hyb->tetra[i+1].v[2] << " "
//                        << mmgMesh_hyb->tetra[i+1].v[3] << std::endl;
            
//            if( fabs(mmgMesh_hyb->tetra[i+1].v[0]-refmesh->Elements[i][0])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[1]-refmesh->Elements[i][1])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[2]-refmesh->Elements[i][2])>tol || fabs(mmgMesh_hyb->tetra[i+1].v[3]-refmesh->Elements[i][3])>tol)
//            {
//                flip = 1;
//            }
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
        
        int t333=0;
        for(int ntt=0;ntt<mmgMesh_hyb->nt;ntt++)
        {
            
            if(mmgMesh_hyb->tria[ntt+1].ref == 3)
            {
                t333++;
            }
            st++;
        }
        
        std::cout << "t333 - after " << t333 << std::endl;
        
       // WriteUS3DGridFromMMG_itN(mmgMesh_hyb, us3d, bnd_face_map);
        
//      std::cout << "tets = " << mmgMesh_TETCOPY->ne << " " << nTets << " " << nPrisms << std::endl;
        
        std::cout<<"Started writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
        //OutputMesh_MMG(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
        std::cout<<"Finished writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
        */
    
    }
    /**/
    

    MPI_Finalize();
    
}

