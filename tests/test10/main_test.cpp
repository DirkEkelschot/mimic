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
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid     = LocElem[i];
        UvariaV     = Uvaria[i];
        Ui_map[gid] = UvariaV;
    }

    std::vector<Vert> Verts = P->getLocalVerts();
    
    std::map<int,double> Uaux = P->CommunicateStateAdjacentElements(Ui_map,comm);

    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }

    t = clock();
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
    
    std::map<int,double > u_vm = P->ReduceFieldToAllVertices(Uaux);
    
    std::map<int,double >::iterator its;

//    for(its=u_vm.begin();its!=u_vm.end();its++)
//    {
//        std::cout << its->first << " > " << its->second << std::endl;
//    }
    
    int size_b = u_vm.size();
    
    P->AddAdjacentVertexDataUS3D(u_vm,comm);
    
    std::map<int,std::vector<int> > scheme_E2V = meshTopo->getScheme_E2V();

    std::map<int,std::vector<int> >::iterator itsch;
    
//    for(itsch=scheme_E2V.begin();itsch!=scheme_E2V.end();itsch++)
//    {
//        std::cout << itsch->first << " N= " << itsch->second.size() << " :: " ;
//        for(int q=0;q<itsch->second.size();q++)
//        {
//            int vid = itsch->second[q];
//            std::cout << u_vm[vid] << " ";
//        }
//        std::cout << std::endl;
//    }
    
    
    
    
    std::map<int,std::vector<double> >::iterator itmm;
    
    Domain* pDom = P->getPartitionDomain();
    
    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_Vrt_US3D(P,Uaux,u_vm,meshTopo,gB,comm);
    
//    std::map<int,Array<double>* > dUdXi = ComputedUdx_MGG(P,Uaux,meshTopo,gB,comm);
    
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
    
    std::map<int,Array<double>* >::iterator grit;
    std::map<int,double> dUidxi_map;
    std::map<int,double> dUidyi_map;
    std::map<int,double> dUidzi_map;

    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        lE2gE->setVal(i,0,grit->first);
        dUidxi->setVal(i,0,grit->second->getVal(0,0));
        dUidyi->setVal(i,0,grit->second->getVal(1,0));
        dUidzi->setVal(i,0,grit->second->getVal(2,0));
        dUidxi_map[grit->first]=grit->second->getVal(0,0);
        dUidyi_map[grit->first]=grit->second->getVal(1,0);
        dUidzi_map[grit->first]=grit->second->getVal(2,0);

        //std::cout << grit->second->getVal(0,0) << std::endl;
        i++;
    }
    
    

    //==================================================================================
//    Domain* pDom = P->getPartitionDomain();
//    std::vector<int> loc_part_verts = pDom->loc_part_verts;
//    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
//    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
//    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    std::map<int,double > dUidxiAux_map  = P->CommunicateStateAdjacentElements(dUidxi_map,comm);
    std::map<int,double > dUidyiAux_map  = P->CommunicateStateAdjacentElements(dUidyi_map,comm);
    std::map<int,double > dUidziAux_map  = P->CommunicateStateAdjacentElements(dUidzi_map,comm);

    std::map<int,double> dudx_vmap = P->ReduceFieldToAllVertices(dUidxiAux_map);
    std::map<int,double> dudy_vmap = P->ReduceFieldToAllVertices(dUidyiAux_map);
    std::map<int,double> dudz_vmap = P->ReduceFieldToAllVertices(dUidziAux_map);

    std::map<int,std::map<int,double> > n2n = P->getNode2NodeMap();
    std::map<int,std::map<int,double> >::iterator n2nit;
    double sum      = 0.0;
    double sum_dist = 0.0;
    double di = 0.0;
    double valx = 0.0;
    double valy = 0.0;
    double valz = 0.0;
    std::map<int,double> dudx_vmap_new;
    std::map<int,double> dudy_vmap_new;
    std::map<int,double> dudz_vmap_new;
    int vid,gvid;
    for(n2nit=n2n.begin();n2nit!=n2n.end();n2nit++)
    {
        gvid     = n2nit->first;
        sum      = 0.0;
        valx     = 0.0;
        valy     = 0.0;
        valz     = 0.0;
        sum_dist = 0.0;
        std::map<int,double>::iterator n2dit;
        for(n2dit=n2nit->second.begin();n2dit!=n2nit->second.end();n2dit++)
        {
            vid      = n2dit->first;
            di       = n2dit->second;
            sum_dist = sum_dist+di;
            valx     = valx+dudx_vmap[vid]*di;
            valy     = valy+dudy_vmap[vid]*di;
            valz     = valz+dudz_vmap[vid]*di;
        }
        
        //std::cout << n2nit->second.size() << std::endl;
        dudx_vmap_new[gvid] = valx/sum_dist;
        dudy_vmap_new[gvid] = valy/sum_dist;
        dudz_vmap_new[gvid] = valz/sum_dist;
    }
    

    std::vector<std::vector<int> > tetras;
    std::vector<std::vector<int> > prisms;
    std::vector<std::vector<int> > hexes;

    std::vector<std::vector<int> > Elements = pDom->Elements;

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
    myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dUdx\", \"dUdy\", \"dUdz\"" << std::endl;
    myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetras.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfilet << Verts[loc_vid].x << " " << Verts[loc_vid].y << " " << Verts[loc_vid].z << " " << dudx_vmap_new[glob_vid] << " " << dudy_vmap_new[glob_vid] << " " << dudz_vmap_new[glob_vid] << std::endl;
    }

    for(int i=0;i<tetras.size();i++)
    {
        myfilet << tetras[i][0]+1 << " " << tetras[i][1]+1 << " "
                << tetras[i][2]+1 << " " << tetras[i][3]+1 <<  std::endl;
    }

    myfilet.close();
    //==================================================================================

    /*
    
    
     
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
    */
    MPI_Finalize();
    
}

