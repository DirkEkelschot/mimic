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
                                 us3d->ifn, us3d->ife, us3d->if_ref,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, comm);
    
    double duration = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Ptime = 0.0;
    MPI_Allreduce(&duration, &Ptime, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
    std::cout << "Timing partitioning: " << duration << std::endl;
    }
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

    std::map<int,double> u_vmap = P->ReduceFieldToVertices(Ui_map);
    Domain* pDom = P->getPartitionDomain();

    std::vector<Vert> Verts = P->getLocalVerts();
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
    
    //Array<double>* xcn_g = ReadDataSetFromFile<double>(fn_grid,"xcn");

    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    
    std::ofstream myfilet;
    myfilet.open("parttet_" + std::to_string(world_rank) + ".dat");
    myfilet << "TITLE=\"new_volume.tec\"" << std::endl;
    myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\"" << std::endl;
    myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetras.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
    

    
    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid  = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfilet << Verts[loc_vid].x << " " << Verts[loc_vid].y << " " << Verts[loc_vid].z << " " << u_vmap[glob_vid] <<  std::endl;
    }
    
//        for(int i=0;i<xcn_g->getNrow();i++)
//        {
//            myfilet << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2) << std::endl;
//        }
    
    for(int i=0;i<tetras.size();i++)
    {
        myfilet << tetras[i][0]+1 << " " << tetras[i][1]+1 << " " << tetras[i][2]+1 << " " << tetras[i][3]+1 << std::endl;
    }
    
//    for(int i=0;i<hexes.size();i++)
//    {
//        myfilet << hexes[i][0]+1 << " " << hexes[i][1]+1 << " " << hexes[i][2]+1 << " " << hexes[i][3]+1 << " " << hexes[i][4]+1 << " " << hexes[i][5]+1 << " " << hexes[i][6]+1 << " " << hexes[i][7]+1 <<  std::endl;
//    }
    
    myfilet.close();
    
//    for(int i=0;i<xcn_g->getNrow();i++)
//    {
//        myfilet << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2) << std::endl;
//    }
//
//    for(int i=0;i<tetras.size();i++)
//    {
//        int g0 = lpartv2gv[tetras[i][0]];
//        int g1 = lpartv2gv[tetras[i][1]];
//        int g2 = lpartv2gv[tetras[i][2]];
//        int g3 = lpartv2gv[tetras[i][3]];
//
//        myfilet << g0+1 << " " << g1+1 << " " << g2+1 << " " << g3+1 << std::endl;
//    }
    
//    if(world_rank == 0)
//    {
//        for(int i=0;i<LocElem.size();i++)
//        {
//            int glob_id = LocElem[i];
//
//            std::cout <<   pDom->LocElem2LocNode->getVal(i,0)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,1)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,2)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,3)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,4)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,5)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,6)+1 << "  " <<
//            pDom->LocElem2LocNode->getVal(i,7)+1 << std::endl;
//        }
//    }
    
    
    std::cout << world_rank << " # prisms = " << prisms.size() << "  " << " # tetras = " << tetras.size() << std::endl;
    
    MPI_Finalize();
}
