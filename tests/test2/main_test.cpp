#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include <iomanip>

std::vector<double> ReadReferenceData()
{
    std::ifstream fin;
    fin.open("GuX.ref");
    
    // Read the file row by row
    double val;
    std::vector<double> Gref;
    int t=0;
    while(fin >> val)
    {
       Gref.push_back(val);
    }
    
    return Gref;
}


int main(int argc, char** argv) {
    
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    const char* fn_grid="../test_mesh/grid.h5";
    const char* fn_conn="../test_mesh/conn.h5";
    const char* fn_data="../test_mesh/data.h5";
    
    US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);

    int Nel_part = us3d->ien->getNrow();

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    for(int i=0;i<Nel_part;i++)
    {
        Ui->setVal(i,0,us3d->interior->getVal(i,varia));
    }
    
    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
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
    std::vector<int> LocElem    = P->getLocElem();
    std::vector<double> Uvaria  = P->getLocElemVaria();

    std::map<int,double> Ui_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid     = LocElem[i];
        UvariaV     = Uvaria[i];
        Ui_map[gid] = UvariaV;

    }

    Domain* pDom = P->getPartitionDomain();

    std::map<int,double> u_vmap = P->ReduceFieldToVertices(Ui_map);

    std::vector<Vert> Verts  = P->getLocalVerts();

    std::string filename = "output_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\"" << std::endl;
    myfile <<"ZONE N = " << u_vmap.size() << ", E = " << LocElem.size() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;

    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::map<int,int> gv2lpartv     = pDom->gv2lpartv;
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    std::map<int,int> gv2lpv        = pDom->gv2lpv;
    i = 0;
    std::map<int,double>::iterator itm;
    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfile << Verts[loc_vid].x << " " << Verts[loc_vid].y << " " << Verts[loc_vid].z << " " << u_vmap[glob_vid] << std::endl;
    }
    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
    for(int i=0;i<LocElem.size();i++)
    {
        int glob_id = LocElem[i];

        myfile <<   pDom->LocElem2LocNode->getVal(i,0)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,1)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,2)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,3)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,4)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,5)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,6)+1 << "  " <<
        pDom->LocElem2LocNode->getVal(i,7)+1 << std::endl;
    }

    myfile.close();


    
    
    
    
    MPI_Finalize();
    
}
