#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"

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
    
    const char* fn_grid="test_mesh/grid.h5";
    const char* fn_conn="test_mesh/conn.h5";
    const char* fn_data="test_mesh/data.h5";
    
    US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);

    int Nel_part = us3d->ien->getNrow();

    Array<double>* Uivar = new Array<double>(Nel_part,1);
    int varia = 4;
    for(int i=0;i<Nel_part;i++)
    {
        Uivar->setVal(i,0,us3d->interior->getVal(i,varia));
    }
    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief,
                                 us3d->ifn, us3d->ife, us3d->if_ref,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Uivar, comm);
    
    
    MPI_Finalize();
    
}
