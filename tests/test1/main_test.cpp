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
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
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
    
    std::map<int,double> Uadj = P->CommunicateAdjacentDataUS3D(Ui_map,comm);
    int* bnd_map;
    int nBnd = 4;

    if(world_rank == 0)
    {
        std::cout << "Started creating mesh topology object... " << std::endl;
    }

    //Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
    
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
   
    if(world_rank == 0)
    {
        std::cout << "Started reconstructing the gradient... " << std::endl;
    }
    
    t = clock();
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,Uadj,gB,comm);
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

    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        lE2gE->setVal(i,0,grit->first);
        dUidxi->setVal(i,0,grit->second->getVal(0,0));
        dUidyi->setVal(i,0,grit->second->getVal(1,0));
        dUidzi->setVal(i,0,grit->second->getVal(2,0));
        i++;
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
    
    MPI_Finalize();
    
}
