#include "adapt_parops.h"

Array<double>* GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, std::map<int,Array<double>*> mv_map, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    MMG_Mesh* mmg = new MMG_Mesh;
    
    std::map<int,Array<double>* >::iterator grit;
    Array<int>* lE2gE = new Array<int>(mv_map.size(),1);
    Array<double>* mv = new Array<double>(mv_map.size(),6);
    i = 0;
    for(grit=mv_map.begin();grit!=mv_map.end();grit++)
    {
        lE2gE->setVal(i,0,grit->first);
        mv->setVal(i,0,grit->second->getVal(0,0));
        mv->setVal(i,1,grit->second->getVal(0,1));
        mv->setVal(i,2,grit->second->getVal(0,2));
        mv->setVal(i,3,grit->second->getVal(1,1));
        mv->setVal(i,4,grit->second->getVal(1,2));
        mv->setVal(i,5,grit->second->getVal(2,2));
        i++;
    }
//
    
    int* lid_nlocs      = new int[world_size];
    int* red_lid_nlocs  = new int[world_size];
    int* lid_offsets    = new int[world_size];

    int* mv_nlocs      = new int[world_size];
    int* red_mv_nlocs  = new int[world_size];
    int* mv_offsets    = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        lid_nlocs[i] = 0;
        mv_nlocs[i] = 0;
        if(i==world_rank)
        {
            lid_nlocs[i] = mv->getNrow();
            mv_nlocs[i]  = mv->getNrow()*6;
        }
        else
        {
            lid_nlocs[i] = 0;
            mv_nlocs[i]  = 0;
        }
    }

    MPI_Allreduce(mv_nlocs, red_mv_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(lid_nlocs, red_lid_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    int offset_lid = 0;
    int offset_mv = 0;
    for(i=0;i<world_size;i++)
    {
        lid_offsets[i] = offset_lid;
        offset_lid = offset_lid+red_lid_nlocs[i];
        
        mv_offsets[i] = offset_mv;
        offset_mv = offset_mv+red_mv_nlocs[i];
    }

    Array<int>*  lE2gE_g;
    Array<double>*  mv_g;

    int n_glob_lid = offset_lid;

    if(world_rank == 0)
    {
        lE2gE_g   = new Array<int>(n_glob_lid,1);
        mv_g      = new Array<double>(n_glob_lid,6);
    }
    else
    {
        lE2gE_g   = new Array<int>(1,1);
        mv_g      = new Array<double>(1,1);
    }
    int ncol_lid = 1;
    MPI_Gatherv(&lE2gE->data[0],
                lE2gE->getNrow()*ncol_lid,
                MPI_INT,
                &lE2gE_g->data[0],
                red_lid_nlocs,
                lid_offsets,
                MPI_INT, 0, comm);

    int ncol_mv = 6;
    MPI_Gatherv(&mv->data[0],
                mv->getNrow()*ncol_mv,
                MPI_DOUBLE,
                &mv_g->data[0],
                red_mv_nlocs,
                mv_offsets,
                MPI_DOUBLE, 0, comm);
    
    
    
    
    Array<double>* xcn_g;
    Array<int>* ien_g;
    int nvg   = us3d->xcn->getNglob();
    int nElem = us3d->ien->getNglob();
    ParallelState* xcn_pstate = P->getXcnParallelState();
    ParallelState* ien_pstate = P->getIenParallelState();
    if(world_rank == 0)
    {
        xcn_g = new Array<double>(nvg,3);
        ien_g = new Array<int>(nElem,8);
    }
    else
    {
        xcn_g = new Array<double>(1,1);
        ien_g = new Array<int>(1,1);
    }
    int* ien_nlocs      = new int[world_size];
    int* ien_offsets    = new int[world_size];
    int* xcn_nlocs      = new int[world_size];
    int* xcn_offsets    = new int[world_size];

    for(int i=0;i<world_size;i++)
    {
        xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]*3;
        xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;

        ien_nlocs[i]   = ien_pstate->getNlocs()[i]*8;
        ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
    }

    MPI_Gatherv(&us3d->xcn->data[0],
                us3d->xcn->getNrow()*3,
                MPI_DOUBLE,
                &xcn_g->data[0],
                xcn_nlocs,
                xcn_offsets,
                MPI_DOUBLE, 0, comm);


    MPI_Gatherv(&us3d->ien->data[0],
                us3d->ien->getNrow()*8,
                MPI_INT,
                &ien_g->data[0],
                ien_nlocs,
                ien_offsets,
                MPI_INT, 0, comm);
    
    
    
    Array<double>* Mg;
    Array<int>* tel;
    if(world_rank == 0)
    {
        Mg  = new Array<double>(us3d->xcn->getNglob(),6);
        for(int u=0;u<n_glob_lid;u++)
        {
            int gvid     = lE2gE_g->getVal(u,0);
            Mg->setVal(gvid,0,mv_g->getVal(u,0));
            Mg->setVal(gvid,1,mv_g->getVal(u,1));
            Mg->setVal(gvid,2,mv_g->getVal(u,2));
            Mg->setVal(gvid,3,mv_g->getVal(u,3));
            Mg->setVal(gvid,4,mv_g->getVal(u,4));
            Mg->setVal(gvid,5,mv_g->getVal(u,5));
        }
        
        
        std::string filename = "MetricRoot.dat";
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
        myfile <<"ZONE N = " << us3d->xcn->getNglob() << ", E = " << us3d->ien->getNglob() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
        
        string filename2 = "metric.dat";
        ofstream myfile2;
        myfile2.open(filename2);

        string filename3 = "elements.dat";
        ofstream myfile3;
        myfile3.open(filename3);

        for(int i=0;i<us3d->xcn->getNglob();i++)
        {
            myfile <<     xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
                    << " " << Mg->getVal(i,0) << " " <<    Mg->getVal(i,1) << " " << Mg->getVal(i,2)
                    << " " << Mg->getVal(i,3) << " " <<    Mg->getVal(i,4) << " " << Mg->getVal(i,5) << std::endl;
            
            myfile2 <<std::setprecision(16)<< Mg->getVal(i,0) << " " <<    Mg->getVal(i,1) << " " << Mg->getVal(i,2)
                    << " " << Mg->getVal(i,3) << " " <<    Mg->getVal(i,4) << " " << Mg->getVal(i,5) << std::endl;
        }
        for(int i=0;i<ien_g->getNrow();i++)
        {
            myfile <<  ien_g->getVal(i,0)+1 << " " <<
                       ien_g->getVal(i,1)+1 << " " <<
                       ien_g->getVal(i,2)+1 << " " <<
                       ien_g->getVal(i,3)+1 << " " <<
                       ien_g->getVal(i,4)+1 << " " <<
                       ien_g->getVal(i,5)+1 << " " <<
                       ien_g->getVal(i,6)+1 << " " <<
                       ien_g->getVal(i,7)+1 << std::endl;
            
            myfile3 << ien_g->getVal(i,0)+1 << " " <<
                       ien_g->getVal(i,1)+1 << " " <<
                       ien_g->getVal(i,2)+1 << " " <<
                       ien_g->getVal(i,3)+1 << " " <<
                       ien_g->getVal(i,4)+1 << " " <<
                       ien_g->getVal(i,5)+1 << " " <<
                       ien_g->getVal(i,6)+1 << " " <<
                       ien_g->getVal(i,7)+1 << std::endl;
        }
        myfile.close();
        myfile2.close();
        myfile3.close();
    }
    else
    {
        Mg = new Array<double>(1,1);
    }
    
    return Mg;
}
