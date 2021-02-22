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


Mesh* ReduceMeshToRoot(ParArray<int>* ien,
                       ParArray<int>* ief,
                       ParArray<double>* xcn,
                       ParArray<int>* ifn,
                       ParArray<int>* ife,
                       ParArray<int>* if_ref,
                       MPI_Comm comm, MPI_Info info)
{
    
    Mesh* us3d_root = new Mesh;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    Array<double>*  xcn_g;
    Array<int>*     ief_g;
    Array<int>*     ien_g;
    Array<int>*     ifn_g;
    Array<int>*     if_ref_g;
    Array<int>*     ife_g;
    
    if(world_rank == 0)
    {
        xcn_g       = new Array<double>(xcn->getNglob(),3);
        ief_g       = new Array<int>(ief->getNglob(),6);
        ien_g       = new Array<int>(ien->getNglob(),8);
        if_ref_g    = new Array<int>(ifn->getNglob(),1);
        ifn_g       = new Array<int>(ifn->getNglob(),4);
        ife_g       = new Array<int>(ifn->getNglob(),2);
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
    int* red_ien_nlocs  = new int[world_size];
    int* ien_offsets    = new int[world_size];
    
    int* ief_nlocs      = new int[world_size];
    int* ief_offsets    = new int[world_size];
    int* red_ief_nlocs  = new int[world_size];

    int* xcn_nlocs      = new int[world_size];
    int* xcn_offsets    = new int[world_size];
    int* red_xcn_nlocs  = new int[world_size];

    int* ifn_nlocs      = new int[world_size];
    int* ifn_offsets    = new int[world_size];
    int* red_ifn_nlocs  = new int[world_size];

    int* if_ref_nlocs   = new int[world_size];
    int* if_ref_offsets = new int[world_size];
    int* red_if_ref_nlocs  = new int[world_size];

    int* ife_nlocs       = new int[world_size];
    int* ife_offsets     = new int[world_size];
    int* red_ife_nlocs   = new int[world_size];
    
    for(int i=0;i<world_size;i++)
    {
        ien_nlocs[i] = 0;
        ief_nlocs[i] = 0;
        xcn_nlocs[i] = 0;
        
        if(i==world_rank)
        {
            ien_nlocs[i] = ien->getNrow()*8;
            ief_nlocs[i] = ief->getNrow()*6;
            xcn_nlocs[i] = xcn->getNrow()*3;
            ifn_nlocs[i] = ifn->getNrow()*4;
            ife_nlocs[i] = ife->getNrow()*2;
            if_ref_nlocs[i] = ifn->getNrow()*1;
        }
        else
        {
            ien_nlocs[i]    = 0;
            ief_nlocs[i]    = 0;
            xcn_nlocs[i]    = 0;
            ifn_nlocs[i]    = 0;
            ife_nlocs[i]    = 0;
            if_ref_nlocs[i] = 0;
        }
    }
    
    MPI_Allreduce(ien_nlocs, red_ien_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ief_nlocs, red_ief_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(xcn_nlocs, red_xcn_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ifn_nlocs, red_ifn_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ife_nlocs, red_ife_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(if_ref_nlocs, red_if_ref_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    
    int ien_o=0;
    int ief_o=0;
    int xcn_o=0;
    int ifn_o=0;
    int ife_o=0;
    int if_ref_o=0;
    
    for(int i=0;i<world_size;i++)
    {
        ien_offsets[i] = ien_o;
        ief_offsets[i] = ief_o;
        xcn_offsets[i] = xcn_o;
        ifn_offsets[i] = ifn_o;
        ife_offsets[i] = ife_o;
        if_ref_offsets[i] = if_ref_o;
        
        ien_o = ien_o + red_ien_nlocs[i];
        ief_o = ief_o + red_ief_nlocs[i];
        xcn_o = xcn_o + red_xcn_nlocs[i];
        ifn_o = ifn_o + red_ifn_nlocs[i];
        ife_o = ife_o + red_ife_nlocs[i];
        if_ref_o = if_ref_o + red_if_ref_nlocs[i];
    }

    MPI_Gatherv(&xcn->data[0],
                xcn->getNrow()*3,
                MPI_DOUBLE,
                &xcn_g->data[0],
                red_xcn_nlocs,
                xcn_offsets,
                MPI_DOUBLE, 0, comm);

    MPI_Gatherv(&ief->data[0],
                ief->getNrow()*6,
                MPI_INT,
                &ief_g->data[0],
                red_ief_nlocs,
                ief_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&ien->data[0],
                ien->getNrow()*8,
                MPI_INT,
                &ien_g->data[0],
                red_ien_nlocs,
                ien_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&ifn->data[0],
                ifn->getNrow()*4,
                MPI_INT,
                &ifn_g->data[0],
                red_ifn_nlocs,
                ifn_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&if_ref->data[0],
                if_ref->getNrow()*1,
                MPI_INT,
                &if_ref_g->data[0],
                red_if_ref_nlocs,
                if_ref_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&ife->data[0],
                ife->getNrow()*2,
                MPI_INT,
                &ife_g->data[0],
                red_ife_nlocs,
                ife_offsets,
                MPI_INT, 0, comm);
    
    us3d_root->xcn      = xcn_g;
    us3d_root->ief      = ief_g;
    us3d_root->ien      = ien_g;
    us3d_root->if_ref   = if_ref_g;
    us3d_root->ifn      = ifn_g;
    us3d_root->ife      = ife_g;
    
    return us3d_root;
}

