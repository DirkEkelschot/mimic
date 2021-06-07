#include "adapt_parops.h"
#include "adapt_operations.h"

Array<double>* GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, std::map<int,Array<double>*> mv_map, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    //MMG_Mesh* mmg = new MMG_Mesh;
    
    Domain* pDom = P->getPartitionDomain();
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;

    std::map<int,Array<double>* >::iterator grit;
    Array<int>* lE2gE = new Array<int>(mv_map.size(),1);
    Array<double>* mv = new Array<double>(mv_map.size(),6);
    i = 0;
    for(grit=mv_map.begin();grit!=mv_map.end();grit++)
    {
        //std::cout << "(rowxcol) =(" << grit->second->getNrow() << " " << grit->second->getNcol() << ") -> " << grit->second->getVal(0,0) << " " << grit->second->getVal(0,1) << " " << grit->second->getVal(0,2) << " " << grit->second->getVal(1,0) << " " << grit->second->getVal(1,1) << " " << grit->second->getVal(1,2) << " " << grit->second->getVal(2,0)  << " " << grit->second->getVal(2,1) << " " << grit->second->getVal(2,2) << std::endl;
        
        lE2gE->setVal(i,0,grit->first);
        mv->setVal(i,0,grit->second->getVal(0,0));
        mv->setVal(i,1,grit->second->getVal(0,1));
        mv->setVal(i,2,grit->second->getVal(0,2));
        mv->setVal(i,3,grit->second->getVal(1,1));
        mv->setVal(i,4,grit->second->getVal(1,2));
        mv->setVal(i,5,grit->second->getVal(2,2));
        i++;
    }
    
    int* lid_nlocs      = new int[world_size];
    int* red_lid_nlocs  = new int[world_size];
    int* lid_offsets    = new int[world_size];

    int* mv_nlocs      = new int[world_size];
    int* red_mv_nlocs  = new int[world_size];
    int* mv_offsets    = new int[world_size];


    
    for(i=0;i<world_size;i++)
    {
        lid_nlocs[i] = 0;
        mv_nlocs[i]  = 0;
        
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
    int offset_mv  = 0;
    
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
    Array<int>* iet_g;


    
    int nvg   = us3d->xcn->getNglob();
    int nElem = us3d->ien->getNglob();
    ParallelState* xcn_pstate = P->getXcnParallelState();
    ParallelState* ien_pstate = P->getIenParallelState();
    if(world_rank == 0)
    {
        xcn_g = new Array<double>(nvg,3);
        ien_g = new Array<int>(nElem,8);
        iet_g = new Array<int>(nElem,1);

    }
    else
    {
        xcn_g = new Array<double>(1,1);
        ien_g = new Array<int>(1,1);
        iet_g = new Array<int>(1,1);
    }
    
    int* iet_nlocs      = new int[world_size];
    int* iet_offsets    = new int[world_size];
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
        
        iet_nlocs[i]   = ien_pstate->getNlocs()[i]*1;
        iet_offsets[i] = ien_pstate->getOffsets()[i]*1;
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
    
    MPI_Gatherv(&us3d->iet->data[0],
                us3d->iet->getNrow()*1,
                MPI_INT,
                &iet_g->data[0],
                iet_nlocs,
                iet_offsets,
                MPI_INT, 0, comm);
    
    int NtetLoc = us3d->elTypes->getVal(0,0);
//    int NprismLoc = us3d->elTypes->getVal(1,0);
//    int NhexLoc = us3d->elTypes->getVal(2,0);
//
//
    int Ntetras = 0;
    int Nprisms = 0;
    int Nhexes = 0;

    
    MPI_Reduce(&NtetLoc,   &Ntetras,   1, MPI_INT, MPI_SUM, 0, comm);
//    MPI_Reduce(&NprismLoc, &Nprisms,   1, MPI_INT, MPI_SUM, 0,  comm);
//    MPI_Reduce(&NhexLoc,   &Nhexes,    1, MPI_INT, MPI_SUM, 0,  comm);
    //std::cout << "RANK = " << NtetLoc << " " << Ntetras << std::endl;

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
        
        
//        std::string filename = "MetricRoot.dat";
//        std::ofstream myfile;
//        myfile.open(filename);
//        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
//        myfile <<"ZONE N = " << us3d->xcn->getNglob() << ", E = " << Ntetras << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
//        
//        string filename2 = "metric.dat";
//        ofstream myfile2;
//        myfile2.open(filename2);
//
//        string filename3 = "elements.dat";
//        ofstream myfile3;
//        myfile3.open(filename3);
//        int tel = 0;
//        for(int i=0;i<us3d->xcn->getNglob();i++)
//        {
//            myfile <<     xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
//                    << " " << Mg->getVal(i,0) << " " <<    Mg->getVal(i,1) << " " << Mg->getVal(i,2)
//                    << " " << Mg->getVal(i,3) << " " <<    Mg->getVal(i,4) << " " << Mg->getVal(i,5) << std::endl;
//            
//            myfile2 <<std::setprecision(16)<< Mg->getVal(i,0) << " " <<    Mg->getVal(i,1) << " " << Mg->getVal(i,2)
//                    << " " << Mg->getVal(i,3) << " " <<    Mg->getVal(i,4) << " " << Mg->getVal(i,5) << std::endl;
//        }
//        for(int i=0;i<ien_g->getNrow();i++)
//        {
//            if(iet_g->getVal(i,0)==2)
//            {
//                myfile <<  ien_g->getVal(i,0)+1 << " " <<
//                           ien_g->getVal(i,1)+1 << " " <<
//                           ien_g->getVal(i,2)+1 << " " <<
//                           ien_g->getVal(i,3)+1 << " " << std::endl;
//                tel++;
//            }
//            
//            
//            myfile3 << ien_g->getVal(i,0)+1 << " " <<
//                       ien_g->getVal(i,1)+1 << " " <<
//                       ien_g->getVal(i,2)+1 << " " <<
//                       ien_g->getVal(i,3)+1 << " " <<
//                       ien_g->getVal(i,4)+1 << " " <<
//                       ien_g->getVal(i,5)+1 << " " <<
//                       ien_g->getVal(i,6)+1 << " " <<
//                       ien_g->getVal(i,7)+1 << std::endl;
//        }
//        myfile.close();
//        myfile2.close();
//        myfile3.close();
        
    }
    else
    {
        Mg = new Array<double>(1,1);
    }

    delete xcn_g;
    delete ien_g;
    delete iet_g;
    //delete pDom;
    delete lE2gE_g;
    delete mv_g;
    delete lE2gE;
    delete mv;
    
    delete[] iet_nlocs;
    delete[] iet_offsets;
    delete[] ien_nlocs;
    delete[] ien_offsets;
    delete[] xcn_nlocs;
    delete[] xcn_offsets;
    delete[] mv_nlocs;
    delete[] red_mv_nlocs;
    delete[] mv_offsets;
    
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





MeshTransfer* GetMeshTransfer(std::map<int,std::vector<int> > Ate, std::vector<int> ufaces_vec, std::vector<int> uverts_vec, MPI_Comm comm, MPI_Info info)
{
    
    int nufaces = ufaces_vec.size();
    Array<int>* ufaces_arr = new Array<int>(ufaces_vec.size(),1);
    for(int i=0;i<ufaces_vec.size();i++)
    {
        ufaces_arr->setVal(i,0,ufaces_vec[i]);
    }
    
    int nuverts = uverts_vec.size();
    Array<int>* uverts_arr = new Array<int>(uverts_vec.size(),1);
    for(int i=0;i<uverts_vec.size();i++)
    {
        uverts_arr->setVal(i,0,uverts_vec[i]);
    }
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    int* ufaces_nlocs      = new int[world_size];
    int* red_ufaces_nlocs  = new int[world_size];
    int* offset_ufaces     = new int[world_size];

    int* uverts_nlocs      = new int[world_size];
    int* red_uverts_nlocs  = new int[world_size];
    int* offset_uverts     = new int[world_size];

//
    for(i=0;i<world_size;i++)
    {
        ufaces_nlocs[i] = 0;
        red_ufaces_nlocs[i] = 0;

        uverts_nlocs[i] = 0;
        red_uverts_nlocs[i] = 0;

        if(i==world_rank)
        {
            ufaces_nlocs[i] = ufaces_vec.size();
            uverts_nlocs[i] = uverts_vec.size();

        }
        else
        {
            ufaces_nlocs[i] = 0;
            red_ufaces_nlocs[i] = 0;

            uverts_nlocs[i] = 0;
            red_uverts_nlocs[i] = 0;
        }
    }
//
    MPI_Allreduce(ufaces_nlocs, red_ufaces_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(uverts_nlocs, red_uverts_nlocs, world_size, MPI_INT, MPI_SUM, comm);

//
    int o_uf = 0;
    int o_uv = 0;

    for(i=0;i<world_size;i++)
    {
        offset_ufaces[i] = o_uf;
        o_uf = o_uf+red_ufaces_nlocs[i];

        offset_uverts[i] = o_uv;
        o_uv = o_uv+red_uverts_nlocs[i];
    }

    int n_glob_uf = o_uf;
    int n_glob_uv = o_uv;
    
    Array<int>*  lE2gE_uf_g;
    Array<int>*  lE2gE_uv_g;
    
    if(world_rank == 0)
    {
        lE2gE_uf_g   = new Array<int>(n_glob_uf,1);
        lE2gE_uv_g   = new Array<int>(n_glob_uv,1);
    }
    else
    {
        lE2gE_uf_g   = new Array<int>(1,1);
        lE2gE_uv_g   = new Array<int>(1,1);
    }

    //std::cout <<  "rank " << world_rank  << " " << ufaces_vec.size() << " " << n_glob_uf << " " << ufaces_vec[0] << std::endl;
    std::vector<int> ufaces_tot(n_glob_uf);
    MPI_Gatherv(&ufaces_arr->data[0],
                ufaces_arr->getNrow()*1,
                MPI_INT,
                &lE2gE_uf_g->data[0],
                red_ufaces_nlocs,
                offset_ufaces,
                MPI_INT, 0, comm);
    
    if(world_rank == 0)
    {
        for(int k=0;k<n_glob_uf;k++)
        {
            if(lE2gE_uf_g->getVal(k,0) == 153214)
            {
                std::cout << "afab " << k << " " << lE2gE_uf_g->getVal(k,0) << std::endl;
            }
        }
    }
    
//    MPI_Allgatherv(&ufaces_vec[0],
//                   ufaces_vec.size(), MPI_INT,
//                   &ufaces_tot[0],
//                   red_ufaces_nlocs,
//                   offset_ufaces,
//                   MPI_INT,comm);
//
    
//    std::vector<int> dupli = FindDuplicates(ufaces_tot);
//    std::cout << "dupli.size() " << dupli.size() << std::endl;
    std::vector<int> dupli2;
    
    ParallelState* eqDistF = new ParallelState(n_glob_uf,comm);
    int* offsetsF = eqDistF->getOffsets();
    int* nlocsF   = eqDistF->getNlocs();

    int* eqDistFace = new int[nlocsF[world_rank]];
    
    std::vector<int> eqDistFaceV(nlocsF[world_rank]);
    
    std::cout << " g vs l " << n_glob_uf << " " << nlocsF[world_rank] << std::endl;
    
    MPI_Scatterv(&lE2gE_uf_g->data[0], nlocsF, offsetsF,
                 MPI_INT, &eqDistFaceV[0],
                 nlocsF[world_rank], MPI_INT, 0, comm);
    
//    MPI_Scatterv(lE2gE_uf_g->data,
//                nloc_face, MPI_INT,
//                eqDistFace, nloc_face,
//                MPI_INT, 0, MPI_COMM_WORLD);
//
//    for(int u=0;u<nloc_face;u++)
//    {
//        eqDistFaceV[u] = eqDistFace[u];
//    }
    
    dupli2 = FindDuplicatesInParallel_VecV2(eqDistFaceV,
                                            nlocsF[world_rank],
                                            n_glob_uf, comm);
    
    //InteriorPartitionEntity* iFaces = FindDuplicatesInParallel_Vec(ufaces_vec, ufaces_vec.size(), n_glob_uf, comm);
    
    //InteriorPartitionEntity* iVerts = FindDuplicatesInParallel_Vec(uverts_vec, uverts_vec.size(), n_glob_uv, comm);
 
    
    //std::cout << "Rank  Faces " << world_rank << " " << " " << ufaces_vec.size() << " " <<  iFaces->glob2loc_PartInterEntity.size() << " "  << " " << iFaces->PartInterEntity.size() <<  std::endl;
    
    //std::cout << "Rank  Verts " << world_rank << " " << " " << uverts_vec.size() << " " <<  iVerts->glob2loc_PartInterEntity.size() << " " << iVerts->InterEntity.size() << " " << iVerts->PartInterEntity.size() <<  std::endl;
    //std::cout << n_glob_uf << " " << recv2.size() << std::endl;
//    if(world_rank==0)
//    {
//        for(int i=0;i<recv2.size();i++)
//        {
//            //std::cout << recv2[i] << std::endl;
//        }
//    }
//    MPI_Gatherv(&uverts_vec[0],
//                uverts_vec.size(),
//                MPI_INT,
//                &lE2gE_uv_g->data[0],
//                red_uverts_nlocs,
//                offset_uverts,
//                MPI_INT, 0, comm);
    
    MeshTransfer* Tr = new MeshTransfer;
    
    //Tr->iFaces = iFaces;
    //Tr->iVerts = iVerts;

    Tr->UniqueFacesVec = dupli2;
    
    return Tr;

}





Array<int>* GatherTetrahedraOnRoot(std::map<int,std::vector<int> > Ate, MPI_Comm comm, MPI_Info info)
{


    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    std::map<int,std::vector<int> >::iterator grit;
    Array<int>* lE2gE_tet = new Array<int>(Ate.size(),1);
    Array<int>* ien_tet = new Array<int>(Ate.size(),4);
    
    std::cout << "ien_tet " << ien_tet->getNrow() << std::endl;
    i = 0;
    
    for(grit=Ate.begin();grit!=Ate.end();grit++)
    {
        lE2gE_tet->setVal(i,0,grit->first);
        ien_tet->setVal(i,0,grit->second[0]);
        ien_tet->setVal(i,1,grit->second[1]);
        ien_tet->setVal(i,2,grit->second[2]);
        ien_tet->setVal(i,3,grit->second[3]);
        i++;
    }
    
    int* lid_tet_nlocs      = new int[world_size];
    int* ien_tet_nlocs      = new int[world_size];
    
    int* red_lid_tet_nlocs  = new int[world_size];
    int* red_ien_tet_nlocs  = new int[world_size];
    
    int* lid_tet_offsets = new int[world_size];
    int* ien_tet_offsets = new int[world_size];
    
    for(i=0;i<world_size;i++)
    {
        lid_tet_nlocs[i] = 0;
        ien_tet_nlocs[i] = 0;
        
        if(i==world_rank)
        {
            lid_tet_nlocs[i] = Ate.size();
            ien_tet_nlocs[i] = Ate.size()*4;
        }
        else
        {
            lid_tet_nlocs[i] = 0;
            ien_tet_nlocs[i] = 0;
        }
    }

    MPI_Allreduce(lid_tet_nlocs, red_lid_tet_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ien_tet_nlocs, red_ien_tet_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    
    int offset_lid_tet = 0;
    int offset_ien_tet  = 0;
    for(i=0;i<world_size;i++)
    {
        lid_tet_offsets[i] = offset_lid_tet;
        ien_tet_offsets[i] = offset_ien_tet;
        
        offset_lid_tet = offset_lid_tet+red_lid_tet_nlocs[i];
        offset_ien_tet = offset_ien_tet+red_ien_tet_nlocs[i];
    }

    Array<int>*  lE2gE_tet_g;
    Array<int>*  ien_tet_g;

    int n_glob_lid_tet = offset_lid_tet;
    int n_glob_ien_tet = offset_ien_tet;
    if(world_rank == 0)
    {
        lE2gE_tet_g   = new Array<int>(n_glob_lid_tet,1);
        ien_tet_g   = new Array<int>(n_glob_ien_tet,4);
    }
    else
    {
        lE2gE_tet_g   = new Array<int>(1,1);
        ien_tet_g   = new Array<int>(1,1);
    }
    
    int ncol_lid = 1;
    MPI_Gatherv(&lE2gE_tet->data[0],
                lE2gE_tet->getNrow()*1,
                MPI_INT,
                &lE2gE_tet_g->data[0],
                red_lid_tet_nlocs,
                lid_tet_offsets,
                MPI_INT, 0, comm);
//
    
//
    MPI_Gatherv(&ien_tet->data[0],
                ien_tet->getNrow()*4,
                MPI_INT,
                &ien_tet_g->data[0],
                red_ien_tet_nlocs,
                ien_tet_offsets,
                MPI_INT, 0, comm);

    delete lE2gE_tet_g;
    delete lE2gE_tet;
    delete ien_tet;
    
    /**/
    
    return ien_tet_g;

}

std::map<int,std::vector<int> > GatherElementsOnRoot(std::map<int,std::vector<int> >Apr, std::map<int,std::vector<int> > Ate, MPI_Comm comm, MPI_Info info)
{
    
    std::map<int,std::vector<int> > elements;

    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    std::map<int,std::vector<int> >::iterator grit;
    Array<int>* lE2gE_tet = new Array<int>(Ate.size(),1);
    Array<int>* lE2gE_pri = new Array<int>(Apr.size(),1);

    Array<int>* ien_tet = new Array<int>(Ate.size(),4);
    Array<int>* ien_pri = new Array<int>(Apr.size(),6);
    i = 0;
    for(grit=Ate.begin();grit!=Ate.end();grit++)
    {
        lE2gE_tet->setVal(i,0,grit->first);
        ien_tet->setVal(i,0,grit->second[0]);
        ien_tet->setVal(i,1,grit->second[1]);
        ien_tet->setVal(i,2,grit->second[2]);
        ien_tet->setVal(i,3,grit->second[3]);
        i++;
    }
    
    i = 0;
    for(grit=Apr.begin();grit!=Apr.end();grit++)
    {
        lE2gE_pri->setVal(i,0,grit->first);
        ien_pri->setVal(i,0,grit->second[0]);
        ien_pri->setVal(i,1,grit->second[1]);
        ien_pri->setVal(i,2,grit->second[2]);
        ien_pri->setVal(i,3,grit->second[3]);
        ien_pri->setVal(i,4,grit->second[4]);
        ien_pri->setVal(i,5,grit->second[5]);
        //std::cout << grit->second[0] << " " << grit->second[1] << " " << grit->second[2] << " " << grit->second[3] << " " << grit->second[4] << " " << grit->second[5] << std::endl;
        i++;
    }
    
    int* lid_tet_nlocs      = new int[world_size];
    int* lid_pri_nlocs      = new int[world_size];
    int* ien_tet_nlocs      = new int[world_size];
    int* ien_pri_nlocs      = new int[world_size];
    
    int* red_lid_tet_nlocs  = new int[world_size];
    int* red_lid_pri_nlocs  = new int[world_size];
    int* red_ien_tet_nlocs  = new int[world_size];
    int* red_ien_pri_nlocs  = new int[world_size];
    
    int* lid_tet_offsets = new int[world_size];
    int* lid_pri_offsets = new int[world_size];
    int* ien_tet_offsets = new int[world_size];
    int* ien_pri_offsets = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        lid_tet_nlocs[i] = 0;
        lid_pri_nlocs[i] = 0;
        ien_tet_nlocs[i] = 0;
        ien_pri_nlocs[i] = 0;

        
        if(i==world_rank)
        {
            lid_tet_nlocs[i] = Ate.size();
            lid_pri_nlocs[i] = Apr.size();
            ien_tet_nlocs[i] = Ate.size()*4;
            ien_pri_nlocs[i] = Apr.size()*6;
        }
        else
        {
            lid_tet_nlocs[i] = 0;
            lid_pri_nlocs[i] = 0;
            ien_tet_nlocs[i] = 0;
            ien_pri_nlocs[i] = 0;
        }
    }

    MPI_Allreduce(lid_tet_nlocs, red_lid_tet_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(lid_pri_nlocs, red_lid_pri_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ien_tet_nlocs, red_ien_tet_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(ien_pri_nlocs, red_ien_pri_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    
    int offset_lid_tet = 0;
    int offset_lid_pri = 0;
    int offset_ien_tet  = 0;
    int offset_ien_pri  = 0;
    for(i=0;i<world_size;i++)
    {
        lid_tet_offsets[i] = offset_lid_tet;
        lid_pri_offsets[i] = offset_lid_pri;
        ien_tet_offsets[i] = offset_ien_tet;
        ien_pri_offsets[i] = offset_ien_pri;
        
        offset_lid_tet = offset_lid_tet+red_lid_tet_nlocs[i];
        offset_lid_pri = offset_lid_pri+red_lid_pri_nlocs[i];
        offset_ien_tet = offset_ien_tet+red_ien_tet_nlocs[i];
        offset_ien_pri = offset_ien_pri+red_ien_pri_nlocs[i];
    }

    Array<int>*  lE2gE_tet_g;
    Array<int>*  lE2gE_pri_g;
    Array<int>*  ien_tet_g;
    Array<int>*  ien_pri_g;

    int n_glob_lid_tet = offset_lid_tet;
    int n_glob_lid_pri = offset_lid_pri;
    int n_glob_ien_tet = offset_ien_tet;
    int n_glob_ien_pri = offset_ien_pri;
    if(world_rank == 0)
    {
        lE2gE_tet_g   = new Array<int>(n_glob_lid_tet,1);
        lE2gE_pri_g   = new Array<int>(n_glob_lid_pri,1);
        ien_tet_g   = new Array<int>(n_glob_ien_tet,4);
        ien_pri_g   = new Array<int>(n_glob_ien_pri,6);
    }
    else
    {
        lE2gE_tet_g   = new Array<int>(1,1);
        lE2gE_pri_g   = new Array<int>(1,1);
        ien_tet_g   = new Array<int>(1,1);
        ien_pri_g   = new Array<int>(1,1);
    }
    
    int ncol_lid = 1;
    MPI_Gatherv(&lE2gE_tet->data[0],
                lE2gE_tet->getNrow()*1,
                MPI_INT,
                &lE2gE_tet_g->data[0],
                red_lid_tet_nlocs,
                lid_tet_offsets,
                MPI_INT, 0, comm);
//
    MPI_Gatherv(&lE2gE_pri->data[0],
                lE2gE_pri->getNrow()*1,
                MPI_INT,
                &lE2gE_pri_g->data[0],
                red_lid_pri_nlocs,
                lid_pri_offsets,
                MPI_INT, 0, comm);
//
    MPI_Gatherv(&ien_tet->data[0],
                ien_tet->getNrow()*4,
                MPI_INT,
                &ien_tet_g->data[0],
                red_ien_tet_nlocs,
                ien_tet_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&ien_pri->data[0],
                ien_pri->getNrow()*6,
                MPI_INT,
                &ien_pri_g->data[0],
                red_ien_pri_nlocs,
                ien_pri_offsets,
                MPI_INT, 0, comm);



    for(int i=0;i<lE2gE_tet_g->getNrow();i++)
    {
        int gEl = lE2gE_tet_g->getVal(i,0);
        std::vector<int> tet(4);
        tet[0] = ien_tet_g->getVal(i,0);
        tet[1] = ien_tet_g->getVal(i,1);
        tet[2] = ien_tet_g->getVal(i,2);
        tet[3] = ien_tet_g->getVal(i,3);
        elements[gEl] = tet;
    }

    for(int i=0;i<lE2gE_pri_g->getNrow();i++)
    {
        int gEl = lE2gE_pri_g->getVal(i,0);
        std::vector<int> pri(6);
        pri[0] = ien_pri_g->getVal(i,0);
        pri[1] = ien_pri_g->getVal(i,1);
        pri[2] = ien_pri_g->getVal(i,2);
        pri[3] = ien_pri_g->getVal(i,3);
        pri[4] = ien_pri_g->getVal(i,4);
        pri[5] = ien_pri_g->getVal(i,5);
        
        elements[gEl] = pri;
    }


    delete lE2gE_tet_g;
    delete lE2gE_pri_g;
    delete ien_tet_g;
    delete ien_pri_g;

    delete lE2gE_tet;
    delete lE2gE_pri;
    delete ien_tet;
    delete ien_pri;
    
    /**/
    
    
    return elements;

}
