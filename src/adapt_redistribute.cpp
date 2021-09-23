#include "adapt_redistribute.h"


RedistributePartitionObject::RedistributePartitionObject(US3D* us3d,
                        std::map<int,std::vector<int> > tetras,
						std::map<int,std::vector<int> > iferank_map,
                        std::map<int,std::vector<int> > ief_part_map,
                        std::map<int,std::vector<int> > ifn_part_map,
                        std::map<int,std::vector<int> > ife_part_map,
                        std::map<int,int > if_ref_part_map,
                        std::map<int,std::vector<int> > ushell,
                        std::map<int,Array<double>* > M_vmap,
                        MPI_Comm comm)
{
    mpi_comm = comm;
    

    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process

    MPI_Comm_rank(comm, &world_rank);
    
    
    RebasePartitionObject(tetras,
						  iferank_map,
                          ief_part_map,
                          ifn_part_map,
                          ife_part_map,
                          if_ref_part_map,
                          ushell,
                          M_vmap);
    
    
    GetNewGlobalPartitioningTetrahedraMesh();
    
    int* newSizesOnRanks_tmp = new int[world_size];
    int* newSizesOnRanks     = new int[world_size];

    for(int i=0;i<world_size;i++)
    {
        newSizesOnRanks_tmp[i] = 0;
        newSizesOnRanks[i] = 0;
    }
    
    int pid;
    for(int i=0;i<m_part->getNrow();i++)
    {
        pid = m_part->getVal(i,0);
        newSizesOnRanks_tmp[pid] = newSizesOnRanks_tmp[pid]+1;
    }
    
    MPI_Allreduce(newSizesOnRanks_tmp, newSizesOnRanks, world_size, MPI_INT, MPI_SUM, mpi_comm);
    
    // int* newSizesOnRanks    = CommunicatePartitionLayout();
    //
    // Grab the optimal number of elements on current rank based on the partitioning array.
    
    int nTetrahedra         = newSizesOnRanks[world_rank];
    int nTetrahedraGlob     = 0;

    int tetoff = 0;
    for(int i=0;i<world_size;i++)
    {
        nTetrahedraGlob=nTetrahedraGlob+newSizesOnRanks[i];
        tetoff = tetoff + newSizesOnRanks[i];
    }

    
    ParallelState* xcn_pstate = new ParallelState(us3d->xcn->getNglob(),mpi_comm);
    ParallelState* ife_pstate = new ParallelState(us3d->ifn->getNglob(),mpi_comm);
    
    //std::cout << "before world " << world_rank << " #elements " << ief_part_hybrid->getNrow() << " " << m_M_vmap.size() << std::endl;

    UpdateTetrahedraOnPartition(nTetrahedraGlob, nTetrahedra,
                                us3d->xcn, xcn_pstate,
                                ife_pstate);
    
    //std::cout << "after world " << world_rank << " #elements " << ief_part_hybrid->getNrow() << " " << m_M_vmap.size() << std::endl;
    
    GetFace2NodeRedistributedMesh(us3d->ifn, 3, ife_pstate,
                                  nTetrahedraGlob,mpi_comm);
    
    GetFace2RankTetrahedraMesh();
    //GetShellVert2RefMap();
    GetShellVert2RefMap_Global();
    /**/
    
}


void RedistributePartitionObject::RebasePartitionObject(
                                                   std::map<int,std::vector<int> > tetras,
												   std::map<int,std::vector<int> > iferank_part_map,
                                                   std::map<int,std::vector<int> > ief_part_map,
                                                   std::map<int,std::vector<int> > ifn_part_map,
                                                   std::map<int,std::vector<int> > ife_part_map,
                                                   std::map<int,int> if_ref_part_map,
                                                   std::map<int,std::vector<int> > ushell,
                                                   std::map<int,Array<double>* > M_vmap)
{
        
    int shfn = 0;
    int shf = 0;
    int i;
    int nTetras = tetras.size();
    int gvid,gfid;
    
    std::map<int,std::vector<int> >::iterator ite;
    int r0,r1,el0,el1,pos,ra;
    int ref;
    int lcv = 0;
    std::map<int,std::vector<int> > ref2bcface;

    int lf  = 0;
    std::map<int,int> sharedFaces;

    //int* ielement_offsets     = new int[world_size];
    //int* ielement_nlocs       = new int[world_size];

    int elTel = 0;
    
    std::set<int> gvid_set;
    int nLocalVerts = 0;
    
//    std::map<int,int> lF2gF_tets;
//    std::map<int,int> gF2lF_tets;
    std::set<int> gfid_set;
    int nLocalFaces = 0;
    std::map<int,int> face2ref;
    std::set<int> ufaces;
    for(ite=tetras.begin();ite!=tetras.end();ite++)
    {
        //key = GID, value global node nmber;
        int gEl = ite->first;
        
        for(int j=0;j<ief_part_map[gEl].size();j++)
        {
            gfid = ief_part_map[gEl][j];
            
            if(ushell.find(gfid)!=ushell.end())
            {
                ref = 13;
            }
            else
            {
                ref = if_ref_part_map[gfid];
            }
            
            
            if(gfid_set.find(gfid)==gfid_set.end())
            {
                gfid_set.insert(gfid);
                face2ref[gfid] = ref;
                nLocalFaces++;
            }
            
            for(int k=0;k<ifn_part_map[gfid].size();k++)
            {
                gvid = ifn_part_map[gfid][k];
                
                if(gvid_set.find(gvid)==gvid_set.end())
                {
                    gvid_set.insert(gvid);
                    nLocalVerts++;
                }
            }
            
            if(ufaces.find(gfid)==ufaces.end())
            {
                ufaces.insert(gfid);

                el0    = ife_part_map[gfid][0];
                el1    = ife_part_map[gfid][1];
                
                r0 	   = iferank_part_map[gfid][0];
                r1 	   = iferank_part_map[gfid][1];
                
                if(ref==2)
                {
                    if(r0==world_rank && r1!=world_rank)
                    {
                        sharedFaces[gfid] = r0;
                        shf++;
                    }
                    if(r0!=world_rank && r1==world_rank)
                    {
                        sharedFaces[gfid] = r1;
                        shf++;
                    }
                }
            }
        }
        elTel++;
    }
    
    int nSharedFaces                          = sharedFaces.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,mpi_comm);
    DistributedParallelState* distLocalVerts  = new DistributedParallelState(nLocalVerts,mpi_comm);
    DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,mpi_comm);

    int Nt_shFaces               = distSharedFaces->getNel();
    int* shFace_offsets          = distSharedFaces->getOffsets();
    int* shFace_nlocs            = distSharedFaces->getNlocs();
    int* shFacesIDs              = new int[nSharedFaces];
    int* shFaces_RankIDs         = new int[nSharedFaces];

    int iter = 0;
    std::set<int> UniqueSharedVertsOnRank_set;
    std::vector<int> UniqueSharedVertsOnRank;
    std::vector<int> UniqueSharedVertsOnRank_RankID;
    
    std::map<int,int>::iterator itsf;
    int lvrtid = 0;
    int tel = shFace_offsets[world_rank];
    for(itsf=sharedFaces.begin();itsf!=sharedFaces.end();itsf++)
    {
        shFacesIDs[iter] = itsf->first;
        shFaces_RankIDs[iter] = itsf->second;
        gfid = itsf->first;

        for(int q=0;q<ifn_part_map[gfid].size();q++)
        {
            gvid   = ifn_part_map[gfid][q];
            
            if(UniqueSharedVertsOnRank_set.find(gvid)==UniqueSharedVertsOnRank_set.end())
            {
                UniqueSharedVertsOnRank_set.insert(gvid);
                UniqueSharedVertsOnRank.push_back(gvid);
                UniqueSharedVertsOnRank_RankID.push_back(world_rank);
                lvrtid++;
            }
        }
        tel++;
        iter++;
    }

    int nSharedVerts = UniqueSharedVertsOnRank.size();

    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,mpi_comm);
    
    int Nt_shVerts               = distSharedVerts->getNel();
    int* shVerts_nlocs           = distSharedVerts->getNlocs();
    int* shVerts_offsets         = distSharedVerts->getOffsets();
    
    int* TotalSharedVerts        = new int[Nt_shVerts];
    int* TotalSharedVerts_RankID = new int[Nt_shVerts];
    int* TotalSharedFaces        = new int[Nt_shFaces];
    int* TotalSharedFaces_RankID = new int[Nt_shFaces];
    
    // Communicate vert map to all ranks.
    MPI_Allgatherv(&UniqueSharedVertsOnRank[0],
                   nSharedVerts,
                   MPI_INT,
                   TotalSharedVerts,
                   shVerts_nlocs,
                   shVerts_offsets,
                   MPI_INT, mpi_comm);
    
    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
                   nSharedVerts,
                   MPI_INT,
                   TotalSharedVerts_RankID,
                   shVerts_nlocs,
                   shVerts_offsets,
                   MPI_INT, mpi_comm);
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(shFacesIDs,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFaces,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, mpi_comm);
    
    MPI_Allgatherv(shFaces_RankIDs,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFaces_RankID,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, mpi_comm);
    
    int tmp;
    std::map<int,int> f2r;
    std::set<int> f2r_s;
    
    
    int* NewGlobVertCountPerRank = new int[world_size];
    int* NewGlobFaceCountPerRank = new int[world_size];

    for(int u=0;u<world_size;u++)
    {
        NewGlobVertCountPerRank[u] = 0;
        NewGlobFaceCountPerRank[u] = 0;
    }
        
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];
        
        if(f2r_s.find(key)==f2r_s.end())
        {
            f2r_s.insert(key);
            f2r[key]=val;
            NewGlobFaceCountPerRank[val]=NewGlobFaceCountPerRank[val]+1;
        }
        else
        {
            tmp = f2r[key];
            if(val<tmp)
            {
                f2r[key]=val;
            }
            if(val>tmp)
            {
                f2r[key]=tmp;
            }
        }
    }
    
    std::map<int,int> v2r;

    for(int i=0;i<Nt_shVerts;i++)
    {
        int key = TotalSharedVerts[i];
        int val = TotalSharedVerts_RankID[i];
        
        if(v2r.find(key)==v2r.end())
        {
            v2r[key]=val;
            NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
        }
        else
        {
            tmp = v2r[key];
            
            if(val<tmp)
            {
                v2r[key]=val;

            }
            if(val>tmp)
            {
                v2r[key]=tmp;
            }
        }
    }

    int* NewGlobVertOffset = new int[world_size];
    int* NewGlobFaceOffset = new int[world_size];


    for(int u=1;u<world_size;u++)
    {
        NewGlobVertOffset[u] = NewGlobVertOffset[u-1]+NewGlobVertCountPerRank[u-1];
        NewGlobFaceOffset[u] = NewGlobFaceOffset[u-1]+NewGlobFaceCountPerRank[u-1];
    }
     
    std::map<int,std::vector<int> >::iterator itv;

    std::map<int,int> sharedVertsGlobal;
    
    for(int u=0;u<world_size;u++)
    {
        NewGlobVertCountPerRank[u] = 0;
        NewGlobFaceCountPerRank[u] = 0;
    }
    
    int iVshared = distLocalVerts->getNel()-v2r.size();

    std::map<int,int >::iterator itvv;
    std::map<int,int> sharedVmap;
    for(itvv=v2r.begin();itvv!=v2r.end();itvv++)
    {
        sharedVmap[itvv->first] = iVshared;
        iVshared++;
    }
    
    std::map<int,int> sharedFmap;
    int iFshared = distLocalFaces->getNel()-f2r.size();

    for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
    {
        sharedFmap[itvv->first] = iFshared;
        iFshared++;
    }
    
//    int nNonSharedVerts          = nLocalVerts-owned_verts[world_rank];
//    int nNonSharedFaces          = nLocalFaces-owned_faces[world_rank];
    
    int nNonSharedVerts          = nLocalVerts-nSharedVerts;
    int nNonSharedFaces          = nLocalFaces-nSharedFaces;

    int* nNonSharedArray         = new int[world_size];
    int* nNonSharedArrayRed      = new int[world_size];
    int* nNonSharedVertsArrayOff = new int[world_size];

    int* nNonSharedFacesArray    = new int[world_size];
    int* nNonSharedFacesArrayRed = new int[world_size];
    int* nNonSharedFacesArrayOff = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        nNonSharedArray[i]      = 0;
        nNonSharedFacesArray[i] = 0;
        if(i==world_rank)
        {
            nNonSharedArray[i]      = nNonSharedVerts;
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }
    
    MPI_Allreduce(nNonSharedArray,
                  nNonSharedArrayRed,
                  world_size,
                  MPI_INT, MPI_SUM, mpi_comm);
    
    MPI_Allreduce(nNonSharedFacesArray,
                  nNonSharedFacesArrayRed,
                  world_size,
                  MPI_INT, MPI_SUM, mpi_comm);
    
    int nonSharedOff          = 0;
    int nonFacesSharedOff     = 0;
    for(i=0;i<world_size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }
    
    //=================================================================================================
    //=================================================================================================
    //=================================================================================================

    int* ini_nEl      = new int[world_size];
    int* red_ini_nEl  = new int[world_size];
    int* ini_offsetEl = new int[world_size];

    

    
    for(int i=0;i<world_size;i++)
    {
        ini_nEl[i]      = 0;
        red_ini_nEl[i]  = 0;
        ini_offsetEl[i] = 0;
        
        if(i==world_rank)
        {
            ini_nEl[i] = nTetras;
        }
    }
    MPI_Allreduce(ini_nEl, red_ini_nEl, world_size, MPI_INT, MPI_SUM, mpi_comm);

    int offsetEl = 0;
    
    for(int i=0;i<world_size;i++)
    {
        ini_offsetEl[i] = offsetEl;
        offsetEl        = offsetEl + red_ini_nEl[i];
    }
    
    int size = world_size;
    int optimalSize = int(offsetEl/size) + ( world_rank < offsetEl%size );
    //double rat = (double)nTetras/(double)optimalSize;
        
    int NtoRecv = 0;
    int NtoSend = 0;
    
    if(nTetras>optimalSize)
    {
        NtoSend = nTetras-optimalSize;
    }
    if(nTetras<optimalSize)
    {
        NtoRecv = optimalSize-nTetras;
    }

   
    int* toR_red2 = new int[world_size];
    int* toS_red2 = new int[world_size];
    int* toS_red = new int[world_size];
    int* toR_red = new int[world_size];
    int* optiSize = new int[world_size];
    int* optiSize_red = new int[world_size];
    int* old_ntets = new int[world_size];
    int* old_ntets_red = new int[world_size];
    int* toS = new int[world_size];
    int* toR = new int[world_size];
    

    for(int i=0;i<world_size;i++)
    {
        toS[i] = 0;
        toR[i] = 0;
        optiSize[i] = 0;
        optiSize_red[i] = 0;
        toS_red[i] = 0;
        toR_red[i] = 0;
        toS_red2[i] = 0;
        toR_red2[i] = 0;
        old_ntets[i] = 0;
        old_ntets_red[i] = 0;
        if(i==world_rank)
        {
            optiSize[i] = optimalSize;
            toS[i] = NtoSend;
            toR[i] = NtoRecv;
            old_ntets[i] = nTetras;

        }
    }
    MPI_Allreduce(optiSize, optiSize_red, world_size, MPI_INT, MPI_SUM, mpi_comm);
    MPI_Allreduce(toS,           toS_red, world_size, MPI_INT, MPI_SUM, mpi_comm);
    MPI_Allreduce(toR,           toR_red, world_size, MPI_INT, MPI_SUM, mpi_comm);
    MPI_Allreduce(&old_ntets[0], &old_ntets_red[0], world_size, MPI_INT, MPI_SUM, mpi_comm);
    int sent;
    int sendUpdate;
    
    int opti_Rank,ntet_Rank;
    int ns=0;
    int nr=0;
    int nc=0;
    for(int i=0;i<world_size;i++)
    {
        opti_Rank = optiSize_red[i];
        ntet_Rank = old_ntets_red[i];
        
        if(ntet_Rank>opti_Rank)
        {
            toS_red2[i] = old_ntets_red[i]-optiSize_red[i];
            toR_red2[i] = 0;
            
            ns++;
        }
        if(opti_Rank>ntet_Rank)
        {
            toS_red2[i] = 0;
            toR_red2[i] = optiSize_red[i]-old_ntets_red[i];
            nr++;
        }
        if(opti_Rank==ntet_Rank)
        {
            toS_red2[i] = 0;
            toR_red2[i] = 0;
            nc++;
        }
    }
    
    int* Psending   = new int[ns];
    int* Nsending   = new int[ns];
    int* Preceiving = new int[nr];
    int* Nreceiving = new int[nr];

    std::set<int> Pstay;
    
    std::map<int,std::vector<int> > recvRa;
    std::map<int,std::vector<int> > recvNe;
    
    std::map<int,std::vector<int> > sendRa;
    std::map<int,std::vector<int> > sendNe;
    
    int r = 0;
    int s = 0;
    
    for(int i=0;i<world_size;i++)
    {
        if(toS_red2[i]>0)
        {
            Psending[s] = i;
            Nsending[s] = toS_red2[i];
            s++;
        }
        if(toR_red2[i]>0)
        {
            Preceiving[r] = i;
            Nreceiving[r] = toR_red2[i];
            r++;
        }
        if(toS_red2[i]==0 && toR_red2[i]==0)
        {
            Pstay.insert(i);
        }
    }
    
    
    
    
    
    int adv = 0;
    int st = 0;
    int residual = 0;
    int Psend;
    while(adv<ns)
    {
        
        int dist  = Nsending[adv];
        
        if(dist!=0)
        {
            Psend = Psending[adv];
        }

        std::vector<int> toRank;
        std::vector<int> NtoRank;
        //std::cout << Psend << " with " << dist << " sends -> ";
        
        while(dist!=0)
        {

            int PtoS = Preceiving[st];

            toRank.push_back(PtoS);

            if(residual != 0)
            {
                if(dist>residual)
                {
                    dist = dist - residual;
                    NtoRank.push_back(residual);
                    residual = 0;
                    st++;
                }
                else
                {
                    NtoRank.push_back(dist);
                    residual = residual - dist;
                    dist     = 0;
                }
            }
            else if(dist>Nreceiving[st])
            {
                dist = dist - Nreceiving[st];
                NtoRank.push_back(Nreceiving[st]);
                st++;
            }
            else
            {
                NtoRank.push_back(dist);
                residual = Nreceiving[st]-dist;
                dist = 0;
            }
        }
        
        //std::cout << std::endl;
//
        sendRa[Psend]=toRank;
        sendNe[Psend]=NtoRank;
        
//        if(dist!=0)
//        {
//
//        }
        
        adv++;
    }
    

    
    std::map<int,std::vector<int> >::iterator its;

    for(its=sendRa.begin();its!=sendRa.end();its++)
    {
        for(int q=0;q<its->second.size();q++)
        {
            recvRa[its->second[q]].push_back(its->first);
            recvNe[its->second[q]].push_back(sendNe[its->first][q]);
        }
    }

    
//    if(world_rank==0)
//    {
//
//        for(int i=0;i<world_size;i++)
//        {
//            std::cout << i << " " << toS_red2[i]<< " " <<toR_red2[i] << " " << optiSize_red[i] << " " << old_ntets_red[i] << std::endl;
//        }
//        std::map<int,std::vector<int> >::iterator its;
//
//        for(its=sendRa.begin();its!=sendRa.end();its++)
//        {
//            std::cout << its->first << " S-> ";
//            for(int q=0;q<its->second.size();q++)
//            {
//                std::cout << its->second[q] << " " << sendNe[its->first][q] << " ";
//            }
//            std::cout << std::endl;
//        }
//
//        for(its=recvRa.begin();its!=recvRa.end();its++)
//        {
//            std::cout << its->first << " R-> ";
//            for(int q=0;q<its->second.size();q++)
//            {
//                std::cout << its->second[q] << " " << recvNe[its->first][q] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    
    //=================================================================================================
    //================================================================================================
    //=================================================================================================
    
    int nTet    = 0;
    int elloc   = 0;
    int lbvid   = nNonSharedVertsArrayOff[world_rank];
    
    int lbfid   = nNonSharedFacesArrayOff[world_rank];
    int ltetvid = 0;
    
    std::set<int> gl_set;
    std::map<int,int> gl_map;
    int lvid2   = 0;
    
    std::set<int> glf_set;
    std::map<int,int> glf_map;
    int lbfids  = 0;
    int lfid2   = 0;
    
    int dd      = 0;
    int Nsend   = 0;
    int nv      = 4;
    int c       = 0;
    int if_ref  = 0;
    
    Array<int>* new_GidEl   = new Array<int>(optimalSize,4);
    Array<int>* new_ien     = new Array<int>(optimalSize,4);
    Array<int>* new_ien_or  = new Array<int>(optimalSize,4);
    Array<int>* new_ief     = new Array<int>(optimalSize,4);
    Array<int>* new_ief_or  = new Array<int>(optimalSize,4);
    Array<int>* new_iefref  = new Array<int>(optimalSize,4);
    int nadded = 0;
    
    if(Pstay.find(world_rank)!=Pstay.end())
    {
        //copy all data into tmesh
        int uloc = 0;
        int u = 0;
        for(ite=tetras.begin();ite!=tetras.end();ite++)
        {
            int gEl       = ite->first;
            //int lEl       = ini_offsetEl[world_rank]+u;
            int* ien      = new int[4];
            int* ien_o    = new int[4];
            int* ief      = new int[4];
            int* ief_o    = new int[4];
            int* iefref   = new int[4];
            double* Met   = new double[6*4];
            
            for(int q=0;q<ite->second.size();q++)
            {
                gvid = ite->second[q];
                
                ien_o[q]   = gvid;
                
                Met[q*6+0] = M_vmap[gvid]->getVal(0,0);
                Met[q*6+1] = M_vmap[gvid]->getVal(0,1);
                Met[q*6+2] = M_vmap[gvid]->getVal(0,2);
                Met[q*6+3] = M_vmap[gvid]->getVal(1,1);
                Met[q*6+4] = M_vmap[gvid]->getVal(1,2);
                Met[q*6+5] = M_vmap[gvid]->getVal(2,2);
                                 
                if(v2r.find(gvid)!=v2r.end())
                {
                    lvid2  = sharedVmap[gvid];
                    ien[q] = lvid2;
                }
                else
                {
                    if(gl_set.find(gvid)==gl_set.end())
                    {
                        gl_set.insert(gvid);
                        gl_map[gvid] = lbvid;
                        ien[q]       = lbvid;
                        lbvid        = lbvid + 1;
                    }
                    else
                    {
                        int lbbvid   = gl_map[gvid];
                        ien[q]       = lbbvid;
                    }
                }
            }

            for(int q=0;q<ief_part_map[gEl].size();q++)
            {
                gfid = ief_part_map[gEl][q];
                
                if(f2r.find(gfid)!=f2r.end())
                {
                    lfid2      = sharedFmap[gfid];
                    ief[q]     = lfid2;
                    ief_o[q]   = gfid;
                    iefref[q]  = face2ref[gfid];
                }
                else
                {
                    if(glf_set.find(gfid)==glf_set.end())
                    {
                        glf_set.insert(gfid);
                        glf_map[gfid]   = lbfid;
                        ief[q]          = lbfid;
                        ief_o[q]        = gfid;
                        iefref[q]       = face2ref[gfid];
                        lbfid           = lbfid + 1;
                    }
                    else
                    {
                        int lbbfid  = glf_map[gfid];
                        ief[q]      = lbbfid;
                        ief_o[q]    = gfid;
                        iefref[q]   = face2ref[gfid];
                    }
                }
            }
            
            new_GidEl->setVal(uloc,0,gEl);
            new_ien->setVal(uloc,0,ien[0]);
            new_ien->setVal(uloc,1,ien[1]);
            new_ien->setVal(uloc,2,ien[2]);
            new_ien->setVal(uloc,3,ien[3]);

            new_ien_or->setVal(uloc,0,ien_o[0]);
            new_ien_or->setVal(uloc,1,ien_o[1]);
            new_ien_or->setVal(uloc,2,ien_o[2]);
            new_ien_or->setVal(uloc,3,ien_o[3]);
                                        
            new_ief->setVal(uloc,0,ief[0]);
            new_ief->setVal(uloc,1,ief[1]);
            new_ief->setVal(uloc,2,ief[2]);
            new_ief->setVal(uloc,3,ief[3]);
            
            new_ief_or->setVal(uloc,0,ief_o[0]);
            new_ief_or->setVal(uloc,1,ief_o[1]);
            new_ief_or->setVal(uloc,2,ief_o[2]);
            new_ief_or->setVal(uloc,3,ief_o[3]);

            new_iefref->setVal(uloc,0,iefref[0]);
            new_iefref->setVal(uloc,1,iefref[1]);
            new_iefref->setVal(uloc,2,iefref[2]);
            new_iefref->setVal(uloc,3,iefref[3]);
            
            
            for(int s=0;s<4;s++)
            {
                if(m_M_vmap.find(ien[s])==m_M_vmap.end())
                {
                    Array<double>* Mtrcs = new Array<double>(6,1);
                    Mtrcs->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                    Mtrcs->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                    Mtrcs->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                    Mtrcs->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                    Mtrcs->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                    Mtrcs->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                    
                    m_M_vmap[ien[s]]=Mtrcs;
                }
            }
            
            delete[] Met;
            delete[] iefref;
            delete[] ief_o;
            delete[] ief;
            delete[] ien_o;
            delete[] ien;
            
            
            uloc++;
            u++;
        }
    }
    
    
    if(sendRa.find(world_rank)!=sendRa.end())
    {
        std::vector<int> toRanks    = sendRa[world_rank];
        std::vector<int> NeltoRanks = sendNe[world_rank];
        
        std::vector<std::vector<int> > elIDs;
        std::vector<std::vector<int> > elNodeIDs;
        std::vector<std::vector<int> > elNodeOriginalIDs;
        std::vector<std::vector<double> > elNodeMetrics;
        std::vector<std::vector<int> > elFaceIDs;
        std::vector<std::vector<int> > elFaceRefs;
        std::vector<std::vector<int> > elFaceOriginalIDs;
        
        for(int i=0;i<toRanks.size();i++)
        {
            int Nel = NeltoRanks[i];
            std::vector<int> rowEl(Nel);
            std::vector<int> rowNode(Nel*4);
            std::vector<double> rowNodeMetric(Nel*4*6);
            std::vector<int> rowFace(Nel*4);
            std::vector<int> rowFaceRef(Nel*4);
            std::vector<int> rowNodeOriginal(Nel*4);
            std::vector<int> rowFaceOriginal(Nel*4);

            elIDs.push_back(rowEl);
            elNodeIDs.push_back(rowNodeOriginal);
            elNodeOriginalIDs.push_back(rowNode);
            elNodeMetrics.push_back(rowNodeMetric);
            elFaceIDs.push_back(rowFace);
            elFaceRefs.push_back(rowFaceRef);
            elFaceOriginalIDs.push_back(rowFaceOriginal);
        }
        
        int cc = 0;
        //int sRank    = toRanks[0];
        
        int offPrank = 0;
        int cntv     = 0;
        
        int t = 0;
        int nuloc    = 0;
        int uloc     = 0;
        int u        = 0;

        for(ite=tetras.begin();ite!=tetras.end();ite++)
        {
            int nelPrank  = NeltoRanks[cc];
            //int sRank     = toRanks[cc];
            int gEl       = ite->first;
            //int lEl       = ini_offsetEl[world_rank]+u;
            int* ien      = new int[4];
            int* ien_o    = new int[4];
            int* ief      = new int[4];
            int* ief_o    = new int[4];
            int* iefref   = new int[4];
            double* Met      = new double[6*4];
            
            for(int q=0;q<4;q++)
            {
                gvid = ite->second[q];
                
                ien_o[q]   = gvid;
                
                Met[q*6+0] = M_vmap[gvid]->getVal(0,0);
                Met[q*6+1] = M_vmap[gvid]->getVal(0,1);
                Met[q*6+2] = M_vmap[gvid]->getVal(0,2);
                Met[q*6+3] = M_vmap[gvid]->getVal(1,1);
                Met[q*6+4] = M_vmap[gvid]->getVal(1,2);
                Met[q*6+5] = M_vmap[gvid]->getVal(2,2);
                                 
                if(v2r.find(gvid)!=v2r.end())
                {
                    lvid2    = sharedVmap[gvid];
                    ien[q]   = lvid2;
                }
                else
                {
                    if(gl_set.find(gvid)==gl_set.end())
                    {
                        gl_set.insert(gvid);
                        gl_map[gvid]        = lbvid;
                        ien[q]              = lbvid;
                        lbvid               = lbvid + 1;
                    }
                    else
                    {
                        int lbbvid = gl_map[gvid];
                        ien[q] = lbbvid;
                    }
                }
            }

            for(int q=0;q<4;q++)
            {
                gfid = ief_part_map[gEl][q];
                
                if(f2r.find(gfid)!=f2r.end())
                {
                    lfid2      = sharedFmap[gfid];
                    ief[q]     = lfid2;
                    ief_o[q]   = gfid;
                    iefref[q]  = face2ref[gfid];
                }
                else
                {
                    if(glf_set.find(gfid)==glf_set.end())
                    {
                        glf_set.insert(gfid);
                        glf_map[gfid]   = lbfid;
                        ief[q]          = lbfid;
                        ief_o[q]        = gfid;
                        iefref[q]       = face2ref[gfid];
                        lbfid           = lbfid + 1;
                    }
                    else
                    {
                        int lbbfid  = glf_map[gfid];
                        ief[q]      = lbbfid;
                        ief_o[q]    = gfid;
                        iefref[q]   = face2ref[gfid];
                    }
                }
            }
            
            if(u<toS_red2[world_rank])
            {
                if(t<(nelPrank))
                {
                    elIDs[cc][t] = gEl;
                    
                    elNodeIDs[cc][4*t+0]=ien[0];
                    elNodeIDs[cc][4*t+1]=ien[1];
                    elNodeIDs[cc][4*t+2]=ien[2];
                    elNodeIDs[cc][4*t+3]=ien[3];
                    
                    elNodeOriginalIDs[cc][4*t+0]=ien_o[0];
                    elNodeOriginalIDs[cc][4*t+1]=ien_o[1];
                    elNodeOriginalIDs[cc][4*t+2]=ien_o[2];
                    elNodeOriginalIDs[cc][4*t+3]=ien_o[3];
                    
                    elNodeMetrics[cc][4*6*t+0] = Met[0*6+0];
                    elNodeMetrics[cc][4*6*t+1] = Met[0*6+1];
                    elNodeMetrics[cc][4*6*t+2] = Met[0*6+2];
                    elNodeMetrics[cc][4*6*t+3] = Met[0*6+3];
                    elNodeMetrics[cc][4*6*t+4] = Met[0*6+4];
                    elNodeMetrics[cc][4*6*t+5] = Met[0*6+5];

                    elNodeMetrics[cc][4*6*t+6]  = Met[1*6+0];
                    elNodeMetrics[cc][4*6*t+7]  = Met[1*6+1];
                    elNodeMetrics[cc][4*6*t+8]  = Met[1*6+2];
                    elNodeMetrics[cc][4*6*t+9]  = Met[1*6+3];
                    elNodeMetrics[cc][4*6*t+10] = Met[1*6+4];
                    elNodeMetrics[cc][4*6*t+11] = Met[1*6+5];

                    elNodeMetrics[cc][4*6*t+12] = Met[2*6+0];
                    elNodeMetrics[cc][4*6*t+13] = Met[2*6+1];
                    elNodeMetrics[cc][4*6*t+14] = Met[2*6+2];
                    elNodeMetrics[cc][4*6*t+15] = Met[2*6+3];
                    elNodeMetrics[cc][4*6*t+16] = Met[2*6+4];
                    elNodeMetrics[cc][4*6*t+17] = Met[2*6+5];

                    elNodeMetrics[cc][4*6*t+18] = Met[3*6+0];
                    elNodeMetrics[cc][4*6*t+19] = Met[3*6+1];
                    elNodeMetrics[cc][4*6*t+20] = Met[3*6+2];
                    elNodeMetrics[cc][4*6*t+21] = Met[3*6+3];
                    elNodeMetrics[cc][4*6*t+22] = Met[3*6+4];
                    elNodeMetrics[cc][4*6*t+23] = Met[3*6+5];
                    
                    
                    elFaceIDs[cc][4*t+0]=ief[0];
                    elFaceIDs[cc][4*t+1]=ief[1];
                    elFaceIDs[cc][4*t+2]=ief[2];
                    elFaceIDs[cc][4*t+3]=ief[3];
                    
                    elFaceOriginalIDs[cc][4*t+0]=ief_o[0];
                    elFaceOriginalIDs[cc][4*t+1]=ief_o[1];
                    elFaceOriginalIDs[cc][4*t+2]=ief_o[2];
                    elFaceOriginalIDs[cc][4*t+3]=ief_o[3];
                
                    elFaceRefs[cc][4*t+0]=iefref[0];
                    elFaceRefs[cc][4*t+1]=iefref[1];
                    elFaceRefs[cc][4*t+2]=iefref[2];
                    elFaceRefs[cc][4*t+3]=iefref[3];
                    
                   
                    
                    if(t==nelPrank-1)
                    {
                        t = 0;
                        cc=cc+1;
                    }
                    else
                    {
                        t=t+1;
                    }
                }
                
                nuloc++;
            }
            else
            {
                new_GidEl->setVal(uloc,0,gEl);
                
                new_ien->setVal(uloc,0,ien[0]);
                new_ien->setVal(uloc,1,ien[1]);
                new_ien->setVal(uloc,2,ien[2]);
                new_ien->setVal(uloc,3,ien[3]);

                new_ien_or->setVal(uloc,0,ien_o[0]);
                new_ien_or->setVal(uloc,1,ien_o[1]);
                new_ien_or->setVal(uloc,2,ien_o[2]);
                new_ien_or->setVal(uloc,3,ien_o[3]);
                                            
                new_ief->setVal(uloc,0,ief[0]);
                new_ief->setVal(uloc,1,ief[1]);
                new_ief->setVal(uloc,2,ief[2]);
                new_ief->setVal(uloc,3,ief[3]);
                
                new_ief_or->setVal(uloc,0,ief_o[0]);
                new_ief_or->setVal(uloc,1,ief_o[1]);
                new_ief_or->setVal(uloc,2,ief_o[2]);
                new_ief_or->setVal(uloc,3,ief_o[3]);

                new_iefref->setVal(uloc,0,iefref[0]);
                new_iefref->setVal(uloc,1,iefref[1]);
                new_iefref->setVal(uloc,2,iefref[2]);
                new_iefref->setVal(uloc,3,iefref[3]);
                
                for(int s=0;s<4;s++)
                {
                    if(m_M_vmap.find(ien[s])==m_M_vmap.end())
                    {
                        Array<double>* Mtrcs = new Array<double>(6,1);
                        Mtrcs->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                        Mtrcs->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                        Mtrcs->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                        Mtrcs->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                        Mtrcs->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                        Mtrcs->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                        
                        m_M_vmap[ien[s]]=Mtrcs;
                    }
                }
                

                uloc++;
            }
            //delete[] Met;
            u++;
        }
        
        int acull = 0;
        for(int i=0;i<toRanks.size();i++)
        {
            int dest     = toRanks[i];
            int n_Ele    = NeltoRanks[i];
            int n_Vrt    = NeltoRanks[i]*4;
            int n_VrtMet = NeltoRanks[i]*4*6;
            
            std::vector<int> ElIDs          = elIDs[i];
            std::vector<int> Elvec          = elNodeIDs[i];
            std::vector<int> ElFvec         = elFaceIDs[i];
            std::vector<int> Elovec         = elNodeOriginalIDs[i];
            std::vector<double> Vmetrics    = elNodeMetrics[i];
            std::vector<int> ElFovec        = elFaceOriginalIDs[i];
            std::vector<int> ElFrefvec      = elFaceRefs[i];
            MPI_Send(&n_Vrt         ,     1,     MPI_INT, dest, dest,       mpi_comm);
            MPI_Send(&Elvec[0]      , n_Vrt,     MPI_INT, dest, dest*100,   mpi_comm);
            MPI_Send(&ElFvec[0]     , n_Vrt,     MPI_INT, dest, dest*500,   mpi_comm);
            MPI_Send(&Elovec[0]     , n_Vrt,     MPI_INT, dest, dest*200,   mpi_comm);
            MPI_Send(&ElFrefvec[0]  , n_Vrt,     MPI_INT, dest, dest*20000, mpi_comm);
            MPI_Send(&ElFovec[0]    , n_Vrt,     MPI_INT, dest, dest*40000, mpi_comm);
            MPI_Send(&n_Ele         ,     1,     MPI_INT, dest, dest*50000, mpi_comm);
            MPI_Send(&ElIDs[0]      , n_Ele,     MPI_INT, dest, dest*60000, mpi_comm);
            MPI_Send(&n_VrtMet      ,     1,     MPI_INT, dest, dest*70000, mpi_comm);
            MPI_Send(&Vmetrics[0]   , n_VrtMet,  MPI_DOUBLE, dest, dest*80000, mpi_comm);

            acull = acull + n_Vrt;
        }
        
        //std::cout << "Stored at first  " << world_rank << " " << M_vmap.size() << std::endl;
    }
    
    if(recvRa.find(world_rank)!=recvRa.end())
    {
        std::vector<int > expFromRank = recvRa[world_rank];
        
        std::map<int,std::vector<int> > collected_ElIds;
        std::map<int,std::vector<int> > collected_NIds;
        std::map<int,std::vector<int> > collected_OriginalNIds;
        std::map<int,std::vector<int> > collected_FIds;
        std::map<int,std::vector<int> > collected_Frefs;
        std::map<int,std::vector<int> > collected_FOriginalNIds;
        std::map<int,std::vector<double> > collected_VrtMetrics;
        for(int i=0;i<expFromRank.size();i++)
        {
            int origin = expFromRank[i];
            int n_Elr;
            MPI_Recv(&n_Elr,   1, MPI_INT, origin, world_rank, mpi_comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvNElVec(n_Elr);
            MPI_Recv(&recvNElVec[0], n_Elr, MPI_INT, origin, world_rank*100, mpi_comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFElVec(n_Elr);
            MPI_Recv(&recvFElVec[0], n_Elr, MPI_INT, origin, world_rank*500, mpi_comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvONElVec(n_Elr);
            MPI_Recv(&recvONElVec[0],n_Elr, MPI_INT, origin, world_rank*200, mpi_comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFrefElVec(n_Elr);
            MPI_Recv(&recvFrefElVec[0], n_Elr, MPI_INT, origin, world_rank*20000, mpi_comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFONElVec(n_Elr);
            MPI_Recv(&recvFONElVec[0], n_Elr, MPI_INT, origin, world_rank*40000, mpi_comm, MPI_STATUS_IGNORE);
            
            int n_EleRecv;
            MPI_Recv(&n_EleRecv,   1, MPI_INT, origin, world_rank*50000, mpi_comm, MPI_STATUS_IGNORE);

            std::vector<int> recvElIDsvec(n_EleRecv);
            MPI_Recv(&recvElIDsvec[0],   n_EleRecv, MPI_INT, origin, world_rank*60000, mpi_comm, MPI_STATUS_IGNORE);
            
            int n_VrtMetRecv;
            MPI_Recv(&n_VrtMetRecv,   1, MPI_INT, origin, world_rank*70000, mpi_comm, MPI_STATUS_IGNORE);

            std::vector<double> VrtMetric(n_VrtMetRecv);
            MPI_Recv(&VrtMetric[0],   n_VrtMetRecv, MPI_DOUBLE, origin, world_rank*80000, mpi_comm, MPI_STATUS_IGNORE);

            collected_ElIds[origin]         = recvElIDsvec;
            collected_NIds[origin]             = recvNElVec;
            collected_OriginalNIds[origin]     = recvONElVec;
            collected_FIds[origin]             = recvFElVec;
            collected_Frefs[origin]         = recvFrefElVec;
            collected_FOriginalNIds[origin] = recvFONElVec;
            collected_VrtMetrics[origin]    = VrtMetric;

        }
        
        int el = 0;
        int u  = 0;

        for(ite=tetras.begin();ite!=tetras.end();ite++)
        {
            int gEl       = ite->first;
            //int lEl       = ini_offsetEl[world_rank]+u;
            int* ien      = new int[nv];
            int* ien_o    = new int[nv];
            int* ief      = new int[nv];
            int* ief_o    = new int[nv];
            int* iefref   = new int[nv];

            for(int q=0;q<4;q++)
            {
                gvid = ite->second[q];
                
                if(v2r.find(gvid)!=v2r.end())
                {
                    lvid2 = sharedVmap[gvid];
                    ien[q] = lvid2;
                    ien_o[q] = gvid;
                }
                else
                {
                    if(gl_set.find(gvid)==gl_set.end())
                    {
                        gl_set.insert(gvid);
                        gl_map[gvid]        = lbvid;
                        ien[q]              = lbvid;
                        ien_o[q]            = gvid;
                        lbvid               = lbvid + 1;
                    }
                    else
                    {
                        int lbbvid = gl_map[gvid];
                        ien[q] = lbbvid;
                        ien_o[q] = gvid;
                    }
                }
            }
            
            for(int q=0;q<4;q++)
            {
                gfid       = ief_part_map[gEl][q];

                if(f2r.find(gfid)!=f2r.end())
                {
                    lfid2       = sharedFmap[gfid];
                    ief[q]      = lfid2;
                    iefref[q]   = face2ref[gfid];
                    ief_o[q]    = gfid;
                }
                else
                {
                    if(glf_set.find(gfid)==glf_set.end())
                    {
                        glf_set.insert(gfid);
                        glf_map[gfid]   = lbfid;
                        ief[q]          = lbfid;
                        iefref[q]       = face2ref[gfid];
                        ief_o[q]        = gfid;
                        lbfid           = lbfid + 1;
                    }
                    else
                    {
                        int lbbfid = glf_map[gfid];
                        iefref[q]  = face2ref[gfid];
                        ief_o[q]   = gfid;
                        ief[q]     = lbbfid;
                    }
                }
            }
            
            new_GidEl->setVal(u,0,gEl);
            
            new_ien->setVal(u,0,ien[0]);
            new_ien->setVal(u,1,ien[1]);
            new_ien->setVal(u,2,ien[2]);
            new_ien->setVal(u,3,ien[3]);
            
            new_ien_or->setVal(u,0,ien_o[0]);
            new_ien_or->setVal(u,1,ien_o[1]);
            new_ien_or->setVal(u,2,ien_o[2]);
            new_ien_or->setVal(u,3,ien_o[3]);
            
            new_ief->setVal(u,0,ief[0]);
            new_ief->setVal(u,1,ief[1]);
            new_ief->setVal(u,2,ief[2]);
            new_ief->setVal(u,3,ief[3]);
            
            new_ief_or->setVal(u,0,ief_o[0]);
            new_ief_or->setVal(u,1,ief_o[1]);
            new_ief_or->setVal(u,2,ief_o[2]);
            new_ief_or->setVal(u,3,ief_o[3]);
            
            new_iefref->setVal(u,0,iefref[0]);
            new_iefref->setVal(u,1,iefref[1]);
            new_iefref->setVal(u,2,iefref[2]);
            new_iefref->setVal(u,3,iefref[3]);
            
            for(int s=0;s<4;s++)
            {
                if(m_M_vmap.find(ien[s])==m_M_vmap.end())
                {
                    Array<double>* Mtrc = new Array<double>(6,1);
                    
                    Mtrc->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                    Mtrc->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                    Mtrc->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                    Mtrc->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                    Mtrc->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                    Mtrc->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                    
                    m_M_vmap[ien[s]] = Mtrc;
                    
                }
            }
            
            delete[] ien;
            delete[] ien_o;
            delete[] ief;
            delete[] ief_o;
            delete[] iefref;
            u++;
        }
        
        std::map<int,std::vector<int> >::iterator collit;
        int ntot = nTetras;
        for(collit=collected_NIds.begin();collit!=collected_NIds.end();collit++)
        {
            int nel =  collit->second.size()/4;
            ntot = ntot + nel;
            
            //int fromRank = collit->first;
            for(int q=0;q<nel;q++)
            {
                new_GidEl->setVal(u,0,collected_ElIds[collit->first][q]);
                
                new_ien->setVal(u,0,collit->second[q*4+0]);
                new_ien->setVal(u,1,collit->second[q*4+1]);
                new_ien->setVal(u,2,collit->second[q*4+2]);
                new_ien->setVal(u,3,collit->second[q*4+3]);
                
                for(int s=0;s<4;s++)
                {
                    if(m_M_vmap.find(collit->second[q*4+s])==m_M_vmap.end())
                    {
                        Array<double>* Mtrc = new Array<double>(6,1);
                        Mtrc->setVal(0,0,collected_VrtMetrics[collit->first][6*4*q+6*s+0]);
                        Mtrc->setVal(1,0,collected_VrtMetrics[collit->first][6*4*q+6*s+1]);
                        Mtrc->setVal(2,0,collected_VrtMetrics[collit->first][6*4*q+6*s+2]);
                        Mtrc->setVal(3,0,collected_VrtMetrics[collit->first][6*4*q+6*s+3]);
                        Mtrc->setVal(4,0,collected_VrtMetrics[collit->first][6*4*q+6*s+4]);
                        Mtrc->setVal(5,0,collected_VrtMetrics[collit->first][6*4*q+6*s+5]);
                        nadded++;
                        m_M_vmap[collit->second[q*4+s]] = Mtrc;
                    }
                }
                     
                new_ien_or->setVal(u,0,collected_OriginalNIds[collit->first][4*q+0]);
                new_ien_or->setVal(u,1,collected_OriginalNIds[collit->first][4*q+1]);
                new_ien_or->setVal(u,2,collected_OriginalNIds[collit->first][4*q+2]);
                new_ien_or->setVal(u,3,collected_OriginalNIds[collit->first][4*q+3]);
                
                new_ief->setVal(u,0,collected_FIds[collit->first][4*q+0]);
                new_ief->setVal(u,1,collected_FIds[collit->first][4*q+1]);
                new_ief->setVal(u,2,collected_FIds[collit->first][4*q+2]);
                new_ief->setVal(u,3,collected_FIds[collit->first][4*q+3]);
                
                new_ief_or->setVal(u,0,collected_FOriginalNIds[collit->first][4*q+0]);
                new_ief_or->setVal(u,1,collected_FOriginalNIds[collit->first][4*q+1]);
                new_ief_or->setVal(u,2,collected_FOriginalNIds[collit->first][4*q+2]);
                new_ief_or->setVal(u,3,collected_FOriginalNIds[collit->first][4*q+3]);
                
                new_iefref->setVal(u,0,collected_Frefs[collit->first][4*q+0]);
                new_iefref->setVal(u,1,collected_Frefs[collit->first][4*q+1]);
                new_iefref->setVal(u,2,collected_Frefs[collit->first][4*q+2]);
                new_iefref->setVal(u,3,collected_Frefs[collit->first][4*q+3]);
                                
                u++;
            }
        }
        
        collected_ElIds.clear();
        collected_NIds.clear();
        collected_OriginalNIds.clear();
        collected_FIds.clear();
        collected_FOriginalNIds.clear();
        collected_Frefs.clear();
        collected_VrtMetrics.clear();
    }
    
    delete[] toR_red2;
    delete[] toS_red2;
    delete[] toS_red;
    delete[] toR_red;
    delete[] optiSize;
    delete[] optiSize_red;
    delete[] old_ntets;
    delete[] old_ntets_red;
    delete[] toS;
    delete[] toR;
    delete[] nNonSharedVertsArrayOff;
    delete[] nNonSharedFacesArrayOff;
    delete[] NewGlobVertCountPerRank;
    delete[] NewGlobFaceCountPerRank;
    delete[] Psending;
    delete[] Nsending;
    delete[] Preceiving;
    delete[] Nreceiving;
    delete[] NewGlobVertOffset;
    delete[] NewGlobFaceOffset;
    delete[] ini_nEl;
    delete[] red_ini_nEl;
    delete[] ini_offsetEl;
    
//    delete distSharedFaces;
//    delete distLocalFaces;
//    delete distLocalVerts;
    // return these guys.
    
    ElGids             = new_GidEl;
    ien_part_tetra     = new_ien;
    ien_part_hybrid    = new_ien_or;
    ief_part_tetra     = new_ief;
    ief_part_hybrid    = new_ief_or;
    iefref_part_tetra  = new_iefref;
}




void RedistributePartitionObject::GetNewGlobalPartitioningTetrahedraMesh()
{
    int i;

//    Array<int>* new_ien     = ien_part_tetra;
//    Array<int>* new_ief     = ief_part_tetra;
//    Array<int>* new_ien_or  = ien_part_hybrid;
    
    int optimalSize         = ien_part_tetra->getNrow();
    int* ielement_nlocs     = new int[world_size];
    int* ielement_offsets   = new int[world_size];

    int* red_ielement_nlocs = new int[world_size];
    int* elmdist            = new int[world_size+1];

    for(i=0;i<world_size;i++)
    {
        ielement_nlocs[i]     = 0;
        red_ielement_nlocs[i] = 0;

        if(i==world_rank)
        {
            ielement_nlocs[i] = optimalSize;
        }
        else
        {
            ielement_nlocs[i] = 0;
        }
    }
    
    MPI_Allreduce(ielement_nlocs,
                  &red_ielement_nlocs[0],
                  world_size, MPI_INT,
                  MPI_SUM, mpi_comm);
    
    int o_ie = 0;
    
    for(i=0;i<world_size;i++)
    {
        ielement_offsets[i] = o_ie;
        ielement_nlocs[i]   = red_ielement_nlocs[i];
        elmdist[i]          = o_ie;
        o_ie                = o_ie+red_ielement_nlocs[i];
    }
    
    elmdist[world_size] = o_ie;

    int* eptr     = new int[optimalSize+1];
    int* eind     = new int[optimalSize*4];

    eptr[0]  = 0;
    for(int i=0;i<optimalSize;i++)
    {
        eptr[i+1] = eptr[i]+4;
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j]    = ien_part_tetra->data[j];
        }
    }

    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {3};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    int edgecut      = 0;
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];
    idx_t wgtflag_[] = {2};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.1};
    real_t *ubvec    = ubvec_;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    int* part_arr = new int[optimalSize];
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
    real_t itr_[]    = {1.05};
    real_t *itr = itr_;
    
    
    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    int *elmwgt = new int[optimalSize];
    for(int i=0;i<optimalSize;i++)
    {
        elmwgt[i] = 1;
    }

    
    ParMETIS_V3_Mesh2Dual(elmdist,
                          eptr,
                          eind,
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&mpi_comm);
    
    ParMETIS_V3_PartKway(elmdist,
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &mpi_comm);

    Array<int>* part_new         = new Array<int>(optimalSize,1);

    part_new->data = part_arr;
    
//    MPI_Allgatherv(&part_new->data[0],
//                   red_ielement_nlocs[world_rank], MPI_INT,
//                   &part_global_new->data[0],
//                   ielement_nlocs,
//                   ielement_offsets,
//                   MPI_INT,mpi_comm);
    
    
    m_part          = part_new;
    //m_part_global   = part_global_new;

    
}





void RedistributePartitionObject::UpdateTetrahedraOnPartition(int nglob,
                                                              int nElexpctd,
                                                              ParArray<double>* xcn,
                                                              ParallelState* xcn_pstate,
                                                              ParallelState* ife_pstate)
{
//    m_locV2globV.clear();
//    m_globV2locV.clear();

    Array<int>* ElGids_update             = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_tetra_update     = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_hybrid_update    = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_tetra_update     = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_hybrid_update    = new Array<int>(nElexpctd,4);
    Array<int>* iefref_part_tetra_update  = new Array<int>(nElexpctd,4);
    //Array<int>* loc2globEL_update         = new Array<int>(nElexpctd,1);
    ParallelState* ien_pstate             = new ParallelState(nglob,mpi_comm);

    std::map<int,Array<double>* > M_vmap_copy_update;
    int floc_tmp=0;
    int vloc_tmp=0;
    int q=0;
    int i=0;
    int size = world_size;
    int rank = world_rank;
    

    int el_id;
    int p_id;
    int v_id;
    int v_id_o;
    Vert V;
    std::vector<Vert> part_verts;
    std::vector<std::vector<int> >  part_elem2verts;
    std::map<int,std::vector<int> > rank2req_face;
    std::map<int,std::vector<int> > rank2req_faceOri;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> part_v;
    std::vector<int> loc_r_elem;
    
    std::set<int> unique_vertIDs_on_rank_set;
    std::set<int> unique_faceIDs_on_rank_set;
    
    int r      = 0;
    int lv_id  = 0;
    int lf_id  = 0;
    int f_id   = 0;
    int f_id_o = 0;
    int xcn_o = xcn->getOffset(rank);
    double varia = 0.0;
    int not_on_rank=0;
    int on_rank = 0;
    int* new_V_offsets = new int[size];
    
    for(i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }
    
    int nvPerEl;
    int nfPerEl;
    
    int fref;
    int tellert = 0;
    int sum = 0;
    int sum3 = 0;
    int newElID;
    //int gEL;
//
//    hybV2tetV.clear();
//    tetV2hybV.clear();
    // face2ref.clear();

    
    for(i=0;i<m_part->getNrow();i++)
    {
        p_id        = m_part->getVal(i,0);
        el_id       = ien_pstate->getOffsets()[rank]+i;
        int gEL     = ElGids->getVal(i,0);
       //loc2globEL->setVal(i,0,el_id);
////        m_TetEl2HybEl[el_id] = gEL;
//        nvPerEl = 4;
//        nfPerEl = 4;
//        sum  = 0;
//        sum3 = 0;
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            m_TetraToSendToRanks[p_id].push_back(el_id); // rank to element map.
            m_HybridToSendToRanks[p_id].push_back(gEL);
            
            //====================Hybrid=======================
            for(int k=0;k<4;k++)
            {
                v_id   = ien_part_tetra->getVal(i,k);
                v_id_o = ien_part_hybrid->getVal(i,k);
                
                m_TetraVertIDToSendToRanks[p_id].push_back(v_id);
                m_HybridVertIDToSendToRanks[p_id].push_back(v_id_o);
                sum3=sum3+v_id;
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            
            for(int k=0;k<4;k++)
            {
                f_id   = ief_part_tetra->getVal(i,k);
                f_id_o = ief_part_hybrid->getVal(i,k);
                fref   = iefref_part_tetra->getVal(i,k);
                
                
                m_TetraFaceIDToSendToRanks[p_id].push_back(f_id);
                m_TetraFaceRefToSendToRanks[p_id].push_back(fref);
                m_HybridFaceIDToSendToRanks[p_id].push_back(f_id_o);
            }
            //====================Hybrid=======================
            not_on_rank++;/**/
        }
        
        
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            
            std::vector<int> elem;

            ElGids_update->setVal(on_rank,0,gEL);

            hybE2tetE[gEL]=el_id;
            tetE2hybE[el_id]=gEL;
            
            //hybE2tetE[] =
            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                v_id   = ien_part_tetra->getVal(i,k);
                v_id_o = ien_part_hybrid->getVal(i,k);
                
                ien_part_tetra_update->setVal(on_rank,k,v_id);
                ien_part_hybrid_update->setVal(on_rank,k,v_id_o);
                
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    
                    hybV2tetV[v_id_o]        = v_id;
                    tetV2hybV[v_id]          = v_id_o;
                    M_vmap_copy_update[v_id] = m_M_vmap[v_id];

                    r = FindRank(new_V_offsets,size,v_id_o);
                    
                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        m_TetraRank2ReqVerts[r].push_back(v_id);
                        // add the vertex id that needs to be requested from rank r.
                        m_HybridRank2ReqVerts[r].push_back(v_id_o);
                    }
                    else
                    {
                        m_TetraVertsOnRank.push_back(v_id);
                        // add the vertex to list that is already available on rank.
                        m_HybridVertsOnRank.push_back(v_id_o);

                        vloc_tmp++;
                    }
                    
                    lv_id++;
                    
                }
            }
            
    
            
            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                f_id    = ief_part_tetra->getVal(i,k);
                f_id_o  = ief_part_hybrid->getVal(i,k);
                fref    = iefref_part_tetra->getVal(i,k);
                
                ief_part_tetra_update->setVal(on_rank,k,f_id);
                ief_part_hybrid_update->setVal(on_rank,k,f_id_o);
                iefref_part_tetra_update->setVal(on_rank,k,fref);
                
                if(fref!=2 && m_face2ref.find(f_id)==m_face2ref.end())
                {
                    m_face2ref[f_id] = fref;
                    m_Boundary_Ref2Face[fref].push_back(f_id);

                }
//                //std::cout << fref << " ";
//                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end() && f_id != -1) // add the required unique vertex for current rank.
//                {
//                    unique_faceIDs_on_rank_set.insert(f_id);
//                    //unique_verts_on_rank_vec.push_back(v_id);
//
//                    r = FindRank(new_F_offsets,size,f_id_o);
//
//                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
//                    {
////                        m_TetraRank2ReqFaces[r].push_back(f_id); // add the vertex id that needs to be requested from rank r.
////                        m_HybridRank2ReqFaces[r].push_back(f_id_o);
//
//                    }
//                    else
//                    {
//                        //faceIDs_on_rank.push_back(f_id);  // add the vertex to list that is already available on rank.
////                      hybF2tetF[f_id_o]=f_id;
////                      tetF2hybF[f_id]=f_id_o;
////                        ref2face[fref].push_back(f_id);
////                        face2ref[f_id]=fref;
//                        floc_tmp++;
//                    }
//                    lf_id++;
//                }
            }
            //std::cout << std::endl;
            
            loc_r_elem.push_back(el_id);
            
            on_rank++;
            
        }/**/
        
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    int Ell,Elg;

    //std::cout << "ONRANK + " << rank << " " << on_rank << std::endl;
    
    //============================================================
    //============================================================
    //============================================================
    std::map<int,std::vector<int> >::iterator itt;
    
    int globvid;
    for(itt=m_TetraVertIDToSendToRanks.begin();itt!=m_TetraVertIDToSendToRanks.end();itt++)
    {
        std::vector<double> metrics(6*itt->second.size());

        for(int q=0;q<itt->second.size();q++)
        {
            globvid = itt->second[q];
                        
            metrics[q*6+0] = m_M_vmap[globvid]->getVal(0,0);
            metrics[q*6+1] = m_M_vmap[globvid]->getVal(1,0);
            metrics[q*6+2] = m_M_vmap[globvid]->getVal(2,0);
            metrics[q*6+3] = m_M_vmap[globvid]->getVal(3,0);
            metrics[q*6+4] = m_M_vmap[globvid]->getVal(4,0);
            metrics[q*6+5] = m_M_vmap[globvid]->getVal(5,0);
            
        }

        metricsToSend[itt->first] = metrics;
    }
    
    //============================================================
    //============================================================
    //============================================================
    
    std::set<int> hybFaces;
    std::set<int> tetFaces;
    for(int i=0;i<on_rank;i++)
    {
        for(int j=0;j<4;j++)
        {
            int hybF = ien_part_hybrid_update->getVal(i,j);
            int tetF = ien_part_tetra_update->getVal(i,j);
            
            if(hybFaces.find(hybF)==hybFaces.end())
            {
                hybFaces.insert(hybF);
            }
            
            if(tetFaces.find(tetF)==tetFaces.end())
            {
                tetFaces.insert(tetF);
            }
        }
    }
    
    
    ScheduleObj* part_schedule_elem = DoScheduling(m_TetraToSendToRanks,mpi_comm);
    std::map<int,std::vector<int> >  part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_ov_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_f_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_fref_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_of_map;
    std::map<int,std::vector<int> >  TotRecvElement_GIDs_map;
    std::map<int,std::vector<double> >  TotRecvVertMtrcs_map;
    std::map<int,std::vector<int> >::iterator it;
        
    int n_req_recv;
    int n_req_recv_v;
    int n_req_recv_f;
    int n_req_recv_M;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = m_TetraToSendToRanks.begin(); it != m_TetraToSendToRanks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = m_TetraVertIDToSendToRanks[it->first].size();
                int n_req_f         = m_TetraFaceIDToSendToRanks[it->first].size();
                int n_req_M         = n_req_v*6;
                
                
                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, mpi_comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, mpi_comm);
                MPI_Send(&n_req_f, 1, MPI_INT, dest, dest*222, mpi_comm);
                MPI_Send(&n_req_M, 1, MPI_INT, dest, dest*33300, mpi_comm);
                
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, dest*66666+5555, mpi_comm);
                MPI_Send(&m_TetraVertIDToSendToRanks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, mpi_comm);
                MPI_Send(&m_HybridVertIDToSendToRanks[it->first][0], n_req_v, MPI_INT, dest, 339000+100+dest*2, mpi_comm);
                MPI_Send(&m_TetraFaceIDToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, mpi_comm);
                MPI_Send(&m_TetraFaceRefToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 449000+100+dest*2, mpi_comm);
                MPI_Send(&m_HybridFaceIDToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 889000+100+dest*2, mpi_comm);
                MPI_Send(&m_HybridToSendToRanks[it->first][0], n_req, MPI_INT, dest, dest*7000000, mpi_comm);
                MPI_Send(&metricsToSend[it->first][0], n_req_M, MPI_DOUBLE, dest, dest*8000000, mpi_comm);

                

                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank,    mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111,    mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222,    mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_M, 1, MPI_INT, q, rank*33300,  mpi_comm, MPI_STATUS_IGNORE);

            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_vrt_id(n_req_recv_v);
            std::vector<int>    part_recv_orivrt_id(n_req_recv_v);
            std::vector<int>    part_recv_face_id(n_req_recv_f);
            std::vector<int>    part_recv_face_ref(n_req_recv_f);
            std::vector<int>    part_recv_orifaces_id(n_req_recv_v);
            std::vector<int>    part_recv_Gel_id(n_req_recv);
            std::vector<double>    part_recv_Mtrcs(n_req_recv_M);


            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0], n_req_recv_v, MPI_INT, q, 9000+100+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orivrt_id[0], n_req_recv_v, MPI_INT, q, 339000+100+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_ref[0], n_req_recv_f, MPI_INT, q, 449000+100+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orifaces_id[0], n_req_recv_f, MPI_INT, q, 889000+100+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_Gel_id[0], n_req_recv, MPI_INT, q, rank*7000000, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_Mtrcs[0], n_req_recv_M, MPI_DOUBLE, q, rank*8000000, mpi_comm, MPI_STATUS_IGNORE);

            TotRecvElement_IDs_ov_map[q]    = part_recv_orivrt_id;
            TotRecvElement_IDs_v_map[q]     = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q]     = part_recv_face_id;
            part_tot_recv_elIDs_map[q]      = part_recv_el_id;
            TotRecvElement_IDs_fref_map[q]  = part_recv_face_ref;
            TotRecvElement_IDs_of_map[q]    = part_recv_orifaces_id;
            TotRecvElement_GIDs_map[q]      = part_recv_Gel_id;
            TotRecvVertMtrcs_map[q]         = part_recv_Mtrcs;
        }
    }
    
    std::vector<int> TotRecvElement_GIDs;
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvVerts_IDs;
    std::vector<int> TotRecvOriVerts_IDs;
    std::vector<int> TotRecvFaces_IDs;
    std::vector<int> TotRecvFaces_Refs;
    std::vector<int> TotRecvOriFaces_IDs;
    //std::vector<double> TotRecvVrtMtrcs;

    std::set<int> uovert;
    std::set<int> uvert;
    
    std::map<int,std::vector<int> >::iterator totrecv;
    //unpack the element IDs and their corresponding variable values.
    int TotNelem_recv = 0;
    for(totrecv=part_tot_recv_elIDs_map.begin();totrecv!=part_tot_recv_elIDs_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElement_IDs.push_back(part_tot_recv_elIDs_map[totrecv->first][r]);
            TotRecvElement_GIDs.push_back(TotRecvElement_GIDs_map[totrecv->first][r]);
        }
        TotNelem_recv = TotNelem_recv + totrecv->second.size();
    }
    //unpack the vertex IDs and their corresponding variable values.
    int tt=0;
    int cc=0;
    for(totrecv=TotRecvElement_IDs_v_map.begin();totrecv!=TotRecvElement_IDs_v_map.end();totrecv++)
    {
        for(int l=0;l<totrecv->second.size()/4;l++)
        {
            for(int h=0;h<4;h++)
            {
                TotRecvVerts_IDs.push_back(TotRecvElement_IDs_v_map[totrecv->first][l*4+h]);
                TotRecvOriVerts_IDs.push_back(TotRecvElement_IDs_ov_map[totrecv->first][l*4+h]);
        
                
                if(M_vmap_copy_update.find(TotRecvElement_IDs_v_map[totrecv->first][l*4+h])==M_vmap_copy_update.end())
                {
                    Array<double>* mtrcs = new Array<double>(6,1);
                    
                    mtrcs->setVal(0,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+0]);
                    mtrcs->setVal(1,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+1]);
                    mtrcs->setVal(2,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+2]);
                    mtrcs->setVal(3,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+3]);
                    mtrcs->setVal(4,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+4]);
                    mtrcs->setVal(5,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+5]);
                    
                    M_vmap_copy_update[TotRecvElement_IDs_v_map[totrecv->first][l*4+h]] = mtrcs;
                }
            }
            tt=cc/4;
            cc++;
        }
        
    }
    int addedFace=0;
    //unpack the face IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_f_map.begin();totrecv!=TotRecvElement_IDs_f_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvFaces_IDs.push_back(TotRecvElement_IDs_f_map[totrecv->first][r]);
            TotRecvOriFaces_IDs.push_back(TotRecvElement_IDs_of_map[totrecv->first][r]);
            TotRecvFaces_Refs.push_back(TotRecvElement_IDs_fref_map[totrecv->first][r]);
            
        addedFace++;
    }
    }
    
    TotRecvElement_IDs_ov_map.clear();
    TotRecvElement_IDs_v_map.clear();
    TotRecvElement_IDs_f_map.clear();
    part_tot_recv_elIDs_map.clear();
    TotRecvElement_IDs_fref_map.clear();
    TotRecvElement_IDs_of_map.clear();
    
    int Nel_extra = TotNelem_recv;
    int cnt_v = 0;
    int cnt_f = 0;
    
    int onrank_before = on_rank;
    int sum2 = 0;
    
    
    
    
    
    
    
    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;
        
        int elID  = TotRecvElement_IDs[i];
        int GelID = TotRecvElement_GIDs[i];
        sum2 = 0;
        
        ElGids_update->setVal(on_rank,0,GelID);
        hybE2tetE[GelID]=elID;
        tetE2hybE[elID]=GelID;
        
        for(int k=0;k<4;k++)
        {
            int v_id_n   = TotRecvVerts_IDs[cnt_v+k];
            int v_id_o_n = TotRecvOriVerts_IDs[cnt_v+k];
            
            ien_part_tetra_update->setVal(on_rank,k,v_id_n);
            ien_part_hybrid_update->setVal(on_rank,k,v_id_o_n);
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                hybV2tetV[v_id_o_n]=v_id_n;
                tetV2hybV[v_id_n]=v_id_o_n;
                
                r = FindRank(new_V_offsets, size, v_id_o_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    m_TetraRank2ReqVerts[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                    m_HybridRank2ReqVerts[r].push_back(v_id_o_n);
                }
                else
                {
                    m_TetraVertsOnRank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    m_HybridVertsOnRank.push_back(v_id_o_n);
                
                    vloc_tmp++;
                }
                lv_id++;
            }
        }
        
        for(int k=0;k<4;k++)// looping over the vertices for element "i".
        {
            int f_id_n   = TotRecvFaces_IDs[cnt_f+k];
            int f_id_o_n = TotRecvOriFaces_IDs[cnt_f+k];
            int fref_n   = TotRecvFaces_Refs[cnt_f+k];
            
            ief_part_tetra_update->setVal(on_rank,k,f_id_n);
            ief_part_hybrid_update->setVal(on_rank,k,f_id_o_n);
            iefref_part_tetra_update->setVal(on_rank,k,fref_n);
                 
            if(fref_n!=2 && m_face2ref.find(f_id_n)==m_face2ref.end())
            {
                m_face2ref[f_id_n] = fref_n;
                m_Boundary_Ref2Face[fref_n].push_back(f_id_n);

            }
//            if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
//            {
//                unique_faceIDs_on_rank_set.insert(f_id_n);
//                //unique_verts_on_rank_vec.push_back(v_id);
//
//                r = FindRank(new_F_offsets,size,f_id_o_n);
//
//                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
//                {
////                    m_TetraRank2ReqFaces[r].push_back(f_id_n); // add the vertex id that needs to be requested from rank r.
////                    m_HybridRank2ReqFaces[r].push_back(f_id_o_n);
//                }
//                else
//                {
////                    face2ref[f_id_n] = fref_n;
////                    ref2face[fref_n].push_back(f_id_n);
//                    floc_tmp++;
//                }
//                lf_id++;
//            }
        }
        
        cnt_v=cnt_v+4;
        cnt_f=cnt_f+4;
        
        on_rank++;
        
    }
    
    TotRecvVerts_IDs.clear();
    TotRecvElement_IDs.clear();
    TotRecvOriVerts_IDs.clear();
    TotRecvOriFaces_IDs.clear();
    TotRecvFaces_IDs.clear();
    TotRecvFaces_Refs.clear();
    
    

    //std::cout << "WORLD_RANK - " << rank << " " << on_rank << " " << ief_part_tetra->getNrow() << std::endl;
    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    

    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This m_TetraRank2ReqVerts map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "m_TetraRank2ReqVerts" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    int m = 0;
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* part_schedule = DoScheduling(m_TetraRank2ReqVerts,mpi_comm);
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Ori_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Metrics_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = m_TetraRank2ReqVerts.begin(); it != m_TetraRank2ReqVerts.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                
                //    MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, mpi_comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, mpi_comm);
                //    MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, mpi_comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, mpi_comm);
                MPI_Send(&m_HybridRank2ReqVerts[it->first][0], n_req, MPI_INT, dest, 2229876*2+dest*2, mpi_comm);

                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, mpi_comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_Ori_ids(n_reqstd_ids);
            std::vector<int> recv_metrics(n_reqstd_ids*6);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_reqstd_Ori_ids[0], n_reqstd_ids, MPI_INT, q, 2229876*2+rank*2, mpi_comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
            reqstd_Ori_ids_per_rank[q] = recv_reqstd_Ori_ids;
        }
    }
    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nverts;
    std::map<int,double* > recv_back_verts;
    std::map<int,int* > recv_back_verts_ids;
    
    int n_recv_back;
    
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_Ori_ids_per_rank.begin(); it != reqstd_Ori_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                double* vert_send = new double[nv_send*3];
                offset_xcn        = xcn_pstate->getOffset(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                    vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                    vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);
                }
                
                int dest = it->first;
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, mpi_comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, mpi_comm);
            
                MPI_Send(&vert_send[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, mpi_comm);
                MPI_Send(&reqstd_ids_per_rank[it->first][0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,mpi_comm);
                
                delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, mpi_comm, MPI_STATUS_IGNORE);
            
            double* recv_back_arr = new double[n_recv_back*3];
            int* recv_back_arr_ids = new int[n_recv_back];
            //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, mpi_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr_ids[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, mpi_comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q]     = n_recv_back;
            recv_back_verts[q]      = recv_back_arr;
            recv_back_verts_ids[q]  = recv_back_arr_ids;
        
         }
    }
    
    int vfor = 0;
    std::map<int,double* >::iterator it_f;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int c   =   0;
        vfor    =   vfor+recv_back_Nverts[it_f->first];
    }

    int gvid=0;
    int lvid=0;
    int gvid_gl = 0;
    
    for(m=0;m<m_TetraVertsOnRank.size();m++)
    {
        gvid    = m_TetraVertsOnRank[m];
        gvid_gl = m_HybridVertsOnRank[m];
        
        Vert* V = new Vert;

        V->x = xcn->getVal(gvid_gl-xcn_o,0);
        V->y = xcn->getVal(gvid_gl-xcn_o,1);
        V->z = xcn->getVal(gvid_gl-xcn_o,2);

        LocalVerts.push_back(V);
        
        if(m_locV2globV.find(lvid)==m_locV2globV.end())
        {
            m_locV2globV[lvid] = gvid;
            m_globV2locV[gvid] = lvid;
        }
        else
        {
            std::cout << "huh not unique ?"<<std::endl;
        }
        
        lvid++;
    }
   
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = m_TetraRank2ReqVerts[it_f->first][u];
            
            Vert* V = new Vert;
            
            V->x = it_f->second[u*3+0];
            V->y = it_f->second[u*3+1];
            V->z = it_f->second[u*3+2];
            
            LocalVerts.push_back(V);
            m_locV2globV[lvid]=gvid;
            m_globV2locV[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }
    

    int nLoc_Verts = LocalVerts.size();
    
    
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    reqstd_ids_per_rank.clear();
    reqstd_Ori_ids_per_rank.clear();
    part_elem2verts.clear();
    
    rank2req_face.clear();
    rank2req_faceOri.clear();
    faceIDs_on_rank.clear();
    
    part_v.clear();
    loc_r_elem.clear();
    unique_vertIDs_on_rank_set.clear();
    unique_faceIDs_on_rank_set.clear();
    delete[] new_V_offsets;
    delete ien_pstate;
    m_TetraToSendToRanks.clear();
    m_HybridToSendToRanks.clear();
    m_TetraVertIDToSendToRanks.clear();
    m_HybridVertIDToSendToRanks.clear();
    m_TetraFaceIDToSendToRanks.clear();
    m_TetraFaceRefToSendToRanks.clear();
    m_HybridFaceIDToSendToRanks.clear();
    m_TetraRank2ReqVerts.clear();
    m_HybridRank2ReqVerts.clear();

    metricsToSend.clear();
    
    std::vector<int> m_TetraVertsOnRank;
    std::vector<int> m_HybridVertsOnRank;

    std::map<int,int> m_TetEl2HybEl;
    
    ElGids               = ElGids_update;
    ien_part_tetra       = ien_part_tetra_update;
    ien_part_hybrid      = ien_part_hybrid_update;
    ief_part_tetra       = ief_part_tetra_update;
    ief_part_hybrid      = ief_part_hybrid_update;
    iefref_part_tetra    = iefref_part_tetra_update;
    m_M_vmap             = M_vmap_copy_update;
    
    
        //std::cout << "Testje =  "<< ien_part_tetra->getVal(5125,0) << " " << ien_part_tetra->getVal(5125,1) << " " << ien_part_tetra->getVal(5125,2) << " " << ien_part_tetra->getVal(5125,3) << " on_rank " << on_rank <<" "<<ien_part_tetra->getNrow()<< std::endl;
    
    
    /**/
    
    
    //return tmesh_ret;
}




void RedistributePartitionObject::GetFace2NodeRedistributedMesh(ParArray<int>* ife, int ncol, ParallelState* ife_pstate, int nGlob, MPI_Comm comm)
{
    
    //std::map<int,std::vector<int> > face2element;

    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    ParallelState* ien_pstate      = new ParallelState(nGlob,comm);
    
    int ife_o = ife_pstate->getOffset(rank);
    int face_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_FacesTet;
    std::map<int,std::vector<int> > rank2req_FacesHyb;
    int* new_offsets = new int[size];
    std::map<int,std::vector<int> > ife_tetra;
    
    int gvid;
    int lvid;
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
//    std::map<int,int> hybF2tetF;
//    std::map<int,int> tetF2hybF;
    
    tetF2hybF.clear();
    hybF2tetF.clear();
    std::set<int> ufacesHyb;
    std::set<int> ufacesTet;

    std::vector<int> ee;
    int el_id;
    std::map<int,int> face2ref_v2;
    std::map<int,std::vector<int> >::iterator itefmap;
    std::map<int,std::vector<int> >ife_part_hyb_map;
    
    
    for(int i=0;i<ief_part_hybrid->getNrow();i++)
    {
        el_id   = ien_pstate->getOffsets()[rank]+i;
        
        for(int q=0;q<4;q++)
        {
            int face_hyb = ief_part_hybrid->getVal(i,q);
            int face_tet = ief_part_tetra->getVal(i,q);
            int face_ref = iefref_part_tetra->getVal(i,q);

            //face2element[face_tet].push_back(el_id);
            
            if(ufacesHyb.find(face_hyb)==ufacesHyb.end())
            {
                ufacesHyb.insert(face_hyb);
                tetF2hybF[face_tet] = face_hyb;
                hybF2tetF[face_hyb] = face_tet;
                face2ref_v2[face_tet] = face_ref;
                r = FindRank(new_offsets,size,face_hyb);
                
                if(r != rank)
                {
                    rank2req_FacesHyb[r].push_back(face_hyb);
                    rank2req_FacesTet[r].push_back(face_tet);
                }
                else
                {
                    for(int j=0;j<ncol;j++)
                    {
                        gvid = ife->getVal(face_hyb-ife_o,j);
                        ife_part_hyb_map[face_hyb].push_back(gvid);
                    }
                }
            }
        }
    }
    
    ScheduleObj* ife_schedule = DoScheduling(rank2req_FacesHyb,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_FacesHyb.begin(); it != rank2req_FacesHyb.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nife;
    std::map<int,int* > recv_back_face_ids;
    std::map<int,double* > recv_back_ife;
    std::map<int,int* > recv_back_ifeRef;
    int n_recv_back;
    
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send             = it->second.size();
                double* ife_send        = new double[nf_send*ncol];
                double* ifeRef_send     = new double[nf_send*1];

                int offset_ife          = ife_pstate->getOffset(rank);
                
                for(int u=0;u<it->second.size();u++)
                {
                    int faceID     = it->second[u];
                    ifeRef_send[u] = ife->getVal(it->second[u]-offset_ife,4);
                    
                    for(int h=0;h<ncol;h++)
                    {
                        ife_send[u*ncol+h] = ife->getVal(it->second[u]-offset_ife,h);
                    }
                }

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&ife_send[0], nf_send*ncol, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);
                
                MPI_Send(&ifeRef_send[0], nf_send*1, MPI_INT, dest, 9876*9999+dest*8888, comm);

                delete[] ife_send;
            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_ife_arr   = new double[n_recv_back*ncol];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            int*  recv_back_ifeRef_arr  = new int[n_recv_back];
             
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_ife_arr[0], n_recv_back*ncol, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);
             MPI_Recv(&recv_back_ifeRef_arr[0], n_recv_back*1, MPI_INT, q, 9876*9999+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nife[q]       = n_recv_back;
            recv_back_face_ids[q]   = recv_back_ids_arr;
            recv_back_ife[q]        = recv_back_ife_arr;
            recv_back_ifeRef[q]     = recv_back_ifeRef_arr;


         }
    }
   
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    int face_ref,fhyb,ftet,hvid,tvid;
    for(iter=recv_back_Nife.begin();iter!=recv_back_Nife.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            face_id  = recv_back_face_ids[iter->first][s];
            face_ref = recv_back_ifeRef[iter->first][s];
            
            if(hybF2tetF.find(face_id)!=hybF2tetF.end())
            {
                ftet = hybF2tetF[iter->first];
                face2ref_v2[ftet] = face_ref;
            }
            
            for(int r=0;r<ncol;r++)
            {
                ife_part_hyb_map[face_id].push_back(recv_back_ife[iter->first][s*ncol+r]);
            }
        }
        ntotal=ntotal+L;
    }

    delete[] new_offsets;
    
    
    
    
    std::map<int,std::vector<int> >::iterator itf;
    int tyel   = 0;
    int fref13 = 0;
    
    
    int iref = 0;
    
    for(itf=ife_part_hyb_map.begin();itf!=ife_part_hyb_map.end();itf++)
    {
        fhyb = itf->first;
        
        if(hybF2tetF.find(itf->first)!=hybF2tetF.end())
        {
            ftet = hybF2tetF[itf->first];
            face_ref = face2ref_v2[ftet];
            face2ref_v2[ftet] = face_ref;
            if(face_ref == 13 && m_shell_tet2hyb_loc.find(ftet)==m_shell_tet2hyb_loc.end())
            {
                m_shell_tet2hyb_loc[ftet] = fhyb;
                //std::set<int> unF;
                std::vector<int> shF(3);
                for(int j=0;j<ncol;j++)
                {
                    hvid    = itf->second[j];
                    tvid    = hybV2tetV[hvid];
                    shF[j]  = tvid;
                    if(m_shellVert2RefMapLocal.find(tvid)==m_shellVert2RefMapLocal.end())
                    {
                        m_shellTagID2TetVID[tvid] = hvid;
                        m_shellVert2RefMapLocal[tvid]=iref;
                        m_VertTag2RefLocal[hvid]=iref;
                        iref++;
                    }
                }
                
                m_ShellFace2Vert[fhyb] = shF;
              
                fref13++;
            }

        }
        else
        {
            std::cout << "OOk geen feesie " << fhyb << " ON RANK " << rank << std::endl;
        }
        int* nodes = new int[ncol];
        
        for(int j=0;j<ncol;j++)
        {
            hvid = itf->second[j];

            if(hybV2tetV.find(hvid)!=hybV2tetV.end())
            {
                tvid = hybV2tetV[hvid];
            }
            else
            {
                std::cout << "Errror hvid doesnt exist " << hvid << " " << " This error occurs for face " << ftet << " " << fhyb << " on rank " << rank  << " " << hybV2tetV.size() << " " << j << std::endl;
            }
            nodes[j]     = tvid;
        }
        
        face2node[ftet] = nodes;
    }
    
    
    delete ien_pstate;
    //std::cout << "m_shellVert2RefMapLocal.size() " << world_rank << " " << m_shellVert2RefMapLocal.size() << std::endl;
    
    
    //std::cout << "alt version " << world_rank << " " << hybV2tetV.size() << std::endl;
    //std::cout << rank << " ==> " << face2ref_v2.size() << " -- " << m_face2ref.size() << std::endl;
    /**/
    
}


void RedistributePartitionObject::GetFace2RankTetrahedraMesh()
{
    std::map<int,int*>::iterator itst;
    int element, fid;
    DistributedParallelState* newElemDist = new DistributedParallelState(ief_part_tetra->getNrow(),mpi_comm);
    int el_offset = newElemDist->getOffsets()[world_rank];
    
    for(int i=0;i<ief_part_tetra->getNrow();i++)
    {
        element = i + el_offset;
        
        for(int j=0;j<4;j++)
        {
            fid = ief_part_tetra->getVal(i,j);
            face2element[fid].push_back(element);
        }
    }
    
    std::map<int,std::vector<int> >::iterator itf;
    std::vector<int> sharedFonRank;
    std::vector<int> interiorFonRank;
    std::vector<Vert*> locVs = LocalVerts;

    for(itf=face2element.begin();itf!=face2element.end();itf++)
    {
        if(itf->second.size()==1)
        {
            sharedFonRank.push_back(itf->first);
        }
        if(itf->second.size()==2)
        {
            interiorFonRank.push_back(itf->first);
        }
    }
    
    int nSharedFonRank = sharedFonRank.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFonRank,mpi_comm);

    int Nt_shFaces               = distSharedFaces->getNel();
    int* shFace_offsets          = distSharedFaces->getOffsets();
    int* shFace_nlocs            = distSharedFaces->getNlocs();
    int* shFacesIDs              = new int[nSharedFonRank];
    int* shFaces_RankIDs         = new int[nSharedFonRank];

    for(int i=0;i<nSharedFonRank;i++)
    {
        shFacesIDs[i]      = sharedFonRank[i];
        shFaces_RankIDs[i] = world_rank;
    }
    
    int* TotalSharedFaces        = new int[Nt_shFaces];
    int* TotalSharedFaces_RankID = new int[Nt_shFaces];
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(shFacesIDs,
                   nSharedFonRank,
                   MPI_INT,
                   TotalSharedFaces,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, mpi_comm);
    
    MPI_Allgatherv(shFaces_RankIDs,
                   nSharedFonRank,
                   MPI_INT,
                   TotalSharedFaces_RankID,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, mpi_comm);
    
    std::map<int,std::vector<int> > face2rank;
    
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];
        face2rank[key].push_back(val);
    }
    
    delete[] TotalSharedFaces;
    delete[] TotalSharedFaces_RankID;
    delete[] shFace_offsets;
    delete[] shFace_nlocs;
    delete[] shFacesIDs;
    delete[] shFaces_RankIDs;
    //delete distSharedFaces;
    
    std::vector<int> faces4parmmg;
    std::map<int,std::vector<int> >::iterator itff;
    int shf = 0;
    int bf  = 0;
    std::map<int,std::vector<int> > Boundary_Ref2Face;
    std::set<int> uSharedVert;
    std::set<int> uBoundVert;
    std::vector<std::vector<int> > pltSfaces;
    std::vector<std::vector<int> > pltBfaces;
    std::map<int,int> g2l_s;
    std::map<int,int> g2l_b;
    std::map<int,int> b2l_g;
    int Bvid=0,lBvid;
    int Svid=0,lSvid;
    std::vector<Vert*> svs;
    std::vector<Vert*> bvs;
    std::map<int,int> v2ref;
    int ref;
    int f = 0;
    std::map<int,int> face2ref;
  
    for(itff=face2rank.begin();itff!=face2rank.end();itff++)
    {
        int faceID = itff->first;
        
        if(itff->second.size()==2)
        {
            if(itff->second[0]==world_rank)
            {
                m_ColorsFaces[itff->second[1]].push_back(faceID);
            }
            else if(itff->second[1]==world_rank)
            {
                m_ColorsFaces[itff->second[0]].push_back(faceID);
            }
            
            std::vector<int> sface(3);
            if(face2node.find(faceID)!=face2node.end())
            {
                m_faces4parmmg.push_back(faceID);

                for(int u=0;u<3;u++)
                {
                    int vid = face2node[faceID][u];

                    if(uSharedVert.find(vid)==uSharedVert.end())
                    {
                        uSharedVert.insert(vid);
                        g2l_s[vid] = Svid;

                        lSvid = m_globV2locV[vid];

                        Vert* sv = new Vert;
                        sv->x = locVs[lSvid]->x;
                        sv->y = locVs[lSvid]->y;
                        sv->z = locVs[lSvid]->z;

                        sface[u]=Svid;
                        svs.push_back(sv);
                        Svid++;
                    }
                    else
                    {
                        int Svid2 = g2l_s[vid];
                        sface[u]=Svid2;
                    }
                }
                pltSfaces.push_back(sface);
                m_globShF2locShF[faceID] = f;
                m_locShF2globShF[f] = faceID;
                f++;
                
                shf++;
            }
            
        }
        
        
        if(itff->second.size()==1)
        {
            std::vector<int> bface(3);
            
            if(face2node.find(faceID)!=face2node.end())
            {
                m_faces4parmmg.push_back(faceID);

                if(face2ref.find(faceID)!=face2ref.end())
                {
                    ref = face2ref[faceID];
                    Boundary_Ref2Face[ref].push_back(faceID);
                }
                
                for(int u=0;u<3;u++)
                {
                    int vid = face2node[faceID][u];
                    
                    if(uBoundVert.find(vid)==uBoundVert.end())
                    {
                        v2ref[vid]   = ref;
                        int hybv     = tetV2hybV[vid];
                        uBoundVert.insert(vid);
                        g2l_b[vid]   = Bvid;
                        b2l_g[Bvid]  = vid;
                        lBvid        = m_globV2locV[vid];
                        Vert* bv     = new Vert;
                        bv->x        = locVs[lBvid]->x;
                        bv->y        = locVs[lBvid]->y;
                        bv->z        = locVs[lBvid]->z;
                        bface[u]     = Bvid;
                        
                        bvs.push_back(bv);
                        Bvid++;
                    }
                    else
                    {
                        int Bvid2    = g2l_b[vid];
                        bface[u]     = Bvid2;
                    }
                }
                
                pltBfaces.push_back(bface);
                m_globShF2locShF[faceID] = f;
                m_locShF2globShF[f] = faceID;
                bf++;
                f++;
                
            }
        }
    }
    
    //std::cout << "face2rank.size  " << face2rank.size() << std::endl;
    m_nPartBoundFaces = m_faces4parmmg.size();
    m_ncomm           = m_ColorsFaces.size();
    m_color_face = (int *) malloc(m_ncomm*sizeof(int));
    m_ntifc = (int *) malloc(m_ncomm*sizeof(int));
    
    m_ifc_tria_loc = (int **)malloc(m_ncomm*sizeof(int *));
    m_ifc_tria_glo = (int **)malloc(m_ncomm*sizeof(int *));

//    int* ifc_tria_glob[m_ncomm];
//    int* ifc_tria_loc[m_ncomm];
    int icomm=0;
    std::map<int,std::vector<int> >::iterator itc;
    for(itc=m_ColorsFaces.begin();itc!=m_ColorsFaces.end();itc++)
    {
        m_color_face[icomm]     = itc->first;
        m_ntifc[icomm]          = itc->second.size();
        m_ifc_tria_loc[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));
        m_ifc_tria_glo[icomm]  = (int *) malloc(itc->second.size()*sizeof(int));
        
        for(int q=0;q<itc->second.size();q++)
        {
            m_ifc_tria_glo[icomm][q] = itc->second[q]+1;
            m_ifc_tria_loc[icomm][q]  = m_globShF2locShF[itc->second[q]]+1;
        }
        icomm++;
    }
    
//    m_ifc_tria_loc  = ifc_tria_loc;
//    m_ifc_tria_glo  = ifc_tria_glob;
    
    
//    if(debug == 1)
//    {
//        //================================================================
//        std::string filename = "SharedPartFaces_" + std::to_string(world_rank) + ".dat";
//        std::ofstream myfile;
//        myfile.open(filename);
//        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//        myfile <<"ZONE N = " << svs.size() << ", E = " << pltSfaces.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
//
//        for(int i=0;i<svs.size();i++)
//        {
//            myfile << svs[i]->x << " " << svs[i]->y << " " << svs[i]->z << std::endl;
//        }
//        int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
//        int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
//        for(int i=0;i<pltSfaces.size();i++)
//        {
//            myfile <<   pltSfaces[i][0]+1 << "  " <<
//                        pltSfaces[i][1]+1 << "  " <<
//                        pltSfaces[i][2]+1 << std::endl;
//        }
//
//
//        myfile.close();
//
//
//        //================================================================
//
//        int nBfaces = pltBfaces.size();
//
//
//        //================================================================
//        std::string filename2 = "BoundaryPartFaces_" + std::to_string(world_rank) + ".dat";
//        std::ofstream myfile2;
//        myfile2.open(filename2);
//        myfile2 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//        myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//        myfile2 <<"ZONE N = " << bvs.size() << ", E = " << nBfaces << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
//        int evv = -1;
//        for(int i=0;i<bvs.size();i++)
//        {
//            myfile2 << bvs[i]->x << " " << bvs[i]->y << " " << bvs[i]->z << std::endl;
//        }
//
//
//        for(int i=0;i<nBfaces;i++)
//        {
//            myfile2 <<   pltBfaces[i][0]+1 << "  " <<
//                         pltBfaces[i][1]+1 << "  " <<
//                         pltBfaces[i][2]+1 << std::endl;
//
//        }
//
//        myfile2.close();
//    }
    
    
    
//    pb->faces4parmgg;
//    pb->globShF2locShF;
//    pb->ifc_tria_glob = ifc_tria_glob;
//    pb->ifc_tria_loc = ifc_tria_loc;
//    pb->m_color_face;
//    pb->ntifc;
}


std::map<int,Array<double>* > RedistributePartitionObject::GetVert2MetricMap()
{
    return m_M_vmap;
}
Array<int>* RedistributePartitionObject::GetElement2NodeMap()
{
    return ien_part_tetra;
}
std::map<int,int* > RedistributePartitionObject::GetFace2NodeMap()
{
    return face2node;
}
std::map<int,std::vector<int> > RedistributePartitionObject::GetFace2ElementMap()
{
    return face2element;
}
std::vector<Vert*> RedistributePartitionObject::GetLocalVertices()
{
    return LocalVerts;
}
int** RedistributePartitionObject::GetFace2LocalNode()
{
    return m_ifc_tria_loc;
}
int** RedistributePartitionObject::GetFace2GlobalNode()
{
    return m_ifc_tria_glo;
}
int RedistributePartitionObject::GetNBoundaryFaces()
{
    return m_nPartBoundFaces;
}
std::vector<int> RedistributePartitionObject::GetFaces4ParMMG()
{
    return m_faces4parmmg;
}
std::map<int,int> RedistributePartitionObject::GetGlobalVert2LocalVertMap()
{
    return m_globV2locV;
}
std::map<int,int> RedistributePartitionObject::GetLocalVert2GlobalVertMap()
{
    return m_locV2globV;
}
int* RedistributePartitionObject::GetColorFace()
{
    return m_color_face;
}
int* RedistributePartitionObject::GetNFacesPerColor()
{
    return m_ntifc;
}
int RedistributePartitionObject::GetNcomm()
{
    return m_ncomm;
}
std::map<int,int> RedistributePartitionObject::GetLocalSharedFace2GlobalSharedFace()
{
    return m_locShF2globShF;
}
std::map<int,int> RedistributePartitionObject::GetFace2RefMap()
{
    return m_face2ref;
}
std::map<int,std::vector<int> > RedistributePartitionObject::GetBndRef2FaceMap()
{
    return m_Boundary_Ref2Face;
}
std::map<int,int> RedistributePartitionObject::GetTetF2HybFMap()
{
    return tetF2hybF;
}
std::map<int,int> RedistributePartitionObject::GetShellTet2HybFaceMap()
{
    int mapSizeLoc = m_shell_tet2hyb_loc.size();
    std::map<int,int>::iterator it;
    
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);

    int mapSizeTot = distrimap->getNel();
    int* ftet_loc = new int[mapSizeLoc];
    int* fhyb_loc = new int[mapSizeLoc];
    int* ftet_tot = new int[mapSizeTot];
    int* fhyb_tot = new int[mapSizeTot];
    int ftet,fhyb;
    int i = 0;
    for(it=m_shell_tet2hyb_loc.begin();it!=m_shell_tet2hyb_loc.end();it++)
    {
        ftet_loc[i] = it->first;
        fhyb_loc[i] = it->second;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    MPI_Allgatherv(ftet_loc,
                   mapSizeLoc,
                   MPI_INT,
                   ftet_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    
    MPI_Allgatherv(fhyb_loc,
                   mapSizeLoc,
                   MPI_INT,
                   fhyb_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    std::map<int,int> m_shell_tet2hyb;
    
    for(int i=0;i<mapSizeTot;i++)
    {
        if(m_shell_tet2hyb.find(ftet_tot[i])==m_shell_tet2hyb.end())
        {
            m_shell_tet2hyb[ftet_tot[i]] = fhyb_tot[i];
        }
    }
    
    return m_shell_tet2hyb;
    
    
}





std::map<int,int > RedistributePartitionObject::GetShellVert2RefMap()
{
    return m_shellVert2RefMapLocal;
}


std::map<int,int> RedistributePartitionObject::GetVertTag2RefLocalMap()
{
    return m_VertTag2RefLocal;
}

std::map<int,int > RedistributePartitionObject::GetShellVert2RefMap_Local()
{
    return m_shellVert2RefMapLocal;
}

std::map<int,int > RedistributePartitionObject::GetShellVert2RefMap_Global()
{
    
    int mapSizeLoc = m_shellVert2RefMapLocal.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();
    int* hvid_loc = new int[mapSizeLoc];
    int* tvid_loc = new int[mapSizeLoc];
    int* iref_loc = new int[mapSizeLoc];
    
    int* hvid_tot = new int[mapSizeTot];
    int* tvid_tot = new int[mapSizeTot];
    int* iref_tot = new int[mapSizeTot];
    int tvid_tmp,iref_tmp;
    int i = 0;
    
    //std::cout << "m_shellVert2RefMapLocal  " << world_rank << " " << m_shellVert2RefMapLocal.size() << std::endl;
    
    std::map<int,int>::iterator itred;
    
    for(itred=m_shellVert2RefMapLocal.begin();itred!=m_shellVert2RefMapLocal.end();itred++)
    {
        tvid_loc[i] = itred->first;
        iref_loc[i] = itred->second;
        hvid_loc[i] = m_shellTagID2TetVID[itred->first];
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    MPI_Allgatherv(hvid_loc,
                   mapSizeLoc,
                   MPI_INT,
                   hvid_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    MPI_Allgatherv(tvid_loc,
                   mapSizeLoc,
                   MPI_INT,
                   tvid_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    
    MPI_Allgatherv(iref_loc,
                   mapSizeLoc,
                   MPI_INT,
                   iref_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    int tmp;
    int irefn = 100;
    int refcount = 0;
    int key,key2;
    for(int i=0;i<mapSizeTot;i++)
    {
        key  = tvid_tot[i];
        key2 = hvid_tot[i];
        if(m_shellVert2RefMapGlobal.find(key)==m_shellVert2RefMapGlobal.end())
        {
            m_shellVert2RefMapGlobal[key]     = irefn;
            m_shellVertTag2RefMapGlobal[key2] = irefn;
            irefn++;
            refcount++;
        }
    }
    
    //std::cout << " m_shellVert2RefMapGlobal.size() " << world_rank << " " << m_shellVert2RefMapGlobal.size() << " refcount--> " << refcount << " " << irefn << " mapSizeTot " << mapSizeTot << std::endl;
    
    
    return m_shellVert2RefMapGlobal;
}

std::map<int,int > RedistributePartitionObject::GetShellVertTag2RefMap_Global()
{
    return m_shellVertTag2RefMapGlobal;
}

std::map<int,int > RedistributePartitionObject::GetTag2TetVertMap()
{
	return hybV2tetV;
}

std::map<int,int > RedistributePartitionObject::GetTet2TagVertMap()
{
	return tetV2hybV;
}


std::map<int,std::set<int> > RedistributePartitionObject::GetShellFace2VertRefMap()
{
    std::map<int,std::vector<int> >::iterator itm;
    
    int mapSizeLoc = m_ShellFace2Vert.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    DistributedParallelState* distriFacemap = new DistributedParallelState(mapSizeLoc*3,mpi_comm);

    int mapSizeTot = distrimap->getNel();
    
    int* fhyb_loc = new int[mapSizeLoc];
    int* triv_loc = new int[mapSizeLoc*3];
    int* fhyb_tot = new int[mapSizeTot];
    int* triv_tot = new int[mapSizeTot*3];
    
    int tvid_tmp,iref_tmp;
    int i = 0;
    
    std::map<int,std::vector<int> >::iterator itred;
    for(itred=m_ShellFace2Vert.begin();itred!=m_ShellFace2Vert.end();itred++)
    {
        fhyb_loc[i] = itred->first;
        
        for(int q=0;q<3;q++)
        {
            triv_loc[i*3+q] = itred->second[q];
        }
        
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    
    MPI_Allgatherv(fhyb_loc,
                   mapSizeLoc,
                   MPI_INT,
                   fhyb_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    int* offsetsFace = distriFacemap->getOffsets();
    int* nlocsFace   = distriFacemap->getNlocs();
    
    MPI_Allgatherv(triv_loc,
                   mapSizeLoc*3,
                   MPI_INT,
                   triv_tot,
                   nlocsFace,
                   offsetsFace,
                   MPI_INT, mpi_comm);
    
    std::map<int,std::vector<int> > m_ShellFace2VertTot;
    for(int i=0;i<mapSizeTot;i++)
    {
        int key = fhyb_tot[i];
        std::vector<int> val(3);
        for(int j=0;j<3;j++)
        {
            val[j] = triv_tot[i*3+j];
        }
        
        m_ShellFace2VertTot[key] = val;
    }

    for(itm = m_ShellFace2VertTot.begin();itm!=m_ShellFace2VertTot.end();itm++)
    {
        int fhyb = itm->first;
        std::set<int> refUpdate;
        for(int i=0;i<itm->second.size();i++)
        {
            int vid  = itm->second[i];
            int iref = m_shellVert2RefMapGlobal[vid];
            refUpdate.insert(iref);
        }
        
        
        m_ShellFace2VertRef[fhyb] = refUpdate;
        
        refUpdate.clear();
    }

    //std::cout << " m_ShellFace2VertRef .size " << m_ShellFace2VertRef.size() << std::endl;
    
    return m_ShellFace2VertRef;
}




