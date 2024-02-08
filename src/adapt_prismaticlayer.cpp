#include "adapt_prismaticlayer.h"



    
PrismaticLayer::PrismaticLayer(std::map<int,std::vector<int> > elements,
                               std::map<int,std::vector<int> > ief_part_map,
                               std::map<int,std::vector<int> > ifn_part_map,
                               std::map<int,std::vector<int> > ife_part_map,
                               std::map<int,int> if_ref_part_map,
                               std::map<int,int> if_Nv_part_map,
							   std::map<int,std::vector<int> > if_Erank_part_map,
                               std::map<int,std::vector<int> > ushell,
                               std::map<int,int> tag2locV,
                               std::vector<std::vector<double> > locVerts,
                               std::map<int,int> shellvert2ref_glob,
							   std::map<int,std::vector<int> > ranges_id, 
                               MPI_Comm comm )
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::map<int,std::vector<int> >::iterator prit;
    std::map<int,int> SharedVertsOwned;
    std::map<int,int> NonSharedVertsOwned;
    
    int lf  = 0;
    int shf = 0;

    std::map<int,int> sharedFaces;
    std::set<int> gfid_set;
    int ref = 0;
    int fzone = -1;
    int nLocalFaces = 0;
    int nLocalVerts = 0;
    std::set<int> ufaces;
    std::set<int> gvid_set;
    int gvid,gfid,r0,r1,el0,el1;
    int elTel = 0;
    std::map<int,int> loc2globV;
    int foundshellvert = 0;
    std::map<int,int> unf;

    std::map<int,int* > ifn_updt;
    int ell = 0;

    for(prit=elements.begin();prit!=elements.end();prit++)
    {
        //key = GID, value global node nmber;
        int gEl = prit->first;
        
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
            
            int nppf = if_Nv_part_map[gfid];
            for(int h=0;h<nppf;h++)
            {
                gvid = ifn_part_map[gfid][h];
                
                if(gvid_set.find(gvid)==gvid_set.end() &&
                   shellvert2ref_glob.find(gvid)==shellvert2ref_glob.end() )
                {
                    gvid_set.insert(gvid);
                }
                
            }
            
            if(gfid_set.find(gfid)==gfid_set.end())
            {
                gfid_set.insert(gfid);
                nLocalFaces++;
            }
            
            if(ufaces.find(gfid)==ufaces.end())
            {
                ufaces.insert(gfid);

//                el0    = ife_part_map[gfid][0];
//                el1    = ife_part_map[gfid][1];

                r0 	   = if_Erank_part_map[gfid][0];
                r1 	   = if_Erank_part_map[gfid][1];
                
                if(ref==2)
                {
//                    r0 = part_global->getVal(el0,0);
//                    r1 = part_global->getVal(el1,0);

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
    
    
    nLocalVerts = gvid_set.size();
    
    int nSharedFaces                          = sharedFaces.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
    DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,comm);

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
        int nppf = if_Nv_part_map[gfid];
        for(int q=0;q<nppf;q++)
        {
            gvid   = ifn_part_map[gfid][q];
            
            if(UniqueSharedVertsOnRank_set.find(gvid)==UniqueSharedVertsOnRank_set.end() &&
               shellvert2ref_glob.find(gvid)==shellvert2ref_glob.end() )
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

    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
    
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
                   MPI_INT, comm);
    
    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
                   nSharedVerts,
                   MPI_INT,
                   TotalSharedVerts_RankID,
                   shVerts_nlocs,
                   shVerts_offsets,
                   MPI_INT, comm);
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(shFacesIDs,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFaces,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    MPI_Allgatherv(shFaces_RankIDs,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFaces_RankID,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    int tmp;
    std::map<int,int> f2r;
    
        
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];
        
        if(f2r.find(key)==f2r.end())
        {
            f2r[key]=val;
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
    
    std::map<int,int>::iterator itmm;
    int nOwnedSharedVerts = 0;
    for(itmm=v2r.begin();itmm!=v2r.end();itmm++)
    {
        if(itmm->second==world_rank)
        {
            nOwnedSharedVerts++;
        }
    }

    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,int> sharedVertsGlobal;
    
    
    
    int nNonSharedVerts          = nLocalVerts-nSharedVerts;
    int nNonSharedFaces          = nLocalFaces-nSharedFaces;

    DistributedParallelState* nonSharedVertDistr = new DistributedParallelState(nNonSharedVerts,comm);
    DistributedParallelState* nonSharedFaceDistr = new DistributedParallelState(nNonSharedFaces,comm);

    
    int nNonSharedVertsTot = nonSharedVertDistr->getNel();
    
    std::map<int,int >::iterator itvv;

    DistributedParallelState* ownedSharedVrtsDist = new DistributedParallelState(nOwnedSharedVerts,comm);
    
    int lshaVrt = 0;
    int tshaVrt = 0;

    
    int* ownedVs = new int[world_size];
    for(int i=0;i<world_size;i++)
    {
        ownedVs[i] = 0;
    }
    
    for(itmm=v2r.begin();itmm!=v2r.end();itmm++)
    {
        int globid = nNonSharedVertsTot+ownedSharedVrtsDist->getOffsets()[itmm->second]+ownedVs[itmm->second]+1;
        
        sharedVmap[itmm->first] = globid;
        
        if(itmm->second==world_rank &&
           shellvert2ref_glob.find(itmm->first)==shellvert2ref_glob.end())
        {
            SharedVertsOwned[itmm->first] = globid;
        }
        
        ownedVs[itmm->second] = ownedVs[itmm->second]+1;
    }

    // std::map<int,int> sharedFmap;
    // int iFshared = distLocalFaces->getNel()-f2r.size();

    // for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
    // {
    //     sharedFmap[itvv->first] = iFshared;
    //     iFshared++;
    // }
    
    DistributedParallelState* ElementDistr = new DistributedParallelState(elements.size(),comm);
    
    int uloc = 0;
    int u    = 0;
    
    std::map<int,std::vector<int> > colorRh;
    std::map<int,std::vector<int> > colorFh;

    iet = new Array<int>(elements.size(),1);
//    std::map<int,std::vector<int> > shellFace2Node_prism;
//    std::map<int,std::vector<int> > BoundaryFace2Node;
//    std::map<int,std::vector<int> > SharedFace2Node;
//    std::map<int,std::vector<int> > InternalFace2Node;
//    std::map<int,int> tagE2gE;
//    std::map<int,int> gE2tagE;
//    std::map<int,int> pface2ref;
//
    std::map<int,int> tag2ElementID;
    
    int fowned = 0;
    int fownedInmap = 0;
    double orient0 = 1.0;
    std::map<int,double> orient_map;
    int ngf = 0;
    std::set<int> facemap;

    
    int notowned = 0;
    std::cout << "gvid " << gvid << std::endl;
    for(prit=elements.begin();prit!=elements.end();prit++)
    {
        int gEl       = prit->first;
        int rank      = world_rank;
        int lEl       = ElementDistr->getOffsets()[rank]+u+1;
        int nVrt	  = prit->second.size();
        iet->setVal(u,0,nVrt);
        
        tagE2gE[gEl]  = lEl;
        gE2tagE[lEl]  = gEl;
       
        int nFaces = 5; // for both pyramids and prisms;
        
        for(int q=0;q<nFaces;q++)
        {
        	int fzone = 3;
            gfid = ief_part_map[gEl][q];
            int nfvrts = if_Nv_part_map[gfid];
            // Setting the reference to 13 for the shell faces;
            if(ushell.find(gfid)!=ushell.end())
            {
                ref                 = 13;
                tag2ElementID[gEl]  = lEl;
            }
            else
            {
                ref                 = if_ref_part_map[gfid];
                
                
                //std::cout<< "fzone " << world_rank << " " << fzone << std::endl;
                
            }
            
            if(facemap.find(gfid)==facemap.end())
            {
                facemap.insert(gfid);
                
                if(f2r.find(gfid)!=f2r.end() && ref != 13)
                {
                    if(f2r[gfid]==world_rank)
                    {
                        fowned++;
                        
                        if(SharedFace2Node.find(gfid)==SharedFace2Node.end())
                        {
                            fownedInmap++;
                            
                            el0                     = ife_part_map[gfid][0];
                            el1                     = ife_part_map[gfid][1];
                            
                            r0 	   = if_Erank_part_map[gfid][0];
                            r1 	   = if_Erank_part_map[gfid][1];
                            
//                            r0                      = part_global->getVal(el0,0);
//                            r1                      = part_global->getVal(el1,0);
                            
                            if(r0==world_rank && r1!=world_rank)
                            {
                                colorRh[r1].push_back(el1);
                                colorFh[r1].push_back(gfid);
                            }
                            else if(r1==world_rank && r0!=world_rank)
                            {
                                colorRh[r0].push_back(el0);
                                colorFh[r0].push_back(gfid);
                            }
                            
                            std::vector<int> fn_tag(nfvrts);
                            for(int n=0;n<nfvrts;n++)
                            {
                                fn_tag[n] = ifn_part_map[gfid][n];
                                
                                if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() 
                                		&& shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                {
                                    if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end())
                                    {
                                        SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
                                    }
                                }
                            }
                            
                            SharedFace2Node[gfid] = fn_tag;
                            lhp[gfid]             = lEl;

                        }
                    }
                }
                else
                {
                    if(ref == 2)
                    {
                        if(InternalFace2Node.find(gfid)==InternalFace2Node.end())
                        {
                            int nfvrts         = if_Nv_part_map[gfid];
                            
                            std::vector<int> fn_tag(nfvrts);

                            for(int n=0;n<nfvrts;n++)
                            {
                                fn_tag[n] = ifn_part_map[gfid][n];
                                if(v2r.find(fn_tag[n])!=v2r.end())
                                {
                                    if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                    {
                                        if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end())
                                        {
                                            SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
                                        }
                                    }
                                }
                                else
                                {
                                    if(NonSharedVertsOwned.find(fn_tag[n])==NonSharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                    {
                                        NonSharedVertsOwned[fn_tag[n]]  = gvid;
                                        gvid++;
                                    }
                                }
                            }
                            
                           
                            InternalFace2Node[gfid] = fn_tag;
                            lhp[gfid]          = lEl;
                        }
                    }
                    
                    if(ref != 2 && ref!=13)
                    {
                        ref2bcface[ref].push_back(gfid);
                        
                        fzone = ProvideBoundaryID(gfid,ranges_id);
                        
                        zone2bcface[fzone].push_back(gfid);
                        
                        if(BoundaryFace2Node.find(gfid)==BoundaryFace2Node.end())
                        {
                            int nfvrts         = if_Nv_part_map[gfid];
                            std::vector<int> fn_tag(nfvrts);
                            for(int n=0;n<nfvrts;n++)
                            {
                                fn_tag[n] = ifn_part_map[gfid][n];

                                if(v2r.find(fn_tag[n])!=v2r.end())
                                {
                                    if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                    {
                                        if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end())
                                        {
                                            SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
                                        }
                                    }
                                }
                                else
                                {
                                    if(NonSharedVertsOwned.find(fn_tag[n])==NonSharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                    {
                                        NonSharedVertsOwned[fn_tag[n]]  = gvid;
                                        gvid++;
                                    }
                                }
                            }
                            

                            
                            BoundaryFace2Node[gfid]     = fn_tag;
                            lhp[gfid]                   = lEl;
                        }
                    }
                }
            }
            else
            {
                rhp[gfid]           = lEl;
            }
        }
        
        uloc++;
        u++;
    
    }
    
    
    ScheduleObj* rh_schedule = DoScheduling(colorRh,comm);
    std::map<int,std::vector<int> > recv_rhElIDs;
    for(int q=0;q<world_size;q++)
    {
        if(world_rank==q)
        {
            int i=0;
            
            for (itv = colorRh.begin(); itv != colorRh.end(); itv++)
            {
                int n_req           = itv->second.size();
                int dest            = itv->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
                MPI_Send(&itv->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
                i++;
            }
        }
        else if (rh_schedule->SendFromRank2Rank[q].find( world_rank ) != rh_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
            recv_rhElIDs[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> > sendEl;
    std::map<int,std::vector<int> >::iterator rcvit;
    for(rcvit=recv_rhElIDs.begin();rcvit!=recv_rhElIDs.end();rcvit++)
    {
        int frank = rcvit->first;
        int nE    = rcvit->second.size();
        
        for(int j=0;j<nE;j++)
        {
            int gEl = tagE2gE[rcvit->second[j]];
            sendEl[frank].push_back(gEl);
        }
    }
    
    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);

    std::map<int,std::vector<int> > adj_ids;
    for(int q=0;q<world_size;q++)
    {
        if(world_rank==q)
        {
            int i=0;
            for (itv = sendEl.begin(); itv != sendEl.end(); itv++)
            {
                int n_req           = itv->second.size();
                int dest            = itv->first;

                MPI_Send(&n_req, 1,
                        MPI_INT, dest,
                        6798+78000*dest, comm);
                
                MPI_Send(&itv->second[0],
                        n_req, MPI_INT,
                        dest, 14876000+dest, comm);

                i++;
            }
        }
        else if (ishBack_schedule->SendFromRank2Rank[q].find( world_rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            
            MPI_Recv(&n_reqstd_ids,
                    1, MPI_INT, q,
                    6798+78000*world_rank,
                    comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0],
                    n_reqstd_ids,
                    MPI_INT, q,
                    14876000+world_rank,
                    comm, MPI_STATUS_IGNORE);
            
            adj_ids[q] = recv_reqstd_ids;

        }
    }
    
    std::map<int,int> shFid2el_rh;
    int adde = 0;
    int adde2 = 0;
    for(rcvit=adj_ids.begin();rcvit!=adj_ids.end();rcvit++)
    {
        int rrank = rcvit->first;
        int nelem = rcvit->second.size();
        
        for(int q=0;q<nelem;q++)
        {
            int fid = colorFh[rrank][q];
            
            if(shFid2el_rh.find(fid)==shFid2el_rh.end())
            {
                shFid2el_rh[fid] = rcvit->second[q];
                adde2++;
                if(rhp.find(fid)==rhp.end())
                {
                    adde++;
                    rhp[fid] = rcvit->second[q];
                }
            }
        }
    }
    
    std::map<int,std::vector<int> > face2node;
    
    int lfid;
    std::map<int,int>::iterator pritm;

    std::set<int> excludeShellV;
    
    //int nUniqueVerts_prisms = NonSharedVertsOwned.size() + SharedVertsOwned.size();
    
    int lvl             = 0;

    std::map<int,int>::iterator itmp;
    std::map<int,int> loc2glob_prism;
    std::map<int,int> loc2glob_prism_realz;

    int mapSizeLoc = tag2ElementID.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,comm);

    int mapSizeTot = distrimap->getNel();
    int* tag_loc = new int[mapSizeLoc];
    int* eid_loc = new int[mapSizeLoc];
    int* tag_tot = new int[mapSizeTot];
    int* eid_tot = new int[mapSizeTot];

    int tvid_tmp,iref_tmp;
    int i = 0;
    std::map<int,int>::iterator itred;

    for(itred=tag2ElementID.begin();itred!=tag2ElementID.end();itred++)
    {
        tag_loc[i] = itred->first;
        eid_loc[i] = itred->second;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    MPI_Allgatherv(tag_loc,
                   mapSizeLoc,
                   MPI_INT,
                   tag_tot,
                   nlocs,
                   offsets,
                   MPI_INT, comm);
    
    
    MPI_Allgatherv(eid_loc,
                   mapSizeLoc,
                   MPI_INT,
                   eid_tot,
                   nlocs,
                   offsets,
                   MPI_INT, comm);
    
    
    
    int key,val;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = tag_tot[i];
        val = eid_tot[i];
        if(tag2element_TetPrismInterface.find(key)==tag2element_TetPrismInterface.end())
        {
            tag2element_TetPrismInterface[key] = val;
        }
    }
    
    int nLocIntVrts         = NonSharedVertsOwned.size();
    int nLocShVrts          = SharedVertsOwned.size();
    DistributedParallelState* distnLocIntVrts = new DistributedParallelState(nLocIntVrts,comm);
    DistributedParallelState* distnLocShVrts  = new DistributedParallelState(nLocShVrts,comm);

    std::map<int,int> globV2tagV;
    int vert = 0;
    std::map<int,int>::iterator iitm;
    double xc,yc,zc;
    
    xcn_int    = new Array<double>(nLocIntVrts,3);
    xcn_shared = new Array<double>(nLocShVrts,3);
   
    for(iitm=NonSharedVertsOwned.begin();iitm!=NonSharedVertsOwned.end();iitm++)
    {
        int tag = iitm->first;
        int lvpart                   = tag2locV[tag];
        xc                           = locVerts[lvpart][0];
        yc                           = locVerts[lvpart][1];
        zc                           = locVerts[lvpart][2];
        
        xcn_int->setVal(vert,0,xc);
        xcn_int->setVal(vert,1,yc);
        xcn_int->setVal(vert,2,zc);
        
        tagV2globV[tag] = vert+distnLocIntVrts->getOffsets()[world_rank]+1;
        globV2tagV[vert+distnLocIntVrts->getOffsets()[world_rank]+1] = tag;
        
        vert++;
    }
    
    int vertsh = 0;
    for(iitm=SharedVertsOwned.begin();iitm!=SharedVertsOwned.end();iitm++)
    {
        int tag = iitm->first;
        int lvpart                   = tag2locV[tag];
        xc                           = locVerts[lvpart][0];
        yc                           = locVerts[lvpart][1];
        zc                           = locVerts[lvpart][2];
        
        xcn_shared->setVal(vertsh,0,xc);
        xcn_shared->setVal(vertsh,1,yc);
        xcn_shared->setVal(vertsh,2,zc);
        
        tagV2globV[tag] = SharedVertsOwned[tag];
        globV2tagV[SharedVertsOwned[tag]] = tag;
        
        vertsh++;
    }
    
    delete distnLocIntVrts;
    delete distnLocShVrts;
    delete nonSharedVertDistr;
    delete nonSharedFaceDistr;
    delete ownedSharedVrtsDist;
    delete[] ownedVs;
    delete ElementDistr;
    delete[] offsets;
    delete[] nlocs;
    

//    nnf->tagE2gE                    = tagE2gE;
//    nnf->gE2tagE                    = gE2tagE;
//
//    nnf->sharedVmap                 = sharedVmap;

    //nnf->tagV2globV             = tagV2globV;
//    nnf->SharedVertsOwned           = SharedVertsOwned;
//    nnf->NonSharedVertsOwned        = NonSharedVertsOwned;
//    nnf->SharedVertsNotOwned        = SharedVertsNotOwned;
    
//    nnf->sharedFace2Node            = SharedFace2Node;
//    nnf->intFace2Node               = InternalFace2Node;
//    nnf->bcFace2Node                = BoundaryFace2Node;
    
//    nnf->face2node                  = face2node;
//    nnf->ien                        = elements;
    //nnf->pbcmap                     = ref2bcface;
    //nnf->tag2ElementID              = tag2element_TetPrismInterface;

    //nnf->iet                        = iet;
    
    /**/
    //return nnf;
}

std::map<int,int> PrismaticLayer::getLeftElementGlobalIDForFaces()
{
    return lhp;
}

std::map<int,int> PrismaticLayer::getRightElementGlobalIDForFaces()
{
    return rhp;
}

std::map<int,std::vector<int> > PrismaticLayer::getOwnedSharedFace2NodeMap()
{
    return SharedFace2Node;
}

std::map<int,std::vector<int> > PrismaticLayer::getInternalFace2NodeMap()
{
    return InternalFace2Node;
}

std::map<int,std::vector<int> > PrismaticLayer::getBoundaryFace2NodeMap()
{
    return BoundaryFace2Node;
}

std::map<int,int> PrismaticLayer::getGlobal2TagElementMap()
{
    return gE2tagE;
}

std::map<int,int> PrismaticLayer::getTag2GlobalElementMap()
{
    return tagE2gE;
}

Array<double>* PrismaticLayer::getInternalCoordinates()
{
    return xcn_int;
}

Array<double>* PrismaticLayer::getSharedCoordinates()
{
    return xcn_shared;
}

std::map<int,int> PrismaticLayer::getNotOwnedSharedVerticesMap()
{
    return SharedVertsNotOwned;
}

std::map<int,int> PrismaticLayer::getVertexTag2GlobalMap()
{
    return tagV2globV;
}

std::map<int,int> PrismaticLayer::getTag2Element4TetPrismInterface()
{
    return tag2element_TetPrismInterface;
}

std::map<int,std::vector<int> > PrismaticLayer::getBoundaryCondition2FaceID()
{
    return ref2bcface;
}

std::map<int,std::vector<int> > PrismaticLayer::getBoundaryConditionZone2FaceID()
{
    return zone2bcface;
}

Array<int>* PrismaticLayer::getElementType()
{
    return iet;
}

std::map<int,int> PrismaticLayer::getSharedVertexMap()
{
    return sharedVmap;
}

PrismaticLayer::~PrismaticLayer()
{
    SharedFace2Node.clear();
    InternalFace2Node.clear();
    BoundaryFace2Node.clear();
    tagE2gE.clear();
    gE2tagE.clear();
    lhp.clear();
    rhp.clear();
    ref2bcface.clear();
    tag2element_TetPrismInterface.clear();
    sharedVmap.clear();
    zone2bcface.clear();
    SharedVertsNotOwned.clear();

    tagV2globV.clear();

    delete xcn_int;
    delete xcn_shared;
    delete iet;
}
