#include "adapt_DefinePrismMesh.h"

newNumberingNodesFaces* ExtractPrismsForHybridMesh(Array<int>* part_global,
                                          std::map<int,std::vector<int> > elements,
                                          std::map<int,std::vector<int> > ief_part_map,
                                          std::map<int,std::vector<int> > ifn_part_map,
                                          std::map<int,std::vector<int> > ife_part_map,
                                          std::map<int,std::vector<int> > if_ref_part_map,
                                          std::map<int,std::vector<int> > if_Nv_part_map,
                                          std::map<int,std::vector<int> > ushell,
                                          std::map<int,int> tag2locV,
                                          std::vector<Vert*> locVerts,
                                          std::map<int,int> shellvert2ref_glob,
                                          MPI_Comm comm )
{
    
    newNumberingNodesFaces* nnf = new newNumberingNodesFaces;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::map<int,std::vector<int> >::iterator prit;

    int testElem = 9;
    std::map<int,std::vector<int> > ref2bcface;
    int lf  = 0;
    int shf = 0;
    std::map<int,int> sharedFaces;
    std::set<int> gfid_set;
    std::map<int,int> gl_map;
    std::set<int> glf_set;
    std::map<int,int> glf_map;
    std::map<int,int> glf_map_inv;
    
    std::map<int,int> glf_map_sha;
    std::map<int,int> glf_map_sha_inv;
    int ref = 0;
    std::map<int,int> face2ref_prism;
    int nLocalFaces = 0;
    int nLocalVerts = 0;
    std::set<int> ufaces;
    std::set<int> gvid_set;
    int gvid,gfid,r0,r1,el0,el1;
    int elTel = 0;
    std::map<int,int> loc2globV;
    int foundshellvert = 0;
    std::map<int,int> unf;

    int changed = 0;
    int unchanged = 0;
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
                ref = if_ref_part_map[gfid][0];
            }
            
            int nppf = if_Nv_part_map[gfid][0];
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
                face2ref_prism[gfid] = ref;
                nLocalFaces++;
            }
            
            if(ufaces.find(gfid)==ufaces.end())
            {
                ufaces.insert(gfid);

                el0    = ife_part_map[gfid][0];
                el1    = ife_part_map[gfid][1];

                if(ref==2)
                {
                    r0 = part_global->getVal(el0,0);
                    r1 = part_global->getVal(el1,0);

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
        int nppf = if_Nv_part_map[gfid][0];
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
    std::map<int,int> SharedVertsOwned;
    std::map<int,int> SharedVertsNotOwned;
    std::map<int,int> sharedVmap;
    
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

    std::map<int,int> sharedFmap;
    int iFshared = distLocalFaces->getNel()-f2r.size();

    for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
    {
        sharedFmap[itvv->first] = iFshared;
        iFshared++;
    }
    
    DistributedParallelState* ElementDistr = new DistributedParallelState(elements.size(),comm);
    
    int uloc = 0;
    int u    = 0;
    
    std::map<int,std::vector<int> > colorRh;
    std::map<int,std::vector<int> > colorFh;

    Array<int>* iet = new Array<int>(elements.size(),1);
    std::map<int,std::vector<int> > shellFace2Node_prism;
    std::map<int,std::vector<int> > bcFace2Node_prism;
    std::map<int,std::vector<int> > sharedFace2Node_prism;
    std::map<int,std::vector<int> > intFace2Node_prism;
    std::map<int,int> tagE2gE;
    std::map<int,int> gE2tagE;
    std::map<int,int> pface2ref;
    std::map<int,int> lhp;
    std::map<int,int> rhp;
    std::map<int,int> tag2ElementID;
    
    int fowned = 0;
    int fownedInmap = 0;
    double orient0 = 1.0;
    std::map<int,double> orient_map;
    int ngf = 0;
    std::set<int> facemap;

    std::map<int,int> NonSharedVertsOwned;
    
    int notowned = 0;
    for(prit=elements.begin();prit!=elements.end();prit++)
    {
        int gEl     = prit->first;
        int rank    = world_rank;
        int lEl     = ElementDistr->getOffsets()[rank]+u+1;
        tagE2gE[gEl]  = lEl;
        gE2tagE[lEl]  = gEl;
        iet->setVal(u,0,6);
        
        for(int q=0;q<5;q++)
        {
            gfid = ief_part_map[gEl][q];
            int nfvrts = if_Nv_part_map[gfid][0];
            // Setting the reference to 13 for the shell faces;
            if(ushell.find(gfid)!=ushell.end())
            {
                ref                 = 13;
                tag2ElementID[gEl]  = lEl;
            }
            else
            {
                ref                 = if_ref_part_map[gfid][0];
            }
            
            if(facemap.find(gfid)==facemap.end())
            {
                facemap.insert(gfid);
                
                if(f2r.find(gfid)!=f2r.end() && ref != 13)
                {
                    if(f2r[gfid]==world_rank)
                    {
                        fowned++;
                        
                        if(sharedFace2Node_prism.find(gfid)==sharedFace2Node_prism.end())
                        {
                            fownedInmap++;
                            
                            el0                     = ife_part_map[gfid][0];
                            el1                     = ife_part_map[gfid][1];
                            
                            r0                      = part_global->getVal(el0,0);
                            r1                      = part_global->getVal(el1,0);
                            
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
                                
                                if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
                                {
                                    if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end())
                                    {
                                        SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
                                    }
                                }
                            }
                            
                            sharedFace2Node_prism[gfid] = fn_tag;
                            lhp[gfid]                   = lEl;

                        }
                    }
                }
                else
                {
                    if(ref == 2)
                    {
                        if(intFace2Node_prism.find(gfid)==intFace2Node_prism.end())
                        {
                            int nfvrts         = if_Nv_part_map[gfid][0];
                            
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
                            
                           
                            intFace2Node_prism[gfid] = fn_tag;
                            lhp[gfid]          = lEl;
                        }
                    }
                    
                    if(ref != 2 && ref!=13)
                    {
                        ref2bcface[ref].push_back(gfid);
                        
                        if(bcFace2Node_prism.find(gfid)==bcFace2Node_prism.end())
                        {
                            int nfvrts         = if_Nv_part_map[gfid][0];
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
                            

                            
                            bcFace2Node_prism[gfid]     = fn_tag;
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
    
    
    
    std::map<int,int> tag2ElementID_tot;
    int key,val;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = tag_tot[i];
        val = eid_tot[i];
        if(tag2ElementID_tot.find(key)==tag2ElementID_tot.end())
        {
            tag2ElementID_tot[key] = val;
        }
    }
    
    int nLocIntVrts         = NonSharedVertsOwned.size();
    int nLocShVrts          = SharedVertsOwned.size();
    DistributedParallelState* distnLocIntVrts = new DistributedParallelState(nLocIntVrts,comm);
    DistributedParallelState* distnLocShVrts  = new DistributedParallelState(nLocShVrts,comm);

    std::map<int,int> tag2glob_prism;
    std::map<int,int> glob2tag_prism;
    int vert = 0;
    std::map<int,int>::iterator iitm;
    double xc,yc,zc;
    
    Array<double>* xcn_parmmg_prism_int    = new Array<double>(nLocIntVrts,3);
    Array<double>* xcn_parmmg_prism_shared = new Array<double>(nLocShVrts,3);
   
    for(iitm=NonSharedVertsOwned.begin();iitm!=NonSharedVertsOwned.end();iitm++)
    {
        int tag = iitm->first;
        int lvpart                   = tag2locV[tag];
        xc                           = locVerts[lvpart]->x;
        yc                           = locVerts[lvpart]->y;
        zc                           = locVerts[lvpart]->z;
        
        xcn_parmmg_prism_int->setVal(vert,0,xc);
        xcn_parmmg_prism_int->setVal(vert,1,yc);
        xcn_parmmg_prism_int->setVal(vert,2,zc);
        
        tag2glob_prism[tag] = vert+distnLocIntVrts->getOffsets()[world_rank]+1;
        glob2tag_prism[vert+distnLocIntVrts->getOffsets()[world_rank]+1] = tag;
        
        vert++;
    }
    
    int vertsh = 0;
    for(iitm=SharedVertsOwned.begin();iitm!=SharedVertsOwned.end();iitm++)
    {
        int tag = iitm->first;
        int lvpart                   = tag2locV[tag];
        xc                           = locVerts[lvpart]->x;
        yc                           = locVerts[lvpart]->y;
        zc                           = locVerts[lvpart]->z;
        
        xcn_parmmg_prism_shared->setVal(vertsh,0,xc);
        xcn_parmmg_prism_shared->setVal(vertsh,1,yc);
        xcn_parmmg_prism_shared->setVal(vertsh,2,zc);
        
        tag2glob_prism[tag] = SharedVertsOwned[tag];
        glob2tag_prism[SharedVertsOwned[tag]] = tag;
        
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
    
    nnf->xcn_int                    = xcn_parmmg_prism_int;
    nnf->xcn_shared                 = xcn_parmmg_prism_shared;
    nnf->lh                         = lhp;
    nnf->rh                         = rhp;
    nnf->tagE2gE                    = tagE2gE;
    nnf->gE2tagE                    = gE2tagE;
    nnf->shellFace2Node             = shellFace2Node_prism;
    nnf->local2globalVertMap        = loc2glob_prism;
    nnf->sharedVmap                 = sharedVmap;

    nnf->tag2glob_prism             = tag2glob_prism;
    nnf->SharedVertsOwned           = SharedVertsOwned;
    nnf->NonSharedVertsOwned        = NonSharedVertsOwned;
    nnf->SharedVertsNotOwned        = SharedVertsNotOwned;
    
    nnf->sharedFace2Node            = sharedFace2Node_prism;
    nnf->intFace2Node               = intFace2Node_prism;
    nnf->bcFace2Node                = bcFace2Node_prism;
    
//    nnf->face2node                  = face2node;
//    nnf->ien                        = elements;
    nnf->pbcmap                     = ref2bcface;
    nnf->tag2ElementID              = tag2ElementID_tot;

    nnf->iet                        = iet;
    
    /**/
    return nnf;
    //===========================================================================
    //===========================================================================
    //===========================================================================
    
}





//newNumberingNodesFaces* DetermineNewNumberingOfElementSubset_Test(Array<int>* part_global,
//                                          std::map<int,std::vector<int> > elements,
//                                          std::map<int,std::vector<int> > ief_part_map,
//                                          std::map<int,std::vector<int> > ifn_part_map,
//                                          std::map<int,std::vector<int> > ife_part_map,
//                                          std::map<int,std::vector<int> > if_ref_part_map,
//                                          std::map<int,std::vector<int> > if_Nv_part_map,
//                                          std::map<int,std::vector<int> > ushell,
//                                          std::map<int,int> tag2locV,
//                                          std::vector<Vert*> locVerts,
//                                          std::map<int,int> shellvert2ref_glob,
//                                          MPI_Comm comm )
//{
//
//    newNumberingNodesFaces* nnf = new newNumberingNodesFaces;
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//
//    std::map<int,std::vector<int> >::iterator prit;
//
//
//    int testElem = 9;
//    std::map<int,std::vector<int> > ref2bcface;
//    int lf  = 0;
//    int shf = 0;
//    std::map<int,int> sharedFaces;
//    std::set<int> gfid_set;
//    std::map<int,int> gl_map;
//    std::set<int> glf_set;
//    std::map<int,int> glf_map;
//    std::map<int,int> glf_map_inv;
//
//    std::map<int,int> glf_map_sha;
//    std::map<int,int> glf_map_sha_inv;
//    int ref = 0;
//    std::map<int,int> face2ref_prism;
//    int nLocalFaces = 0;
//    int nLocalVerts = 0;
//    std::set<int> ufaces;
//    std::set<int> gvid_set;
//    int gvid,gfid,r0,r1,el0,el1;
//    int elTel = 0;
//    std::map<int,int> loc2globV;
//    int foundshellvert = 0;
//    std::map<int,int> unf;
//
//    int changed = 0;
//    int unchanged = 0;
//    std::map<int,int* > ifn_updt;
//    int ell = 0;
//
//    std::vector<std::set<int> > prism_faces;
//    std::set<int> prow_0;
//    prow_0.insert(0);
//    prow_0.insert(1);
//    prow_0.insert(2);
//    prism_faces.push_back(prow_0);
//    std::set<int> prow_1;
//    prow_1.insert(3);
//    prow_1.insert(5);
//    prow_1.insert(4);
//    prism_faces.push_back(prow_1);
//    std::set<int> prow_2;
//    prow_2.insert(0);
//    prow_2.insert(2);
//    prow_2.insert(4);
//    prow_2.insert(3);
//    prism_faces.push_back(prow_2);
//    std::set<int> prow_3;
//    prow_3.insert(1);
//    prow_3.insert(5);
//    prow_3.insert(4);
//    prow_3.insert(2);
//    prism_faces.push_back(prow_3);
//    std::set<int> prow_4;
//    prow_4.insert(0);
//    prow_4.insert(3);
//    prow_4.insert(5);
//    prow_4.insert(1);
//    prism_faces.push_back(prow_4);
//
//    for(prit=elements.begin();prit!=elements.end();prit++)
//    {
//        //key = GID, value global node nmber;
//        int gEl = prit->first;
//
////        if(world_rank == 0)
////        {
////            std::cout << "=================================" << std::endl;
////            std::cout << "              "<< ief_part_map[gEl].size() << "          --> ";
////        }
//        for(int j=0;j<ief_part_map[gEl].size();j++)
//        {
//            gfid = ief_part_map[gEl][j];
//
//            if(ushell.find(gfid)!=ushell.end())
//            {
//                ref = 13;
//            }
//            else
//            {
//                ref = if_ref_part_map[gfid][0];
//            }
//
//            int nppf = if_Nv_part_map[gfid][0];
////            if(world_rank == 0)
////            {
////                std::cout << nppf << ", ";
////            }
//            for(int h=0;h<nppf;h++)
//            {
//                gvid = ifn_part_map[gfid][h];
//
//                if(gvid_set.find(gvid)==gvid_set.end())
//                {
//                    gvid_set.insert(gvid);
//                }
//
//            }
//
//            if(gfid_set.find(gfid)==gfid_set.end())
//            {
//                gfid_set.insert(gfid);
//                face2ref_prism[gfid] = ref;
//                nLocalFaces++;
//            }
//
//            if(ufaces.find(gfid)==ufaces.end())
//            {
//                ufaces.insert(gfid);
//
//                el0    = ife_part_map[gfid][0];
//                el1    = ife_part_map[gfid][1];
//
//                if(ref==2)
//                {
//                    r0 = part_global->getVal(el0,0);
//                    r1 = part_global->getVal(el1,0);
//
//                    if(r0==world_rank && r1!=world_rank)
//                    {
//                        sharedFaces[gfid] = r0;
//                        shf++;
//                    }
//                    if(r0!=world_rank && r1==world_rank)
//                    {
//                        sharedFaces[gfid] = r1;
//                        shf++;
//                    }
//                }
//
//                if(ref!=2 && ref!=13)
//                {
//                    //ref2bcface[ref].push_back(gfid);
//                    lf++;
//                }
//            }
//        }
////        if(world_rank == 0)
////        {
////            std::cout << std::endl;
////        }
////        if(world_rank == 0)
////        {
////            std::cout << "=================================" << std::endl;
////
////        }
//        elTel++;
//    }
//
//
//    nLocalVerts = gvid_set.size();
//
//    int nSharedFaces                          = sharedFaces.size();
//    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
//    DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,comm);
//
//    int Nt_shFaces               = distSharedFaces->getNel();
//    int* shFace_offsets          = distSharedFaces->getOffsets();
//    int* shFace_nlocs            = distSharedFaces->getNlocs();
//    int* shFacesIDs              = new int[nSharedFaces];
//    int* shFaces_RankIDs         = new int[nSharedFaces];
//
//    int iter = 0;
//    std::set<int> UniqueSharedVertsOnRank_set;
//    std::vector<int> UniqueSharedVertsOnRank;
//    std::vector<int> UniqueSharedVertsOnRank_RankID;
//
//    std::map<int,int>::iterator itsf;
//    int lvrtid = 0;
//    int tel = shFace_offsets[world_rank];
//    for(itsf=sharedFaces.begin();itsf!=sharedFaces.end();itsf++)
//    {
//        shFacesIDs[iter] = itsf->first;
//        shFaces_RankIDs[iter] = itsf->second;
//        gfid = itsf->first;
//        int nppf = if_Nv_part_map[gfid][0];
//        for(int q=0;q<nppf;q++)
//        {
//            gvid   = ifn_part_map[gfid][q];
//
//            if(UniqueSharedVertsOnRank_set.find(gvid)==UniqueSharedVertsOnRank_set.end())
//            {
//                UniqueSharedVertsOnRank_set.insert(gvid);
//                UniqueSharedVertsOnRank.push_back(gvid);
//                UniqueSharedVertsOnRank_RankID.push_back(world_rank);
//                lvrtid++;
//            }
//        }
//        tel++;
//        iter++;
//    }
//
//    int nSharedVerts = UniqueSharedVertsOnRank.size();
//
//    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
//
//    int Nt_shVerts               = distSharedVerts->getNel();
//    int* shVerts_nlocs           = distSharedVerts->getNlocs();
//    int* shVerts_offsets         = distSharedVerts->getOffsets();
//
//    int* TotalSharedVerts        = new int[Nt_shVerts];
//    int* TotalSharedVerts_RankID = new int[Nt_shVerts];
//    int* TotalSharedFaces        = new int[Nt_shFaces];
//    int* TotalSharedFaces_RankID = new int[Nt_shFaces];
//
//    // Communicate vert map to all ranks.
//    MPI_Allgatherv(&UniqueSharedVertsOnRank[0],
//                   nSharedVerts,
//                   MPI_INT,
//                   TotalSharedVerts,
//                   shVerts_nlocs,
//                   shVerts_offsets,
//                   MPI_INT, comm);
//
//    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
//                   nSharedVerts,
//                   MPI_INT,
//                   TotalSharedVerts_RankID,
//                   shVerts_nlocs,
//                   shVerts_offsets,
//                   MPI_INT, comm);
//
//    // Communicate face map to all ranks.
//    MPI_Allgatherv(shFacesIDs,
//                   nSharedFaces,
//                   MPI_INT,
//                   TotalSharedFaces,
//                   shFace_nlocs,
//                   shFace_offsets,
//                   MPI_INT, comm);
//
//    MPI_Allgatherv(shFaces_RankIDs,
//                   nSharedFaces,
//                   MPI_INT,
//                   TotalSharedFaces_RankID,
//                   shFace_nlocs,
//                   shFace_offsets,
//                   MPI_INT, comm);
//
//    int tmp;
//    std::map<int,int> f2r;
//    std::set<int> f2r_s;
//
//
//    int* NewGlobVertCountPerRank = new int[world_size];
//    int* NewGlobFaceCountPerRank = new int[world_size];
//
//    for(int u=0;u<world_size;u++)
//    {
//        NewGlobVertCountPerRank[u] = 0;
//        NewGlobFaceCountPerRank[u] = 0;
//    }
//
//    for(int i=0;i<Nt_shFaces;i++)
//    {
//        int key = TotalSharedFaces[i];
//        int val = TotalSharedFaces_RankID[i];
//
//        if(f2r_s.find(key)==f2r_s.end())
//        {
//            f2r_s.insert(key);
//            f2r[key]=val;
//            NewGlobFaceCountPerRank[val]=NewGlobFaceCountPerRank[val]+1;
//        }
//        else
//        {
//            tmp = f2r[key];
//            if(val<tmp)
//            {
//                f2r[key]=val;
//            }
//            if(val>tmp)
//            {
//                f2r[key]=tmp;
//            }
//        }
//    }
//
//    std::map<int,int> v2r;
//    std::set<int> v2r_s;
//    int* NewGlobVertOffset = new int[world_size];
//    int* NewGlobFaceOffset = new int[world_size];
//    int* owned_verts       = new int[world_size];
//    for(int u=0;u<world_size;u++)
//    {
//        owned_verts[u]       = 0;
//        NewGlobVertOffset[u] = 0;
//        NewGlobFaceOffset[u] = 0;
//    }
//
//    for(int i=0;i<Nt_shVerts;i++)
//    {
//        int key = TotalSharedVerts[i];
//        int val = TotalSharedVerts_RankID[i];
//
//        if(v2r_s.find(key)==v2r_s.end())
//        {
//            v2r_s.insert(key);
//            v2r[key]=val;
//            NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
//        }
//        else
//        {
//            tmp = v2r[key];
//
//            if(val<tmp)
//            {
//                v2r[key]=val;
//
//            }
//            if(val>tmp)
//            {
//                v2r[key]=tmp;
//            }
//        }
//    }
//
//    std::map<int,int>::iterator itmm;
//    int nOwnedSharedVerts = 0;
//    for(itmm=v2r.begin();itmm!=v2r.end();itmm++)
//    {
//        if(itmm->second==world_rank)
//        {
//            nOwnedSharedVerts++;
//        }
//    }
//
//
//
//
//
//
//    for(int u=1;u<world_size;u++)
//    {
//        NewGlobVertOffset[u] = NewGlobVertOffset[u-1]+NewGlobVertCountPerRank[u-1];
//        NewGlobFaceOffset[u] = NewGlobFaceOffset[u-1]+NewGlobFaceCountPerRank[u-1];
//    }
//
//    std::map<int,std::vector<int> >::iterator itv;
//
//    std::map<int,int> sharedVertsGlobal;
//
//    for(int u=0;u<world_size;u++)
//    {
//        NewGlobVertCountPerRank[u] = 0;
//        NewGlobFaceCountPerRank[u] = 0;
//    }
//
////    int nNonSharedVerts          = nLocalVerts-owned_verts[world_rank];
////    int nNonSharedFaces          = nLocalFaces-owned_faces[world_rank];
//
//    int nNonSharedVerts          = nLocalVerts-nSharedVerts;
//    int nNonSharedFaces          = nLocalFaces-nSharedFaces;
//
//    //std::cout << "LOCALVERTS " << world_rank << "   -- = " << nNonSharedVerts <<" "<< nLocalVerts << " " << owned_verts[world_rank] << " " << nSharedVerts << std::endl;
//
//    DistributedParallelState* nonSharedVertDistr = new DistributedParallelState(nNonSharedVerts,comm);
//    DistributedParallelState* nonSharedFaceDistr = new DistributedParallelState(nNonSharedFaces,comm);
//
//
//    int nNonSharedVertsTot = nonSharedVertDistr->getNel();
//    int iVshared          = nonSharedVertDistr->getNel()+1;//-shellvert2ref_glob.size();
//
//    int iVshared_bef       = iVshared;
//    std::map<int,int >::iterator itvv;
//
//    DistributedParallelState* ownedSharedVrtsDist = new DistributedParallelState(nOwnedSharedVerts,comm);
//
//    int* offsetsOwnedSharedVerts = ownedSharedVrtsDist->getOffsets();
//    int* locsOwnedSharedVerts = ownedSharedVrtsDist->getNlocs();
//    int lshaVrt = 0;
//    int tshaVrt = 0;
//    std::map<int,int> SharedVertsOwned;
//    std::map<int,int> SharedVertsNotOwned;
//    std::map<int,int> sharedVmap;
//
//    int* ownedVs = new int[world_size];
//    for(int i=0;i<world_size;i++)
//    {
//        ownedVs[i] = 0;
//    }
//
//    for(itmm=v2r.begin();itmm!=v2r.end();itmm++)
//    {
//        int globid = nNonSharedVertsTot+ownedSharedVrtsDist->getOffsets()[itmm->second]+ownedVs[itmm->second]+1;
//
//        sharedVmap[itmm->first] = globid;
//
//        if(itmm->second==world_rank &&
//           shellvert2ref_glob.find(itmm->first)==shellvert2ref_glob.end())
//        {
//            SharedVertsOwned[itmm->first] = globid;
//        }
//
//        ownedVs[itmm->second] = ownedVs[itmm->second]+1;
//    }
//
//
//    if(world_rank == 0)
//    {
//        for(int i=0;i<world_size;i++)
//        {
//            std::cout << "owned " << i << " " << ownedVs[i] << std::endl;
//        }
//    }
//
//
//
//    std::map<int,int> sharedFmap;
//    int iFshared = distLocalFaces->getNel()-f2r.size();
//
//    for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
//    {
//        sharedFmap[itvv->first] = iFshared;
//        iFshared++;
//    }
//
//    int lbvid         = nonSharedVertDistr->getOffsets()[world_rank];
//    int lbfid         = nonSharedFaceDistr->getOffsets()[world_rank];
//
//    //std::cout << "===checking values = " << " "<< " world_rank  " << world_rank  << " " << nLocalVerts << " " << nSharedVerts << " " << ownedsharedVmap << " " << nNonSharedVerts << " " << unf.size() <<  std::endl;
//
//    DistributedParallelState* ElementDistr = new DistributedParallelState(elements.size(),comm);
//
//    int offsets_element = ElementDistr->getOffsets()[world_rank];
//    std::map<int,std::vector<int> > ifref;
//    std::map<int,std::vector<int> > ifn;
//    int uloc = 0;
//    int u    = 0;
//    std::map<int,std::vector<int> > ief_new;
//    std::map<int,int> pface2ref;
//    std::map<int,int> lhp;
//    std::map<int,int> rhp;
//    std::map<int,int> tag2ElementID;
//    std::map<int,std::vector<int> > colorRh;
//    std::map<int,std::vector<int> > colorFh;
//    std::map<int,int> tagE2gE;
//    std::map<int,int> gE2tagE;
//    std::map<int,int> shF2El;
//    std::map<int,int> lhp_sha;
//    std::map<int,int> rhp_sha;
//    Array<int>* iet = new Array<int>(elements.size(),1);
//
//
//    std::map<int,int> owned_glov2tagv_NEW;
//    std::map<int,int> owned_tagv2glov_NEW;
//    std::map<int,int> owned_glov2tagv;
//    std::map<int,int> owned_tagv2glov;
//    std::map<int,Vert*> tag2coords;
//
//    std::map<int,std::vector<int> > shellFace2Node_prism;
//    std::map<int,std::vector<int> > bcFace2Node_prism;
//    std::map<int,std::vector<int> > sharedFace2Node_prism;
//    std::map<int,std::vector<int> > intFace2Node_prism;
//
//    int fowned = 0;
//    int fownedInmap = 0;
//    double orient0 = 1.0;
//    std::map<int,double> orient_map;
//    int ngf = 0;
//    std::set<int> facemap;
//
//    std::map<int,int> NonSharedVertsOwned;
//
//    int notowned = 0;
//    for(prit=elements.begin();prit!=elements.end();prit++)
//    {
//        int gEl     = prit->first;
//        int rank    = world_rank;
//        int lEl     = ElementDistr->getOffsets()[rank]+u+1;
//        tagE2gE[gEl]  = lEl;
//        gE2tagE[lEl]  = gEl;
//        iet->setVal(u,0,6);
//
//        // Compute the center of the element in order;
//        Vert* Vijk = new Vert;
//        Vijk->x = 0.0;
//        Vijk->y = 0.0;
//        Vijk->z = 0.0;
//        // compute element center;
//        int nvp = prit->second.size();
//        for(int q=0;q<nvp;q++)
//        {
//            int tag  = prit->second[q];
//            int lvp  = tag2locV[tag];
//
//            Vijk->x = Vijk->x + locVerts[lvp]->x;
//            Vijk->y = Vijk->y + locVerts[lvp]->y;
//            Vijk->z = Vijk->z + locVerts[lvp]->z;
//        }
//
//        Vijk->x = Vijk->x/nvp;
//        Vijk->y = Vijk->y/nvp;
//        Vijk->z = Vijk->z/nvp;
//
//        std::vector<int> newFids(ief_part_map[gEl].size());
//
////        if(world_rank == 2)
////        {
////            std::cout <<"=======" << std::endl;
////
////            std::cout <<ief_part_map[gEl].size() << " ---> ";
////        }
//        for(int q=0;q<5;q++)
//        {
//            gfid = ief_part_map[gEl][q];
//            int nfvrts = if_Nv_part_map[gfid][0];
//            // Setting the reference to 13 for the shell faces;
//            if(ushell.find(gfid)!=ushell.end())
//            {
//                ref                 = 13;
//                tag2ElementID[gEl]  = lEl;
//            }
//            else
//            {
//                ref                 = if_ref_part_map[gfid][0];
//            }
//
//            if(facemap.find(gfid)==facemap.end())
//            {
//                facemap.insert(gfid);
//                if(f2r.find(gfid)!=f2r.end() && ref != 13)
//                {
//                    if(f2r[gfid]==world_rank)
//                    {
//                        fowned++;
//
//                        if(sharedFace2Node_prism.find(gfid)==sharedFace2Node_prism.end())
//                        {
//                            fownedInmap++;
//                            int lfid2               = sharedFmap[gfid];
//
//                            el0                     = ife_part_map[gfid][0];
//                            el1                     = ife_part_map[gfid][1];
//
//                            r0                      = part_global->getVal(el0,0);
//                            r1                      = part_global->getVal(el1,0);
//
//                            if(r0==world_rank && r1!=world_rank)
//                            {
//                                colorRh[r1].push_back(el1);
//                                colorFh[r1].push_back(gfid);
//                            }
//                            else if(r1==world_rank && r0!=world_rank)
//                            {
//                                colorRh[r0].push_back(el0);
//                                colorFh[r0].push_back(gfid);
//                            }
//
//                            int nfvrts         = if_Nv_part_map[gfid][0];
//
//                            std::vector<int> fn_tag(nfvrts);
////                            std::vector<int> fn_tag1(nfvrts);
////                            Vert* VcF = new Vert;
////                            std::vector<Vert*> Vfaces;
//                            for(int n=0;n<nfvrts;n++)
//                            {
//                                fn_tag[n] = ifn_part_map[gfid][n];
////                                fn_tag1[n] = ifn_part_map[gfid][n];
////
////                                int lvp  = tag2locV[fn_tag[n]];
////                                Vert* Vf = new Vert;
////                                Vf->x = locVerts[lvp]->x;
////                                Vf->y = locVerts[lvp]->y;
////                                Vf->z = locVerts[lvp]->z;
////
////                                VcF->x = VcF->x + Vf->x;
////                                VcF->y = VcF->y + Vf->y;
////                                VcF->z = VcF->z + Vf->z;
////
////                                Vfaces.push_back(Vf);
//
//                                if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                {
//                                    if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end())
//                                    {
//                                        SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
//                                    }
//                                }
//                            }
//
////                            VcF->x = VcF->x/nfvrts;
////                            VcF->y = VcF->y/nfvrts;
////                            VcF->z = VcF->z/nfvrts;
////
////                            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
////
////                            if(orient0<0.0)
////                            {
////                                if(nfvrts==4)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[3];
////                                    fn_tag[2]=fn_tag1[2];
////                                    fn_tag[3]=fn_tag1[1];
////                                }
////                                if(nfvrts==3)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[2];
////                                    fn_tag[2]=fn_tag1[1];
////                                }
////                            }
//
//                            sharedFace2Node_prism[gfid] = fn_tag;
//                            lhp[gfid]                   = lEl;
//
//                        }
//                    }
//                }
//                else
//                {
//                    if(ref == 2)
//                    {
//                        if(intFace2Node_prism.find(gfid)==intFace2Node_prism.end())
//                        {
//                            int nfvrts         = if_Nv_part_map[gfid][0];
//
//                            std::vector<int> fn_glo(nfvrts);
//                            std::vector<int> fn_tag(nfvrts);
////                            std::vector<int> fn_tag1(nfvrts);
////                            Vert* VcF = new Vert;
////                            std::vector<Vert*> Vfaces;
//                            for(int n=0;n<nfvrts;n++)
//                            {
//                                fn_tag[n] = ifn_part_map[gfid][n];
////                                fn_tag1[n] = ifn_part_map[gfid][n];
////
////                                int lvp  = tag2locV[fn_tag[n]];
////                                Vert* Vf = new Vert;
////                                Vf->x = locVerts[lvp]->x;
////                                Vf->y = locVerts[lvp]->y;
////                                Vf->z = locVerts[lvp]->z;
////
////                                VcF->x = VcF->x + Vf->x;
////                                VcF->y = VcF->y + Vf->y;
////                                VcF->z = VcF->z + Vf->z;
////
////                                Vfaces.push_back(Vf);
//
//                                if(v2r.find(fn_tag[n])!=v2r.end())
//                                {
//                                    if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                    {
//                                        if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                        {
//                                            SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
//                                        }
//                                    }
//                                }
//                                else
//                                {
//                                    if(NonSharedVertsOwned.find(fn_tag[n])==NonSharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                    {
//                                        NonSharedVertsOwned[fn_tag[n]]  = gvid;
//                                        gvid++;
//                                    }
//                                }
//                            }
//
////                            VcF->x = VcF->x/nfvrts;
////                            VcF->y = VcF->y/nfvrts;
////                            VcF->z = VcF->z/nfvrts;
////
////                            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
////
////                            if(orient0<0.0)
////                            {
////                                if(nfvrts==4)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[3];
////                                    fn_tag[2]=fn_tag1[2];
////                                    fn_tag[3]=fn_tag1[1];
////                                }
////                                if(nfvrts==3)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[2];
////                                    fn_tag[2]=fn_tag1[1];
////                                }
////                            }
//
//
//                            intFace2Node_prism[gfid] = fn_tag;
//                            lhp[gfid]          = lEl;
//                            lbfid              = lbfid + 1;
//                        }
//                    }
//
//                    if(ref != 2)
//                    {
//                        ref2bcface[ref].push_back(gfid);
//
//                        if(bcFace2Node_prism.find(gfid)==bcFace2Node_prism.end())
//                        {
//                            int nfvrts         = if_Nv_part_map[gfid][0];
//
//                            //std::vector<int> fn_glo(nfvrts);
//                            std::vector<int> fn_tag(nfvrts);
//                            //std::vector<int> fn_tag1(nfvrts);
//
//                            Vert* VcF = new Vert;
//                            std::vector<Vert*> Vfaces;
//                            for(int n=0;n<nfvrts;n++)
//                            {
//                                fn_tag[n] = ifn_part_map[gfid][n];
////                                fn_tag1[n] = ifn_part_map[gfid][n];
////                                int lvp  = tag2locV[fn_tag[n]];
////                                Vert* Vf = new Vert;
////                                Vf->x = locVerts[lvp]->x;
////                                Vf->y = locVerts[lvp]->y;
////                                Vf->z = locVerts[lvp]->z;
////
////                                VcF->x = VcF->x + Vf->x;
////                                VcF->y = VcF->y + Vf->y;
////                                VcF->z = VcF->z + Vf->z;
////
////                                Vfaces.push_back(Vf);
//
//                                if(v2r.find(fn_tag[n])!=v2r.end())
//                                {
//                                    if(SharedVertsOwned.find(fn_tag[n])==SharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                    {
//                                        if(SharedVertsNotOwned.find(fn_tag[n])==SharedVertsNotOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                        {
//                                            SharedVertsNotOwned[fn_tag[n]] = sharedVmap[fn_tag[n]];
//                                        }
//                                    }
//                                }
//                                else
//                                {
//                                    if(NonSharedVertsOwned.find(fn_tag[n])==NonSharedVertsOwned.end() && shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                    {
//                                        NonSharedVertsOwned[fn_tag[n]]  = gvid;
//                                        gvid++;
//                                    }
//                                }
//                            }
//
////                            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
////
////                            if(orient0<0.0)
////                            {
////                                if(nfvrts==4)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[3];
////                                    fn_tag[2]=fn_tag1[2];
////                                    fn_tag[3]=fn_tag1[1];
////                                }
////                                if(nfvrts==3)
////                                {
////                                    fn_tag[0]=fn_tag1[0];
////                                    fn_tag[1]=fn_tag1[2];
////                                    fn_tag[2]=fn_tag1[1];
////                                }
////                            }
//
//                            bcFace2Node_prism[gfid]     = fn_tag;
//                            lhp[gfid]                   = lEl;
//                            lbfid                       = lbfid + 1;
//                        }
//                    }
//                }
//            }
//            else
//            {
//                rhp[gfid]           = lEl;
//            }
//        }
//
//        ief_new[gEl] = newFids;
//        uloc++;
//        u++;
//
//    }
//
//    //std::cout << "ngf " << ngf << " ON WORKRANK " << world_rank << " " << bcFace2Node_prism.size() << " " << intFace2Node_prism.size() << " " << sharedFace2Node_prism.size() << " " << NonSharedVertsOwned.size() << " " << SharedVertsOwned.size() << " " << SharedVertsNotOwned.size() << " " << notowned <<  std::endl;
//    //std::cout << "facemap sizes " << world_rank << " " << intFace2Node_prism.size() << " " << bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << std::endl;
//    //std::cout << " before " << world_rank << " :: " << lhp.size() << " " << rhp.size() + bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << std::endl;
//
//    ScheduleObj* rh_schedule = DoScheduling(colorRh,comm);
//    std::map<int,std::vector<int> > recv_rhElIDs;
//    for(int q=0;q<world_size;q++)
//    {
//        if(world_rank==q)
//        {
//            int i=0;
//
//            for (itv = colorRh.begin(); itv != colorRh.end(); itv++)
//            {
//                int n_req           = itv->second.size();
//                int dest            = itv->first;
//
//                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
//                MPI_Send(&itv->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
//                i++;
//            }
//        }
//        else if (rh_schedule->SendFromRank2Rank[q].find( world_rank ) != rh_schedule->SendFromRank2Rank[q].end())
//        {
//            int n_reqstd_ids;
//            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);
//
//            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
//
//            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
//            recv_rhElIDs[q] = recv_reqstd_ids;
//        }
//    }
//
//    std::map<int,std::vector<int> > sendEl;
//    std::map<int,std::vector<int> >::iterator rcvit;
//    for(rcvit=recv_rhElIDs.begin();rcvit!=recv_rhElIDs.end();rcvit++)
//    {
//        int frank = rcvit->first;
//        int nE    = rcvit->second.size();
//
//        for(int j=0;j<nE;j++)
//        {
//            int gEl = tagE2gE[rcvit->second[j]];
//            sendEl[frank].push_back(gEl);
//        }
//    }
//
//    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);
//
//    std::map<int,std::vector<int> > adj_ids;
//    for(int q=0;q<world_size;q++)
//    {
//        if(world_rank==q)
//        {
//            int i=0;
//            for (itv = sendEl.begin(); itv != sendEl.end(); itv++)
//            {
//                int n_req           = itv->second.size();
//                int dest            = itv->first;
//
//                MPI_Send(&n_req, 1,
//                        MPI_INT, dest,
//                        6798+78000*dest, comm);
//
//                MPI_Send(&itv->second[0],
//                        n_req, MPI_INT,
//                        dest, 14876000+dest, comm);
//
//                i++;
//            }
//        }
//        else if (ishBack_schedule->SendFromRank2Rank[q].find( world_rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
//        {
//            int n_reqstd_ids;
//
//            MPI_Recv(&n_reqstd_ids,
//                    1, MPI_INT, q,
//                    6798+78000*world_rank,
//                    comm, MPI_STATUS_IGNORE);
//
//            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
//
//            MPI_Recv(&recv_reqstd_ids[0],
//                    n_reqstd_ids,
//                    MPI_INT, q,
//                    14876000+world_rank,
//                    comm, MPI_STATUS_IGNORE);
//
//            adj_ids[q] = recv_reqstd_ids;
//
//        }
//    }
//
//    std::cout << " lhp vs rhp sizes BEFORE " << bcFace2Node_prism.size() << " " << lhp.size()-rhp.size() << std::endl;
//
//
//    std::map<int,int> shFid2el_rh;
//    int adde = 0;
//    int adde2 = 0;
//    for(rcvit=adj_ids.begin();rcvit!=adj_ids.end();rcvit++)
//    {
//        int rrank = rcvit->first;
//        int nelem = rcvit->second.size();
//
//        //std::cout << "wr " << world_rank << " sending " << rrank << " " << nelem << std::endl;
//
//        for(int q=0;q<nelem;q++)
//        {
//            int tag = colorRh[rrank][q];
//            int fid = colorFh[rrank][q];
//
//            if(shFid2el_rh.find(fid)==shFid2el_rh.end())
//            {
//                shFid2el_rh[fid] = rcvit->second[q];
//                adde2++;
//                if(rhp.find(fid)==rhp.end())
//                {
//                    adde++;
//                    rhp[fid] = rcvit->second[q];
//                }
//            }
//        }
//    }
//
//    std::cout << " lhp vs rhp sizes AFTER " << bcFace2Node_prism.size() << " " << lhp.size()-rhp.size() << std::endl;
//
//    //std::cout << " after " << world_rank << " :: " << lhp.size() << " " << rhp.size() + bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << " adde " << adde << " " << adde2 << std::endl;
//
//    std::map<int,std::vector<int> > face2node;
//
//    int lfid;
//    std::map<int,int>::iterator pritm;
//    int ff13 = 0;
//
//    std::set<int> excludeShellV;
//
//    int nUniqueVerts_prisms = NonSharedVertsOwned.size() + SharedVertsOwned.size();
//
//    DistributedParallelState* distri_prism_verts  = new DistributedParallelState(nUniqueVerts_prisms,comm);
//
//    int NtotVerts_prism = distri_prism_verts->getNel();
//
//    //std::cout << "Number of total unique " << world_rank << " " << NtotVerts_prism << " face sizes " << intFace2Node_prism.size() << " " << bcFace2Node_prism.size() << " " << sharedFace2Node_prism.size() << " " << shellFace2Node_prism.size() << " " << lhp.size() << " " << rhp.size() << " " << lhp_sha.size() << " " << rhp_sha.size() << " " << adde << std::endl;
//
//    int lvg             = distri_prism_verts->getOffsets()[world_rank]+1;
//    int lvl             = 0;
//
//    std::map<int,int>::iterator itmp;
//    std::map<int,int> loc2glob_prism;
//    std::map<int,int> loc2glob_prism_realz;
//    //std::cout << "world Rank " << world_rank << " " << gl_map_i.size() << " " << gl_map.size() << " vakjes " << vak1.size() << " " << vak2.size() << std::endl;
//
//    //owned_tagv2glov_NEW
//
////    for(itmp=owned_tagv2glov_NEW.begin();itmp!=owned_tagv2glov_NEW.end();itmp++)
////    {
////        int tag     = itmp->first;
////
////        double xc = tag2coords[tag]->x;
////        double yc = tag2coords[tag]->y;
////        double zc = tag2coords[tag]->z;
////
////        xcn_parmmg_prism->setVal(lvl,0,xc);
////        xcn_parmmg_prism->setVal(lvl,1,yc);
////        xcn_parmmg_prism->setVal(lvl,2,zc);
////
////        loc2glob_prism[tag]       = lvg;
////        loc2glob_prism_realz[tag] = lvl;
////
////        lvl++;
////        lvg++;
////    }
//
//    int mapSizeLoc = tag2ElementID.size();
//    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,comm);
//
//    int mapSizeTot = distrimap->getNel();
//    int* tag_loc = new int[mapSizeLoc];
//    int* eid_loc = new int[mapSizeLoc];
//    int* tag_tot = new int[mapSizeTot];
//    int* eid_tot = new int[mapSizeTot];
//
//    int tvid_tmp,iref_tmp;
//    int i = 0;
//    std::map<int,int>::iterator itred;
//
//    for(itred=tag2ElementID.begin();itred!=tag2ElementID.end();itred++)
//    {
//        tag_loc[i] = itred->first;
//        eid_loc[i] = itred->second;
//        i++;
//    }
//
//    int* offsets = distrimap->getOffsets();
//    int* nlocs   = distrimap->getNlocs();
//
//
//    MPI_Allgatherv(tag_loc,
//                   mapSizeLoc,
//                   MPI_INT,
//                   tag_tot,
//                   nlocs,
//                   offsets,
//                   MPI_INT, comm);
//
//
//    MPI_Allgatherv(eid_loc,
//                   mapSizeLoc,
//                   MPI_INT,
//                   eid_tot,
//                   nlocs,
//                   offsets,
//                   MPI_INT, comm);
//
//
//
//    std::map<int,int> tag2ElementID_tot;
//    int key,val;
//    for(int i=0;i<mapSizeTot;i++)
//    {
//        key = tag_tot[i];
//        val = eid_tot[i];
//        if(tag2ElementID_tot.find(key)==tag2ElementID_tot.end())
//        {
//            tag2ElementID_tot[key] = val;
//        }
//    }
//
//    int lvert = 0;
//
//    // testing site!
//
//
//
////    for(prit=elements.begin();prit!=elements.end();prit++)
////    {
////        int gEl     = prit->first;
////        int rank    = world_rank;
////        int lEl     = ElementDistr->getOffsets()[rank]+u+1;
////        iet->setVal(u,0,6);
////        tagE2gE[gEl]  = lEl;
////
////        std::vector<int> newFids(ief_part_map[gEl].size());
////
////        for(int q=0;q<ief_part_map[gEl].size();q++)
////        {
////            gfid = ief_part_map[gEl][q];
////
////
////
////
////        }
////    }
////
//    int nLocIntVrts         = NonSharedVertsOwned.size();
//    int nLocShVrts          = SharedVertsOwned.size();
//    int nLocTotVrts         = nLocIntVrts+nLocShVrts;
//    DistributedParallelState* distnLocIntVrts = new DistributedParallelState(nLocIntVrts,comm);
//    DistributedParallelState* distnLocShVrts  = new DistributedParallelState(nLocShVrts,comm);
//
//    int nIntVrts = distnLocIntVrts->getNel();
//    std::map<int,int> tag2glob_prism;
//    std::map<int,int> glob2tag_prism;
//    int vert = 0;
//    std::map<int,int>::iterator iitm;
//    double xc,yc,zc;
//    Array<double>* xcn_parmmg_prism_int    = new Array<double>(nLocIntVrts,3);
//    Array<double>* xcn_parmmg_prism_shared = new Array<double>(nLocShVrts,3);
//
//
//    for(iitm=NonSharedVertsOwned.begin();iitm!=NonSharedVertsOwned.end();iitm++)
//    {
//        int tag = iitm->first;
//        int lvpart                   = tag2locV[tag];
//        xc                           = locVerts[lvpart]->x;
//        yc                           = locVerts[lvpart]->y;
//        zc                           = locVerts[lvpart]->z;
//
//        xcn_parmmg_prism_int->setVal(vert,0,xc);
//        xcn_parmmg_prism_int->setVal(vert,1,yc);
//        xcn_parmmg_prism_int->setVal(vert,2,zc);
//
//        tag2glob_prism[tag] = vert+distnLocIntVrts->getOffsets()[world_rank]+1;
//        glob2tag_prism[vert+distnLocIntVrts->getOffsets()[world_rank]+1] = tag;
//
//        vert++;
//    }
//
//    int vertsh = 0;
//    for(iitm=SharedVertsOwned.begin();iitm!=SharedVertsOwned.end();iitm++)
//    {
//        int tag = iitm->first;
//        int lvpart                   = tag2locV[tag];
//        xc                           = locVerts[lvpart]->x;
//        yc                           = locVerts[lvpart]->y;
//        zc                           = locVerts[lvpart]->z;
//
//        xcn_parmmg_prism_shared->setVal(vertsh,0,xc);
//        xcn_parmmg_prism_shared->setVal(vertsh,1,yc);
//        xcn_parmmg_prism_shared->setVal(vertsh,2,zc);
//
//        tag2glob_prism[tag] = SharedVertsOwned[tag];
//        glob2tag_prism[SharedVertsOwned[tag]] = tag;
//
//        vertsh++;
//    }
//
//    std::cout << "check lengths Verts = " << nUniqueVerts_prisms << " " <<NonSharedVertsOwned.size() << " " << SharedVertsOwned.size() << " " << tag2glob_prism.size() << std::endl;
//
//
//
//    nnf->xcn_int                    = xcn_parmmg_prism_int;
//    nnf->xcn_shared                 = xcn_parmmg_prism_shared;
//    nnf->tag2coords                 = tag2coords;
//    nnf->lh                         = lhp;
//    nnf->rh                         = rhp;
//    nnf->tagE2gE                    = tagE2gE;
//    nnf->gE2tagE                    = gE2tagE;
//    nnf->lh_sha                     = lhp_sha;
//    nnf->rh_sha                     = rhp_sha;
//    nnf->nTotUniqueNonShellVerts    = iVshared-1;
//    nnf->shellFace2Node             = shellFace2Node_prism;
//    nnf->sharedvert2rank            = v2r;
//    nnf->local2globalVertMap        = loc2glob_prism;
//    nnf->local2globalVertMapReal    = loc2glob_prism_realz;
//    nnf->sharedVmap                 = sharedVmap;
//
//    nnf->tag2glob_prism             = tag2glob_prism;
//    nnf->glob2tag_prism             = glob2tag_prism;
//    nnf->SharedVertsOwned           = SharedVertsOwned;
//    nnf->NonSharedVertsOwned        = NonSharedVertsOwned;
//    nnf->SharedVertsNotOwned        = SharedVertsNotOwned;
//    nnf->sharedFace2Node            = sharedFace2Node_prism;
//    nnf->intFace2Node               = intFace2Node_prism;
//    nnf->bcFace2Node                = bcFace2Node_prism;
//
//    nnf->glob2locF                  = glf_map;
//    nnf->loc2globF                  = glf_map_inv;
//    nnf->glob2locF_sha              = glf_map_sha;
//    nnf->loc2globF_sha              = glf_map_sha_inv;
//    nnf->face2node                  = face2node;
//    nnf->ifref                      = pface2ref;
//    nnf->ien                        = elements;
//    nnf->ief                        = ief_new;
//    nnf->dist                       = ElementDistr;
//    nnf->pbcmap                     = ref2bcface;
//    nnf->tag2ElementID              = tag2ElementID_tot;
//    nnf->localV2tagV                = owned_glov2tagv_NEW;
//    nnf->tagV2localV                = owned_tagv2glov_NEW;
////    nnf->localV2tagV_all            = glov2tagv_all;
////    nnf->tagV2localV_all            = tagv2glov_all;
//    nnf->iet                       = iet;
//    /**/
//    return nnf;
//    //===========================================================================
//    //===========================================================================
//    //===========================================================================
//
//}
//
//
//
//
//
//
//
//newNumberingNodesFaces* DetermineNewNumberingOfElementSubset(Array<int>* part_global,
//                                          std::map<int,std::vector<int> > elements,
//                                          std::map<int,std::vector<int> > ief_part_map,
//                                          std::map<int,std::vector<int> > ifn_part_map,
//                                          std::map<int,std::vector<int> > ife_part_map,
//                                          std::map<int,std::vector<int> > if_ref_part_map,
//                                          std::map<int,std::vector<int> > if_Nv_part_map,
//                                          std::map<int,std::vector<int> > ushell,
//                                          std::map<int,int> tag2locV,
//                                          std::vector<Vert*> locVerts,
//                                          std::map<int,int> shellvert2ref_glob,
//                                          MPI_Comm comm )
//{
//
//    newNumberingNodesFaces* nnf = new newNumberingNodesFaces;
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//    std::map<int,std::vector<int> >::iterator prit;
//
//
//
//    int testElem = 9;
//    std::map<int,std::vector<int> > ref2bcface;
//    int lf  = 0;
//    int shf = 0;
//    std::map<int,int> sharedFaces;
//    std::set<int> gfid_set;
//    std::map<int,int> gl_map;
//    std::set<int> glf_set;
//    std::map<int,int> glf_map;
//    std::map<int,int> glf_map_inv;
//
//    std::map<int,int> glf_map_sha;
//    std::map<int,int> glf_map_sha_inv;
//    int ref = 0;
//    std::map<int,int> face2ref_prism;
//    int nLocalFaces = 0;
//    int nLocalVerts = 0;
//    std::set<int> ufaces;
//    std::set<int> gvid_set;
//    int gvid,gfid,r0,r1,el0,el1;
//    int elTel = 0;
//    std::map<int,int> loc2globV;
//    int foundshellvert = 0;
//    std::map<int,int> unf;
//
//    int changed = 0;
//    int unchanged = 0;
//    std::map<int,int* > ifn_updt;
//    int ell = 0;
//
//    for(prit=elements.begin();prit!=elements.end();prit++)
//    {
//        //key = GID, value global node nmber;
//        int gEl = prit->first;
//
//        for(int j=0;j<ief_part_map[gEl].size();j++)
//        {
//            gfid = ief_part_map[gEl][j];
//
//            if(ushell.find(gfid)!=ushell.end())
//            {
//                ref = 13;
//            }
//            else
//            {
//                ref = if_ref_part_map[gfid][0];
//            }
//
//            int nppf = if_Nv_part_map[gfid][0];
//
//            for(int h=0;h<nppf;h++)
//            {
//                gvid = ifn_part_map[gfid][h];
//
//                if(gvid_set.find(gvid)==gvid_set.end())
//                {
//                    gvid_set.insert(gvid);
//                }
//
//            }
//
//            if(gfid_set.find(gfid)==gfid_set.end())
//            {
//                gfid_set.insert(gfid);
//                face2ref_prism[gfid] = ref;
//                nLocalFaces++;
//            }
//
//            if(ufaces.find(gfid)==ufaces.end())
//            {
//                ufaces.insert(gfid);
//
//                el0    = ife_part_map[gfid][0];
//                el1    = ife_part_map[gfid][1];
//
//                if(ref==2)
//                {
//                    r0 = part_global->getVal(el0,0);
//                    r1 = part_global->getVal(el1,0);
//
//                    if(r0==world_rank && r1!=world_rank)
//                    {
//                        sharedFaces[gfid] = r0;
//                        shf++;
//                    }
//                    if(r0!=world_rank && r1==world_rank)
//                    {
//                        sharedFaces[gfid] = r1;
//                        shf++;
//                    }
//                }
//            }
//        }
//        elTel++;
//    }
//
//
//    nLocalVerts = gvid_set.size();
//
//    int nSharedFaces                          = sharedFaces.size();
//    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
//
//    DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,comm);
//
//    int Nt_shFaces               = distSharedFaces->getNel();
//    int* shFace_offsets          = distSharedFaces->getOffsets();
//    int* shFace_nlocs            = distSharedFaces->getNlocs();
//    int* shFacesIDs              = new int[nSharedFaces];
//    int* shFaces_RankIDs         = new int[nSharedFaces];
//
//    int iter = 0;
//    std::set<int> UniqueSharedVertsOnRank_set;
//    std::vector<int> UniqueSharedVertsOnRank;
//    std::vector<int> UniqueSharedVertsOnRank_RankID;
//
//    std::map<int,int>::iterator itsf;
//    int lvrtid = 0;
//    int tel = shFace_offsets[world_rank];
//    for(itsf=sharedFaces.begin();itsf!=sharedFaces.end();itsf++)
//    {
//        shFacesIDs[iter] = itsf->first;
//        shFaces_RankIDs[iter] = itsf->second;
//        gfid = itsf->first;
//
//        for(int q=0;q<ifn_part_map[gfid].size();q++)
//        {
//            gvid   = ifn_part_map[gfid][q];
//
//            if(UniqueSharedVertsOnRank_set.find(gvid)==UniqueSharedVertsOnRank_set.end())
//            {
//                UniqueSharedVertsOnRank_set.insert(gvid);
//                UniqueSharedVertsOnRank.push_back(gvid);
//                UniqueSharedVertsOnRank_RankID.push_back(world_rank);
//                lvrtid++;
//            }
//        }
//        tel++;
//        iter++;
//    }
//
//    int nSharedVerts = UniqueSharedVertsOnRank.size();
//
//    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
//
//    int Nt_shVerts               = distSharedVerts->getNel();
//    int* shVerts_nlocs           = distSharedVerts->getNlocs();
//    int* shVerts_offsets         = distSharedVerts->getOffsets();
//
//    int* TotalSharedVerts        = new int[Nt_shVerts];
//    int* TotalSharedVerts_RankID = new int[Nt_shVerts];
//    int* TotalSharedFaces        = new int[Nt_shFaces];
//    int* TotalSharedFaces_RankID = new int[Nt_shFaces];
//
//    // Communicate vert map to all ranks.
//    MPI_Allgatherv(&UniqueSharedVertsOnRank[0],
//                   nSharedVerts,
//                   MPI_INT,
//                   TotalSharedVerts,
//                   shVerts_nlocs,
//                   shVerts_offsets,
//                   MPI_INT, comm);
//
//    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
//                   nSharedVerts,
//                   MPI_INT,
//                   TotalSharedVerts_RankID,
//                   shVerts_nlocs,
//                   shVerts_offsets,
//                   MPI_INT, comm);
//
//    // Communicate face map to all ranks.
//    MPI_Allgatherv(shFacesIDs,
//                   nSharedFaces,
//                   MPI_INT,
//                   TotalSharedFaces,
//                   shFace_nlocs,
//                   shFace_offsets,
//                   MPI_INT, comm);
//
//    MPI_Allgatherv(shFaces_RankIDs,
//                   nSharedFaces,
//                   MPI_INT,
//                   TotalSharedFaces_RankID,
//                   shFace_nlocs,
//                   shFace_offsets,
//                   MPI_INT, comm);
//
//    int tmp;
//    std::map<int,int> f2r;
//    std::set<int> f2r_s;
//
//
//    int* NewGlobVertCountPerRank = new int[world_size];
//    int* NewGlobFaceCountPerRank = new int[world_size];
//
//    for(int u=0;u<world_size;u++)
//    {
//        NewGlobVertCountPerRank[u] = 0;
//        NewGlobFaceCountPerRank[u] = 0;
//    }
//
//    for(int i=0;i<Nt_shFaces;i++)
//    {
//        int key = TotalSharedFaces[i];
//        int val = TotalSharedFaces_RankID[i];
//
//        if(f2r_s.find(key)==f2r_s.end())
//        {
//            f2r_s.insert(key);
//            f2r[key]=val;
//            NewGlobFaceCountPerRank[val]=NewGlobFaceCountPerRank[val]+1;
//        }
//        else
//        {
//            tmp = f2r[key];
//            if(val<tmp)
//            {
//                f2r[key]=val;
//            }
//            if(val>tmp)
//            {
//                f2r[key]=tmp;
//            }
//        }
//    }
//
//    std::map<int,int> v2r;
//    std::set<int> v2r_s;
//    int* NewGlobVertOffset = new int[world_size];
//    int* NewGlobFaceOffset = new int[world_size];
//    int* owned_verts       = new int[world_size];
//    for(int u=0;u<world_size;u++)
//    {
//        owned_verts[u]       = 0;
//        NewGlobVertOffset[u] = 0;
//        NewGlobFaceOffset[u] = 0;
//    }
//    for(int i=0;i<Nt_shVerts;i++)
//    {
//        int key = TotalSharedVerts[i];
//        int val = TotalSharedVerts_RankID[i];
//
//        if(v2r_s.find(key)==v2r_s.end())
//        {
//            v2r_s.insert(key);
//            v2r[key]=val;
//            NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
//        }
//        else
//        {
//            tmp = v2r[key];
//
//            if(val<tmp)
//            {
//                v2r[key]=val;
//                owned_verts[val]=owned_verts[val]+1;
//
//            }
//            if(val>tmp)
//            {
//                v2r[key]=tmp;
//                owned_verts[tmp]=owned_verts[tmp]+1;
//            }
//        }
//    }
//
//
//
//    for(int u=1;u<world_size;u++)
//    {
//        NewGlobVertOffset[u] = NewGlobVertOffset[u-1]+NewGlobVertCountPerRank[u-1];
//        NewGlobFaceOffset[u] = NewGlobFaceOffset[u-1]+NewGlobFaceCountPerRank[u-1];
//    }
//
//    std::map<int,std::vector<int> >::iterator itv;
//
//    std::map<int,int> sharedVertsGlobal;
//
//    for(int u=0;u<world_size;u++)
//    {
//        NewGlobVertCountPerRank[u] = 0;
//        NewGlobFaceCountPerRank[u] = 0;
//    }
//
////    int nNonSharedVerts          = nLocalVerts-owned_verts[world_rank];
////    int nNonSharedFaces          = nLocalFaces-owned_faces[world_rank];
//
//    int nNonSharedVerts          = nLocalVerts-nSharedVerts;
//    int nNonSharedFaces          = nLocalFaces-nSharedFaces;
//
//    //std::cout << "LOCALVERTS " << world_rank << "   -- = " << nNonSharedVerts <<" "<< nLocalVerts << " " << owned_verts[world_rank] << " " << nSharedVerts << std::endl;
//
//    DistributedParallelState* nonSharedVertDistr = new DistributedParallelState(nNonSharedVerts,comm);
//    DistributedParallelState* nonSharedFaceDistr = new DistributedParallelState(nNonSharedFaces,comm);
//
//    int nNonSharedVertsTot = nonSharedVertDistr->getNel();
//
//    int iVshared           = nonSharedVertDistr->getNel()-shellvert2ref_glob.size();
//    int iVshared_bef = iVshared;
//    std::map<int,int >::iterator itvv;
//    std::map<int,int> sharedVmap;
//    int ownedsharedVmap = 0;
//    std::vector<int> v2r_owned;
//    for(itvv=v2r.begin();itvv!=v2r.end();itvv++)
//    {
//        sharedVmap[itvv->first] = iVshared;
//        if(itvv->second==world_rank)
//        {
//            v2r_owned.push_back(itvv->first);
//            ownedsharedVmap++;
//        }
//        iVshared++;
//    }
//
//    std::cout << "v2r_owned.size() " << v2r_owned.size() << std::endl;
////    int nVertTot_loc = 0;
////    int nVertTot_tot = 0;
////
////    if(world_rank == world_size-1)
////    {
////        nVertTot_tot =  iVshared;
////    }
////    MPI_Allreduce(&nVertTot_loc, &nVertTot_tot, 1, MPI_INT, MPI_MAX, comm);
//
//    //std::cout << "World "  << nVertTot_tot << " " << world_rank << std::endl;
//
//    std::map<int,int> sharedFmap;
//    int iFshared = distLocalFaces->getNel()-f2r.size();
//
//    for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
//    {
//        sharedFmap[itvv->first] = iFshared;
//        iFshared++;
//    }
//
//    int lbvid         = nonSharedVertDistr->getOffsets()[world_rank];
//    int lbfid         = nonSharedFaceDistr->getOffsets()[world_rank];
//
//    //std::cout << "===checking values = " << " "<< " world_rank  " << world_rank  << " " << nLocalVerts << " " << nSharedVerts << " " << ownedsharedVmap << " " << nNonSharedVerts << " " << unf.size() <<  std::endl;
//
//    DistributedParallelState* ElementDistr = new DistributedParallelState(elements.size(),comm);
//
//    int offsets_element = ElementDistr->getOffsets()[world_rank];
//    std::map<int,std::vector<int> > ifref;
//    std::map<int,std::vector<int> > ifn;
//    int uloc = 0;
//    int u    = 0;
//    std::map<int,std::vector<int> > ief_new;
//    std::map<int,int> pface2ref;
//    std::map<int,int> lhp;
//    std::map<int,int> rhp;
//    std::map<int,int> tag2ElementID;
//    std::map<int,std::vector<int> > colorRh;
//    std::map<int,std::vector<int> > colorFh;
//    std::map<int,int> tagE2gE;
//    std::map<int,int> gE2tagE;
//    std::map<int,int> shF2El;
//    std::map<int,int> lhp_sha;
//    std::map<int,int> rhp_sha;
//    Array<int>* iet = new Array<int>(elements.size(),1);
//
//
//    std::map<int,int> owned_glov2tagv_NEW;
//    std::map<int,int> owned_tagv2glov_NEW;
//    std::map<int,int> owned_glov2tagv;
//    std::map<int,int> owned_tagv2glov;
//    std::map<int,Vert*> tag2coords;
//
//    std::map<int,std::vector<int> > shellFace2Node_prism;
//    std::map<int,std::vector<int> > bcFace2Node_prism;
//    std::map<int,std::vector<int> > sharedFace2Node_prism;
//    std::map<int,std::vector<int> > intFace2Node_prism;
//
//    int fowned = 0;
//    int fownedInmap = 0;
//    double orient0 = 1.0;
//    std::map<int,double> orient_map;
//    int ngf = 0;
//
//    for(prit=elements.begin();prit!=elements.end();prit++)
//    {
//        int gEl     = prit->first;
//        int rank    = world_rank;
//        int lEl     = ElementDistr->getOffsets()[rank]+u+1;
//        tagE2gE[gEl]  = lEl;
//        gE2tagE[lEl]  = gEl;
//        iet->setVal(u,0,6);
//
//        // Compute the center of the element in order;
//        Vert* Vijk = new Vert;
//        Vijk->x = 0.0;
//        Vijk->y = 0.0;
//        Vijk->z = 0.0;
//        // compute element center;
//        int nvp = prit->second.size();
//        for(int q=0;q<nvp;q++)
//        {
//            int tag  = prit->second[q];
//            int lvp  = tag2locV[tag];
//
//            Vijk->x = Vijk->x + locVerts[lvp]->x;
//            Vijk->y = Vijk->y + locVerts[lvp]->y;
//            Vijk->z = Vijk->z + locVerts[lvp]->z;
//        }
//
//        Vijk->x = Vijk->x/nvp;
//        Vijk->y = Vijk->y/nvp;
//        Vijk->z = Vijk->z/nvp;
//
//        std::vector<int> newFids(ief_part_map[gEl].size());
//
////        if(world_rank == 2)
////        {
////            std::cout <<"=======" << std::endl;
////
////            std::cout <<ief_part_map[gEl].size() << " ---> ";
////        }
//
//        for(int q=0;q<ief_part_map[gEl].size();q++)
//        {
//            gfid = ief_part_map[gEl][q];
//            int nfvrts = if_Nv_part_map[gfid][0];
//            // Setting the reference to 13 for the shell faces;
//            if(ushell.find(gfid)!=ushell.end())
//            {
//                ref                 = 13;
//                tag2ElementID[gEl]  = lEl;
//            }
//            else
//            {
//                ref                 = if_ref_part_map[gfid][0];
//            }
//
//            if(f2r.find(gfid)!=f2r.end() && ref != 13)
//            {
//                if(f2r[gfid]==world_rank)
//                {
//                    fowned++;
//
//                    if(glf_map.find(gfid)==glf_map.end())
//                    {
//                        fownedInmap++;
//                        int lfid2               = sharedFmap[gfid];
//                        glf_map[gfid]           = lfid2;
//                        glf_map_inv[lfid2]      = gfid;
//                        lhp[gfid]               = lEl;
//
//                        el0                     = ife_part_map[gfid][0];
//                        el1                     = ife_part_map[gfid][1];
//
//                        r0                      = part_global->getVal(el0,0);
//                        r1                      = part_global->getVal(el1,0);
//
//                        if(r0==world_rank && r1!=world_rank)
//                        {
//                            colorRh[r1].push_back(el1);
//                            colorFh[r1].push_back(gfid);
//                        }
//                        else if(r1==world_rank && r0!=world_rank)
//                        {
//                            colorRh[r0].push_back(el0);
//                            colorFh[r0].push_back(gfid);
//                        }
//
//                        std::vector<int> fn_glo(nfvrts);
//                        std::vector<int> fn_tag(nfvrts);
//
//                        Vert* VcF = new Vert;
//                        std::vector<Vert*> Vfaces;
//
//                        for(int n=0;n<nfvrts;n++)
//                        {
//                            fn_tag[n] = ifn_part_map[gfid][n];
//                            fn_glo[n] = sharedVmap[fn_tag[n]];
//
//                            int lvp  = tag2locV[fn_tag[n]];
//                            Vert* Vf = new Vert;
//                            Vf->x = locVerts[lvp]->x;
//                            Vf->y = locVerts[lvp]->y;
//                            Vf->z = locVerts[lvp]->z;
//
//                            VcF->x = VcF->x + Vf->x;
//                            VcF->y = VcF->y + Vf->y;
//                            VcF->z = VcF->z + Vf->z;
//
//                            Vfaces.push_back(Vf);
//
//                            if(v2r.find(fn_tag[n])!=v2r.end())
//                            {
//                                if(v2r[fn_tag[n]]==world_rank)
//                                {
//                                    if(owned_tagv2glov_NEW.find(fn_tag[n])==owned_tagv2glov_NEW.end() &&
//                                       shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                    {
//                                        owned_tagv2glov_NEW[fn_tag[n]]  = fn_glo[n];
//                                        owned_glov2tagv_NEW[fn_glo[n]]  = fn_tag[n];
//
//                                        int lvpart                      = tag2locV[fn_tag[n]];
//                                        Vert* vv                        = new Vert;
//                                        vv->x                           = locVerts[lvpart]->x;
//                                        vv->y                           = locVerts[lvpart]->y;
//                                        vv->z                           = locVerts[lvpart]->z;
//                                        tag2coords[fn_tag[n]]           = vv;
//                                    }
//                                }
//                            }
//                        }
//
//                        // check orientation!!
//
//                        orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
//                        int sw = 0;
//
//                        orient_map[gfid] = orient0;
//
//                        if(orient0<0.0)
//                        {
//                            //std::cout << "negative " << world_rank << " " << nfvrts << std::endl;
//                            ngf++;
//
//                            if(nfvrts == 3)
//                            {
//                                fn_tag[0] = ifn_part_map[gfid][0];
//                                fn_tag[1] = ifn_part_map[gfid][2];
//                                fn_tag[2] = ifn_part_map[gfid][1];
//                            }
//                            if(nfvrts == 4)
//                            {
//                                fn_tag[0] = ifn_part_map[gfid][0];
//                                fn_tag[1] = ifn_part_map[gfid][3];
//                                fn_tag[2] = ifn_part_map[gfid][2];
//                                fn_tag[3] = ifn_part_map[gfid][1];
//
//                                sw=1;
//                            }
//
//                            sw=1;
//                        }
//
//
//                        sharedFace2Node_prism[gfid] = fn_tag;
//                    }
//                    else
//                    {
//                        rhp[gfid] = lEl;
//                    }
//                }
//            }
//            else
//            {
//                if(glf_map.find(gfid)==glf_map.end())
//                {
//                    glf_map[gfid]      = lbfid;
//                    glf_map_inv[lbfid] = gfid;
//                    int nfvrts         = if_Nv_part_map[gfid][0];
//
//                    std::vector<int> fn_tag(nfvrts);
//                    for(int n=0;n<nfvrts;n++)
//                    {
//                        fn_tag[n] = ifn_part_map[gfid][n];
//
//                        int lvp  = tag2locV[fn_tag[n]];
//
//                        tag2coords[fn_tag[n]] = Vf;
//
//                        if(v2r.find(fn_tag[n])!=v2r.end())
//                        {
//                            fn_glo[n] = sharedVmap[fn_tag[n]];
//
//                            if(v2r[fn_tag[n]]==world_rank)
//                            {
//                                if(owned_tagv2glov_NEW.find(fn_tag[n])==owned_tagv2glov_NEW.end() &&
//                                   shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                {
//                                    owned_tagv2glov_NEW[fn_tag[n]]  = fn_glo[n];
//                                    owned_glov2tagv_NEW[fn_glo[n]]  = fn_tag[n];
//                                }
//                                else if(owned_tagv2glov_NEW.find(fn_tag[n])!=owned_tagv2glov_NEW.end() &&
//                                        shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                                {
//                                    int gvid2 = owned_tagv2glov_NEW[fn_tag[n]];
//                                    fn_glo[n] = gvid2;
//                                }
//                            }
//                        }
//                        else
//                        {
//                            if(owned_tagv2glov_NEW.find(fn_tag[n])==owned_tagv2glov_NEW.end() &&
//                               shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                            {
//                                fn_glo[n]                       = gvid;
//                                owned_tagv2glov_NEW[fn_tag[n]]  = fn_glo[n];
//                                owned_glov2tagv_NEW[fn_glo[n]]  = fn_tag[n];
//
//                                gvid++;
//                            }
//                            else if(owned_tagv2glov_NEW.find(fn_tag[n])!=owned_tagv2glov_NEW.end() &&
//                                    shellvert2ref_glob.find(fn_tag[n])==shellvert2ref_glob.end())
//                            {
//                                int gvid2 = owned_tagv2glov_NEW[fn_tag[n]];
//                                fn_glo[n] = gvid2;
//                            }
//                            else if(shellvert2ref_glob.find(fn_tag[n])!=shellvert2ref_glob.end())
//                            {
//                                fn_glo[n] = -1;
//                            }
//                        }
//                    }
//
//
//
//
//
//
//
//
//                    if(ref == 2)
//                    {
//                        intFace2Node_prism[gfid]    = fn_tag;
//                    }
//                    if(ref != 2 && ref != 13)
//                    {
//                        ref2bcface[ref].push_back(gfid);
//
//                        bcFace2Node_prism[gfid]     = fn_tag;
//                    }
//                    if(ref == 13)
//                    {
//                        shellFace2Node_prism[gfid]  = fn_tag;
//                    }
//
//                    lhp[gfid]          = lEl;
//
//                    if(lEl == testElem)
//                    {
//                        std::cout << "lhp[gfid] " << gfid << " " << ref << " " << orient0 << std::endl;
//                    }
//
//                    lbfid              = lbfid + 1;
//                }
//                else
//                {
//                    //int lbbfid  = glf_map[gfid];
//                    rhp[gfid]           = lEl;
//                    double orient = orient_map[gfid];
//                    if(lEl == testElem)
//                    {
//                        std::cout << "rhp[gfid] " << gfid << " " << ref << " " << orient << std::endl;
//                    }
//                }
//            }
//        }
////        if(world_rank == 2)
////        {
////            std::cout <<"=======" << std::endl;
////        }
//        //ief_new[gEl] = newFids;
//        uloc++;
//        u++;
//
//    }
//
//    std::cout << "ngf " << ngf << " ON WORKRANK " << world_rank << std::endl;
//    //std::cout << "facemap sizes " << world_rank << " " << intFace2Node_prism.size() << " " << bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << std::endl;
//    //std::cout << " before " << world_rank << " :: " << lhp.size() << " " << rhp.size() + bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << std::endl;
//
//    ScheduleObj* rh_schedule = DoScheduling(colorRh,comm);
//    std::map<int,std::vector<int> > recv_rhElIDs;
//    for(int q=0;q<world_size;q++)
//    {
//        if(world_rank==q)
//        {
//            int i=0;
//
//            for (itv = colorRh.begin(); itv != colorRh.end(); itv++)
//            {
//                int n_req           = itv->second.size();
//                int dest            = itv->first;
//
//                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
//                MPI_Send(&itv->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
//                i++;
//            }
//        }
//        else if (rh_schedule->SendFromRank2Rank[q].find( world_rank ) != rh_schedule->SendFromRank2Rank[q].end())
//        {
//            int n_reqstd_ids;
//            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);
//
//            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
//
//            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
//            recv_rhElIDs[q] = recv_reqstd_ids;
//        }
//    }
//
//    std::map<int,std::vector<int> > sendEl;
//    std::map<int,std::vector<int> >::iterator rcvit;
//    for(rcvit=recv_rhElIDs.begin();rcvit!=recv_rhElIDs.end();rcvit++)
//    {
//        int frank = rcvit->first;
//        int nE    = rcvit->second.size();
//
//        for(int j=0;j<nE;j++)
//        {
//            int gEl = tagE2gE[rcvit->second[j]];
//            sendEl[frank].push_back(gEl);
//        }
//    }
//
//    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);
//
//    std::map<int,std::vector<int> > adj_ids;
//    for(int q=0;q<world_size;q++)
//    {
//        if(world_rank==q)
//        {
//            int i=0;
//            for (itv = sendEl.begin(); itv != sendEl.end(); itv++)
//            {
//                int n_req           = itv->second.size();
//                int dest            = itv->first;
//
//                MPI_Send(&n_req, 1,
//                        MPI_INT, dest,
//                        6798+78000*dest, comm);
//
//                MPI_Send(&itv->second[0],
//                        n_req, MPI_INT,
//                        dest, 14876000+dest, comm);
//
//                i++;
//            }
//        }
//        else if (ishBack_schedule->SendFromRank2Rank[q].find( world_rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
//        {
//            int n_reqstd_ids;
//
//            MPI_Recv(&n_reqstd_ids,
//                    1, MPI_INT, q,
//                    6798+78000*world_rank,
//                    comm, MPI_STATUS_IGNORE);
//
//            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
//
//            MPI_Recv(&recv_reqstd_ids[0],
//                    n_reqstd_ids,
//                    MPI_INT, q,
//                    14876000+world_rank,
//                    comm, MPI_STATUS_IGNORE);
//
//            adj_ids[q] = recv_reqstd_ids;
//
//        }
//    }
//    std::map<int,int> shFid2el_rh;
//    int adde = 0;
//    int adde2 = 0;
//    for(rcvit=adj_ids.begin();rcvit!=adj_ids.end();rcvit++)
//    {
//        int rrank = rcvit->first;
//        int nelem = rcvit->second.size();
//
//        //std::cout << "wr " << world_rank << " sending " << rrank << " " << nelem << std::endl;
//
//        for(int q=0;q<nelem;q++)
//        {
//            int tag = colorRh[rrank][q];
//            int fid = colorFh[rrank][q];
//
//            if(shFid2el_rh.find(fid)==shFid2el_rh.end())
//            {
//                shFid2el_rh[fid] = rcvit->second[q];
//                adde2++;
//                if(rhp.find(fid)==rhp.end())
//                {
//                    adde++;
//                    rhp[fid] = rcvit->second[q];
//                }
//            }
//        }
//    }
//
//    //std::cout << " after " << world_rank << " :: " << lhp.size() << " " << rhp.size() + bcFace2Node_prism.size() << " " << shellFace2Node_prism.size() << " adde " << adde << " " << adde2 << std::endl;
//
//    std::map<int,std::vector<int> > face2node;
//
//    int lfid;
//    std::map<int,int>::iterator pritm;
//    int ff13 = 0;
//
//    std::set<int> excludeShellV;
//
//    int nUniqueVerts_prisms = owned_tagv2glov_NEW.size();
//
//    DistributedParallelState* distri_prism_verts  = new DistributedParallelState(nUniqueVerts_prisms,comm);
//
//    int NtotVerts_prism = distri_prism_verts->getNel();
//
//    //std::cout << "Number of total unique " << world_rank << " " << NtotVerts_prism << " face sizes " << intFace2Node_prism.size() << " " << bcFace2Node_prism.size() << " " << sharedFace2Node_prism.size() << " " << shellFace2Node_prism.size() << " " << lhp.size() << " " << rhp.size() << " " << lhp_sha.size() << " " << rhp_sha.size() << " " << adde << std::endl;
//
//    int lvg             = distri_prism_verts->getOffsets()[world_rank]+1;
//    int lvl             = 0;
//
//    std::map<int,int>::iterator itmp;
//    std::map<int,int> loc2glob_prism;
//    std::map<int,int> loc2glob_prism_realz;
//    Array<double>* xcn_parmmg_prism = new Array<double>(nUniqueVerts_prisms,3);
//    //std::cout << "world Rank " << world_rank << " " << gl_map_i.size() << " " << gl_map.size() << " vakjes " << vak1.size() << " " << vak2.size() << std::endl;
//
//    //owned_tagv2glov_NEW
//
//    for(itmp=owned_tagv2glov_NEW.begin();itmp!=owned_tagv2glov_NEW.end();itmp++)
//    {
//        int tag     = itmp->first;
//
//        double xc = tag2coords[tag]->x;
//        double yc = tag2coords[tag]->y;
//        double zc = tag2coords[tag]->z;
//
//        xcn_parmmg_prism->setVal(lvl,0,xc);
//        xcn_parmmg_prism->setVal(lvl,1,yc);
//        xcn_parmmg_prism->setVal(lvl,2,zc);
//
//        loc2glob_prism[tag]       = lvg;
//        loc2glob_prism_realz[tag] = lvl;
//
//        lvl++;
//        lvg++;
//    }
//
//
//
//    //std::cout << " Glboal nTotUniqueNonShellVerts = " << world_rank << " PPPP "  << nonSharedVertDistr->getNel()+v2r.size() << " " << gl_map.size() << " " << gl_map_i.size() << " gl_map_int " << gl_map_int.size() << " " <<distrimap_gl_map->getNel() << " " << lbvid-1 << " " << tag2coords.size() << std::endl;
//    //std::cout << "Worl d " << distLocalVerts->getNel() <<  " " << nonSharedVertDistr->getNel() <<  " " << v2r.size() << " iVs " << iVshared << " " << ownedsharedVmap << " " << loc2glob_prism.size() <<  std::endl;
//
//    int nTotF = sharedFace2Node_prism.size()+intFace2Node_prism.size()+bcFace2Node_prism.size()+shellFace2Node_prism.size();
//
////    std::cout << "WORLD " << world_rank << " :: lh = " << lhp.size() << " " << rhp.size() << " face2nodes = "
////            << face2node.size() << " =? " << nTotF << " sharedFaces  " << sharedFmap.size() << " vs " << sharedFace2Node_prism.size() << " "
////            << intFace2Node_prism.size() << " " << bcFace2Node_prism.size() << " :: " << glf_map.size() << " "
////            << shellFace2Node_prism.size() << " " << f2r.size()  << " " << shF2El.size() << std::endl;
////
////
////    std::cout << "lhsh vs rhsh " << world_rank << " " << lhp.size() << " " << rhp.size() << " " << bcFace2Node_prism.size() << " results in " << lhp.size()-rhp.size()-bcFace2Node_prism.size() << " ---- " << shellFace2Node_prism.size() << std::endl;
//
//    //std::cout << tag2ElementID.size() << " tag2ElementID.size() " << elements.size() << std::endl;
//
//
//    int mapSizeLoc = tag2ElementID.size();
//    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,comm);
//
//    int mapSizeTot = distrimap->getNel();
//    int* tag_loc = new int[mapSizeLoc];
//    int* eid_loc = new int[mapSizeLoc];
//    int* tag_tot = new int[mapSizeTot];
//    int* eid_tot = new int[mapSizeTot];
//
//    int tvid_tmp,iref_tmp;
//    int i = 0;
//    std::map<int,int>::iterator itred;
//
//    for(itred=tag2ElementID.begin();itred!=tag2ElementID.end();itred++)
//    {
//        tag_loc[i] = itred->first;
//        eid_loc[i] = itred->second;
//        i++;
//    }
//
//    int* offsets = distrimap->getOffsets();
//    int* nlocs   = distrimap->getNlocs();
//
//
//    MPI_Allgatherv(tag_loc,
//                   mapSizeLoc,
//                   MPI_INT,
//                   tag_tot,
//                   nlocs,
//                   offsets,
//                   MPI_INT, comm);
//
//
//    MPI_Allgatherv(eid_loc,
//                   mapSizeLoc,
//                   MPI_INT,
//                   eid_tot,
//                   nlocs,
//                   offsets,
//                   MPI_INT, comm);
//
//
//
//    std::map<int,int> tag2ElementID_tot;
//    int key,val;
//    for(int i=0;i<mapSizeTot;i++)
//    {
//        key = tag_tot[i];
//        val = eid_tot[i];
//        if(tag2ElementID_tot.find(key)==tag2ElementID_tot.end())
//        {
//            tag2ElementID_tot[key] = val;
//        }
//    }
//
//    int lvert = 0;
//
//    // testing site!
//
//
//
////    for(prit=elements.begin();prit!=elements.end();prit++)
////    {
////        int gEl     = prit->first;
////        int rank    = world_rank;
////        int lEl     = ElementDistr->getOffsets()[rank]+u+1;
////        iet->setVal(u,0,6);
////        tagE2gE[gEl]  = lEl;
////
////        std::vector<int> newFids(ief_part_map[gEl].size());
////
////        for(int q=0;q<ief_part_map[gEl].size();q++)
////        {
////            gfid = ief_part_map[gEl][q];
////
////
////
////
////        }
////    }
////
//
//    nnf->xcn                        = xcn_parmmg_prism;
//    nnf->tag2coords                 = tag2coords;
//    nnf->lh                         = lhp;
//    nnf->rh                         = rhp;
//    nnf->tagE2gE                    = tagE2gE;
//    nnf->gE2tagE                    = gE2tagE;
//    nnf->lh_sha                     = lhp_sha;
//    nnf->rh_sha                     = rhp_sha;
//    nnf->nTotUniqueNonShellVerts    = iVshared-1;
//    nnf->shellFace2Node             = shellFace2Node_prism;
//    nnf->sharedvert2rank            = v2r;
//    nnf->local2globalVertMap        = loc2glob_prism;
//    nnf->local2globalVertMapReal    = loc2glob_prism_realz;
//    nnf->sharedVmap                 = sharedVmap;
//
//    nnf->sharedFace2Node            = sharedFace2Node_prism;
//    nnf->intFace2Node               = intFace2Node_prism;
//    nnf->bcFace2Node                = bcFace2Node_prism;
//
//    nnf->glob2locF                  = glf_map;
//    nnf->loc2globF                  = glf_map_inv;
//    nnf->glob2locF_sha              = glf_map_sha;
//    nnf->loc2globF_sha              = glf_map_sha_inv;
//    nnf->face2node                  = face2node;
//    nnf->ifref                      = pface2ref;
//    nnf->ien                        = elements;
//    nnf->ief                        = ief_new;
//    nnf->dist                       = ElementDistr;
//    nnf->pbcmap                     = ref2bcface;
//    nnf->tag2ElementID              = tag2ElementID_tot;
//    nnf->localV2tagV                = owned_glov2tagv_NEW;
//    nnf->tagV2localV                = owned_tagv2glov_NEW;
////    nnf->localV2tagV_all            = glov2tagv_all;
////    nnf->tagV2localV_all            = tagv2glov_all;
//    nnf->iet                       = iet;
//    /**/
//    return nnf;
//    //===========================================================================
//    //===========================================================================
//    //===========================================================================
//
//}
