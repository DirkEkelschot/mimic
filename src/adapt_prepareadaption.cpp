#include "adapt_prepareadaption.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"

PrepareAdaption::PrepareAdaption(PartObject* part, 
                                MPI_Comm comm,
                                int                                             neglob,
                                std::map<int,std::vector<double> >&&            LocalVertsMap,
                                std::map<int,int>&&                             localV2globalV,
                                std::map<int,std::vector<int> >&&               Elem2Face,
                                std::map<int,std::vector<int> >&&               Elem2Vert,
                                std::map<int,std::vector<int> >&&               Face2Vert,
                                std::map<int,std::vector<int> >&&               Face2Elem,
                                std::map<int,int>&&                             Elem2Rank,
                                std::set<int>                                   TraceVertsOnRank,
                                std::set<int>                                   TraceFacesOnRank,
                                std::set<int>&&                                 ElemSet,
                                std::map<int,std::vector<int> >                 ranges_id,
                                std::map<int,std::vector<int> >                 ranges_ref)
                                :   m_ElemSet(std::move(ElemSet)), 
                                    m_Elem2Face(std::move(Elem2Face)),
                                    m_Elem2Vert(std::move(Elem2Vert)),
                                    m_Face2Vert(std::move(Face2Vert)),
                                    m_Face2Elem(std::move(Face2Elem)),
                                    m_Elem2Rank(std::move(Elem2Rank)),
                                    m_locV2gloV(std::move(localV2globalV)),
                                    m_LocalVertsMap(std::move(LocalVertsMap))
        {
    //int faceidtrace = 2353748;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    Ne_glob = neglob;


    buildUpdatedVertexAndFaceNumbering(comm,
                                        TraceVertsOnRank,
                                        TraceFacesOnRank,
                                        ranges_id,
                                        ranges_ref);

    buildInteriorSharedAndBoundaryFaceMaps(comm,
                                        TraceVertsOnRank,
                                        TraceFacesOnRank,
                                        ranges_id,
                                        ranges_ref);

}



void PrepareAdaption::buildUpdatedVertexAndFaceNumbering(MPI_Comm                   comm, 
                                                    std::set<int>                   TraceVertsOnRank,
                                                    std::set<int>                   TraceFacesOnRank,
                                                    std::map<int,std::vector<int> > ranges_id,
                                                    std::map<int,std::vector<int> > ranges_ref)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::set<int> TraceVertsOnRank_glob = AllGatherSet(TraceVertsOnRank,comm);

    int Nel_loc  = m_ElemSet.size();
    int ref;
    
    std::set<int> gvid_set;
    int inf = 0;
    int shf = 0;
    int tracef = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::map<int,int> sharedFace2Rank;
    std::set<int> VertsOnRank;
    int bface = 0;

    std::set<int> refsc;

    std::map<int,int> sharedFaces;
    std::map<int,int> interiorFaces;
    std::map<int,int> boundaryFaces;

    std::set<int>::iterator its;
    int ntrace = 0;
    int teller = 0;
    int teller2 = 0;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        int elid = *its;
        int Nf = m_Elem2Face[elid].size();

        for(int j=0;j<Nf;j++)
        {
            int gfid = m_Elem2Face[elid][j];
            
            if(TraceFacesOnRank.find(gfid)!=TraceFacesOnRank.end())
            {
                ref = 13;
                ntrace++;
            }
            else
            {
                
                if(m_Face2Vert.find(gfid)!=m_Face2Vert.end())
                {
                    ref         = ProvideBoundaryRef(gfid,ranges_ref);
                    
                    int e0      = m_Face2Elem[gfid][0];
                    int e1      = m_Face2Elem[gfid][1];
                    int Nv      = m_Face2Vert[gfid].size();
                    
                    // copy the correct amount of nodes into a new face vector.
                    std::vector<int> face(Nv,0);
                    for(int k=0;k<Nv;k++)
                    {
                        int vid = m_Face2Vert[gfid][k];
                        face[k] = vid;

                        if(VertsOnRank.find(vid)==VertsOnRank.end() &&
                        TraceVertsOnRank_glob.find(vid)==TraceVertsOnRank_glob.end())
                        {
                            VertsOnRank.insert(vid);
                        }
                    }
                    if(e0<Ne_glob && e1<Ne_glob)
                    {
                        
                        int r0; // = part_global[e0]; // rank of first adjacent element.
                        int r1; // = part_global[e1]; // rank of second adjacent element.

                        if(m_Elem2Rank.find(e0)!=m_Elem2Rank.end())
                        {
                            r0 = m_Elem2Rank[e0];
                        }
                        if(m_Elem2Rank.find(e1)!=m_Elem2Rank.end())
                        {
                            r1 = m_Elem2Rank[e1];
                        }

                        // int r0 = part_global[e0]; // rank of first adjacent element.
                        // int r1 = part_global[e1]; // rank of second adjacent element.

                        if(ref==2)// Internal and Shared faces are here.
                        {
                            //std::cout << "r0rank  " << r0 << " " << r1 << " :: " << rank << std::endl;
                            if(r0==rank && r1!=rank)
                            {
                                sharedFace2Rank[gfid] = r0;
                            }
                            if(r0!=rank && r1==rank)
                            {
                                sharedFace2Rank[gfid] = r1;
                            }
                            teller2++;
                        }
                        else
                        {
                            teller++;
                        }
                    }

                    
                }
            }
        }
    }

    // std::cout << "TraceFacesOnRank " << TraceFacesOnRank.size() << " " << rank << " " << sharedFace2Rank.size() << " " << teller << " " << teller2 << " " << m_Face2Elem.size()<< std::endl; 

    // Take care of shared faces and vertices
    std::set<int> UniqueSharedVertsOnRank_set;
    std::vector<int> UniqueSharedVertsOnRank_RankID;

    std::map<int,int>::iterator itmii;
    for(itmii=sharedFace2Rank.begin();itmii!=sharedFace2Rank.end();itmii++)
    {
        int gfid = itmii->first;
        int Nv   = m_Face2Vert[gfid].size();

        for(int q=0;q<Nv;q++)
        {
            int vid = m_Face2Vert[gfid][q];

            if(UniqueSharedVertsOnRank_set.find(vid)==UniqueSharedVertsOnRank_set.end() &&
                   TraceVertsOnRank_glob.find(vid)==TraceVertsOnRank_glob.end()
            )
            {
                UniqueSharedVertsOnRank_set.insert(vid);
                UniqueSharedVertsOnRank_RankID.push_back(rank);
            }

           
        }
    }

    int nLocalVertsTot = VertsOnRank.size();
    std::vector<int> nLocalVertsTot_loc_test(size,0);
    nLocalVertsTot_loc_test[rank] = nLocalVertsTot;
    std::vector<int> nLocalVertsTot_red_test(size,0);

    MPI_Allreduce(nLocalVertsTot_loc_test.data(), 
                  nLocalVertsTot_red_test.data(), 
                  size, MPI_INT, MPI_SUM, comm);

    int N_un_shareV = UniqueSharedVertsOnRank_set.size();
    int summednLocalVertsTot = 0;

    for(int i=0;i<size;i++)
    {
        summednLocalVertsTot = summednLocalVertsTot+nLocalVertsTot_red_test[i];
    }
    int testNv = summednLocalVertsTot - N_un_shareV;
    
    //=====================================================================================
    //================================SHARED FACES=========================================
    //=====================================================================================
    int nSharedFaces = sharedFace2Rank.size();
    int N_localFaces = m_Face2Elem.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
    DistributedParallelState* distLocalFaces  = new DistributedParallelState(N_localFaces,comm);

    int Nt_shFaces   = distSharedFaces->getNel();
    std::vector<int> sharedFaceIDs_tmp(sharedFace2Rank.size(),0);
    std::vector<int> sharedFaceRankIDs_tmp(sharedFace2Rank.size(),0);
    
    int c = 0;

    // Copy the face ID in the sharedface2mode map into temporary vector in order to communicate 
    
    for(itmii=sharedFace2Rank.begin();itmii!=sharedFace2Rank.end();itmii++)
    {
        sharedFaceIDs_tmp[c]     = itmii->first;
        sharedFaceRankIDs_tmp[c] = itmii->second;
        c++;
    }

    std::vector<int> TotalSharedFaces(Nt_shFaces,0);
    std::vector<int> TotalSharedFaces_RankID(Nt_shFaces,0);
    // Communicate face map to all ranks.
    MPI_Allgatherv(&sharedFaceIDs_tmp.data()[0],
                   nSharedFaces,
                   MPI_INT,
                   &TotalSharedFaces.data()[0],
                   distSharedFaces->getNlocs(),
                   distSharedFaces->getOffsets(),
                   MPI_INT, comm);
    
    MPI_Allgatherv(&sharedFaceRankIDs_tmp.data()[0],
                   nSharedFaces,
                   MPI_INT,
                  &TotalSharedFaces_RankID.data()[0],
                   distSharedFaces->getNlocs(),
                   distSharedFaces->getOffsets(),
                   MPI_INT, comm);

    int tmp;
    
    
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];

        m_SharedFace2Rank[key].push_back(val);

        if(m_sharedFace2RankMap.find(key)==m_sharedFace2RankMap.end())
        {
            m_sharedFace2RankMap[key]=val;
        }
        else
        {
            tmp = m_sharedFace2RankMap[key];
            if(val<tmp)
            {
                m_sharedFace2RankMap[key]=val;
            }
            if(val>tmp)
            {
                m_sharedFace2RankMap[key]=tmp;
            }
        }
    }

    
    //std::cout << "m_sharedFace2RankMap " << m_sharedFace2RankMap.size() << " " << rank << std::endl;   


    std::map<int,int> sharedFmap;
    int iFshared = distLocalFaces->getNel()-m_sharedFace2RankMap.size();
    std::map<int,int>::iterator itii;
    
    //=====================================================================================
    //===============================SHARED VERTEX=========================================
    //=====================================================================================
    int nSharedVerts = UniqueSharedVertsOnRank_set.size();
    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
    
    int Nt_shVerts = distSharedVerts->getNel();
    std::vector<int> TotalSharedVerts(Nt_shVerts,0);
    std::vector<int> TotalSharedVerts_RankID(Nt_shVerts,0);

    std::vector<int> UniqueSharedVertsOnRank_vec(UniqueSharedVertsOnRank_set.begin(),
                                                 UniqueSharedVertsOnRank_set.end());


    MPI_Allgatherv(&UniqueSharedVertsOnRank_vec[0],
                    nSharedVerts,
                    MPI_INT,
                   &TotalSharedVerts[0],
                    distSharedVerts->getNlocs(),
                    distSharedVerts->getOffsets(),
                    MPI_INT, comm);

    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
                    nSharedVerts,
                    MPI_INT,
                   &TotalSharedVerts_RankID[0],
                    distSharedVerts->getNlocs(),
                    distSharedVerts->getOffsets(),
                    MPI_INT, comm);


     
    for(int i=0;i<Nt_shVerts;i++)
    {
        int key = TotalSharedVerts[i];
        int val = TotalSharedVerts_RankID[i];
        
        if(m_sharedVertex2RankMap.find(key)==m_sharedVertex2RankMap.end())
        {
            m_sharedVertex2RankMap[key]=val;
        }
        else
        {
            tmp = m_sharedVertex2RankMap[key];
            
            if(val<tmp)
            {
                m_sharedVertex2RankMap[key]=val;
            }
            if(val>tmp)
            {
                m_sharedVertex2RankMap[key]=tmp;
            }
        }
    }

    std::map<int,int>::iterator itmm;
    int nOwnedSharedVerts = 0;
    for(itmm=m_sharedVertex2RankMap.begin();itmm!=m_sharedVertex2RankMap.end();itmm++)
    {
        if(itmm->second==rank)
        {
            nOwnedSharedVerts++;
        }
    }
    
    DistributedParallelState* ownedSharedVrtsDist   = new DistributedParallelState(nOwnedSharedVerts,comm);
    int nNonSharedVerts                             = nLocalVertsTot-nSharedVerts;
    DistributedParallelState* distLocalVerts        = new DistributedParallelState(nLocalVertsTot,comm);
    int nNonSharedFaces                             = N_localFaces-nSharedFaces;
    DistributedParallelState* nonSharedVertDistr    = new DistributedParallelState(nNonSharedVerts,comm);
    DistributedParallelState* nonSharedFaceDistr    = new DistributedParallelState(nNonSharedFaces,comm);


    int iVshared = distLocalVerts->getNel()-m_sharedVertex2RankMap.size();

    std::map<int,int >::iterator itvv;
    std::map<int,int> sharedVmap;
    for(itvv=m_sharedVertex2RankMap.begin();itvv!=m_sharedVertex2RankMap.end();itvv++)
    {
        sharedVmap[itvv->first] = iVshared;
        iVshared++;
    }


    std::vector<int> ownedVs(size,0);

    for(int i=0;i<size;i++)
    {
        ownedVs[i] = 0;
    }

    int nNonSharedVertsTot = nonSharedVertDistr->getNel();

    
    for(itmm=m_sharedVertex2RankMap.begin();itmm!=m_sharedVertex2RankMap.end();itmm++)
    {
        int globid = nNonSharedVertsTot+ownedSharedVrtsDist->getOffsets()[itmm->second]+ownedVs[itmm->second]+1;
        
        m_sharedVertexMapUpdatedGlobalID[itmm->first] = globid;
        
        if(itmm->second==rank &&
           TraceVertsOnRank_glob.find(itmm->first)==TraceVertsOnRank_glob.end())
        {
            m_SharedVertsOwned[itmm->first] = globid;
        }
        
        ownedVs[itmm->second] = ownedVs[itmm->second]+1;
    }

    
    int Fid_shared = distLocalFaces->getNel()-m_sharedFace2RankMap.size();

    for(itii=m_sharedFace2RankMap.begin();itii!=m_sharedFace2RankMap.end();itii++)
    {
        m_sharedFaceMapUpdatedGlobalID[itii->first] = Fid_shared;
        Fid_shared++;
    } 
    
    
    //std::cout << "m_sharedFace2RankMap " << m_sharedFace2RankMap.size() << " " << rank << std::endl;   

    DistributedParallelState* ElementDistr  = new DistributedParallelState(m_ElemSet.size(),comm);
    int u                                   = 0;
    int gvidd                               = 0;

    std::map<int,int> tag2ElementID;
    std::set<int> FaceOnRank;
    std::vector<int> nNonSharedFacesArray(size,0);
    std::vector<int> nNonSharedFacesArrayRed(size,0);
    std::vector<int> nNonSharedArrayRed(size,0);
    std::vector<int> nNonSharedVertsArrayOff(size,0);
    for(int i=0;i<size;i++)
    {
        nNonSharedFacesArray[i] = 0;
        if(i==rank)
        {
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }

    MPI_Allreduce(nNonSharedFacesArray.data(),
                nNonSharedFacesArrayRed.data(),
                size,
                MPI_INT, MPI_SUM, comm);

    std::vector<int> nNonSharedFacesArrayOff(size,0);
    int nonFacesSharedOff     = 0;
    int nonSharedOff = 0;
    for(int i=0;i<size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }

    int lfid   = nNonSharedFacesArrayOff[rank];
    int gloVid = nNonSharedFacesArrayOff[rank];
    int gloFid = nNonSharedFacesArrayOff[rank];
    int locVid = 0;
    int locFid = 0;

    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        int gelid   = *its;     
        int gEl     = ElementDistr->getOffsets()[rank]+u+1;
        int Nf      = m_Elem2Face[gelid].size();
        int Nv      = m_Elem2Vert[gelid].size();
        
        for(int p=0;p<Nv;p++)
        {
            int TagVid = m_Elem2Vert[gelid][p];

            //std::cout << "TagVid " << TagVid << " ";

            if(m_sharedVertexMapUpdatedGlobalID.find(TagVid)!=m_sharedVertexMapUpdatedGlobalID.end())
            {
                int GlobVID             = m_sharedVertexMapUpdatedGlobalID[TagVid];

                if(m_tag2globV.find(TagVid)==m_tag2globV.end())
                {
                    m_tag2globV[TagVid]     = GlobVID;
                    m_glob2locV[GlobVID]    = locVid;
                    locVid = locVid + 1;
                }         
            }
            else
            {
                if(m_tag2globV.find(TagVid)==m_tag2globV.end())
                {
                    m_tag2globV[TagVid] = gloVid;
                    m_glob2locV[gloVid] = locVid;
                    gloVid = gloVid + 1;
                    locVid = locVid + 1;
                } 
            }         
        }
        
        //std::cout << std::endl;

        std::vector<int> faces(Nf,0);
        for(int q=0;q<Nf;q++)
        {
            int TagFid = m_Elem2Face[gelid][q];

            if(m_sharedFaceMapUpdatedGlobalID.find(TagFid)!=m_sharedFaceMapUpdatedGlobalID.end())
            {
                int GlobFID = m_sharedFaceMapUpdatedGlobalID[TagFid]; 
                faces[q]    = GlobFID; 
                m_Face2Elem_New[GlobFID].push_back(gEl);

                if(m_tag2globF.find(TagFid)==m_tag2globF.end())
                {
                    m_tag2globF[TagFid] = GlobFID;
                    m_glob2tagF[GlobFID] = TagFid;
                }
            }
            else
            {
                if(m_tag2globF.find(TagFid)==m_tag2globF.end())
                {
                    m_tag2globF[TagFid] = gloFid;
                    m_glob2tagF[gloFid] = TagFid;
                    faces[q]            = gloFid;
                    m_Face2Elem_New[gloFid].push_back(gEl);
                    gloFid              = gloFid + 1;
                }
                else
                {
                    int gloFid_tmp = m_tag2globF[TagFid];
                    faces[q]       = gloFid_tmp;
                    m_Face2Elem_New[gloFid_tmp].push_back(gEl);
                }

                lfid++;
            }
        }

        m_Elem2Face_New[gEl] = faces;

        u = u + 1;

    }

    // std::cout << "m_sharedVertexMapUpdatedGlobalID " << m_sharedVertexMapUpdatedGlobalID.size() << " " << rank << std::endl;
    
    for(itmiv=m_Elem2Face_New.begin();itmiv!=m_Elem2Face_New.end();itmiv++)
    {
        int eid = itmiv->first;
        int nf = itmiv->second.size();

        for(int i=0;i<nf;i++)
        {
            int fid    = m_Elem2Face_New[eid][i];
            int tagfid = m_glob2tagF[fid];

            if(m_Face2Vert.find(tagfid)!=m_Face2Vert.end())
            {
                int nv     = m_Face2Vert[tagfid].size();
                std::vector<int> verts(nv,0);
                for(int j=0;j<nv;j++)
                {
                    int tagvid  =  m_Face2Vert[tagfid][j];
                    verts[j]    =  m_tag2globV[tagvid];
                }
                if(m_Face2Vert_New.find(fid)==m_Face2Vert_New.end())
                {
                    m_Face2Vert_New[fid] = verts;
                }
            }
        }
    }

    TraceVertsOnRank_glob.clear();

    int element, fid;
    
    std::map<int,std::vector<int> >::iterator itf;
    std::vector<int> sharedFonRank;
    std::vector<int> interiorFonRank;

    for(itf=m_Face2Elem_New.begin();itf!=m_Face2Elem_New.end();itf++)
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
    
    int nSharedFonRank                           = sharedFonRank.size();
    DistributedParallelState* distSharedFacesNew = new DistributedParallelState(nSharedFonRank,comm);
    int Nt_shFaces_New                           = distSharedFacesNew->getNel();
    int* shFace_offsets_New                      = distSharedFacesNew->getOffsets();
    int* shFace_nlocs_New                        = distSharedFacesNew->getNlocs();

    std::vector<int> shFacesIDs(nSharedFonRank,0);
    std::vector<int> shFaces_RankIDs(nSharedFonRank,0);

    for(int i=0;i<nSharedFonRank;i++)
    {
        shFacesIDs[i]      = sharedFonRank[i];
        shFaces_RankIDs[i] = rank;
    }

    std::vector<int> TotalSharedFacesNew(Nt_shFaces_New,0);
    std::vector<int> TotalSharedFacesNew_RankID(Nt_shFaces_New,0);
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(&shFacesIDs.data()[0],
                   nSharedFonRank,
                   MPI_INT,
                   &TotalSharedFacesNew.data()[0],
                   shFace_nlocs_New,
                   shFace_offsets_New,
                   MPI_INT, comm);
    
    MPI_Allgatherv(&shFaces_RankIDs.data()[0],
                   nSharedFonRank,
                   MPI_INT,
                   &TotalSharedFacesNew_RankID.data()[0],
                   shFace_nlocs_New,
                   shFace_offsets_New,
                   MPI_INT, comm);
    
    std::map<int,std::vector<int> > face2rank;
    
    for(int i=0;i<Nt_shFaces_New;i++)
    {
        int key = TotalSharedFacesNew[i];
        int val = TotalSharedFacesNew_RankID[i];
        
        face2rank[key].push_back(val);
    }
    
    // std::cout << "face2rank " << face2rank.size() << " Nt_shFaces_New " << Nt_shFaces_New << " " << rank << " " << nSharedFonRank << " " << m_Face2Elem.size() << " " << m_tag2globF.size() << std::endl;
    
    //delete[] TotalSharedFaces;
    //delete[] TotalSharedFaces_RankID;
    delete[] shFace_offsets_New;
    delete[] shFace_nlocs_New;
    //delete[] shFacesIDs;
    //delete[] shFaces_RankIDs;
    //delete distSharedFaces;
    std::map<int,std::vector<int> >::iterator itff;
    shf = 0;
    int bf  = 0;
    std::map<int,std::vector<int> > Boundary_Ref2Face;
    std::set<int> uSharedVert;
    std::set<int> uBoundVert;
    int f = 0;

    // std::map<int,int> o_globShF2locShF;
    int once = 0;
    int twice = 0;
    for(itff=face2rank.begin();itff!=face2rank.end();itff++)
    {
        int faceID = itff->first;
        
        if(itff->second.size()==2) // Actual shared faces.
        {
            if(itff->second[0]==rank)
            {
                m_ColorsFaces[itff->second[1]].push_back(faceID);
            }
            else if(itff->second[1]==rank)
            {
                m_ColorsFaces[itff->second[0]].push_back(faceID);
            }

            if(m_Face2Vert_New.find(faceID)!=m_Face2Vert_New.end())
            {
                m_faces4parmmg.push_back(faceID);
                m_globShF2locShF[faceID] = f;
                m_locShF2globShF[f] = faceID;
                f++;
                
                shf++;
            }

            twice++;
            
            
        }
        
        if(itff->second.size()==1)// Boundary faces.
        {         
             if(m_Face2Vert_New.find(faceID)!=m_Face2Vert_New.end())
            {
                m_faces4parmmg.push_back(faceID);
                m_globShF2locShF[faceID] = f;
                m_locShF2globShF[f] = faceID;
                f++;
            } 
            once++;  
        }
    }
    
    m_ncomm             = m_ColorsFaces.size();
    m_color_face        = (int *) malloc(m_ncomm*sizeof(int));
    m_ntifc             = (int *) malloc(m_ncomm*sizeof(int));
    
    m_ifc_tria_loc      = (int **)malloc(m_ncomm*sizeof(int *));
    m_ifc_tria_glo      = (int **)malloc(m_ncomm*sizeof(int *));
    

    int icomm=0;
    std::map<int,std::vector<int> >::iterator itc;
    
    
    for(itc=m_ColorsFaces.begin();itc!=m_ColorsFaces.end();itc++)
    {
        m_color_face[icomm]     = itc->first;
        m_ntifc[icomm]          = itc->second.size();
        m_ifc_tria_loc[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));
        m_ifc_tria_glo[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));

        for(int q=0;q<itc->second.size();q++)
        {
            m_ifc_tria_glo[icomm][q] = itc->second[q]+1;
            m_ifc_tria_loc[icomm][q]  = m_globShF2locShF[itc->second[q]]+1;

        }

        icomm++;
    }

    /**/
}



void PrepareAdaption::buildInteriorSharedAndBoundaryFaceMaps(MPI_Comm               comm, 
                                                    std::set<int>                   TraceVertsOnRank,
                                                    std::set<int>                   TraceFacesOnRank,
                                                    std::map<int,std::vector<int> > ranges_id,
                                                    std::map<int,std::vector<int> > ranges_ref)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> >::iterator itmiv;
    //==============================
    // BEGIN building m_tagE2gE
    // BEGIN building m_gE2tagE
    // BEGIN building LeftHandRightHandElementVec lhp/rhp
    //==============================

    // buildInteriorSharedAndBoundaryFaceMaps()
    // std::map<int,std::map<int,int> > trace_elem     = trace->GetTrace();
    // //std::map<int,std::vector<int> > trace_verts     = trace->GetTraceVerts();
    // std::map<int,int> uniqure_trace_verts2ref       = trace->GetUniqueTraceVerts2RefMap();
    
    DistributedParallelState* ElementDistr = new DistributedParallelState(m_ElemSet.size(),comm);

    std::set<int> TraceVertsOnRank_glob = AllGatherSet(TraceVertsOnRank, comm);

    std::map<int,int> tag2ElementID;
    std::set<int> FaceOnRank;
    int nSharedFaces = m_sharedFace2Nodes.size();
    //std::cout << "nSharedFaces from long " << nSharedFaces << std::endl;
    int N_localFaces = m_Face2Elem.size();
    int nNonSharedFaces = N_localFaces-nSharedFaces;
    std::vector<int> nNonSharedFacesArray(size,0);
    std::vector<int> nNonSharedFacesArrayRed(size,0);
    std::vector<int> nNonSharedArrayRed(size,0);
    std::vector<int> nNonSharedVertsArrayOff(size,0);
    for(int i=0;i<size;i++)
    {
        nNonSharedFacesArray[i] = 0;

        if(i==rank)
        {
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }

    MPI_Allreduce(nNonSharedFacesArray.data(),
                nNonSharedFacesArrayRed.data(),
                size,
                MPI_INT, MPI_SUM, comm);

    std::vector<int> nNonSharedFacesArrayOff(size,0);
    int nonFacesSharedOff       = 0;
    int nonSharedOff            = 0;
    for(int i=0;i<size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff        = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff   = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }

    int lvid  = nNonSharedFacesArrayOff[rank];
    int lfid  = nNonSharedFacesArrayOff[rank];
    // lfid = 0;
    std::map<int,int> tagF2locFID;
    int fownedInmap = 0;
    std::set<int>::iterator its;
    int ref = -1;

    int u = 0;
    int gvidd = 0;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        int gelid           = *its;

        int lEl             = ElementDistr->getOffsets()[rank]+u+1;
        int Nf              = m_Elem2Face[gelid].size();
        int Nv              = m_Elem2Vert[gelid].size();

        m_tagE2gE[gelid]    = lEl;
        m_gE2tagE[lEl]      = gelid;

        for(int q=0;q<Nf;q++)
        {
            int gfid = m_Elem2Face[gelid][q];

            if(m_Face2Vert.find(gfid)!=m_Face2Vert.end())
            {
                int Nv   = m_Face2Vert[gfid].size();

                if(TraceFacesOnRank.find(gfid)!=TraceFacesOnRank.end())
                {
                    ref  = 13;
                }
                else
                {
                    ref  = ProvideBoundaryRef(gfid,ranges_ref);         
                }
                
                if(FaceOnRank.find(gfid)==FaceOnRank.end())
                {
                    FaceOnRank.insert(gfid);
                
                    if(m_sharedFace2RankMap.find(gfid)!=m_sharedFace2RankMap.end() && ref!=13)
                    {

                        if(m_sharedFace2RankMap[gfid] == rank)
                        {                                
                            if(m_sharedFace2Nodes.find(gfid)==m_sharedFace2Nodes.end())
                            {
                                int e0   = m_Face2Elem[gfid][0];
                                int e1   = m_Face2Elem[gfid][1];
                                
                                int r0 = rank; // = part_global[e0]; // rank of first adjacent element.
                                int r1 = rank; // = part_global[e1]; // rank of second adjacent element.
                                if(m_Elem2Rank.find(e0)!=m_Elem2Rank.end())
                                {
                                    r0 = m_Elem2Rank[e0];
                                }
                                if(m_Elem2Rank.find(e1)!=m_Elem2Rank.end())
                                {
                                    r1 = m_Elem2Rank[e1];
                                }

                                // int r0 = part_global[e0]; // rank of first adjacent element.
                                // int r1 = part_global[e1]; // rank of second adjacent element.
                                
                                if(r0==rank && r1!=rank)
                                {
                                    m_colorRh[r1].push_back(e1);
                                    m_colorFh[r1].push_back(gfid);
                                }
                                else if(r1==rank && r0!=rank)
                                {
                                    m_colorRh[r0].push_back(e0);
                                    m_colorFh[r0].push_back(gfid);
                                }
                                
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = m_Face2Vert[gfid][n];

                                    if(m_SharedVertsOwned.find(fn_tag[n])==m_SharedVertsOwned.end() 
                                            && TraceVertsOnRank_glob.find(fn_tag[n])==TraceVertsOnRank_glob.end())
                                    {
                                        if(m_SharedVertsNotOwned.find(fn_tag[n])==m_SharedVertsNotOwned.end())
                                        {
                                            m_SharedVertsNotOwned[fn_tag[n]] = m_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                        }
                                    }
                                }

                                m_sharedFace2Nodes[gfid] = fn_tag;
                                
                                m_lhp[gfid] = lEl;

                            }
                        }
                    }
                    else
                    {
                        if(ref == 2)
                        {
                            if(m_interiorFace2Nodes.find(gfid)==m_interiorFace2Nodes.end())
                            {
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = m_Face2Vert[gfid][n];
                                    
                                    if(m_sharedVertex2RankMap.find(fn_tag[n])!=m_sharedVertex2RankMap.end())
                                    {
                                        if(m_SharedVertsOwned.find(fn_tag[n])==m_SharedVertsOwned.end() 
                                                && TraceVertsOnRank_glob.find(fn_tag[n])==TraceVertsOnRank_glob.end())
                                        {
                                            if(m_SharedVertsNotOwned.find(fn_tag[n])==m_SharedVertsNotOwned.end())
                                            {
                                                m_SharedVertsNotOwned[fn_tag[n]] = m_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if(m_NonSharedVertsOwned.find(fn_tag[n])==m_NonSharedVertsOwned.end() 
                                        && TraceVertsOnRank_glob.find(fn_tag[n])==TraceVertsOnRank_glob.end())
                                        {
                                            m_NonSharedVertsOwned[fn_tag[n]]  = gvidd;
                                            gvidd++;
                                        }
                                    }
                                }

                                m_lhp[gfid] = lEl;
                                m_interiorFace2Nodes[gfid] = fn_tag;

                            }
                        }
                        if(ref != 2 && ref!=13)
                        {

                            int fzone = ProvideBoundaryID(gfid,ranges_id);
                            m_zone2bcface[fzone].push_back(gfid);

                            if(m_boundaryFace2Nodes.find(gfid)==m_boundaryFace2Nodes.end())
                            {
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = m_Face2Vert[gfid][n];

                                    if(m_sharedVertex2RankMap.find(fn_tag[n])!=m_sharedVertex2RankMap.end())
                                    {
                                        if(m_SharedVertsOwned.find(fn_tag[n])==m_SharedVertsOwned.end() 
                                                && TraceVertsOnRank_glob.find(fn_tag[n])==TraceVertsOnRank_glob.end())
                                        {
                                            if(m_SharedVertsNotOwned.find(fn_tag[n])==m_SharedVertsNotOwned.end())
                                            {
                                                m_SharedVertsNotOwned[fn_tag[n]] = m_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if(m_NonSharedVertsOwned.find(fn_tag[n])==m_NonSharedVertsOwned.end() 
                                        && TraceVertsOnRank_glob.find(fn_tag[n])==TraceVertsOnRank_glob.end())
                                        {
                                            m_NonSharedVertsOwned[fn_tag[n]]  = gvidd;
                                            gvidd++;
                                        }
                                    }
                                    
                                }

                                m_lhp[gfid] = lEl;
                                m_boundaryFace2Nodes[gfid] = fn_tag;
                                
                            }
                        }
                    }
                }
                else
                {
                    tag2ElementID[gelid] = lEl;
                    m_rhp[gfid]          = lEl;
                }
            }
        }
        u++;    
    }


    std::map<int,std::vector<int> >::iterator itv;
    ScheduleObj* rh_schedule = DoScheduling(m_colorRh,comm);
    std::map<int,std::vector<int> > recv_rhElIDs;
    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            
            for (itv = m_colorRh.begin(); itv != m_colorRh.end(); itv++)
            {
                int n_req           = itv->second.size();
                int dest            = itv->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
                MPI_Send(&itv->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
                i++;
            }
        }
        else if (rh_schedule->SendFromRank2Rank[q].find( rank ) != rh_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 
                    14876+rank, comm, MPI_STATUS_IGNORE);

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
            int gEl = m_tagE2gE[rcvit->second[j]];
            sendEl[frank].push_back(gEl);
        }
    }
    
    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);

    std::map<int,std::vector<int> > adj_ids;
    for(int q=0;q<size;q++)
    {
        if(rank==q)
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
        else if (ishBack_schedule->SendFromRank2Rank[q].find( rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            
            MPI_Recv(&n_reqstd_ids,
                    1, MPI_INT, q,
                    6798+78000*rank,
                    comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0],
                    n_reqstd_ids,
                    MPI_INT, q,
                    14876000+rank,
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
            int fid = m_colorFh[rrank][q];
            
            if(shFid2el_rh.find(fid)==shFid2el_rh.end())
            {
                shFid2el_rh[fid] = rcvit->second[q];
                adde2++;
                if(m_rhp.find(fid)==m_rhp.end())
                {
                    adde++;
                    m_rhp[fid] = rcvit->second[q];
                }
            }
        }
    }
    TraceVertsOnRank_glob.clear();
}


std::vector<int> PrepareAdaption::getFace4ParMMG()
{
    return m_faces4parmmg;
}
int** PrepareAdaption::getParMMGCommFace2GlobalVertMap()
{
    return m_ifc_tria_glo;
}

int** PrepareAdaption::getParMMGCommFace2LocalVertMap()
{
    return m_ifc_tria_loc;
}

int* PrepareAdaption::getParMMGCommColorFace()
{
    return m_color_face;
}

int* PrepareAdaption::getParMMGCommNFacesPerColor()
{
    return m_ntifc;
}
int PrepareAdaption::getParMMGNComm()
{
    return m_ncomm;
}

std::map<int,int> PrepareAdaption::getGlobal2TagFMap()
{
    return m_glob2tagF;
}

std::map<int,int> PrepareAdaption::getNewGlobalVert2LocalVertMap()
{
    return m_glob2locV;
}

std::map<int,int> PrepareAdaption::getTagVert2GlobalVertMap()
{
    return m_tag2globV;
}

std::set<int> PrepareAdaption::getElemSet()
{
    return m_ElemSet;
}

std::map<int,std::vector<double> > PrepareAdaption::getLocalVertsMap()
{
    return m_LocalVertsMap;
}

std::map<int,std::vector<int> > PrepareAdaption::getFace2VertMap()
{
    return m_Face2Vert;
}

std::map<int,std::vector<int> > PrepareAdaption::getElem2VertMap()
{
    return m_Elem2Vert;
}

std::map<int,int> PrepareAdaption::getLocalVert2GlobalVert()
{
    return m_locV2gloV;
}

std::map<int,std::vector<int> > PrepareAdaption::getFace2VertNewMap()
{
    return m_Face2Vert_New;
}
std::map<int,int> PrepareAdaption::getNonSharedVertsOwned()
{
    return m_NonSharedVertsOwned;
}
std::map<int,int> PrepareAdaption::getSharedVertsOwned()
{
    return m_SharedVertsOwned;
}
std::map<int,std::vector<int> > PrepareAdaption::getInteriorFaceMap()
{
    return m_interiorFace2Nodes;
}
std::map<int,std::vector<int> > PrepareAdaption::getBoundaryFaceMap()
{
    return m_boundaryFace2Nodes;
}
std::map<int,std::vector<int> > PrepareAdaption::getSharedFaceMap()
{
    return m_sharedFace2Nodes;
}
std::map<int,std::vector<int> > PrepareAdaption::getZone2boundaryFaceID()
{
    return m_zone2bcface;
}
std::map<int,int> PrepareAdaption::GetLeftHandFaceElementMap()
{
    return m_lhp;
}
std::map<int,int> PrepareAdaption::GetRightHandFaceElementMap()
{
    return m_rhp;
}
std::map<int,int> PrepareAdaption::getLocalSharedFace2GlobalSharedFace()
{
    return m_locShF2globShF;
}
std::map<int,int> PrepareAdaption::getSharedVertsNotOwned()
{
    return m_SharedVertsNotOwned;
}
std::map<int,int> PrepareAdaption::getGlobalElement2ElementTag()
{
    return m_gE2tagE;
}
std::map<int,int> PrepareAdaption::getElementTag2GlobalElement()
{
    return m_tagE2gE;
}
