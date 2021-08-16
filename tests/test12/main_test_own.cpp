#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted



struct PartitionInfo
{
    Array<int>* part;
    Array<int>* part_global;
};

struct Hyb2TetraMeshTransfer
{
    std::map<int,std::vector<int> > m_TetraToSendToRanks;
    std::map<int,std::vector<int> > m_HybridToSendToRanks;
    std::map<int,std::vector<int> > m_TetraVertIDToSendToRanks;
    std::map<int,std::vector<int> > m_HybridVertIDToSendToRanks;
    std::map<int,std::vector<int> > m_TetraFaceIDToSendToRanks;
    std::map<int,std::vector<int> > m_TetraFaceRefToSendToRanks;
    std::map<int,std::vector<int> > m_HybridFaceIDToSendToRanks;
    std::map<int,std::vector<int> > m_TetraRank2ReqVerts;
    std::map<int,std::vector<int> > m_HybridRank2ReqVerts;
//    std::map<int,std::vector<int> > m_TetraRank2ReqFaces;
//    std::map<int,std::vector<int> > m_HybridRank2ReqFaces;
    std::map<int,std::vector<double> > metricsToSend;
    
    std::vector<int> m_TetraVertsOnRank;
    std::vector<int> m_HybridVertsOnRank;

    std::map<int,int> m_TetEl2HybEl;
};

struct TetrahedraMesh
{
    Array<int>* ElGids;
    Array<int>* ien_part_tetra;
    Array<int>* ien_part_hybrid;
    std::map<int,std::vector<int> > m_TetraToSendToRanks;
    Array<int>* ief_part_tetra;
    Array<int>* ief_part_hybrid;
    
    Array<int>* iefref_part_tetra;
    
    std::vector<Vert*> LocalVerts;
    std::map<int,int> locV2globV;
    std::map<int,int> globV2locV;
    
    std::vector<int*> LocalFaces;
    std::map<int,int> locF2globF;
    std::map<int,int> globF2locF;
    
    std::map<int,int> hybF2tetF;
    std::map<int,int> tetF2hybF;
    
    std::map<int,int> hybV2tetV;
    std::map<int,int> tetV2hybV;
    
    std::map<int,int> hybE2tetE;
    std::map<int,int> tetE2hybE;
    
    std::map<int,int> face2ref;
    std::map<int,std::vector<int> > ref2face;
    std::map<int,Array<double>* > M_vmap;
};


struct PartitionBoundary
{
    int nPartBoundFaces;
    std::vector<int> faces4parmmg;
    std::map<int,int> globShF2locShF;
    std::map<int,int> locShF2globShF;
    int* ifc_tria_glob;
    int** ifc_tria_loc;
    int* color_face;
    int* ntifc;
    int ncomm;
    std::map<int,std::vector<int> > ColorsFaces;
    std::map<int,int> face2ref;
};





PartitionBoundary* ExtractPartitionBoundary(TetrahedraMesh* tmesh, std::map<int,std::vector<int> > face2rank, std::map<int,int*> face2node, MPI_Comm comm)
{
    PartitionBoundary* pb = new PartitionBoundary;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    std::vector<Vert*> locVs     = tmesh->LocalVerts;

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
    std::map<int,std::vector<int> > ColorsFaces;
    std::map<int,int> globShF2locShF;
  
    for(itff=face2rank.begin();itff!=face2rank.end();itff++)
    {
        int faceID = itff->first;
        
        if(itff->second.size()==2)
        {
            if(itff->second[0]==world_rank)
            {
                ColorsFaces[itff->second[1]].push_back(faceID);
            }
            else if(itff->second[1]==world_rank)
            {
                ColorsFaces[itff->second[0]].push_back(faceID);
            }
            
            std::vector<int> sface(3);
            if(face2node.find(faceID)!=face2node.end())
            {
                pb->faces4parmmg.push_back(faceID);

                for(int u=0;u<3;u++)
                {
                    int vid = face2node[faceID][u];

                    if(uSharedVert.find(vid)==uSharedVert.end())
                    {
                        uSharedVert.insert(vid);
                        g2l_s[vid] = Svid;

                        lSvid = tmesh->globV2locV[vid];

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
                pb->globShF2locShF[faceID] = f;
                pb->locShF2globShF[f] = faceID;
                f++;
                
                shf++;
            }
            
        }
        
        
        if(itff->second.size()==1)
        {
            std::vector<int> bface(3);
            
            if(face2node.find(faceID)!=face2node.end())
            {
                pb->faces4parmmg.push_back(faceID);

                if(tmesh->face2ref.find(faceID)!=tmesh->face2ref.end())
                {
                    ref = tmesh->face2ref[faceID];
                    Boundary_Ref2Face[ref].push_back(faceID);
                }
                
                for(int u=0;u<3;u++)
                {
                    int vid = face2node[faceID][u];
                    
                    if(uBoundVert.find(vid)==uBoundVert.end())
                    {
                        v2ref[vid]   = ref;
                        int hybv     = tmesh->tetV2hybV[vid];
                        uBoundVert.insert(vid);
                        g2l_b[vid]   = Bvid;
                        b2l_g[Bvid]  = vid;
                        lBvid        = tmesh->globV2locV[vid];
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
                pb->globShF2locShF[faceID] = f;
                pb->locShF2globShF[f] = faceID;
                bf++;
                f++;
                
            }
            
        }
    }
    
    pb->nPartBoundFaces = pb->faces4parmmg.size();
    pb->ncomm = ColorsFaces.size();
    
    
    //================================================================
    std::string filename = "SharedPartFaces_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << svs.size() << ", E = " << pltSfaces.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;

    for(int i=0;i<svs.size();i++)
    {
        myfile << svs[i]->x << " " << svs[i]->y << " " << svs[i]->z << std::endl;
    }
    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
    for(int i=0;i<pltSfaces.size();i++)
    {
        myfile <<   pltSfaces[i][0]+1 << "  " <<
                    pltSfaces[i][1]+1 << "  " <<
                    pltSfaces[i][2]+1 << std::endl;
    }


    myfile.close();
    
    
    //================================================================
    
    int nBfaces = pltBfaces.size();
    
    
    //================================================================
    std::string filename2 = "BoundaryPartFaces_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile2;
    myfile2.open(filename2);
    myfile2 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile2 <<"ZONE N = " << bvs.size() << ", E = " << nBfaces << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
    int evv = -1;
    for(int i=0;i<bvs.size();i++)
    {
        myfile2 << bvs[i]->x << " " << bvs[i]->y << " " << bvs[i]->z << std::endl;
    }

    
    for(int i=0;i<nBfaces;i++)
    {
        myfile2 <<   pltBfaces[i][0]+1 << "  " <<
                     pltBfaces[i][1]+1 << "  " <<
                     pltBfaces[i][2]+1 << std::endl;
   
    }

    myfile2.close();
    
    
//    pb->faces4parmgg;
//    pb->globShF2locShF;
    pb->ColorsFaces=ColorsFaces;
//    pb->ifc_tria_glob = ifc_tria_glob;
//    pb->ifc_tria_loc = ifc_tria_loc;
//    pb->color_face;
//    pb->ntifc;
    
    
    
    return pb;
    //================================================================
}

std::map<int,std::vector<int> > GetFace2RankTetrahedraMesh(TetrahedraMesh* tmesh, ParallelState* ien_pstate_tet, MPI_Comm comm)
{
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    std::map<int,int*>::iterator itst;
    int element, fid;
    std::map<int,std::vector<int> > face2element;
    for(int i=0;i<tmesh->ief_part_tetra->getNrow();i++)
    {
        element = i + ien_pstate_tet->getOffsets()[world_rank];
        
        for(int j=0;j<4;j++)
        {
            fid = tmesh->ief_part_tetra->getVal(i,j);
            face2element[fid].push_back(element);
        }
    }
    
    std::map<int,std::vector<int> >::iterator itf;
    std::vector<int> sharedFonRank;
    std::vector<int> interiorFonRank;
    std::vector<Vert*> locVs = tmesh->LocalVerts;

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
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFonRank,comm);

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
                   MPI_INT, comm);
    
    MPI_Allgatherv(shFaces_RankIDs,
                   nSharedFonRank,
                   MPI_INT,
                   TotalSharedFaces_RankID,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
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
    delete distSharedFaces;
    
    return face2rank;
}

std::map<int,int*> GetFace2EntityTetrahedraMesh(TetrahedraMesh* tmesh, ParArray<int>* ife, int ncol, ParallelState* ife_pstate, int nGlob, MPI_Comm comm)
{
    std::map<int,int*> face2node;
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
    
    tmesh->tetF2hybF.clear();
    tmesh->hybF2tetF.clear();
    std::set<int> ufacesHyb;
    std::set<int> ufacesTet;

    std::vector<int> ee;
    int el_id;
    std::map<int,int> face2ref;
    std::map<int,std::vector<int> >::iterator itefmap;
    std::map<int,std::vector<int> >ife_part_hyb_map;
    for(int i=0;i<tmesh->ief_part_hybrid->getNrow();i++)
    {
        el_id   = ien_pstate->getOffsets()[rank]+i;
        
        for(int q=0;q<4;q++)
        {
            int face_hyb = tmesh->ief_part_hybrid->getVal(i,q);
            int face_tet = tmesh->ief_part_tetra->getVal(i,q);
            int face_ref = tmesh->iefref_part_tetra->getVal(i,q);

            //face2element[face_tet].push_back(el_id);
            
            if(ufacesHyb.find(face_hyb)==ufacesHyb.end())
            {
                ufacesHyb.insert(face_hyb);
                tmesh->tetF2hybF[face_tet] = face_hyb;
                tmesh->hybF2tetF[face_hyb] = face_tet;
                face2ref[face_tet] = face_ref;
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
            
            if(tmesh->hybF2tetF.find(face_id)!=tmesh->hybF2tetF.end())
            {
                ftet = tmesh->hybF2tetF[iter->first];
                face2ref[ftet] = face_ref;
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
    for(itf=ife_part_hyb_map.begin();itf!=ife_part_hyb_map.end();itf++)
    {
        fhyb = itf->first;
        
        if(tmesh->hybF2tetF.find(itf->first)!=tmesh->hybF2tetF.end())
        {
            ftet = tmesh->hybF2tetF[itf->first];
            face_ref =
            face2ref[ftet] = face_ref;

            if(face2ref.find(ftet)!=face2ref.end())
            {
                
            }
            else
            {
                std::cout << "Not here " << std::endl;
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

            if(tmesh->hybV2tetV.find(hvid)!=tmesh->hybV2tetV.end())
            {
                tvid = tmesh->hybV2tetV[hvid];
            }
            else
            {
                std::cout << "Errror hvid doesnt exist " << hvid << " " << " This error occurs for face " << ftet << " " << fhyb << " on rank " << rank << std::endl;
            }
            
            nodes[j]     = tvid;
        }
        
        face2node[ftet] = nodes;
    }
    
    //std::cout << rank << " ==> " << face2ref.size() << " -- " << tmesh->face2ref.size() << std::endl;
    /**/
    return face2node;
    
    
    
}


void OutputTetrahedralMeshOnPartition(TetrahedraMesh* tmesh, MPI_Comm comm)
{
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::vector<int> lverts;
    std::map<int,int> lpartv2gv_v2;
    std::map<int,int> gv2lpv2;

    std::set<int> gv_set;
    int lcv2 = 0;
    Array<int>* ien_part_tetra     = tmesh->ien_part_tetra;
    Array<int>* ien_part_hybrid    = tmesh->ien_part_hybrid;
    std::vector<Vert*> locVs       = tmesh->LocalVerts;
    int nElonRank = ien_part_tetra->getNrow();

    Array<int>* locelem2locnode= new Array<int>(nElonRank,4);
    
    std::vector<Vert*> printVs;
    
    for(int i=0;i<ien_part_tetra->getNrow();i++)
    {
        for(int q=0;q<ien_part_tetra->getNcol();q++)
        {
            int gv = ien_part_tetra->getVal(i,q);
            int lvv = tmesh->globV2locV[gv];
            
            if(gv_set.find(gv)==gv_set.end())
            {
                gv_set.insert(gv);
                lverts.push_back(lvv);
                lpartv2gv_v2[lvv]=gv;
                gv2lpv2[gv]=lcv2;
                locelem2locnode->setVal(i,q,lcv2);
                
                printVs.push_back(locVs[lvv]);
                
                lcv2=lcv2+1;
            }
            else
            {
                int lcv_u = gv2lpv2[gv];
                locelem2locnode->setVal(i,q,lcv_u);
            }
        }
    }
        
    std::vector<Vert*> lv = tmesh->LocalVerts;
    std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << printVs.size() << ", E = " << nElonRank << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<printVs.size();i++)
    {
        myfile << printVs[i]->x << " " << printVs[i]->y << " " << printVs[i]->z << std::endl;
    }
    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
    for(int i=0;i<ien_part_hybrid->getNrow();i++)
    {
        myfile <<   locelem2locnode->getVal(i,0)+1 << "  " <<
        locelem2locnode->getVal(i,1)+1 << "  " <<
        locelem2locnode->getVal(i,2)+1 << "  " <<
        locelem2locnode->getVal(i,3)+1 << "  " << std::endl;
    }


    myfile.close();
}

int* CommunicatePartitionLayout(Array<int>* part, MPI_Comm comm)
{

    int pid;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int* newSizesOnRanks_red = new int[world_size];
    int* newSizesOnRanks     = new int[world_size];

    for(int i=0;i<world_size;i++)
    {
        newSizesOnRanks[i] = 0;
        newSizesOnRanks_red[i] = 0;
    }
    
    for(int i=0;i<part->getNrow();i++)
    {
        pid = part->getVal(i,0);
        newSizesOnRanks[pid] = newSizesOnRanks[pid]+1;
    }
    
    MPI_Allreduce(newSizesOnRanks, newSizesOnRanks_red, world_size, MPI_INT, MPI_SUM, comm);

    return newSizesOnRanks_red;
}






void UpdateTetrahedraOnPartition(int nglob, int nElexpctd, Array<int>* part,
                                 TetrahedraMesh* &tmesh,
                                 ParArray<double>* xcn,
                                 ParallelState* xcn_pstate,
                                 ParallelState* ife_pstate,
                                 MPI_Comm comm)
{
//    tmesh->locV2globV.clear();
//    tmesh->globV2locV.clear();

    Array<int>* ElGids             = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_tetra     = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_hybrid    = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_tetra     = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_hybrid    = new Array<int>(nElexpctd,4);
    Array<int>* iefref_part_tetra  = new Array<int>(nElexpctd,4);
    Array<int>* loc2globEL         = new Array<int>(nElexpctd,1);
    ParallelState* ien_pstate      = new ParallelState(nglob,comm);

    std::map<int,Array<double>* > M_vmap_copy;
    int floc_tmp=0;
    int vloc_tmp=0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

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
    int* new_F_offsets = new int[size];
    
    for(i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
        new_F_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    int nvPerEl;
    int nfPerEl;
    
    int fref;
    int tellert = 0;
    int sum = 0;
    int sum3 = 0;
    int newElID;
    //int gEL;
    Hyb2TetraMeshTransfer* tmesh_trans = new Hyb2TetraMeshTransfer;

    tmesh->hybV2tetV.clear();
    tmesh->tetV2hybV.clear();
    tmesh->face2ref.clear();

    
    for(i=0;i<part->getNrow();i++)
    {
        p_id        = part->getVal(i,0);
        el_id       = ien_pstate->getOffsets()[rank]+i;
        int gEL     = tmesh->ElGids->getVal(i,0);
	   //loc2globEL->setVal(i,0,el_id);
////        tmesh_trans->m_TetEl2HybEl[el_id] = gEL;
//        nvPerEl = 4;
//        nfPerEl = 4;
//        sum  = 0;
//        sum3 = 0;
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            tmesh_trans->m_TetraToSendToRanks[p_id].push_back(el_id); // rank to element map.
            tmesh_trans->m_HybridToSendToRanks[p_id].push_back(gEL);
            
            //====================Hybrid=======================
            for(int k=0;k<4;k++)
            {
                v_id   = tmesh->ien_part_tetra->getVal(i,k);
                v_id_o = tmesh->ien_part_hybrid->getVal(i,k);
                
                tmesh_trans->m_TetraVertIDToSendToRanks[p_id].push_back(v_id);
                tmesh_trans->m_HybridVertIDToSendToRanks[p_id].push_back(v_id_o);
                sum3=sum3+v_id;
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            
            for(int k=0;k<4;k++)
            {
                f_id   = tmesh->ief_part_tetra->getVal(i,k);
                f_id_o = tmesh->ief_part_hybrid->getVal(i,k);
                fref   = tmesh->iefref_part_tetra->getVal(i,k);
                
                
                tmesh_trans->m_TetraFaceIDToSendToRanks[p_id].push_back(f_id);
                tmesh_trans->m_TetraFaceRefToSendToRanks[p_id].push_back(fref);
                tmesh_trans->m_HybridFaceIDToSendToRanks[p_id].push_back(f_id_o);
            }
            //====================Hybrid=======================
            not_on_rank++;/**/
        }
        
        
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            
            std::vector<int> elem;

            ElGids->setVal(on_rank,0,gEL);

            tmesh->hybE2tetE[gEL]=el_id;
            tmesh->tetE2hybE[el_id]=gEL;
            
            //tmesh->hybE2tetE[] =
            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                v_id   = tmesh->ien_part_tetra->getVal(i,k);
                v_id_o = tmesh->ien_part_hybrid->getVal(i,k);
                
                ien_part_tetra->setVal(on_rank,k,v_id);
                ien_part_hybrid->setVal(on_rank,k,v_id_o);
                
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    
                    tmesh->hybV2tetV[v_id_o] = v_id;
                    tmesh->tetV2hybV[v_id]   = v_id_o;
                    M_vmap_copy[v_id]        = tmesh->M_vmap[v_id];

                    r = FindRank(new_V_offsets,size,v_id_o);
                    
                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        tmesh_trans->m_TetraRank2ReqVerts[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                        tmesh_trans->m_HybridRank2ReqVerts[r].push_back(v_id_o);
                    }
                    else
                    {
                        tmesh_trans->m_TetraVertsOnRank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        tmesh_trans->m_HybridVertsOnRank.push_back(v_id_o);

                        vloc_tmp++;
                    }
                    
                    lv_id++;
                    
                }
            }
            
    
            
            for(int k=0;k<4;k++)// looping over the vertices for element "i".
            {
                f_id    = tmesh->ief_part_tetra->getVal(i,k);
                f_id_o  = tmesh->ief_part_hybrid->getVal(i,k);
                fref    = tmesh->iefref_part_tetra->getVal(i,k);
                
                ief_part_tetra->setVal(on_rank,k,f_id);
                ief_part_hybrid->setVal(on_rank,k,f_id_o);
                iefref_part_tetra->setVal(on_rank,k,fref);
                
                if(tmesh->face2ref.find(f_id)==tmesh->face2ref.end())
                {
                    tmesh->face2ref[f_id] = fref;
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
////                        tmesh_trans->m_TetraRank2ReqFaces[r].push_back(f_id); // add the vertex id that needs to be requested from rank r.
////                        tmesh_trans->m_HybridRank2ReqFaces[r].push_back(f_id_o);
//
//                    }
//                    else
//                    {
//                        //faceIDs_on_rank.push_back(f_id);  // add the vertex to list that is already available on rank.
////                      tmesh->hybF2tetF[f_id_o]=f_id;
////                      tmesh->tetF2hybF[f_id]=f_id_o;
////                        tmesh->ref2face[fref].push_back(f_id);
////                        tmesh->face2ref[f_id]=fref;
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
    for(itt=tmesh_trans->m_TetraVertIDToSendToRanks.begin();itt!=tmesh_trans->m_TetraVertIDToSendToRanks.end();itt++)
    {
        std::vector<double> metrics(6*itt->second.size());

        for(int q=0;q<itt->second.size();q++)
        {
            globvid = itt->second[q];

            metrics[q*6+0] = tmesh->M_vmap[globvid]->getVal(0,0);
            metrics[q*6+1] = tmesh->M_vmap[globvid]->getVal(1,0);
            metrics[q*6+2] = tmesh->M_vmap[globvid]->getVal(2,0);
            metrics[q*6+3] = tmesh->M_vmap[globvid]->getVal(3,0);
            metrics[q*6+4] = tmesh->M_vmap[globvid]->getVal(4,0);
            metrics[q*6+5] = tmesh->M_vmap[globvid]->getVal(5,0);
            
        }

        tmesh_trans->metricsToSend[itt->first] = metrics;
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
            int hybF = ien_part_hybrid->getVal(i,j);
            int tetF = ien_part_tetra->getVal(i,j);
            
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
    
    
    ScheduleObj* part_schedule_elem = DoScheduling(tmesh_trans->m_TetraToSendToRanks,comm);
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
            for (it = tmesh_trans->m_TetraToSendToRanks.begin(); it != tmesh_trans->m_TetraToSendToRanks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = tmesh_trans->m_TetraVertIDToSendToRanks[it->first].size();
                int n_req_f         = tmesh_trans->m_TetraFaceIDToSendToRanks[it->first].size();
                int n_req_M         = n_req_v*6;
                
                
                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, comm);
                MPI_Send(&n_req_f, 1, MPI_INT, dest, dest*222, comm);
                MPI_Send(&n_req_M, 1, MPI_INT, dest, dest*33300, comm);
                
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&tmesh_trans->m_TetraVertIDToSendToRanks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&tmesh_trans->m_HybridVertIDToSendToRanks[it->first][0], n_req_v, MPI_INT, dest, 339000+100+dest*2, comm);
                MPI_Send(&tmesh_trans->m_TetraFaceIDToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                MPI_Send(&tmesh_trans->m_TetraFaceRefToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 449000+100+dest*2, comm);
                MPI_Send(&tmesh_trans->m_HybridFaceIDToSendToRanks[it->first][0], n_req_f, MPI_INT, dest, 889000+100+dest*2, comm);
                MPI_Send(&tmesh_trans->m_HybridToSendToRanks[it->first][0], n_req, MPI_INT, dest, dest*7000000, comm);
                MPI_Send(&tmesh_trans->metricsToSend[it->first][0], n_req_M, MPI_DOUBLE, dest, dest*8000000, comm);

                

                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank,    comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111,    comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222,    comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_M, 1, MPI_INT, q, rank*33300,  comm, MPI_STATUS_IGNORE);

            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_vrt_id(n_req_recv_v);
            std::vector<int>    part_recv_orivrt_id(n_req_recv_v);
            std::vector<int>    part_recv_face_id(n_req_recv_f);
            std::vector<int>    part_recv_face_ref(n_req_recv_f);
            std::vector<int>    part_recv_orifaces_id(n_req_recv_v);
            std::vector<int>    part_recv_Gel_id(n_req_recv);
            std::vector<double>    part_recv_Mtrcs(n_req_recv_M);


            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0], n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orivrt_id[0], n_req_recv_v, MPI_INT, q, 339000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_ref[0], n_req_recv_f, MPI_INT, q, 449000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orifaces_id[0], n_req_recv_f, MPI_INT, q, 889000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_Gel_id[0], n_req_recv, MPI_INT, q, rank*7000000, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_Mtrcs[0], n_req_recv_M, MPI_DOUBLE, q, rank*8000000, comm, MPI_STATUS_IGNORE);

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
        
                
                if(M_vmap_copy.find(TotRecvElement_IDs_v_map[totrecv->first][l*4+h])==M_vmap_copy.end())
                {
                    Array<double>* mtrcs = new Array<double>(6,1);
                    
                    mtrcs->setVal(0,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+0]);
                    mtrcs->setVal(1,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+1]);
                    mtrcs->setVal(2,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+2]);
                    mtrcs->setVal(3,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+3]);
                    mtrcs->setVal(4,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+4]);
                    mtrcs->setVal(5,0,TotRecvVertMtrcs_map[totrecv->first][6*4*l+6*h+5]);
                    
                    M_vmap_copy[TotRecvElement_IDs_v_map[totrecv->first][l*4+h]] = mtrcs;
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
        
        ElGids->setVal(on_rank,0,GelID);
        tmesh->hybE2tetE[GelID]=elID;
        tmesh->tetE2hybE[elID]=GelID;
        
        for(int k=0;k<4;k++)
        {
            int v_id_n   = TotRecvVerts_IDs[cnt_v+k];
            int v_id_o_n = TotRecvOriVerts_IDs[cnt_v+k];
            
            ien_part_tetra->setVal(on_rank,k,v_id_n);
            ien_part_hybrid->setVal(on_rank,k,v_id_o_n);
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                tmesh->hybV2tetV[v_id_o_n]=v_id_n;
                tmesh->tetV2hybV[v_id_n]=v_id_o_n;
                
                r = FindRank(new_V_offsets, size, v_id_o_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    tmesh_trans->m_TetraRank2ReqVerts[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                    tmesh_trans->m_HybridRank2ReqVerts[r].push_back(v_id_o_n);
                }
                else
                {
                    tmesh_trans->m_TetraVertsOnRank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    tmesh_trans->m_HybridVertsOnRank.push_back(v_id_o_n);
                
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
            
            ief_part_tetra->setVal(on_rank,k,f_id_n);
            ief_part_hybrid->setVal(on_rank,k,f_id_o_n);
            iefref_part_tetra->setVal(on_rank,k,fref_n);
     	        
            if(tmesh->face2ref.find(f_id_n)==tmesh->face2ref.end())
            {
                tmesh->face2ref[f_id_n] = fref_n;
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
////                    tmesh_trans->m_TetraRank2ReqFaces[r].push_back(f_id_n); // add the vertex id that needs to be requested from rank r.
////                    tmesh_trans->m_HybridRank2ReqFaces[r].push_back(f_id_o_n);
//                }
//                else
//                {
////                    tmesh->face2ref[f_id_n] = fref_n;
////                    tmesh->ref2face[fref_n].push_back(f_id_n);
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
    
    ScheduleObj* part_schedule = DoScheduling(tmesh_trans->m_TetraRank2ReqVerts,comm);
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Ori_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Metrics_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = tmesh_trans->m_TetraRank2ReqVerts.begin(); it != tmesh_trans->m_TetraRank2ReqVerts.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                
                //	MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //	MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                MPI_Send(&tmesh_trans->m_HybridRank2ReqVerts[it->first][0], n_req, MPI_INT, dest, 2229876*2+dest*2, comm);

                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_Ori_ids(n_reqstd_ids);
            std::vector<int> recv_metrics(n_reqstd_ids*6);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_reqstd_Ori_ids[0], n_reqstd_ids, MPI_INT, q, 2229876*2+rank*2, comm, MPI_STATUS_IGNORE);
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
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&vert_send[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, comm);
                MPI_Send(&reqstd_ids_per_rank[it->first][0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
                delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, comm, MPI_STATUS_IGNORE);
            
            double* recv_back_arr = new double[n_recv_back*3];
            int* recv_back_arr_ids = new int[n_recv_back];
            //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr_ids[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, comm, MPI_STATUS_IGNORE);

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
    
    for(m=0;m<tmesh_trans->m_TetraVertsOnRank.size();m++)
    {
        gvid    = tmesh_trans->m_TetraVertsOnRank[m];
        gvid_gl = tmesh_trans->m_HybridVertsOnRank[m];
        
        Vert* V = new Vert;

        V->x = xcn->getVal(gvid_gl-xcn_o,0);
        V->y = xcn->getVal(gvid_gl-xcn_o,1);
        V->z = xcn->getVal(gvid_gl-xcn_o,2);

        tmesh->LocalVerts.push_back(V);
        
        if(tmesh->locV2globV.find(lvid)==tmesh->locV2globV.end())
        {
            tmesh->locV2globV[lvid] = gvid;
            tmesh->globV2locV[gvid] = lvid;
        }
        
        lvid++;
    }
   
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = tmesh_trans->m_TetraRank2ReqVerts[it_f->first][u];
            
            Vert* V = new Vert;
            
            V->x = it_f->second[u*3+0];
            V->y = it_f->second[u*3+1];
            V->z = it_f->second[u*3+2];
            
            tmesh->LocalVerts.push_back(V);
            tmesh->locV2globV[lvid]=gvid;
            tmesh->globV2locV[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }
    

    int nLoc_Verts = tmesh->LocalVerts.size();
    
    
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
    
    tmesh_trans->m_TetraToSendToRanks.clear();
    tmesh_trans->m_HybridToSendToRanks.clear();
    tmesh_trans->m_TetraVertIDToSendToRanks.clear();
    tmesh_trans->m_HybridVertIDToSendToRanks.clear();
    tmesh_trans->m_TetraFaceIDToSendToRanks.clear();
    tmesh_trans->m_TetraFaceRefToSendToRanks.clear();
    tmesh_trans->m_HybridFaceIDToSendToRanks.clear();
    tmesh_trans->m_TetraRank2ReqVerts.clear();
    tmesh_trans->m_HybridRank2ReqVerts.clear();
    //tmesh_trans->m_TetraRank2ReqFaces.clear();
    //tmesh_trans->m_HybridRank2ReqFaces.clear();
    tmesh_trans->metricsToSend.clear();
    
    std::vector<int> m_TetraVertsOnRank;
    std::vector<int> m_HybridVertsOnRank;

    std::map<int,int> m_TetEl2HybEl;
    
    tmesh->ElGids               = ElGids;
    tmesh->ien_part_tetra       = ien_part_tetra;
    tmesh->ien_part_hybrid      = ien_part_hybrid;
    tmesh->ief_part_tetra       = ief_part_tetra;
    tmesh->ief_part_hybrid      = ief_part_hybrid;
    tmesh->iefref_part_tetra    = iefref_part_tetra;
    tmesh->M_vmap               = M_vmap_copy;
    
    
        //std::cout << "Testje =  "<< ien_part_tetra->getVal(5125,0) << " " << ien_part_tetra->getVal(5125,1) << " " << ien_part_tetra->getVal(5125,2) << " " << ien_part_tetra->getVal(5125,3) << " on_rank " << on_rank <<" "<<ien_part_tetra->getNrow()<< std::endl;
    
    
    /**/
    
    
    //return tmesh_ret;
}


//void ExtractSharedFaces(Array<int>* part_global,
//                        Array<int>* ien,Array<int>* ief,
//                        std::map<int,int*> face2node,
//                        std::map<int,int*> face2element,
//                        Array<int>* if_ref,
//                        std::set<int> ushell, MPI_Comm comm)
//{
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//
//    std::map<int,int> sharedFaces;
//
//    int* ielement_offsets     = new int[world_size];
//    int* ielement_nlocs       = new int[world_size];
//
//    int elTel = 0;
//
////    std::map<int,int> lV2gV_tets;
////    std::map<int,int> gV2lV_tets;
//    std::set<int> gvid_set;
//    int nLocalVerts = 0;
//
////    std::map<int,int> lF2gF_tets;
////    std::map<int,int> gF2lF_tets;
//    std::set<int> gfid_set;
//    int nLocalFaces = 0;
//    std::map<int,int> face2ref;
//    int ref,gfid,gvid,el0,el1,r0,r1;
//    int shf = 0;
//    int nbcfaces = 0;
//    std::set<int> ufaces;
//
//    for(int i=0;i<ien->getNrow();i++)
//    {
//        //key = GID, value global node nmber;
//        //int gEl = ite->first;
//
//        for(int j=0;j<4;j++)
//        {
//            gfid = ief->getVal(i,j);
//
//            if(ushell.find(gfid)!=ushell.end())
//            {
//                ref = 13;
//            }
//            else
//            {
//                ref = if_ref->getVal(i,j);
//            }
//
//            if(gfid_set.find(gfid)==gfid_set.end())
//            {
//                gfid_set.insert(gfid);
//                face2ref[gfid] = ref;
//
//                nLocalFaces++;
//            }
//
//            for(int k=0;k<3;k++)
//            {
//                gvid = face2node[gfid][k];
//
//                if(gvid_set.find(gvid)==gvid_set.end())
//                {
//                    gvid_set.insert(gvid);
//                    nLocalVerts++;
//                }
//            }
//
//
//
//
//            if(ufaces.find(gfid)==ufaces.end())
//            {
//                ufaces.insert(gfid);
//
//                el0 = face2element[gfid][0];
//                el1 = face2element[gfid][1];
//
//                std::cout << el0 << " " << el1 << std::endl;
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
//                if(ref!=2)
//                {
//                    nbcfaces++;
//                    //lf++;
//                }
//            }
//        }
//        elTel++;
//    }
//
//    int nSharedFaces   = sharedFaces.size();
//
//    std::cout << "nSharedFaces = " << nSharedFaces << std::endl;
//    int nTetras = ien->getNrow();
//    DistributedParallelState* distTetra       = new DistributedParallelState(nTetras,comm);
//    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
//    DistributedParallelState* distLocalVerts  = new DistributedParallelState(nLocalVerts,comm);
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
//        for(int q=0;q<3;q++)
//        {
//            gvid   = face2node[gfid][q];
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
//
//}



TetrahedraMesh* ExtractTetrahedralMesh(Array<int>* part_global,
                                       std::map<int,std::vector<int> > tetras,
                                       std::map<int,std::vector<int> > ief_part_map,
                                       std::map<int,std::vector<int> > ifn_part_map,
                                       std::map<int,std::vector<int> > ife_part_map,
                                       std::map<int,std::vector<int> > if_ref_part_map,
                                       std::set<int> ushell,
                                       std::map<int,Array<double>* > M_vmap,
                                       MPI_Comm comm)
{
    
	TetrahedraMesh* tmesh = new TetrahedraMesh;
	int world_size;
	MPI_Comm_size(comm, &world_size);
	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(comm, &world_rank);
	    
    int shfn = 0;
    int shf = 0;
    int i;
    int nTetras = tetras.size();
	int gvid,gfid;
	
	std::map<int,std::vector<int> >::iterator ite;
	int r0,r1,el0,el1,pos,ra;
	int ref,nbcfaces;
	int lcv = 0;
	std::map<int,std::vector<int> > ref2bcface;

	int lf  = 0;
	std::map<int,int> sharedFaces;

//  int* red_ielement_nlocs   = new int[world_size];
	int* ielement_offsets     = new int[world_size];
	int* ielement_nlocs       = new int[world_size];

	int elTel = 0;
	
//	std::map<int,int> lV2gV_tets;
//	std::map<int,int> gV2lV_tets;
	std::set<int> gvid_set;
	int nLocalVerts = 0;
	
//	std::map<int,int> lF2gF_tets;
//	std::map<int,int> gF2lF_tets;
	std::set<int> gfid_set;
	int nLocalFaces = 0;
    std::map<int,int> face2ref;
    std::set<int> ufaces;
	for(ite=tetras.begin();ite!=tetras.end();ite++)
	{
		//key = GID, value global node nmber;
		int gEl = ite->first;
		
		for(int j=0;j<4;j++)
		{
			gfid = ief_part_map[gEl][j];
            
            if(ushell.find(gfid)!=ushell.end())
            {
                ref = 13;
            }
            else
            {
                ref = if_ref_part_map[gfid][0];
                //std::cout << ref << std::endl;
            }
            
            
			if(gfid_set.find(gfid)==gfid_set.end())
			{
				gfid_set.insert(gfid);
//				gF2lF_tets[gvid] = lfid;
//				lF2gF_tets[lvid] = gfid;
                face2ref[gfid] = ref;

				nLocalFaces++;
			}
			
			for(int k=0;k<3;k++)
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

				if(ref==2)
				{
					r0 = part_global->getVal(el0,0);
					r1 = part_global->getVal(el1,0);
					
//					if(r0==world_rank && r1==world_rank)
//					{
//						lf++;
//					}
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
				
				if(ref!=2)
				{
					//ref2bcface[ref].push_back(gfid);
					nbcfaces++;
					lf++;
				}
			}
		}
		elTel++;
	}
	
	int nSharedFaces   = sharedFaces.size();
	DistributedParallelState* distTetra 	  = new DistributedParallelState(nTetras,comm);
	DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
	DistributedParallelState* distLocalVerts  = new DistributedParallelState(nLocalVerts,comm);
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

		for(int q=0;q<3;q++)
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
	std::set<int> f2r_s;
	int* owned_faces = new int[world_size];
	for(int u=0;u<world_size;u++)
	{
		owned_faces[u] = 0;
	}
	
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
				owned_faces[val]=owned_faces[val]+1;
			}
			if(val>tmp)
			{
				f2r[key]=tmp;
				owned_faces[tmp]=owned_faces[tmp]+1;
			}
		}
	}
	
	std::map<int,int> v2r;
	std::set<int> v2r_s;
	int* owned_verts = new int[world_size];
	for(int u=0;u<world_size;u++)
	{
		owned_verts[u] = 0;
	}
	for(int i=0;i<Nt_shVerts;i++)
	{
		int key = TotalSharedVerts[i];
		int val = TotalSharedVerts_RankID[i];
		
		if(v2r_s.find(key)==v2r_s.end())
		{
			v2r_s.insert(key);
			v2r[key]=val;
			NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
		}
		else
		{
			tmp = v2r[key];
            
			if(val<tmp)
			{
				v2r[key]=val;
				owned_verts[val]=owned_verts[val]+1;

			}
			if(val>tmp)
			{
				v2r[key]=tmp;
				owned_verts[tmp]=owned_verts[tmp]+1;
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
//	int nNonSharedVerts 		 = nLocalVerts-owned_verts[world_rank];
//	int nNonSharedFaces 		 = nLocalFaces-owned_faces[world_rank];
    
    
    int nNonSharedVerts          = nLocalVerts-nSharedVerts;
    int nNonSharedFaces          = nLocalFaces-nSharedFaces;
    
    //std::cout << "LOCALVERTS " << world_rank << "   -- = " << nNonSharedVerts <<" "<< nLocalVerts << " " << owned_verts[world_rank] << " " << nSharedVerts << std::endl;

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
			nNonSharedArray[i] 		= nNonSharedVerts;
			nNonSharedFacesArray[i] = nNonSharedFaces;
		}
	}
	
	MPI_Allreduce(nNonSharedArray,
				  nNonSharedArrayRed,
				  world_size,
				  MPI_INT, MPI_SUM, comm);
	
	MPI_Allreduce(nNonSharedFacesArray,
				  nNonSharedFacesArrayRed,
				  world_size,
				  MPI_INT, MPI_SUM, comm);
	
	int nonSharedOff 		= 0;
	int nonFacesSharedOff 	= 0;
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
	MPI_Allreduce(ini_nEl, red_ini_nEl, world_size, MPI_INT, MPI_SUM, comm);

	int offsetEl = 0;
	
	for(int i=0;i<world_size;i++)
	{
		ini_offsetEl[i] = offsetEl;
		offsetEl        = offsetEl + red_ini_nEl[i];
	}
	
	int size = world_size;
	int optimalSize = int(offsetEl/size) + ( world_rank < offsetEl%size );
	double rat = (double)nTetras/(double)optimalSize;
		
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
    int* to_Send_copy = new int[world_size];
    int* to_Recv_copy = new int[world_size];
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
    MPI_Allreduce(optiSize, optiSize_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(toS,           toS_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(toR,           toR_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&old_ntets[0], &old_ntets_red[0], world_size, MPI_INT, MPI_SUM, comm);
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
            int lEl       = ini_offsetEl[world_rank]+u;
            int* ien      = new int[4];
            int* ien_o    = new int[4];
            int* ief      = new int[4];
            int* ief_o    = new int[4];
            int* iefref   = new int[4];
            double* Met   = new double[6*4];
            
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
                if(tmesh->M_vmap.find(ien[s])==tmesh->M_vmap.end())
                {
                    Array<double>* Mtrcs = new Array<double>(6,1);
                    Mtrcs->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                    Mtrcs->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                    Mtrcs->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                    Mtrcs->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                    Mtrcs->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                    Mtrcs->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                    
                    tmesh->M_vmap[ien[s]]=Mtrcs;
                }
            }
            delete[] Met;
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
		int sRank    = toRanks[0];
		
		int offPrank = 0;
		int cntv     = 0;
		
		int t = 0;
		int nuloc    = 0;
		int uloc     = 0;
        int u        = 0;

		for(ite=tetras.begin();ite!=tetras.end();ite++)
		{
			int nelPrank  = NeltoRanks[cc];
			int sRank     = toRanks[cc];
			int gEl       = ite->first;
			int lEl       = ini_offsetEl[world_rank]+u;
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
                    
                    //std::cout << iefref[0] << " " << iefref[1] << " " << iefref[2] << " " << iefref[3] << std::endl;
                    
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
                    if(tmesh->M_vmap.find(ien[s])==tmesh->M_vmap.end())
                    {
                        Array<double>* Mtrcs = new Array<double>(6,1);
                        Mtrcs->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                        Mtrcs->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                        Mtrcs->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                        Mtrcs->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                        Mtrcs->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                        Mtrcs->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                        
                        tmesh->M_vmap[ien[s]]=Mtrcs;
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
			MPI_Send(&n_Vrt         ,     1,     MPI_INT, dest, dest,         comm);
			MPI_Send(&Elvec[0]      , n_Vrt,     MPI_INT, dest, dest*100,   comm);
			MPI_Send(&ElFvec[0]     , n_Vrt,     MPI_INT, dest, dest*500,   comm);
			MPI_Send(&Elovec[0]     , n_Vrt,     MPI_INT, dest, dest*200,   comm);
            MPI_Send(&ElFrefvec[0]  , n_Vrt,     MPI_INT, dest, dest*20000, comm);
            MPI_Send(&ElFovec[0]    , n_Vrt,     MPI_INT, dest, dest*40000, comm);
            MPI_Send(&n_Ele         ,     1,     MPI_INT, dest, dest*50000, comm);
            MPI_Send(&ElIDs[0]      , n_Ele,     MPI_INT, dest, dest*60000, comm);
            MPI_Send(&n_VrtMet      ,     1,     MPI_INT, dest, dest*70000, comm);
            MPI_Send(&Vmetrics[0]   , n_VrtMet,  MPI_DOUBLE, dest, dest*80000, comm);

			acull = acull + n_Vrt;
		}
        
        //std::cout << "Stored at first  " << world_rank << " " << tmesh->M_vmap.size() << std::endl;
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
			MPI_Recv(&n_Elr,   1, MPI_INT, origin, world_rank, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvNElVec(n_Elr);
			MPI_Recv(&recvNElVec[0], n_Elr, MPI_INT, origin, world_rank*100, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvFElVec(n_Elr);
			MPI_Recv(&recvFElVec[0], n_Elr, MPI_INT, origin, world_rank*500, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvONElVec(n_Elr);
			MPI_Recv(&recvONElVec[0],n_Elr, MPI_INT, origin, world_rank*200, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFrefElVec(n_Elr);
            MPI_Recv(&recvFrefElVec[0], n_Elr, MPI_INT, origin, world_rank*20000, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFONElVec(n_Elr);
            MPI_Recv(&recvFONElVec[0], n_Elr, MPI_INT, origin, world_rank*40000, comm, MPI_STATUS_IGNORE);
			
            int n_EleRecv;
            MPI_Recv(&n_EleRecv,   1, MPI_INT, origin, world_rank*50000, comm, MPI_STATUS_IGNORE);

            std::vector<int> recvElIDsvec(n_EleRecv);
            MPI_Recv(&recvElIDsvec[0],   n_EleRecv, MPI_INT, origin, world_rank*60000, comm, MPI_STATUS_IGNORE);
            
            int n_VrtMetRecv;
            MPI_Recv(&n_VrtMetRecv,   1, MPI_INT, origin, world_rank*70000, comm, MPI_STATUS_IGNORE);

            std::vector<double> VrtMetric(n_VrtMetRecv);
            MPI_Recv(&VrtMetric[0],   n_VrtMetRecv, MPI_DOUBLE, origin, world_rank*80000, comm, MPI_STATUS_IGNORE);

            collected_ElIds[origin]         = recvElIDsvec;
			collected_NIds[origin] 			= recvNElVec;
			collected_OriginalNIds[origin] 	= recvONElVec;
			collected_FIds[origin] 			= recvFElVec;
            collected_Frefs[origin]         = recvFrefElVec;
            collected_FOriginalNIds[origin] = recvFONElVec;
            collected_VrtMetrics[origin]    = VrtMetric;

		}
		
		int el = 0;
        int u  = 0;

		for(ite=tetras.begin();ite!=tetras.end();ite++)
		{	
			int gEl       = ite->first;
			int lEl       = ini_offsetEl[world_rank]+u;
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
            new_iefref->setVal(u,0,iefref[1]);
            new_iefref->setVal(u,0,iefref[2]);
            new_iefref->setVal(u,0,iefref[3]);
            
            //std::cout << iefref[0] << " " << iefref[1] << " " << iefref[2] << " " << iefref[3] << std::endl;
            
            for(int s=0;s<4;s++)
            {
                if(tmesh->M_vmap.find(ien[s])==tmesh->M_vmap.end())
                {
                    Array<double>* Mtrc = new Array<double>(6,1);
                    
                    Mtrc->setVal(0,0,M_vmap[ien_o[s]]->getVal(0,0));
                    Mtrc->setVal(1,0,M_vmap[ien_o[s]]->getVal(0,1));
                    Mtrc->setVal(2,0,M_vmap[ien_o[s]]->getVal(0,2));
                    Mtrc->setVal(3,0,M_vmap[ien_o[s]]->getVal(1,1));
                    Mtrc->setVal(4,0,M_vmap[ien_o[s]]->getVal(1,2));
                    Mtrc->setVal(5,0,M_vmap[ien_o[s]]->getVal(2,2));
                    
                    tmesh->M_vmap[ien[s]] = Mtrc;
                    
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
            
            int fromRank = collit->first;
			for(int q=0;q<nel;q++)
			{
                new_GidEl->setVal(u,0,collected_ElIds[collit->first][q]);
                
				new_ien->setVal(u,0,collit->second[q*4+0]);
				new_ien->setVal(u,1,collit->second[q*4+1]);
				new_ien->setVal(u,2,collit->second[q*4+2]);
				new_ien->setVal(u,3,collit->second[q*4+3]);
                
                for(int s=0;s<4;s++)
                {
                    if(tmesh->M_vmap.find(collit->second[q*4+s])==tmesh->M_vmap.end())
                    {
                        Array<double>* Mtrc = new Array<double>(6,1);
                        Mtrc->setVal(0,0,collected_VrtMetrics[collit->first][6*4*q+6*s+0]);
                        Mtrc->setVal(1,0,collected_VrtMetrics[collit->first][6*4*q+6*s+1]);
                        Mtrc->setVal(2,0,collected_VrtMetrics[collit->first][6*4*q+6*s+2]);
                        Mtrc->setVal(3,0,collected_VrtMetrics[collit->first][6*4*q+6*s+3]);
                        Mtrc->setVal(4,0,collected_VrtMetrics[collit->first][6*4*q+6*s+4]);
                        Mtrc->setVal(5,0,collected_VrtMetrics[collit->first][6*4*q+6*s+5]);
                        nadded++;
                        tmesh->M_vmap[collit->second[q*4+s]] = Mtrc;
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
				
//                std::cout << "aff " << collected_Frefs[collit->first][4*q+0] << " " << collected_Frefs[collit->first][4*q+1] << " " << collected_Frefs[collit->first][4*q+2] << " " << collected_Frefs[collit->first][4*q+3] << std::endl;
                
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
    
    
    

    
	// return these guys.
    tmesh->ElGids             = new_GidEl;
	tmesh->ien_part_tetra     = new_ien;
	tmesh->ien_part_hybrid    = new_ien_or;
	tmesh->ief_part_tetra     = new_ief;
    tmesh->ief_part_hybrid    = new_ief_or;
    tmesh->iefref_part_tetra  = new_iefref;
    
    
     
    return tmesh;
	
}


PartitionInfo* GetNewGlobalPartitioningTetrahedraMesh(TetrahedraMesh* tmesh, MPI_Comm comm)
{
    int i;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    Array<int>* new_ien     = tmesh->ien_part_tetra;
    Array<int>* new_ief     = tmesh->ief_part_tetra;
    Array<int>* new_ien_or  = tmesh->ien_part_hybrid;
    //std::cout << "WORLD_RANK = " << world_rank << " " << lbvid << std::endl;
    
    int optimalSize = new_ien->getNrow();
    int* ielement_nlocs = new int[world_size];
    int* ielement_offsets = new int[world_size];

    int* red_ielement_nlocs = new int[world_size];

    int* elmdist = new int[world_size+1];

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
    
    //std::cout << " optimalSize " << world_rank << " " << optimalSize << std::endl;
    
    MPI_Allreduce(ielement_nlocs,
                  &red_ielement_nlocs[0],
                  world_size, MPI_INT,
                  MPI_SUM, comm);
    
    int o_ie = 0;
    
    for(i=0;i<world_size;i++)
    {
        ielement_offsets[i] = o_ie;
        ielement_nlocs[i]   = red_ielement_nlocs[i];
        elmdist[i]          = o_ie;
        o_ie                = o_ie+red_ielement_nlocs[i];
//        if(world_rank==0)
//        {
//            std::cout << "elmdist[i]  " << elmdist[i] << std::endl;
//
//        }
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
            eind[j]    = new_ien->data[j];
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
                          &xadj_par,&adjncy_par,&comm);
    
    ParMETIS_V3_PartKway(elmdist,
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);

    Array<int>* part_global_new  = new Array<int>(o_ie,1);
    Array<int>* part_new         = new Array<int>(optimalSize,1);

    part_new->data = part_arr;
    
    MPI_Allgatherv(&part_new->data[0],
                   red_ielement_nlocs[world_rank], MPI_INT,
                   &part_global_new->data[0],
                   ielement_nlocs,
                   ielement_offsets,
                   MPI_INT,comm);
    
    PartitionInfo* pinfo = new PartitionInfo;
    
    pinfo->part = part_new;
    pinfo->part_global = part_global_new;
    return pinfo;
}

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    
    int ier,opt;
    
    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    const char* fn_metric = "metric.inp";
    
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
    int ReadFromStats = 0;
    if(metric_inputs.size()==6)
    {
        ReadFromStats=metric_inputs[5];
    }
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve = us3d->xcn->getNglob();
    
    int Nel_part = us3d->ien->getNrow();
    
    //Array<double>* xcn_ref = ReadDataSetFromFile<double>(fn_grid,"xcn");
    //Array<int>* ien_ref    = ReadDataSetFromFile<int>(fn_conn,"ien");

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    for(int i=0;i<Nel_part;i++)
    {
        rhoState = us3d->interior->getVal(i,0);
        uState   = us3d->interior->getVal(i,1);
        vState   = us3d->interior->getVal(i,2);
        wState   = us3d->interior->getVal(i,3);
        TState   = us3d->interior->getVal(i,4);
        VtotState = sqrt(uState*uState+vState*vState+wState*wState);
        aState   = sqrt(1.4*287.05*TState);
        MState = VtotState/aState;
        Ui->setVal(i,0,MState);
    }
    
    delete us3d->interior;

    
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
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);

    std::vector<int> LocElem = P->getLocElem();
    std::vector<double> Uvaria  = P->getLocElemVaria();
    std::map<int,Array<double>*> Uvaria_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid = LocElem[i];
        UvariaV   = Uvaria[i];
        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map[gid] = Uarr;
    }
    
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
    
    std::map<int,double> Volumes = meshTopo->getVol();
    P->AddStateVecForAdjacentElements(Uvaria_map,1,comm);
    
    std::map<int,Array<double>* > var_vmap = P->ReduceStateVecToAllVertices(Uvaria_map,1);
    std::map<int,Array<double>* >::iterator vm;
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }
    
    delete us3d->ghost;
    
    
    //std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,Uvaria_map,gB,comm);
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_Vrt_US3D(P,Uvaria_map,var_vmap,meshTopo,gB,comm);
    
    P->AddStateVecForAdjacentElements(dUdXi,3,comm);
    
    std::map<int,Array<double>* >::iterator grit;
    std::map<int,Array<double>* > dUidxi_map;
    std::map<int,Array<double>* > dUidyi_map;
    std::map<int,Array<double>* > dUidzi_map;
    std::map<int,Array<double>* > dUidXi_map;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        
        Array<double>* dUdx_E = new Array<double>(1,1);
        dUdx_E->setVal(0,0,grit->second->getVal(0,0));
        Array<double>* dUdy_E = new Array<double>(1,1);
        dUdy_E->setVal(0,0,grit->second->getVal(1,0));
        Array<double>* dUdz_E = new Array<double>(1,1);
        dUdz_E->setVal(0,0,grit->second->getVal(2,0));
        
        dUidxi_map[grit->first]=dUdx_E;
        dUidyi_map[grit->first]=dUdy_E;
        dUidzi_map[grit->first]=dUdz_E;
        
        delete grit->second;
    }

    std::map<int,Array<double>* > dudx_vmap = P->ReduceStateVecToAllVertices(dUidxi_map,1);
    std::map<int,Array<double>* > dudy_vmap = P->ReduceStateVecToAllVertices(dUidyi_map,1);
    std::map<int,Array<double>* > dudz_vmap = P->ReduceStateVecToAllVertices(dUidzi_map,1);


    //std::cout << "second gradient "<<std::endl;
//    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUidxi_map,gB,comm);
//    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUidyi_map,gB,comm);
//    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUidzi_map,gB,comm);

    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidxi_map,dudx_vmap,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidyi_map,dudy_vmap,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidzi_map,dudz_vmap,meshTopo,gB,comm);

//      Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,gB,comm);
            
    std::map<int,Array<double>*> Hess_map;
    //std::cout << "second gradient 2"<<std::endl;

    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;
    
    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
        
        Array<double>* Hess = new Array<double>(6,1);
        
        Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
        Hess->setVal(1,0,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(2,0,dU2dXi2[gid]->getVal(2,0));

        Hess->setVal(3,0,dU2dYi2[gid]->getVal(1,0));
        Hess->setVal(4,0,dU2dYi2[gid]->getVal(2,0));
        Hess->setVal(5,0,dU2dZi2[gid]->getVal(2,0));
        
        Hess_map[gid] = Hess;
        
        delete dU2dXi2[gid];
        delete dU2dYi2[gid];
        delete dU2dZi2[gid];
        
        t++;
    }
    
    dU2dXi2.clear();
    dU2dYi2.clear();
    dU2dZi2.clear();
    
    
    double* Hessie = new double[9];
    double * WRn = new double[3];
    Array<double>* DR  = new Array<double>(3,3);
    Array<double>* UR  = new Array<double>(3,3);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++
    double hmin         = metric_inputs[1];
    double hmax         = metric_inputs[2];
    double f            = metric_inputs[3];
    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double cmplxty = 0.0;
    double cmplxty2 = 0.0;
    double Volu=0.0,cmplxty_red=0.0,cmplxty_tmp=0.0,cmplxty_tmp2=0.0;
    double po = 6.0;
    for(int j=0;j<3;j++)
    {
        for(int k=0;k<3;k++)
        {
            DR->setVal(j,k,0.0);
        }
    }
    for(itgg=Hess_map.begin();itgg!=Hess_map.end();itgg++)
    {
        Hessie[0] = itgg->second->getVal(0,0);
        Hessie[1] = itgg->second->getVal(1,0);
        Hessie[2] = itgg->second->getVal(2,0);
        
        Hessie[3] = itgg->second->getVal(1,0);
        Hessie[4] = itgg->second->getVal(3,0);
        Hessie[5] = itgg->second->getVal(4,0);
        
        Hessie[6] = itgg->second->getVal(2,0);
        Hessie[7] = itgg->second->getVal(4,0);
        Hessie[8] = itgg->second->getVal(5,0);
        
        Eig* eig = ComputeEigenDecomp(3, Hessie);
        
        WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[2] = std::min(std::max(f*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);
        
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        double detRf = sqrt( Rf->getVal(0,0)*(Rf->getVal(1,1)*Rf->getVal(2,2)-Rf->getVal(2,1)*Rf->getVal(1,2))
                            -Rf->getVal(0,1)*(Rf->getVal(1,0)*Rf->getVal(2,2)-Rf->getVal(2,0)*Rf->getVal(1,2))
                            +Rf->getVal(0,2)*(Rf->getVal(1,0)*Rf->getVal(2,1)-Rf->getVal(2,0)*Rf->getVal(1,1)));
        
        Volu = Volumes[itgg->first];
        cmplxty_tmp = detRf;
        cmplxty_tmp2 = Volu;
        //std::pow(cmplxty_tmp,(po+1.0)/(2.0*po+3.0));
        cmplxty = cmplxty + cmplxty_tmp;
        cmplxty2 = cmplxty2 + cmplxty_tmp2;
        delete iVR;
        delete Rs;
        delete Rf;
    }
    

    MPI_Allreduce(&cmplxty, &cmplxty_red, 1, MPI_DOUBLE, MPI_SUM, comm);
    //cmplxty_red = std::pow(cmplxty_red,(po+1)/(2.0*po+3.0));

    P->AddStateVecForAdjacentElements(Hess_map,6,comm);

    std::map<int,Array<double>* > hess_vmap = P->ReduceStateVecToAllVertices(Hess_map,6);

    cmplxty_red = Nve/cmplxty_red;
    
    ComputeMetric(P, metric_inputs, comm, hess_vmap, 1.0, po);
    
    //std::cout << "understand bitshift " << world_rank << "  " << (world_rank & (~(1 << 3)))<< std::endl;
    //std::vector<int> LocElem     = P->getLocElem();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Domain* pDom = P->getPartitionDomain();
    std::map<int,std::vector<int> > tetras     = pDom->GTetras;
    std::map<int,std::vector<int> > prisms     = pDom->GPrisms;
    std::set<int> ushell                       = pDom->ushell;
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* ien_part_map                   = P->getIENpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    Array<int>* part_global                    = P->getGlobalPartition();
    //std::map<int,Array<double>* > M_vmap;
    // Based on a partitioned hybrid mesh, we are extracting the tetra and redistribute them
    // uniformly. The hybrid mesh is partitioned without any weights which means that it might happen that the tetrahedra are initial distributed non-uniformly where a significant number of processors end up with any tetrahedra. This function distributed them equally.
    //std::cout << "Starting extracting..." << tetras.size() <<std::endl;
    
    TetrahedraMesh* tmesh = ExtractTetrahedralMesh(part_global, tetras,
                                                   ief_part_map->i_map,
                                                   ifn_part_map->i_map,
                                                   ife_part_map->i_map,
                                                   if_ref_part_map->i_map,
                                                   ushell,
                                                   hess_vmap,
                                                   comm);

    //std::cout << "Done extracting..."<<std::endl;

    // Once we have a uniform distribution of the tetrahedra, we determine a new OPTIMAL partitioning using the PARMETIS routine.
    
    PartitionInfo* partInfo = GetNewGlobalPartitioningTetrahedraMesh(tmesh,comm);

    // Determining the OPTIMAL number of elements on each rank based on the PARMETIS partitioning.

    //std::cout << "partInfo->part " << world_rank << " " << partInfo->part->getNrow() << std::endl;

    int* newSizesOnRanks    = CommunicatePartitionLayout(partInfo->part,comm);
//
    // Grab the optimal number of elements on current rank based on the partitioning array.
    int nTetrahedra         = newSizesOnRanks[world_rank];
    int nTetrahedraGlob     = 0;
    int* nTetOffsets        = new int[world_size];

    int tetoff = 0;
    for(int i=0;i<world_size;i++)
    {
        nTetrahedraGlob=nTetrahedraGlob+newSizesOnRanks[i];
        nTetOffsets[i] = tetoff;
        tetoff = tetoff + newSizesOnRanks[i];
    }


    int nLocTetOffset = nTetOffsets[world_rank];

    // Grab the global number of elements in the mesh.
    // int GlobalNel           = part_global->getNrow();


    // Based on the PARMETIS partitioning (Array<int>* part_new) and the uniform tetrahedra distribution (tmesh) and the unformly read vertices (xcn) and the face maps (ifn), we communicate the appropriate data such that each processor can compute on tetrahedra while minimizing the number of shared faces between partitions.

    //std::cout << "FIRST + " << world_rank << " " << nTetrahedraGlob << " " << nTetrahedra << std::endl;

    std::map<int,Array<double>* >::iterator itA;

    for(itA=tmesh->M_vmap.begin();itA!=tmesh->M_vmap.end();itA++)
    {
        if(itA->second->getNrow() != 6)
        {
            std::cout << itA->first << " " << itA->second->getNrow() <<std::endl;

        }
    }

    int gg;
    int ll = 0;
    std::map<int,int> g2lv;
    for(int i=0;i<tmesh->ien_part_tetra->getNrow();i++)
    {
        for(int j=0;j<tmesh->ien_part_tetra->getNcol();j++)
        {
            int gg = tmesh->ien_part_tetra->getVal(i,j);
            if(g2lv.find(gg)==g2lv.end())
            {
                g2lv[gg]=ll;
                ll++;
            }
        }
    }
    
    

    if(g2lv.size()!=tmesh->M_vmap.size())
    {
        std::cout <<"Mvert size " <<  world_rank << " -> " << g2lv.size() << "  " << tmesh->M_vmap.size() << std::endl;
    }


    UpdateTetrahedraOnPartition(nTetrahedraGlob, nTetrahedra,
                                partInfo->part, tmesh,
                                us3d->xcn, xcn_pstate,
                                ife_pstate,
                                comm);
    
    
    
     
     
     
     
    ParallelState* ien_pstate_tet      = new ParallelState(nTetrahedraGlob,comm);

    
    std::map<int,int*> face2node = GetFace2EntityTetrahedraMesh(tmesh, us3d->ifn, 3, ife_pstate,nTetrahedraGlob,comm);
    
    std::map<int,std::vector<int> > face2rank = GetFace2RankTetrahedraMesh(tmesh,ien_pstate_tet,comm);
    
    PartitionBoundary* pb = ExtractPartitionBoundary(tmesh,face2rank,face2node,comm);
    
    std::vector<Vert*> locVs = tmesh->LocalVerts;
    
    int* color_face = (int *) malloc(pb->ncomm*sizeof(int));
    int* ntifc = (int *) malloc(pb->ncomm*sizeof(int));
    int* ifc_tria_glob[pb->ncomm];
    int* ifc_tria_loc[pb->ncomm];
    
    int icomm=0;
    std::map<int,std::vector<int> >::iterator itc;
    std::set<int> uniShVrts;
    std::set<int> uniShVrts_owned;
    for(itc=pb->ColorsFaces.begin();itc!=pb->ColorsFaces.end();itc++)
    {
        color_face[icomm]     = itc->first;
        ntifc[icomm]          = itc->second.size();
        ifc_tria_loc[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));
        ifc_tria_glob[icomm]  = (int *) malloc(itc->second.size()*sizeof(int));
        
        for(int q=0;q<itc->second.size();q++)
        {
            ifc_tria_glob[icomm][q] = itc->second[q]+1;
            ifc_tria_loc[icomm][q]  = pb->globShF2locShF[itc->second[q]]+1;
            
            if(world_rank<color_face[icomm])
            {
                for(int l=0;l<3;l++)
                {
                    int vtu = face2node[itc->second[q]][l];
                    
                    if(uniShVrts_owned.find(vtu)==uniShVrts_owned.end())
                    {
                        uniShVrts_owned.insert(vtu);
                        uniShVrts.insert(vtu);
                    }
                }
            }
            else
            {
                for(int l=0;l<3;l++)
                {
                    int vtu = face2node[itc->second[q]][l];
                    
                    if(uniShVrts.find(vtu)==uniShVrts.end())
                    {
                        uniShVrts.insert(vtu);
                    }
                }
            }
        }
        icomm++;
    }
    
    std::set<int> intVrts;
    for(int i=0;i<tmesh->ien_part_tetra->getNrow();i++)
    {
        for(int j=0;j<4;j++)
        {
            int vid = tmesh->ien_part_tetra->getVal(i,j);
            if(intVrts.find(vid)==intVrts.end() && uniShVrts.find(vid)==uniShVrts.end())
            {
                intVrts.insert(vid);
            }
        }
    }
    
    
    std::cout << "fa " << world_rank <<" checking " << uniShVrts_owned.size() << " " << uniShVrts.size() << " " << intVrts.size() << std::endl;
    
    //Based on the new local tetrahedra mesh, we output a tecplot file per processor that has the geometry of the computational domain that is owned by world_rank.
    
    OutputTetrahedralMeshOnPartition(tmesh,comm);
    
    Array<int>* ien_part_tetra   = tmesh->ien_part_tetra;
    Array<int>* ien_part_hybrid  = tmesh->ien_part_hybrid;
    Array<int>* ief_part_tetra   = tmesh->ief_part_tetra;
    Array<int>* ief_part_hybrid  = tmesh->ief_part_hybrid;
    int nVertices                = locVs.size();
    int nTriangles               = pb->nPartBoundFaces;
    int nEdges                   = 0;
    int nPrisms                  = 0;
    int nQuadrilaterals          = 0;
//
//
//
//
    PMMG_pParMesh   parmesh;
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);
//
//
//
//
    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                              nQuadrilaterals,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
//
//
//
//
    //PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol)
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
//    const int nSolsAtVertices = 1; // 3 solutions per vertex
//    int solType[1] = {MMG5_Tensor};
//
//    if ( PMMG_Set_solsAtVerticesSize(parmesh,nSolsAtVertices,nVertices,solType) != 1 ) {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    }
    
    for ( k=0; k<nVertices; ++k )
    {
          double vx = locVs[k]->x;
          double vy = locVs[k]->y;
          double vz = locVs[k]->z;
            
          if ( PMMG_Set_vertex(parmesh,vx,vy,vz, 1.0, k+1) != 1 )
          {
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }

          double* tensor = new double[6];
          int vert = tmesh->locV2globV[k];
        
          tensor[0] = tmesh->M_vmap[vert]->getVal(0,0);
          tensor[1] = tmesh->M_vmap[vert]->getVal(1,0);
          tensor[2] = tmesh->M_vmap[vert]->getVal(2,0);
          tensor[3] = tmesh->M_vmap[vert]->getVal(3,0);
          tensor[4] = tmesh->M_vmap[vert]->getVal(4,0);
          tensor[5] = tmesh->M_vmap[vert]->getVal(5,0);
          
        //if(k==0)
        //{
          //  std::cout << world_rank <<" :: "<<vert << " " <<tensor[0] << " " << tensor[1] << " " << tensor[2] << " " << tensor[3] << " " << tensor[4] << " " << tensor[5] << std::endl;
        //}
          
          if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],k+1)!=1)
          {
             MPI_Finalize();
             exit(EXIT_FAILURE);
          }
    }
   
    int v0,v1,v2,v3;
    int v0l,v1l,v2l,v3l;
    int teller = 0;
    int refer  = 0;
    int cref36 = 0;
    for ( k=0; k<nTriangles; ++k )
    {
        int faceID = pb->faces4parmmg[k];

        v0      = face2node[faceID][0];
        v1      = face2node[faceID][1];
        v2      = face2node[faceID][2];
        
        v0l     = tmesh->globV2locV[v0];
        v1l     = tmesh->globV2locV[v1];
        v2l     = tmesh->globV2locV[v2];
        
        refer = tmesh->face2ref[faceID];

        if(refer == 36)
        {
            cref36++;
        }
        if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,refer,k+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    std::cout << "cref36 before "  << cref36 << " -> " << world_rank << std::endl;
    
    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, pb->ncomm);
    
    for(int icomm=0; icomm<pb->ncomm; icomm++ ) {

      // Set nb. of entities on interface and rank of the outward proc
     
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

      // Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               ifc_tria_loc[icomm],
                                               ifc_tria_glob[icomm], 1);
    }
    
    
    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    for ( int t = 0;t < nTetrahedra; t++  )
    {
        v0 = ien_part_tetra->getVal(t,0);
        v1 = ien_part_tetra->getVal(t,1);
        v2 = ien_part_tetra->getVal(t,2);
        v3 = ien_part_tetra->getVal(t,3);
        
        v0l = tmesh->globV2locV[v0];
        v1l = tmesh->globV2locV[v1];
        v2l = tmesh->globV2locV[v2];
        v3l = tmesh->globV2locV[v3];
        
        double* P = new double[4*3];
        
        P[0*3+0]=locVs[v0l]->x;   P[0*3+1]=locVs[v0l]->y;    P[0*3+2]=locVs[v0l]->z;
        P[1*3+0]=locVs[v1l]->x;   P[1*3+1]=locVs[v1l]->y;    P[1*3+2]=locVs[v1l]->z;
        P[2*3+0]=locVs[v2l]->x;   P[2*3+1]=locVs[v2l]->y;    P[2*3+2]=locVs[v2l]->z;
        P[3*3+0]=locVs[v3l]->x;   P[3*3+1]=locVs[v3l]->y;    P[3*3+2]=locVs[v3l]->z;

        double Vtet = GetQualityTetrahedra(P);
        
        if ( PMMG_Set_tetrahedron(parmesh,v0l+1,v1l+1,v2l+1,v3l+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    
    int niter = 1;
    
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, 2.0) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if( !PMMG_Set_dparameter( parmesh,  PMMG_DPARAM_hgradreq , -1.0 ) ){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    std::string filename1 = "NotAdapted_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile1;
    myfile1.open(filename1);
    myfile1 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile1 <<"VARIABLES = \"X\", \"Y\", \"Z\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
    myfile1 <<"ZONE N = " << locVs.size() << ", E = " << nTetrahedra << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<locVs.size();i++)
    {
        int vert = tmesh->locV2globV[i];
      
        double M00 = tmesh->M_vmap[vert]->getVal(0,0);
        double M01 = tmesh->M_vmap[vert]->getVal(1,0);
        double M02 = tmesh->M_vmap[vert]->getVal(2,0);
        double M11 = tmesh->M_vmap[vert]->getVal(3,0);
        double M12 = tmesh->M_vmap[vert]->getVal(4,0);
        double M22 = tmesh->M_vmap[vert]->getVal(5,0);
        
        myfile1 << locVs[i]->x << " " << locVs[i]->y << " " << locVs[i]->z << " "
        << M00 << " "
        << M01 << " "
        << M02 << " "
        << M11 << " "
        << M12 << " "
        << M22 << " " << std::endl;
    }

    
    for(int i=0;i<nTetrahedra;i++)
    {
        v0 = ien_part_tetra->getVal(i,0);
        v1 = ien_part_tetra->getVal(i,1);
        v2 = ien_part_tetra->getVal(i,2);
        v3 = ien_part_tetra->getVal(i,3);
        
        v0l = tmesh->globV2locV[v0];
        v1l = tmesh->globV2locV[v1];
        v2l = tmesh->globV2locV[v2];
        v3l = tmesh->globV2locV[v3];
        
        myfile1 <<   v0l+1 << " "
                <<   v1l+1 << " "
                <<   v2l+1 << " "
                <<   v3l+1 << std::endl;
    }


    myfile1.close();
    
    int nVerticesIN   = 0;
    int nTetrahedraIN = 0;
    int nTrianglesIN  = 0;
    int nEdgesIN      = 0;
    if ( PMMG_Get_meshSize(parmesh,&nVerticesIN,&nTetrahedraIN,NULL,&nTrianglesIN,NULL,
                           &nEdgesIN) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
  
    // remeshing function //
    /* Compute output nodes and triangles global numbering */
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    if(ierlib==0 && world_rank == 0)
    {
        std::cout << "SUCCESFULLY adapted the mesh in parallel!" << std::endl;
    }
    else if(ierlib==1 && world_rank == 0)
    {
        std::cout << "FAILED to adapt the mesh in parallel! "<< std::endl;
    }
    
    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;
    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                           &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    int             nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;
    std::map<int,int> glob2locVid;
    for( k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k-1]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k-1;
    }
    
    
    
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,nTetrahedraOUT,nTrianglesOUT,nEdgesOUT),sizeof(int));
    int *ref = (int*)calloc(MAX4(nVerticesOUT,nTetrahedraOUT,nTrianglesOUT,nEdgesOUT),sizeof(int));
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));
    int pos;
    
    std::vector<std::vector<int> > outT;
    
    double *vertOUT = (double*)calloc((nVerticesOUT)*3,sizeof(double));
    
    for ( k=0; k<nVerticesOUT; k++ ) {
          pos = 3*k;
          if ( PMMG_Get_vertex(parmesh,&(vertOUT[pos]),&(vertOUT[pos+1]),&(vertOUT[pos+2]),
                               &(ref[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
          }
    }
    
    
    

    
    
    
    
    
    
    
    
    
    //int debug = 1;
    if(!niter)
    {
        
        std::cout << "Check the input and outputted shard faces." << std::endl;
        
        int **out_tria_loc;
        int *nitem_face_comm;
        int next_face_comm;
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        int *color_node_out,*color_face_out;
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        // Check matching of input interface nodes with the set ones
        
        // Get input triangle nodes
        int** faceNodes2 = (int **) malloc(pb->ncomm*sizeof(int *));
        for( icomm = 0; icomm < pb->ncomm; icomm++ ) {
          faceNodes2[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
          for( i = 0; i < ntifc[icomm]; i++ ) {
              
            int faceID = ifc_tria_loc[icomm][i]-1;
            int faceID2 = pb->locShF2globShF[faceID];

            v0 = face2node[faceID2][0];
            v1 = face2node[faceID2][1];
            v2 = face2node[faceID2][2];
    
            v0l = tmesh->globV2locV[v0];
            v1l = tmesh->globV2locV[v1];
            v2l = tmesh->globV2locV[v2];

            //pos = ifc_tria_loc[icomm][i];
            faceNodes2[icomm][3*i]     = v0l+1; // tria_vert[3*(pos-1)];
            faceNodes2[icomm][3*i+1]   = v1l+1; // tria_vert[3*(pos-1)+1];
            faceNodes2[icomm][3*i+2]   = v2l+1; // tria_vert[3*(pos-1)+2];
          }
        }

        // Check matching of input interface triangles with the set ones
        if( !PMMG_Check_Set_FaceCommunicators(parmesh,pb->ncomm,ntifc,
                                              color_face,faceNodes2) ) {
          printf("### FAILED:: Wrong set face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: Set the correct face communicators\n");
        }
        
        int *ref2       = (int*)calloc(nTriangles,sizeof(int));
        int *required2  = (int*)calloc(nTriangles,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTriangles,sizeof(int));

        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
          ier = PMMG_STRONGFAILURE;
        }
        else
        {
            printf("### SUCCES:: retrieved all mesh triangles\n");
        }

        int** faceNodes_out = (int **) malloc(next_face_comm*sizeof(int *));
        for( icomm = 0; icomm < next_face_comm; icomm++ )
        {
              faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
              for( i = 0; i < nitem_face_comm[icomm]; i++ )
              {
                  int pos = out_tria_loc[icomm][i];
                  faceNodes_out[icomm][3*i]   = triaNodes2[3*(pos-1)];
                  faceNodes_out[icomm][3*i+1] = triaNodes2[3*(pos-1)+1];
                  faceNodes_out[icomm][3*i+2] = triaNodes2[3*(pos-1)+2];
              }
        }

        
        // Check matching of input interface triangles with the output ones
        if( !PMMG_Check_Get_FaceCommunicators(parmesh,pb->ncomm,ntifc,
                                              color_face,faceNodes2,
                                              next_face_comm,nitem_face_comm,
                                              color_face_out,faceNodes_out) ) {
          printf("### FAILED:: Wrong retrieved face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: retrieved the correct face communicators\n");
        }
        
        int lvt = 0;
        int lvt_o = 0;
        std::map<int,int> sharedVrts_Owned;
        std::map<int,int> sharedVert;
        for( icomm = 0; icomm < next_face_comm; icomm++ )
        {
            faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
            int nPartFace                 = nPartFace + nitem_face_comm[icomm];
            //rank2icomm[color_face_out[icomm]] = icomm;

            if(world_rank < color_face_out[icomm])
            {
                int nTshared_owned = nTshared_owned + nitem_face_comm[icomm];

                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft   = out_tria_loc[icomm][i];
                    int reff = ref2[ft-1];

                    for(int k=0;k<3;k++)
                    {
                        int vt = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];

                        if(sharedVrts_Owned.find(vt)==sharedVrts_Owned.end())
                        {
                            sharedVrts_Owned[vt] = lvt_o;
                            sharedVert[vt]       = lvt;

                            lvt++;
                            lvt_o++;
                        }
                    }
                }
            }
            else
            {
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];

                    std::set<int> faceSh;

                    for(int k=0;k<3;k++)
                    {
                        int vt3 = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k]   = triaNodes2[3*(ft-1)+k];
                        faceSh.insert(vt3);
                        if(sharedVert.find(vt3)==sharedVert.end())
                        {
                            sharedVert[vt3] = lvt;
                            lvt++;
                        }
                    }
                }
            }
        }
        
        std::cout << "f;dlgja;a;a " << world_rank <<" checking " << uniShVrts_owned.size() << " " << uniShVrts.size() << " " << intVrts.size() << " " << sharedVert.size() << " " << sharedVrts_Owned.size() << std::endl;
        
    }
    else
    {
        //std::cout << "Outputting the new triangles..." << std::endl;
        
        //===================================================================
        int **out_tria_loc, **out_node_loc;
        int *nitem_face_comm,*nitem_node_comm;
        int next_face_comm, next_node_comm;
        int *color_node_out,*color_face_out;
        
        ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);
        
        color_node_out  = (int *) malloc(next_node_comm*sizeof(int));
        nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
        for( icomm=0; icomm<next_node_comm; icomm++ )
          ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                                 &color_node_out[icomm],
                                                 &nitem_node_comm[icomm]);
        
        // Get IDs of nodes on each interface //
        out_node_loc = (int **) malloc(next_node_comm*sizeof(int *));
        for( icomm=0; icomm<next_node_comm; icomm++ )
          out_node_loc[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
        ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);
        
        
        
        //===================================================================
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        
        PMMG_pGrp      listgrp,grpI;
        listgrp         = parmesh->listgrp;
        grpI            = &listgrp[0];
        int            poi_id_int,poi_id_glo,imsh,k;
        
        std::cout << "world _ R " << world_rank << " "  << grpI->nitem_int_node_comm << std::endl;
        
        std::set<int> un_poi_id_int;
        std::set<int> un_poi_id_glo;

        std::map<int,int> locSh2globSh;
        std::map<int,int> globSh2locSh;
        
        std::map<int,int> globSh2Rank;
        for ( k = 0; k < grpI->nitem_int_node_comm; ++k )
        {
            poi_id_int = grpI->node2int_node_comm_index1[k];
            poi_id_glo = grpI->node2int_node_comm_index2[k];
            locSh2globSh[poi_id_int]=poi_id_glo;
            globSh2locSh[poi_id_int]=poi_id_int;
            globSh2Rank[poi_id_glo]=world_rank;
            
            //std::cout << "world_rank " << world_rank << " " << poi_id_int << " " << poi_id_glo << std::endl;
            
            if(un_poi_id_int.find(poi_id_int)==un_poi_id_int.end())
            {
                un_poi_id_int.insert(poi_id_int);
            }
            
            if(un_poi_id_glo.find(poi_id_glo)==un_poi_id_glo.end())
            {
                un_poi_id_glo.insert(poi_id_glo);
            }
        }
        
        
        
        
        
        std::cout << "W_R " << world_rank<< " " << un_poi_id_int.size() << " " << un_poi_id_glo.size() << std::endl;
        
        for( icomm=0; icomm<next_node_comm; icomm++ )
        {
            //std::cout << "world rank NODES " << world_rank << " " << color_node_out[icomm] << " " << nitem_node_comm[icomm] << " " << nVerticesOUT<< std::endl;
//            for(int j=0;j<nitem_node_comm[icomm];j++)
//            {
//                if(out_node_loc[icomm][j]>=nVerticesOUT)
//                {
//                    std::cout << "test " << j << " " << out_node_loc[icomm][j] << std::endl;
//                }
//
//            }
        }
          

        for( icomm=0; icomm<next_face_comm; icomm++ )
        {
            //std::cout << "world rank FACE " << world_rank << " " << color_face_out[icomm] << " " << nitem_face_comm[icomm] << std::endl;
        }
          
        
        
        
        
        
        // Check matching of input interface nodes with the set ones
        int *ref2       = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *required2  = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTrianglesOUT,sizeof(int));
        
//        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
//          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
//          ier = PMMG_STRONGFAILURE;
//        }
        int nreq2   = 0;
        int pos     = 0;
        int cnt_bnd = 0;
        int cnt_th  = 0;
        int cnt_int = 0;
        std::map<std::set<int>, int> PMMG_Face2Ref;
        std::set<std::set<int> > TotalFaces;
        std::set<std::set<int> > TotalFaces2;
        std::set<std::set<int> > TotalFaces3;

        
        int c36 = 0;
        for ( k=0; k<nTrianglesOUT; k++ )
        {
            pos = 3*k;
            if ( PMMG_Get_triangle(parmesh,&(triaNodes2[pos]),
                                         &(triaNodes2[pos+1]),
                                         &(triaNodes2[pos+2]),
                                         &(ref2[k]),
                                         &(required2[k])) != 1 )
            {
            fprintf(inm,"Unable to get mesh triangle %d \n",k);
            ier = PMMG_STRONGFAILURE;
            }

            std::set<int> faceSh;
            
            faceSh.insert(triaNodes2[pos]);
            faceSh.insert(triaNodes2[pos+1]);
            faceSh.insert(triaNodes2[pos+2]);
            TotalFaces.insert(faceSh);
            
            if(ref2[k]==0)
            {
               TotalFaces2.insert(faceSh);
            }

            if(ref2[k]!=0)
            {
              PMMG_Face2Ref[faceSh] = ref2[k];
              TotalFaces3.insert(faceSh);
            }
         
            faceSh.clear();
    
            if ( required2 && required2[k] )  nreq2++;
        }
        
        
        
        
        int itt2            = 0;
        int nTshared_owned  = 0;

        std::map<int,int> sharedVert;
        std::map<int,int> sharedVrts_Owned;
        std::set<int> sharedVert_NotOwned;

        std::map<int,int> sharedFacesIDs;
        std::map<int,int> sharedFace_Owned;
        
        std::map<int,std::vector<int> > Color_SharedOwned;
        std::map<int,int> sharedFOld2New;

        int vt,ft;
        int lvt_o = 0;
        int lft_o = 0;
        int lvt   = 0;
        int lft   = 0;
        std::map<int,std::vector<double> > SharedCoords;
        
        std::set<std::set<int> > PMMG_SharedFaces;
        std::set<std::set<int> > PMMG_SharedFacesOwned;
        int nPartFace = 0;
        std::map<int,int> rank2icomm;
        
        std::map<int,std::vector<int> > LocateOppositeFace;
        std::map<int,int> vertonrank;
        int locID_NotShVrt = 0;
        for( icomm = 0; icomm < next_face_comm; icomm++ )
        {
            nPartFace 				          = nPartFace + nitem_face_comm[icomm];
            rank2icomm[color_face_out[icomm]] = icomm;
                    
            if(world_rank < color_face_out[icomm])
            {
                nTshared_owned = nTshared_owned + nitem_face_comm[icomm];
                
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft       = out_tria_loc[icomm][i];
                    int face_ref = ref2[ft-1];
                    Color_SharedOwned[color_face_out[icomm]].push_back(i);
                    
                    std::set<int> faceSh;
                    for(int k=0;k<3;k++)
                    {
                        int vt2 = triaNodes2[3*(ft-1)+k];
                        faceSh.insert(vt2);
                    }
                    
                    PMMG_SharedFaces.insert(faceSh);
                    PMMG_SharedFacesOwned.insert(faceSh);
                    
                    for(int k=0;k<3;k++)
                    {
                        vt = triaNodes2[3*(ft-1)+k];
                    
                        if(sharedVrts_Owned.find(vt)==sharedVrts_Owned.end())
                        {
                            sharedVrts_Owned[vt] = lvt_o;
                            sharedVert[vt]       = lvt;
                            
                            SharedCoords[vt].push_back(vertOUT[(vt-1)*3+0]);
                            SharedCoords[vt].push_back(vertOUT[(vt-1)*3+1]);
                            SharedCoords[vt].push_back(vertOUT[(vt-1)*3+2]);
                        
                            lvt++;
                            lvt_o++;
                        }
                    }
                    
                    faceSh.clear();
                    lft_o++;
                    lft++;

                }
            }
            else
            {
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];
                    
                    std::set<int> faceSh;
                    
                    if(sharedFacesIDs.find(ft)==sharedFacesIDs.end())
                    {
                        for(int k=0;k<3;k++)
                        {
                            int vt = triaNodes2[3*(ft-1)+k];
                            
                            faceSh.insert(vt);
                            if(sharedVert_NotOwned.find(vt)==sharedVert_NotOwned.end())
                            {
                                sharedVert_NotOwned.insert(vt);
                                LocateOppositeFace[vt].push_back(color_face_out[icomm]);
                                LocateOppositeFace[vt].push_back(i);
                            }
                            
                            if(sharedVert.find(vt)==sharedVert.end())
                            {
                                sharedVert[vt]      = lvt;
                                
                                SharedCoords[vt].push_back(vertOUT[(vt-1)*3+0]);
                                SharedCoords[vt].push_back(vertOUT[(vt-1)*3+1]);
                                SharedCoords[vt].push_back(vertOUT[(vt-1)*3+2]);
                                
                                lvt++;
                            }
                        }
                        
                        sharedFacesIDs[ft] = lft;
                        PMMG_SharedFaces.insert(faceSh);

                        lft++;
                        
                    }
                    faceSh.clear();
                }
            }
        }
        
        // Check for duplicates:
        int nu = 0;
        std::vector<int> dupl_sharedVrts;
        std::map<int,int>::iterator itsh;

        std::set<int> dupl_sharedVrts_set;
        for(itsh=sharedVrts_Owned.begin();itsh!=sharedVrts_Owned.end();itsh++)
        {
            int ggvidd = itsh->first;
            if(dupl_sharedVrts_set.find(ggvidd)==dupl_sharedVrts_set.end() &&
               sharedVert_NotOwned.find(ggvidd)!=sharedVert_NotOwned.end())
            {
                dupl_sharedVrts_set.insert(ggvidd);
                dupl_sharedVrts.push_back(ggvidd);
            }
        }
        
        //============================================================================
        //============================================================================
        //============================================================================
        //============================================================================
        //============================================================================
        
        std::map<int,std::vector<int> > req_locFaceID;
        std::map<int,std::vector<double> > req_locVrtOfCrds;
        
        
        for(int u=0;u<dupl_sharedVrts.size();u++)
        {
            if(sharedVert_NotOwned.find(dupl_sharedVrts[u])!=sharedVert_NotOwned.end())
            {
                int alsoOnRank = LocateOppositeFace[dupl_sharedVrts[u]][0];
                int alsoFid    = LocateOppositeFace[dupl_sharedVrts[u]][1];
                                
                req_locFaceID[alsoOnRank].push_back(alsoFid);
                
                req_locVrtOfCrds[alsoOnRank].push_back(SharedCoords[dupl_sharedVrts[u]][0]);
                req_locVrtOfCrds[alsoOnRank].push_back(SharedCoords[dupl_sharedVrts[u]][1]);
                req_locVrtOfCrds[alsoOnRank].push_back(SharedCoords[dupl_sharedVrts[u]][2]);
                
            }
        }

        ScheduleObj* global_vrt_schedule = DoScheduling(req_locFaceID,comm);

        std::map<int,std::vector<int> > recv_Fids;
        std::map<int,std::vector<double> > recv_Crds;
        std::map<int,std::vector<int> >::iterator itsched;

        for(int q=0;q<world_size;q++)
        {
            if(world_rank==q)
            {
                int i=0;
                for (itsched = req_locFaceID.begin(); itsched != req_locFaceID.end(); itsched++)
                {
                    int n_req           = itsched->second.size();
                    int dest            = itsched->first;
                    
                    MPI_Send(&n_req, 1, MPI_INT, dest, 16798+78*dest, comm);
                    MPI_Send(&itsched->second[0], n_req, MPI_INT, dest, 114876+dest, comm);
                    MPI_Send(&req_locVrtOfCrds[itsched->first][0], n_req*3, MPI_DOUBLE, dest, 214876*5+dest, comm);

                    i++;
                }
            }
            else if (global_vrt_schedule->SendFromRank2Rank[q].find( world_rank ) != global_vrt_schedule->SendFromRank2Rank[q].end())
            {
                int n_reqstd_ids;
                MPI_Recv(&n_reqstd_ids,
                         1, MPI_INT, q,
                         16798+78*world_rank, comm, MPI_STATUS_IGNORE);

                std::vector<int> recv_reqstd_Fids(n_reqstd_ids);
                std::vector<double> recv_reqstd_Crds(n_reqstd_ids*3);
                
                MPI_Recv(&recv_reqstd_Fids[0],
                         n_reqstd_ids, MPI_INT, q,
                         114876+world_rank, comm, MPI_STATUS_IGNORE);
                
                MPI_Recv(&recv_reqstd_Crds[0],
                         n_reqstd_ids*3, MPI_DOUBLE, q,
                         214876*5+world_rank, comm, MPI_STATUS_IGNORE);

                recv_Fids[q] = recv_reqstd_Fids;
                recv_Crds[q] = recv_reqstd_Crds;
            }
        }


        std::map<int,std::vector<int> >::iterator itje;
        double eps = 1.0e-08;
        int mfound;
        int vtfound;
        
        std::set<int> dupliVrtsFound;
        int dpl = 0;
        
        for(itje = recv_Fids.begin();itje!=recv_Fids.end();itje++)
        {
            int rank_recvd = itje->first;
            int icomm_cor = rank2icomm[itje->first];
            
            for(int i=0;i<itje->second.size();i++)
            {
                int fid_cons = itje->second[i];
                int ft       = out_tria_loc[icomm_cor][fid_cons];
                
                double xtest = recv_Crds[itje->first][i*3+0];
                double ytest = recv_Crds[itje->first][i*3+1];
                double ztest = recv_Crds[itje->first][i*3+2];
                
                for(int m=0;m<3;m++)
                {
                    int vt = triaNodes2[3*(ft-1)+m];
                    double resx   = fabs(xtest-SharedCoords[vt][0]);
                    double resy   = fabs(ytest-SharedCoords[vt][1]);
                    double resz   = fabs(ztest-SharedCoords[vt][2]);
                    double restot = resx+resy+resz;
                    if(restot<eps)
                    {
                        mfound  = m;
                        vtfound = vt;
                        dupliVrtsFound.insert(vt);
                        sharedVrts_Owned.erase(vtfound);
                        dpl++;
                    }
                }
            }
        }

        //========================================================================
        //========================================================================
        //========================================================================
        //========================================================================
        //========================================================================
        
        int nSharedVrts          = sharedVert.size();
        int nSharedVrts_own      = sharedVrts_Owned.size();
        int nSharedFaces         = sharedFacesIDs.size();
        int nSharedFaces_own     = sharedFace_Owned.size();
        int nShVrtsCrdsOwn       = SharedCoords.size();
        DistributedParallelState* SharedVrtIds  = new DistributedParallelState(nSharedVrts,comm);
        DistributedParallelState* SharedVrtCrdsOwned = new DistributedParallelState(nShVrtsCrdsOwn,comm);
        DistributedParallelState* SharedVrtIdsOwned  = new DistributedParallelState(nSharedVrts_own,comm);
        int* uvoffset = SharedVrtCrdsOwned->getOffsets();
        std::map<int,int> local2glob_sharedVert;
        std::map<int,int>::iterator unvit;
        int voff = uvoffset[world_rank];
        int guvid;
        for(unvit=sharedVrts_Owned.begin();unvit!=sharedVrts_Owned.end();unvit++)
        {
            int gvid = unvit->first;
            //sharedVrts_Owned[lvid] = guvid;
            
            guvid=voff+1;
        }
        
        int nTotSharedVrts_unique = SharedVrtIdsOwned->getNel();
        int nTotSharedVrts = SharedVrtIds->getNel();

        std::cout << "nTotSharedVrts_unique " << " " << world_rank << " " << nTotSharedVrts_unique << " " << nSharedVrts_own << " " << nTotSharedVrts << " " << nSharedVrts << std::endl;
        
        /*
        //===========================================================================
        //===========================================================================
        //===========================================================================
        
//        std::string filename4 = "AdaptedFaces_" + std::to_string(world_rank) + ".dat";
//        std::ofstream myfile4;
//        myfile4.open(filename4);
//        myfile4 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//        myfile4 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//        myfile4 <<"ZONE N = " << nVerticesOUT << ", E = " << nTrianglesOUT << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
//
//        for(int i=0;i<nVerticesOUT;i++)
//        {
//            pos = 3*i;
//            myfile4 << vertOUT[pos  ] << " " << vertOUT[pos+1  ] << " " << vertOUT[pos+2  ] << std::endl;
//        }
//
//        for ( k=0; k<nTrianglesOUT; k++ )
//        {
//            pos = 3*k;
//
//            myfile4 << triaNodes2[pos] << " "
//                    << triaNodes2[pos+1] << " "
//                    << triaNodes2[pos+2] << std::endl;
//
//        }
//
//        myfile4.close();
        
        //Extract the Local Face and vertex definition for the elements on each partition.
        //===========================================================================
        //===========================================================================
        //===========================================================================
        
        int* tetraOUT = new int[nTetrahedraOUT*4];
        std::vector<std::vector<int> > tetrasOUT;
        std::set<int> NonSharedVrts;
        std::vector<int> NonSharedVrts_vec;

        std::map<int,std::vector<double> > coords_int;

        
        
        std::map<int,int> lh;
        std::map<int,int> rh;
        
        DistributedParallelState* distTetraOut = new DistributedParallelState(nTetrahedraOUT,comm);
        int* TetraOUT_offsets = distTetraOut->getOffsets();
        int curElID = TetraOUT_offsets[world_rank]+1;
        
        std::map<int,std::vector<int> > ienOUT;
        std::map<int,std::vector<int> > fmBnd;
        std::map<int,std::vector<int> > fmInt;
        std::map<int,std::vector<int> > fmSh;
        std::map<int,std::vector<int> > fm;

        
        std::map<std::set<int>, int> facemap;
        std::map<int,std::set<int> > facemap_inv;
        std::map<int,std::vector<int> > bcmap;
        
        
        Array<int>* parmmg_iet = new Array<int>(nTetrahedraOUT,1);
        int q;
        
        std::map<int,int> gvid2shid;
        int fid  = 0;
        int lshf = 0;
        int tetra_faces[4][3] = {{1,2,3},{0,2,3},{0,3,1},{0,2,1}};
        int NoNShFaces = 0;
        
        std::map<int,int> locShF2globShF;
        std::map<int,int> locBF2globBF;
        
        for ( k=0; k<nTetrahedraOUT; k++ )
        {
            if ( PMMG_Get_tetrahedron(parmesh,
                                        &(tetraOUT[pos  ]),&(tetraOUT[pos+1]),
                                        &(tetraOUT[pos+2]),&(tetraOUT[pos+3]),
                                        &(ref[k]),&(required[k])) != 1 )
            {
                fprintf(inm,"Unable to get mesh tetra %d \n",k);
                ier = PMMG_STRONGFAILURE;
            }
            
            std::vector<int> Elvrts(4);
            Elvrts[0] = tetraOUT[pos];
            Elvrts[1] = tetraOUT[pos+1];
            Elvrts[2] = tetraOUT[pos+2];
            Elvrts[3] = tetraOUT[pos+3];
            
            for(int u=0;u<4;u++)
            {
                std::set<int> Face;
                
                for(int w=0;w<3;w++)
                {
                    int lv = tetra_faces[u][w];
                    Face.insert(Elvrts[lv]);
                    
                    int gvidd = Elvrts[lv];
                    
                    //================================================================
                    //================================================================
                    //================================================================
                    // SCENARIO 1!
                    if(NonSharedVrts.find(gvidd)==NonSharedVrts.end() &&
                       sharedVert.find(gvidd)==sharedVert.end())
                    {
                        NonSharedVrts.insert(gvidd);
                        NonSharedVrts_vec.push_back(gvidd);
                    }
                    
                    //================================================================
                    //================================================================
                    //================================================================
                }
                
                if(facemap.find(Face)==facemap.end())
                {
                    facemap[Face]=fid;
                    std::vector<int> fce(3);
                    
                    fce[0]  = Elvrts[tetra_faces[u][0]];
                    fce[1]  = Elvrts[tetra_faces[u][1]];
                    fce[2]  = Elvrts[tetra_faces[u][2]];
                    
                    fm[fid] = fce;
                    
                    if(PMMG_SharedFaces.find(Face) == PMMG_SharedFaces.end()
                       && PMMG_Face2Ref.find(Face) == PMMG_Face2Ref.end())
                    {
                        fmInt[fid]  = fce;
                        lh[fid]     = curElID;
                    }
                    
                    if(PMMG_SharedFacesOwned.find(Face) != PMMG_SharedFacesOwned.end())
                    {
                        fmSh[fid]               = fce;
                        locShF2globShF[lshf]    = fid;
                        lh[fid] = curElID;
                        lshf++;
                    }
                    
                    if(PMMG_Face2Ref.find(Face) != PMMG_Face2Ref.end())
                    {
                        fmBnd[fid]          = fce;
                        int FaceRef         = PMMG_Face2Ref[Face];
                        locBF2globBF[lshf]  = fid;
                        lh[fid]             = curElID;
                        bcmap[FaceRef].push_back(fid);
                    }
                    fid++;
                }
                else
                {
                    int fid_n     = facemap[Face];
                    rh[fid_n]     = curElID;
                }
                
                Face.clear();
            }
            
            double* P = new double[4*3];
            
            P[0*3+0]=vertOUT[(v0-1)*3];
            P[0*3+1]=vertOUT[(v0-1)*3+1];
            P[0*3+2]=vertOUT[(v0-1)*3+2];
            
            P[1*3+0]=vertOUT[(v1-1)*3];
            P[1*3+1]=vertOUT[(v1-1)*3+1];
            P[1*3+2]=vertOUT[(v1-1)*3+2];
            
            P[2*3+0]=vertOUT[(v2-1)*3];
            P[2*3+1]=vertOUT[(v2-1)*3+1];
            P[2*3+2]=vertOUT[(v2-1)*3+2];
            
            P[3*3+0]=vertOUT[(v3-1)*3];
            P[3*3+1]=vertOUT[(v3-1)*3+1];
            P[3*3+2]=vertOUT[(v3-1)*3+2];
            
            std::vector<int> tetra(4);
            tetra[0] = v0;
            tetra[1] = v1;
            tetra[2] = v2;
            tetra[3] = v3;
            
            ienOUT[curElID] = tetra;
            tetrasOUT.push_back(tetra);
            double Vtet = GetQualityTetrahedra(P);
            
            if(Vtet<0.0)
            {
                std::cout << "Error " << Vtet << std::endl;
            }
            
            curElID++;

            delete[] P;
        }
        
        
        
        // GENERATE THE GLOBAL VERTEX IDS AND THE CORRESPONDING MAPS

        //==================================================================
        //==================================================================
        //==================================================================
        //==================================================================
        //==================================================================
        int totalLocalVrts = NonSharedVrts_vec.size()+sharedVrts_Owned.size();
        DistributedParallelState* distUniqueTotVerts     = new DistributedParallelState(totalLocalVrts,comm);
        DistributedParallelState* distUniqueVerts        = new DistributedParallelState(NonSharedVrts_vec.size(),comm);
        DistributedParallelState* distUniqueSharedVerts  = new DistributedParallelState(sharedVrts_Owned.size(),comm);

        int nTotVrtsUnique      = distUniqueTotVerts->getNel();
        int nTotVrts_unique     = distUniqueVerts->getNel();
        
    	std::map<int,int> locVrt2globVrt;
    	std::map<int,int> globVrt2locVrt;

        int ToTVrts_offset      = distUniqueTotVerts->getOffsets()[world_rank];
        int ToTVrts             = distUniqueTotVerts->getNel();
        int gv                  = ToTVrts_offset;
        int gv_b                = ToTVrts_offset;
        int vrtid               = 0;

        Array<double>* xcn_parmmg = new Array<double>(totalLocalVrts,3);
        
    	for(int i=0;i<NonSharedVrts_vec.size();i++)
    	{
            int lcid = NonSharedVrts_vec[i];
            
    		locVrt2globVrt[NonSharedVrts_vec[i]] = gv;
    		globVrt2locVrt[gv] = NonSharedVrts_vec[i];
            
            xcn_parmmg->setVal(vrtid,0,vertOUT[(lcid-1)*3+0]);
            xcn_parmmg->setVal(vrtid,1,vertOUT[(lcid-1)*3+1]);
            xcn_parmmg->setVal(vrtid,2,vertOUT[(lcid-1)*3+2]);
            
            vrtid++;
    		gv++;
    	}
        
        std::map<int,int>::iterator sh_vrt;
        for(sh_vrt=sharedVrts_Owned.begin();sh_vrt!=sharedVrts_Owned.end();sh_vrt++)
        {
            int lcid = sh_vrt->first;
            
            locVrt2globVrt[lcid] = gv;
            globVrt2locVrt[gv] = lcid;
            
            xcn_parmmg->setVal(vrtid,0,vertOUT[(lcid-1)*3+0]);
            xcn_parmmg->setVal(vrtid,1,vertOUT[(lcid-1)*3+1]);
            xcn_parmmg->setVal(vrtid,2,vertOUT[(lcid-1)*3+2]);
            
            vrtid++;
            gv++;
        }
        
        std::cout << "nTotVrtsUnique " << nTotVrtsUnique << " "  << totalLocalVrts << std::endl;
        
        int nTotIntpShFaces = fmInt.size() + fmSh.size();
        DistributedParallelState* distUniqueIntpShFaces   = new DistributedParallelState(nTotIntpShFaces,comm);
        DistributedParallelState* distUniqueInteriorFace  = new DistributedParallelState(fmInt.size(),comm);
        DistributedParallelState* distUniqueSharedFace    = new DistributedParallelState(fmSh.size(),comm);
        DistributedParallelState* distUniqueBndFace       = new DistributedParallelState(fmBnd.size(),comm);
        int nTotBndFaces_unique                           = distUniqueBndFace->getNel();
        int nTotSharedFaces_unique                        = distUniqueSharedFace->getNel();
        int nTotIntpShFaces_unique                        = distUniqueIntpShFaces->getNel();
        int nTotIntFaces_unique                           = distUniqueInteriorFace->getNel();

        int nLocIntpShFace_offset = distUniqueIntpShFaces->getOffsets()[world_rank];

        
        
        // Ordering is InteriorFaces, SharedFaces, BoundaryFaces;

        //=================================================================================
        //=================================================================================
        //=================================================================================
        //=================================================================================
        //=================================================================================
        //=================================================================================
        //=================================================================================
        //=================================================================================
        
        //Add send/receive of the adjacent elements of the shared faces;
        
        
        //Color_SharedOwned
        
        ScheduleObj* ish_schedule = DoScheduling(Color_SharedOwned,comm);
  
        std::map<int,std::vector<int> > recv_ids;
        
        std::map<int,std::vector<int> >::iterator it;
          
        for(int q=0;q<world_size;q++)
		{
			if(world_rank==q)
			{
				int i=0;
				for (it = Color_SharedOwned.begin(); it != Color_SharedOwned.end(); it++)
				{
					int n_req           = it->second.size();
					int dest            = it->first;

					MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
					MPI_Send(&it->second[0], n_req, MPI_INT, dest, 14876+dest, comm);

					i++;
				}
			}
			else if (ish_schedule->SendFromRank2Rank[q].find( world_rank ) != ish_schedule->SendFromRank2Rank[q].end())
			{
				int n_reqstd_ids;
				MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);

				std::vector<int> recv_reqstd_ids(n_reqstd_ids);
				
				MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
				recv_ids[q] = recv_reqstd_ids;
			}
		}
        
        
		std::map<int,std::vector<int> > sendEl;
        
        std::map<int,std::vector<int> >::iterator rcvit;
        for(rcvit=recv_ids.begin();rcvit!=recv_ids.end();rcvit++)
        {
        	int frank = rcvit->first;
        	
        	int icomm = rank2icomm[frank];
        	int nF    = rcvit->second.size();  
        	
			for(int j=0;j<nF;j++)
			{
				int fidInt = rcvit->second[j];
				int ofid   = out_tria_loc[icomm][fidInt];
				int nfid   = locShF2globShF[ofid];
				
				sendEl[frank].push_back(lh[nfid]);
			}
        }
        
        ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);


        std::map<int,std::vector<int> > adj_ids;
        for(int q=0;q<world_size;q++)
		{
			if(world_rank==q)
			{
				int i=0;
				for (it = sendEl.begin(); it != sendEl.end(); it++)
				{
					int n_req           = it->second.size();
					int dest            = it->first;

					MPI_Send(&n_req, 1, 
							MPI_INT, dest, 
							6798+78000*dest, comm);
					
					MPI_Send(&it->second[0], 
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
        
        
        //int nBndFaces	 = nTrianglesOUT-sharedFaces.size();
        int nShaFaces	 = sharedFace_Owned.size();
        int nIntFaces	 = rh.size();
		int nLocFaceNBnd = rh.size()+nShaFaces;
        Array<int>* ifnOUT = new Array<int>(nLocFaceNBnd,8);

        std::map<int,int>::iterator lhit;
        
        int ftot = 0;
		std::map<int,int> updateID;
        std::map<int,int>::iterator itm;
		int gv0,gv1,gv2;
        for(itm=rh.begin();itm!=rh.end();itm++)
        {
        	fid = itm->first;
        	updateID[fid] = ftot;
			ifnOUT->setVal(ftot,0,3);
            if(fm.find(fid)!=fm.end())
            {
                if(locVrt2globVrt.find(fm[fid][0])!=locVrt2globVrt.end()
                   && locVrt2globVrt.find(fm[fid][1])!=locVrt2globVrt.end()
                   && locVrt2globVrt.find(fm[fid][2])!=locVrt2globVrt.end())
                {
                    gv0 = locVrt2globVrt[fm[fid][0]]+1;
                    gv1 = locVrt2globVrt[fm[fid][1]]+1;
                    gv2 = locVrt2globVrt[fm[fid][2]]+1;
                }
            }
            
			ifnOUT->setVal(ftot,1,gv0);
			ifnOUT->setVal(ftot,2,gv1);
			ifnOUT->setVal(ftot,3,gv2);
			ifnOUT->setVal(ftot,4,0);
			ifnOUT->setVal(ftot,5,rh[fid]);
			ifnOUT->setVal(ftot,6,lh[fid]);
			ifnOUT->setVal(ftot,7,2);
            
            
        	ftot++;
        }
        
        std::cout << ftot << std::endl;
        
        
        
		int cntrr=0;
		//		
		int NotFound = 0;
		for(cit=Color_SharedOwned.begin();cit!=Color_SharedOwned.end();cit++)
		{
			int icomm   =  cit->first;
			int srank 	= color_face_out[icomm];
			int nF 		= cit->second.size();
			
			for(int j=0;j<nF;j++)
			{
				int fidInt   = cit->second[j];
				int oldfid   = out_tria_loc[icomm][fidInt];
				int newfid   = sharedFOld2New[oldfid];
				
				reqEl[srank].push_back(fidInt);
				
				int elLh = lh[newfid];
				int elRh = adj_ids[srank][j];
				
				if(rh.find(newfid)==rh.end())
				{
					rh[newfid] = elRh;
					updateID[newfid] = ftot;
					ifnOUT->setVal(ftot,0,3);
					
					gv0 = locVrt2globVrt[fm[newfid][0]]+1;
					gv1 = locVrt2globVrt[fm[newfid][1]]+1;
					gv2 = locVrt2globVrt[fm[newfid][2]]+1;
                    
					ifnOUT->setVal(ftot,1,gv0);
					ifnOUT->setVal(ftot,2,gv1);
					ifnOUT->setVal(ftot,3,gv2);
                    
					ifnOUT->setVal(ftot,4,0);
					ifnOUT->setVal(ftot,5,rh[newfid]);
					ifnOUT->setVal(ftot,6,lh[newfid]);
					ifnOUT->setVal(ftot,7,2);
					ftot++;
				}
				cntrr++;
			}
		}
		
        
        
        
        
        
        
        
        
        
        
        
        
        
        DistributedParallelState* distInteriorFaces = new DistributedParallelState(nLocFaceNBnd,comm);

		int nTotInteriorFaces    = distInteriorFaces->getNel();
        int* TotIntFaces_offsets = distInteriorFaces->getOffsets();

        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        

        
        int nLocFacesInt = rh.size();
        int nLocFacesTot = lh.size();
        
        DistributedParallelState* distIntFace  	= new DistributedParallelState(nLocFacesInt,comm);
        int Nt_IntFaces              			= distIntFace->getNel();
        int* IntFaces_offsets        			= distIntFace->getOffsets();
        int* IntFaces_nlocs          			= distIntFace->getNlocs();
        int* IntFacesIDs             			= new int[nLocFacesInt];
        int* IntFaces_RankIDs        			= new int[nLocFacesInt];
        
        DistributedParallelState* distTotFace  	= new DistributedParallelState(nLocFacesTot,comm);
        int Nt_TotFaces              			= distTotFace->getNel();
        int* TotFaces_offsets        			= distTotFace->getOffsets();
        int* TotFaces_nlocs          			= distTotFace->getNlocs();
        int* TotFacesIDs             			= new int[nLocFacesTot];
        int* TotFaces_RankIDs       	 		= new int[nLocFacesTot];
        
        std::map<int,int> rh_new;
        int iface = 0;
        
        
        std::map<int,std::vector<int> >::iterator bit;
        std::set<int> sorted_BCid;
        std::map<int,int> sorted_BCid_map;
        std::map<int,int> sorted_NBCid_map;
        int ii = 0;
        int Nbf;
        int nloc_bcs = bcmap.size();
        int* Lbcs = new int[nloc_bcs];
        for(bit=bcmap.begin();bit!=bcmap.end();bit++)
        {
            //Nbf         = bit->second.size();
            Lbcs[ii]      = bit->first;
            ii++;
        }
        
        
        DistributedParallelState* distLocBCs = new DistributedParallelState(nloc_bcs,comm);
        
        int Nt_BCs                 = distLocBCs->getNel();
        int* BCs_offsets           = distLocBCs->getOffsets();
        int* BCs_nlocs             = distLocBCs->getNlocs();
        int* BCs_arr               = new int[Nt_BCs];
        
        MPI_Allgatherv(&Lbcs[0],
                       nloc_bcs,
                       MPI_INT,
                       BCs_arr,
                       BCs_nlocs,
                       BCs_offsets,
                       MPI_INT, comm);

        std::set<int> bcsToT;
        std::map<int,int> bcentry;
        for(int i=0;i<Nt_BCs;i++)
        {
            
            if(bcsToT.find(BCs_arr[i])==bcsToT.end())
            {
                
                bcsToT.insert(BCs_arr[i]);
//                bcentry[bcsToT] = entry;
//                entry++;
            }
        }
        
        int* bcid = new int[bcsToT.size()];
        int* nlbc = new int[bcsToT.size()];
        std::set<int>::iterator entry;
        int cnt = 0;
        
        int lac;
        q = 0;
        for(entry=bcsToT.begin();entry!=bcsToT.end();entry++)
        {
            int ee = *entry;
            
            if(bcmap.find(ee)!=bcmap.end())
            {
                Nbf      = bcmap[ee].size();
                bcid[q]  = ee;
                nlbc[q]  = Nbf;
            }
            else
            {
                bcid[q]  = ee;
                nlbc[q]  = 0;
            }
            q++;
        }
        
        std::vector<Array<int>* > bcArrays;
        std::map<int,int> bcsizing;
        std::vector<int> bci_offsets;
        std::vector<int> bciTot_offsets;
        int nTotBCFaces = 0;
        int nTotBCFaces_offset = 0;
        for(int i=0;i<bcsToT.size();i++)
        {
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);
            
            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();

            Array<int>* ifn_bc_i = new Array<int>(NelLoc_bci,8);

            int offsetbci = distBCi->getOffsets()[world_rank];

            if(bcmap.find(bc_id)!=bcmap.end())
            {
                Nbf = bcmap[bc_id].size();
                
                int fbc = 0;

                for(int q=0;q<Nbf;q++)
                {
                    int bcface = bcmap[bc_id][q];
                    
                    gv0 = locVrt2globVrt[fm[bcface][0]]+1;
                    gv1 = locVrt2globVrt[fm[bcface][1]]+1;
                    gv2 = locVrt2globVrt[fm[bcface][2]]+1;
                    
                    ifn_bc_i->setVal(fbc,0,3);
                    ifn_bc_i->setVal(fbc,1,gv0);
                    ifn_bc_i->setVal(fbc,2,gv1);
                    ifn_bc_i->setVal(fbc,3,gv2);
                    
                    ifn_bc_i->setVal(fbc,4,0);
                    ifn_bc_i->setVal(fbc,5,0);
                    ifn_bc_i->setVal(fbc,6,lh[bcface]);
                    
                    ifn_bc_i->setVal(fbc,7,bc_id);
                    
                    fbc++;
                }
            }
            
            bcsizing[bc_id] = NelTot_bci;
            bci_offsets.push_back(offsetbci);
            bciTot_offsets.push_back(nTotBCFaces_offset);
            bcArrays.push_back(ifn_bc_i);
            
            nTotBCFaces_offset = nTotBCFaces_offset + NelTot_bci;
            nTotBCFaces        = nTotBCFaces + NelTot_bci;
        }
        
        DistributedParallelState* distTetra = new DistributedParallelState(nTetrahedraOUT,comm);

        int ToTElements         = distTetra->getNel();
        int ToTElements_offset  = distTetra->getOffsets()[world_rank];
        
        int nbo = bcArrays.size();
        std::cout << "-- Constructing the zdefs array..."<<std::endl;
        Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
        // Collect node data (10) . Starting index-ending index Nodes
        adapt_zdefs->setVal(0,0,10);
        adapt_zdefs->setVal(0,1,-1);
        adapt_zdefs->setVal(0,2,1);
        adapt_zdefs->setVal(0,3,1);
        adapt_zdefs->setVal(0,4,ToTVrts);
        adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
        adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
        // Collect element data (12) . Starting index-ending index Element
        adapt_zdefs->setVal(1,0,12);
        adapt_zdefs->setVal(1,1,-1);
        adapt_zdefs->setVal(1,2,2);
        adapt_zdefs->setVal(1,3,1);
        adapt_zdefs->setVal(1,4,ToTElements);
        adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
        adapt_zdefs->setVal(1,6,2);
        // Collect internal face data (13) . Starting index-ending index internal face.
        adapt_zdefs->setVal(2,0,13);
        adapt_zdefs->setVal(2,1,-1);
        adapt_zdefs->setVal(2,2, 3);
        adapt_zdefs->setVal(2,3, 1);
        adapt_zdefs->setVal(2,4,nTotInteriorFaces);
        adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
        adapt_zdefs->setVal(2,6,2);
        
        int qq  = 1;
        int nb = 0;
        int face_start = nTotInteriorFaces+1;
        int face_end;
        std::map<int,int>::iterator itr;
        for(itr=bcsizing.begin();itr!=bcsizing.end();itr++)
        {
            int bnd_ref  = itr->first;
            int bnd_size = itr->second;
            face_end = face_start+bnd_size-1;
            adapt_zdefs->setVal(3+nb,0,13);
            adapt_zdefs->setVal(3+nb,1,-1);
            adapt_zdefs->setVal(3+nb,2,3+qq);
            adapt_zdefs->setVal(3+nb,3,face_start);
            adapt_zdefs->setVal(3+nb,4,face_end);
            adapt_zdefs->setVal(3+nb,5,bnd_ref);
            adapt_zdefs->setVal(3+nb,6,2);
            //std::cout << "us3d->zdefs->getVal(3+nb,5) " << us3d->zdefs->getVal(3+nb,5) << std::endl;
            face_start = face_end+1;
            
            if(world_rank == 0)
            {
                std::cout << "bc id  " << bnd_ref << " " << bnd_size << std::endl;
            }
            //std::cout << "nb  = " << nb << " " << ref2bface.size() << " " << ref2bqface.size() << std::endl;
            nb++;
            qq++;
        }
        
        int nTotFaces = nTotInteriorFaces + nTotBCFaces;
        
        if(world_rank == 0)
        {
            
            std::cout << "Total faces = " << nTotFaces << std::endl;
            
            PlotBoundaryData(us3d->znames,adapt_zdefs);
        }

        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        
        hid_t ret;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        hid_t file_id = H5Fcreate("par_grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        hid_t    dset_id;
        hid_t    filespace;
        hid_t    memspace;
        hid_t    status;
        hsize_t     dimsf[2];
        hsize_t     countH5[2];
        hsize_t     offsetH5[2];
        
        hsize_t dimsf_att = 1;
        hid_t att_space = H5Screate_simple(1, &dimsf_att, NULL);
        hid_t type =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type, 14);
        ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
        hid_t attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri[] = "US3D Grid File";
        status = H5Awrite(attr_id, type, &stri);
        H5Aclose(attr_id);
        
        hsize_t dimsf_att2 = 1;
         att_space = H5Screate_simple(1, &dimsf_att2, NULL);
        hid_t type2 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type2, 5);
        ret = H5Tset_strpad(type2,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (file_id, "filevers", type2, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri2[] = "1.1.8";
        status = H5Awrite(attr_id, type2, &stri2);
        H5Aclose(attr_id);
        
        hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

        hsize_t dimsf_att3 = 1;
        att_space = H5Screate_simple(1, &dimsf_att3, NULL);
        hid_t type3 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type3, 10);
        ret = H5Tset_strpad(type3,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (group_info_id, "date", type3, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri3[] = "27-05-1987";
        status = H5Awrite(attr_id, type3, &stri3);
        H5Aclose(attr_id);
        
        hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        int value = ToTElements;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotBCFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = ToTVrts;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        
        //====================================================================================
        // Add iet map to the grid.h5 file
        //====================================================================================
        dimsf[0] = ToTElements;
        dimsf[1] = parmmg_iet->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);

        dset_id = H5Dcreate(file_id, "iet",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = parmmg_iet->getNrow();
        countH5[1]  = parmmg_iet->getNcol();
        
        offsetH5[0] = ToTElements_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet->data);
        delete parmmg_iet;
        //====================================================================================
        // Add xcn map to the grid.h5 file
        //====================================================================================
        dimsf[0] = ToTVrts;
        dimsf[1] = xcn_parmmg->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);

        dset_id = H5Dcreate(file_id, "xcn",
                            H5T_NATIVE_DOUBLE, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = xcn_parmmg->getNrow();
        countH5[1]  = xcn_parmmg->getNcol();
        
        offsetH5[0] = ToTVrts_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_parmmg->data);
        delete xcn_parmmg;
        //===================================================================================
        
        
        dimsf[0]  = nTotInteriorFaces+nTotBCFaces;
        dimsf[1]  = ifnOUT->getNcol();
        
        filespace = H5Screate_simple(2, dimsf, NULL);
        dset_id = H5Dcreate(file_id, "ifn",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//
        countH5[0]  = ifnOUT->getNrow();
        countH5[1]  = dimsf[1];
        
        offsetH5[0] = TotIntFaces_offsets[world_rank ];
        offsetH5[1] = 0;
        
        memspace     = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, ifnOUT->data);
        
        
        for(int i=0;i<bcsToT.size();i++)
        {
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);

            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();

            Array<int>* ifn_bc_i = bcArrays[i];

            countH5[0]  = ifn_bc_i->getNrow();
            countH5[1]  = dimsf[1];

            offsetH5[0] = nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i];
            offsetH5[1] = 0;
            memspace     = H5Screate_simple(2, countH5, NULL);
            filespace     = H5Dget_space(dset_id);
            
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);

            plist_id     = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, ifn_bc_i->data);

        }
        
        
        
        
        // Create group;
        //====================================================================================
        hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        // Add attribute to group:
        //====================================================================================
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = 3+nbo;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //====================================================================================
        dimsf[0] = adapt_zdefs->getNrow();
        dimsf[1] = adapt_zdefs->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);
        hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);
        
        countH5[0]  = dimsf[0];
        countH5[1]  = dimsf[1];
        offsetH5[0] = 0;
        offsetH5[1] = 0;
        memspace  = H5Screate_simple(2, countH5, NULL);
        filespace = H5Dget_space(dset_zdefs_id);
        
        //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
        //====================================================================================
        
        dimsf_att = us3d->znames->getNrow();
        filespace = H5Screate_simple(1, &dimsf_att, NULL);
        type =  H5Tcopy (H5T_C_S1);
        ret  = H5Tset_size (type, 20);
        ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
        hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sclose(filespace);

        hsize_t cnt2 = us3d->znames->getNrow();
        memspace  = H5Screate_simple(1, &cnt2, NULL);
        filespace = H5Dget_space(dset_znames_id);


        status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);
    

        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        
                
        if(tetrasOUT.size()!=0)
        {
            std::string filename2 = "AdaptedMin_" + std::to_string(world_rank) + ".dat";
            std::ofstream myfile2;
            myfile2.open(filename2);
            myfile2 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
            myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
            myfile2 <<"ZONE N = " << nVerticesOUT << ", E = " << tetrasOUT.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

            for(int i=0;i<nVerticesOUT;i++)
            {
                pos = 3*i;
                myfile2 << vertOUT[pos  ] << " " << vertOUT[pos+1  ] << " " << vertOUT[pos+2  ] << std::endl;
            }

            
            for(int i=0;i<tetrasOUT.size();i++)
            {
                pos = 4*i;
                
                myfile2 <<   tetrasOUT[i][0] << " "
                        <<   tetrasOUT[i][1] << " "
                        <<   tetrasOUT[i][2] << " "
                        <<   tetrasOUT[i][3] << std::endl;
            }

            myfile2.close();

        }
         
        
    }*/
    
     
    
    
    
    

    

    
    
    MPI_Finalize();
    
}
