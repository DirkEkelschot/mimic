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

struct TetrahedraMesh
{
    Array<int>* ElGids;
    Array<int>* ien_part_tetra;
    Array<int>* ien_part_hybrid;
    
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
};


struct PartitionBoundary
{
    int nPartBoundFaces;
    std::vector<int> faces4parmmg;
    std::map<int,int> globShF2locShF;
    int* ifc_tria_glob;
    int** ifc_tria_loc;
    int* color_face;
    int* ntifc;
    int ncomm;
    std::map<int,std::vector<int> > ColorsFaces;
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
                f++;
            }
            shf++;
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
                        
                        v2ref[vid] = ref;
                        
                        uBoundVert.insert(vid);
                        g2l_b[vid] = Bvid;
                        lBvid = tmesh->globV2locV[vid];
                        Vert* bv = new Vert;
                        bv->x = locVs[lBvid]->x;
                        bv->y = locVs[lBvid]->y;
                        bv->z = locVs[lBvid]->z;
                        
                        bface[u]=Bvid;
                        bvs.push_back(bv);
                        
                        Bvid++;
                    }
                    else
                    {
                        int Bvid2 = g2l_b[vid];
                        bface[u]=Bvid2;
                    }
                }
                
                pltBfaces.push_back(bface);
                pb->globShF2locShF[faceID] = f;
                f++;
            }
            bf++;
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
    
    std::map<int,std::vector<int> >::iterator itefmap;
    std::map<int,std::vector<int> >ife_part_hyb_map;
    for(int i=0;i<tmesh->ief_part_hybrid->getNrow();i++)
    {
        el_id   = ien_pstate->getOffsets()[rank]+i;
        
        for(int q=0;q<4;q++)
        {
            int face_hyb = tmesh->ief_part_hybrid->getVal(i,q);
            int face_tet = tmesh->ief_part_tetra->getVal(i,q);
            
            //face2element[face_tet].push_back(el_id);
            
            if(ufacesHyb.find(face_hyb)==ufacesHyb.end())
            {
                ufacesHyb.insert(face_hyb);
                tmesh->tetF2hybF[face_tet] = face_hyb;
                tmesh->hybF2tetF[face_hyb] = face_tet;
                
                r = FindRank(new_offsets,size,face_hyb);
                if(face_id == 865052)
                {
                    std::cout << " Sends it " << tmesh->ief_part_hybrid->getVal(i,q) << " " << r << std::endl;
                }
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
    int n_recv_back;
    
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send             = it->second.size();
                double* ife_send        = new double[nf_send*ncol];
                int offset_ife          = ife_pstate->getOffset(rank);
                
                for(int u=0;u<it->second.size();u++)
                {
                    int faceID = it->second[u];
                    
                    for(int h=0;h<ncol;h++)
                    {
                        ife_send[u*ncol+h] = ife->getVal(it->second[u]-offset_ife,h);
                    }
                }

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&ife_send[0], nf_send*ncol, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                delete[] ife_send;
            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_ife_arr   = new double[n_recv_back*ncol];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_ife_arr[0], n_recv_back*ncol, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nife[q]       = n_recv_back;
            recv_back_face_ids[q]   = recv_back_ids_arr;
            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }
   
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Nife.begin();iter!=recv_back_Nife.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            face_id = recv_back_face_ids[iter->first][s];
            
            for(int r=0;r<ncol;r++)
            {
                ife_part_hyb_map[face_id].push_back(recv_back_ife[iter->first][s*ncol+r]);
            }
        }
        ntotal=ntotal+L;
    }

    delete[] new_offsets;
    
    
    int fhyb,ftet,hvid,tvid;
    
    std::map<int,std::vector<int> >::iterator itf;
    for(itf=ife_part_hyb_map.begin();itf!=ife_part_hyb_map.end();itf++)
    {
        fhyb       = itf->first;
        ftet       = tmesh->hybF2tetF[itf->first];
        int* nodes = new int[ncol];
        
        for(int j=0;j<ncol;j++)
        {
            hvid         = itf->second[j];
            tvid         = tmesh->hybV2tetV[hvid];
            
            
            nodes[j]     = tvid;
        }
        
        face2node[ftet] = nodes;
    }
    
    
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
    
    Array<int>* ElGids             = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_tetra     = new Array<int>(nElexpctd,4);
    Array<int>* ien_part_hybrid    = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_tetra     = new Array<int>(nElexpctd,4);
    Array<int>* ief_part_hybrid    = new Array<int>(nElexpctd,4);
    Array<int>* iefref_part_tetra  = new Array<int>(nElexpctd,4);
    Array<int>* loc2globEL         = new Array<int>(nElexpctd,1);
    ParallelState* ien_pstate      = new ParallelState(nglob,comm);

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
    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > Gelms_to_send_to_ranks;
    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > OriVertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > OriFaceIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceRefs_to_send_to_ranks;
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_vertOri;
    std::map<int,std::vector<int> > rank2req_face;
    std::map<int,std::vector<int> > rank2req_faceOri;
    
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::vector<int> vertIDs_on_rank_Ori;
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
    int gEL;
    for(i=0;i<part->getNrow();i++)
    {
        p_id    = part->getVal(i,0);
        el_id   = ien_pstate->getOffsets()[rank]+i;
        gEL     = tmesh->ElGids->getVal(i,0);
        
        loc2globEL->setVal(i,0,el_id);
        
        nvPerEl = 4;
        nfPerEl = 4;
        sum =0;
        sum3 = 0;
        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            Gelms_to_send_to_ranks[p_id].push_back(gEL); // rank to element map.
            //====================Hybrid=======================
            for(int k=0;k<4;k++)
            {
                
                v_id   = tmesh->ien_part_tetra->getVal(i,k);
                v_id_o = tmesh->ien_part_hybrid->getVal(i,k);
                
                vertIDs_to_send_to_ranks[p_id].push_back(v_id);
                OriVertIDs_to_send_to_ranks[p_id].push_back(v_id_o);
                
                sum3=sum3+v_id;
                
                
            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            
            for(int k=0;k<4;k++)
            {
                f_id   = tmesh->ief_part_tetra->getVal(i,k);
                f_id_o = tmesh->ief_part_hybrid->getVal(i,k);
                fref   = tmesh->iefref_part_tetra->getVal(i,k);
                
                faceIDs_to_send_to_ranks[p_id].push_back(f_id);
                faceRefs_to_send_to_ranks[p_id].push_back(fref);
                OriFaceIDs_to_send_to_ranks[p_id].push_back(f_id_o);

            }

            
            //====================Hybrid=======================
            not_on_rank++;
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
                sum = sum + v_id;
                
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    tmesh->hybV2tetV[v_id_o]=v_id;
                    tmesh->tetV2hybV[v_id]=v_id_o;
                    
                    
                    r = FindRank(new_V_offsets,size,v_id_o);
                    
                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                        rank2req_vertOri[r].push_back(v_id_o);
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vertIDs_on_rank_Ori.push_back(v_id_o);
                        
                        
                        
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

                if(unique_faceIDs_on_rank_set.find( f_id ) == unique_faceIDs_on_rank_set.end() && f_id != -1) // add the required unique vertex for current rank.
                {
                    unique_faceIDs_on_rank_set.insert(f_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
                    r = FindRank(new_F_offsets,size,f_id_o);

                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_face[r].push_back(f_id); // add the vertex id that needs to be requested from rank r.
                        rank2req_faceOri[r].push_back(f_id_o);

                    }
                    else
                    {
                        faceIDs_on_rank.push_back(f_id);  // add the vertex to list that is already available on rank.
                        
//                      tmesh->hybF2tetF[f_id_o]=f_id;
//                      tmesh->tetF2hybF[f_id]=f_id_o;
                        
                        tmesh->ref2face[fref].push_back(f_id);
                        
                        if(rank==0 && f_id == 160591)
                        {
                            std::cout << "RANK " << rank << " fref = " << f_id << " " << fref << std::endl;
                        }
                        floc_tmp++;
                    }
                    lf_id++;
                }
            }
            
            loc_r_elem.push_back(el_id);
            
            on_rank++;
        }
    }
    
    
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
    
    
    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);
    std::map<int,std::vector<int> >  part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_ov_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_f_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_fref_map;
    std::map<int,std::vector<int> >  TotRecvElement_IDs_of_map;
    std::map<int,std::vector<int> >  TotRecvElement_GIDs_map;

    std::map<int,std::vector<int> >::iterator it;
        
    int n_req_recv;
    int n_req_recv_v;
    int n_req_recv_f;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = elms_to_send_to_ranks.begin(); it != elms_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = vertIDs_to_send_to_ranks[it->first].size();
                
                int n_req_f         = faceIDs_to_send_to_ranks[it->first].size();
                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, comm);
                MPI_Send(&n_req_f, 1, MPI_INT, dest, dest*222, comm);
                
                //std::cout << "vertIDs_to_send_to_ranks[it->first] " << vertIDs_to_send_to_ranks[it->first].size() << " " << vertIDs_to_send_to_ranks[it->first].size() << std::endl;
                
                //MPI_Send(&it->second[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                
//                if(dest==1)
//                {
//                    std::cout << rank << " communicating " << vertIDs_to_send_to_ranks[it->first].size()/4 << " " << OriVertIDs_to_send_to_ranks[it->first].size()/4 << std::endl;
//                    if(rank==3)
//                    {
//                        for(int l=0;l<vertIDs_to_send_to_ranks[it->first].size()/4;l++)
//                        {
//                            std::cout << l << " commu :: ";
//                            for(int h=0;h<4;h++)
//                            {
//                                std::cout << OriVertIDs_to_send_to_ranks[it->first][l*4+h] << " ";
//                            }
//                            std::cout << std::endl;
//                        }
//                    }
//                }
                
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&OriVertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 339000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                MPI_Send(&faceRefs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 449000+100+dest*2, comm);
                MPI_Send(&OriFaceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 889000+100+dest*2, comm);
                MPI_Send(&Gelms_to_send_to_ranks[0], n_req, MPI_INT, dest, dest*7000000, comm);

                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222, comm, MPI_STATUS_IGNORE);

            std::vector<int>    part_recv_el_id(n_req_recv);
            std::vector<int>    part_recv_vrt_id(n_req_recv_v);
            std::vector<int>    part_recv_orivrt_id(n_req_recv_v);
            std::vector<int>    part_recv_face_id(n_req_recv_f);
            std::vector<int>    part_recv_face_ref(n_req_recv_f);
            std::vector<int>    part_recv_orifaces_id(n_req_recv_v);
            std::vector<int>    part_recv_Gel_id(n_req_recv);

            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0], n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orivrt_id[0], n_req_recv_v, MPI_INT, q, 339000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_ref[0], n_req_recv_f, MPI_INT, q, 449000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_orifaces_id[0], n_req_recv_f, MPI_INT, q, 889000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_Gel_id[0], n_req_recv, MPI_INT, q, rank*7000000, comm, MPI_STATUS_IGNORE);
            
            TotRecvElement_IDs_ov_map[q]    = part_recv_orivrt_id;
            TotRecvElement_IDs_v_map[q]     = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q]     = part_recv_face_id;
            part_tot_recv_elIDs_map[q]      = part_recv_el_id;
            TotRecvElement_IDs_fref_map[q]  = part_recv_face_ref;
            TotRecvElement_IDs_of_map[q]    = part_recv_orifaces_id;
            TotRecvElement_GIDs_map[q]      = part_recv_Gel_id;
        }
    }
    std::vector<int> TotRecvElement_GIDs;
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvVerts_IDs;
    std::vector<int> TotRecvOriVerts_IDs;
    std::vector<int> TotRecvFaces_IDs;
    std::vector<int> TotRecvFaces_Refs;
    std::vector<int> TotRecvOriFaces_IDs;

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
            }
            
            tt=cc/4;
        
            cc++;
            
        }
        
    }
    
    //unpack the face IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_f_map.begin();totrecv!=TotRecvElement_IDs_f_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvFaces_IDs.push_back(TotRecvElement_IDs_f_map[totrecv->first][r]);
            TotRecvOriFaces_IDs.push_back(TotRecvElement_IDs_of_map[totrecv->first][r]);
            TotRecvFaces_Refs.push_back(TotRecvElement_IDs_fref_map[totrecv->first][r]);
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
            sum2 = sum2+v_id_n;
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                tmesh->hybV2tetV[v_id_o_n]=v_id_n;
                tmesh->tetV2hybV[v_id_n]=v_id_o_n;
                
                r = FindRank(new_V_offsets, size, v_id_o_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_vert[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                    rank2req_vertOri[r].push_back(v_id_o_n);
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    vertIDs_on_rank_Ori.push_back(v_id_o_n);
                    
                    
                    
                    vloc_tmp++;
                }
                lv_id++;
            }
        }
        
        if(sum2==0)
        {
            std::cout << "help " << rank << " " << i << " " << TotNelem_recv << " " << on_rank << std::endl;
        }
        
        
        for(int k=0;k<4;k++)// looping over the vertices for element "i".
        {
            int f_id_n   = TotRecvFaces_IDs[cnt_f+k];
            int f_id_o_n = TotRecvOriFaces_IDs[cnt_f+k];
            int fref_n = TotRecvFaces_Refs[cnt_f+k];
            
            ief_part_tetra->setVal(on_rank,k,f_id_n);
            ief_part_hybrid->setVal(on_rank,k,f_id_o_n);
            iefref_part_tetra->setVal(on_rank,k,fref_n);
            
            if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_faceIDs_on_rank_set.insert(f_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_F_offsets,size,f_id_o_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_face[r].push_back(f_id_n); // add the vertex id that needs to be requested from rank r.
                    rank2req_faceOri[r].push_back(f_id_o_n);

                }
                else
                {
                    faceIDs_on_rank.push_back(f_id_n);  // add the vertex to list that is already available on rank.
//                  tmesh->hybF2tetF[f_id_o_n]=f_id_n;
//                  tmesh->tetF2hybF[f_id_n]=f_id_o_n;
                    tmesh->face2ref[f_id_n] = fref_n;
                    tmesh->ref2face[fref_n].push_back(f_id_n);
                    
                    if(rank==0 && f_id == 160591)
                    {
                        std::cout << "RANK " << rank << " fref = " << f_id_n << " " << fref_n << std::endl;
                    }
                    
                    floc_tmp++;
                }
                lf_id++;
            }
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
    
    
    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    

    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // =================================================================================
    // =================================================================================
    // =================================================================================
    
    int m = 0;
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;
    std::map<int,std::vector<int> >  reqstd_Ori_ids_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_vert.begin(); it != rank2req_vert.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                
                //	MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //	MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                MPI_Send(&rank2req_vertOri[it->first][0], n_req, MPI_INT, dest, 2229876*2+dest*2, comm);

                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            std::vector<int> recv_reqstd_Ori_ids(n_reqstd_ids);
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
                    
//                    if(it->second[u]==59944)
//                    {
//                        std::cout << xcn->getVal(it->second[u]-offset_xcn,0) << " "<< xcn->getVal(it->second[u]-offset_xcn,1) << " "<< xcn->getVal(it->second[u]-offset_xcn,2) << std::endl;
//                    }
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
        int c  = 0;
        vfor=vfor+recv_back_Nverts[it_f->first];
    }

    int gvid=0;
    int lvid=0;
    int gvid_gl = 0;

    for(m=0;m<vloc_tmp;m++)
    {
        gvid    = vertIDs_on_rank[m];
        gvid_gl = vertIDs_on_rank_Ori[m];
        
        Vert* V = new Vert;
        
        V->x = xcn->getVal(gvid_gl-xcn_o,0);
        V->y = xcn->getVal(gvid_gl-xcn_o,1);
        V->z = xcn->getVal(gvid_gl-xcn_o,2);
            
        tmesh->LocalVerts.push_back(V);
        tmesh->locV2globV[lvid] = gvid;
        tmesh->globV2locV[gvid] = lvid;
        lvid++;
    }
    
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = rank2req_vert[it_f->first][u];
            
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
    elms_to_send_to_ranks.clear();
    vertIDs_to_send_to_ranks.clear();
    OriVertIDs_to_send_to_ranks.clear();
    OriFaceIDs_to_send_to_ranks.clear();
    faceIDs_to_send_to_ranks.clear();
    rank2req_vert.clear();
    rank2req_vertOri.clear();
    rank2req_face.clear();
    rank2req_faceOri.clear();
    faceIDs_on_rank.clear();
    vertIDs_on_rank.clear();
    vertIDs_on_rank_Ori.clear();
    part_v.clear();
    loc_r_elem.clear();
    unique_vertIDs_on_rank_set.clear();
    unique_faceIDs_on_rank_set.clear();
    
    tmesh->ElGids               = ElGids;
    tmesh->ien_part_tetra       = ien_part_tetra;
    tmesh->ien_part_hybrid      = ien_part_hybrid;
    tmesh->ief_part_tetra       = ief_part_tetra;
    tmesh->ief_part_hybrid      = ief_part_hybrid;
    tmesh->iefref_part_tetra    = iefref_part_tetra;
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
	
	int nNonSharedVerts 		 = nLocalVerts-owned_verts[world_rank];
	int nNonSharedFaces 		 = nLocalFaces-owned_faces[world_rank];

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

   
	int* toS_red = new int[world_size];
	int* toR_red = new int[world_size];
	int* to_Send_copy = new int[world_size];
	int* to_Recv_copy = new int[world_size];
	int* optiSize = new int[world_size];
	int* optiSize_red = new int[world_size];

	int* toS = new int[world_size];
	int* toR = new int[world_size];
	for(int i=0;i<world_size;i++)
	{
		toS[i] = 0;
		toR[i] = 0;
		optiSize[i] = 0;
		optiSize_red[i] = 0;
		if(i==world_rank)
		{
			optiSize[i] = optimalSize;

			toS[i] = NtoSend;
			toR[i] = NtoRecv;
		}
	}
	MPI_Allreduce(optiSize, optiSize_red, world_size, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(toS,           toS_red, world_size, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(toR,           toR_red, world_size, MPI_INT, MPI_SUM, comm);

	int sent;
	int sendUpdate;
	
	std::map<int,std::vector<int> > recvRa;
	std::map<int,std::vector<int> > recvNe;
	
	std::map<int,std::vector<int> > sendRa;
	std::map<int,std::vector<int> > sendNe;
	
	for(int i=0;i<world_size;i++)
	{
		to_Send_copy[i] =  toS_red[i];
		to_Recv_copy[i] =  toR_red[i];
        
		if(toR_red[i]!=0)
		{
			for(int j=0;j<world_size;j++)
			{
				if(toS_red[j]>=toR_red[i])
				{
					sendUpdate = toS_red[j]-toR_red[i];
					sent       = toR_red[i];
				}
				
				if(toS_red[j]<toR_red[i])
				{
					sendUpdate = 0;
					sent       = toS_red[j];
				}
				
				toS_red[j] = sendUpdate;
				toR_red[i] = toR_red[i]-sent;
				
				if(sent>0)
				{
					recvRa[i].push_back(j);
					recvNe[i].push_back(sent);
					
					sendRa[j].push_back(i);
					sendNe[j].push_back(sent);
				}
			}
		}
	}
	
//    if(world_rank==0)
//    {
//        std::map<int,std::vector<int> >::iterator its;
//
//        for(its=recvRa.begin();its!=recvRa.end();its++)
//        {
//            std::cout << its->first << " -> ";
//            for(int q=0;q<its->second.size();q++)
//            {
//                std::cout << recvNe[its->first][q] << " ";
//            }
//        }
//    }
    
    
	//=================================================================================================
	//================================================================================================
	//=================================================================================================
	
	
	int nTet    = 0;
	int elloc   = 0;
	int lbvid   = nNonSharedVertsArrayOff[world_rank];
    //std::cout << nNonSharedVertsArrayOff[world_rank] << " " << nNonSharedVertsArrayOff[world_rank] << std::endl;
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
	
	int c = 0;
    int if_ref = 0;
    Array<int>* new_GidEl   = new Array<int>(optimalSize,4);
	Array<int>* new_ien     = new Array<int>(optimalSize,4);
    Array<int>* new_ien_or  = new Array<int>(optimalSize,4);
	Array<int>* new_ief     = new Array<int>(optimalSize,4);
    Array<int>* new_ief_or  = new Array<int>(optimalSize,4);
    Array<int>* new_iefref  = new Array<int>(optimalSize,4);
    
	if(sendRa.find(world_rank)!=sendRa.end())
	{
		std::vector<int> toRanks    = sendRa[world_rank];
		std::vector<int> NeltoRanks = sendNe[world_rank];
		
        std::vector<std::vector<int> > elIDs;
		std::vector<std::vector<int> > elNodeIDs;
		std::vector<std::vector<int> > elNodeOriginalIDs;
		std::vector<std::vector<int> > elFaceIDs;
        std::vector<std::vector<int> > elFaceRefs;
        std::vector<std::vector<int> > elFaceOriginalIDs;
        
		for(int i=0;i<toRanks.size();i++)
		{
			int Nel = NeltoRanks[i];
            std::vector<int> rowEl(Nel);
			std::vector<int> rowNode(Nel*4);
			std::vector<int> rowFace(Nel*4);
            std::vector<int> rowFaceRef(Nel*4);
			std::vector<int> rowNodeOriginal(Nel*4);
            std::vector<int> rowFaceOriginal(Nel*4);

            elIDs.push_back(rowEl);
			elNodeIDs.push_back(rowNodeOriginal);
			elNodeOriginalIDs.push_back(rowNode);
			elFaceIDs.push_back(rowFace);
            elFaceRefs.push_back(rowFaceRef);
            elFaceOriginalIDs.push_back(rowFaceOriginal);
		}
		
		int cc = 0;
		int sRank = toRanks[0];
		
		int offPrank = 0;
		int cntv     = 0;
		
		int t = 0;
		int nuloc = 0;
		int uloc = 0;
        int u = 0;

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
            
			for(int q=0;q<4;q++)
			{
				gvid = ite->second[q];
                
                ien_o[q]        = gvid;

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
						gl_map[gvid]    = lbvid;
						ien[q]          = lbvid;
						lbvid           = lbvid + 1;
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
			
            
            
			if(u<to_Send_copy[world_rank])
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

				uloc++;
			}
			
			u++;
		}
		
		int acull = 0;
		for(int i=0;i<toRanks.size();i++)
		{
			int dest = toRanks[i];
            int n_Ele = NeltoRanks[i];
			int n_Vrt = NeltoRanks[i]*4;
            std::vector<int> ElIDs      = elIDs[i];
			std::vector<int> Elvec      = elNodeIDs[i];
			std::vector<int> ElFvec     = elFaceIDs[i];
			std::vector<int> Elovec     = elNodeOriginalIDs[i];
            std::vector<int> ElFovec    = elFaceOriginalIDs[i];
            std::vector<int> ElFrefvec  = elFaceRefs[i];

			MPI_Send(&n_Vrt        ,     1, MPI_INT, dest, dest,         comm);
			MPI_Send(&Elvec[0]     , n_Vrt, MPI_INT, dest, dest*10000,   comm);
			MPI_Send(&ElFvec[0]    , n_Vrt, MPI_INT, dest, dest*50000,   comm);
			MPI_Send(&Elovec[0]    , n_Vrt, MPI_INT, dest, dest*20000,   comm);
            MPI_Send(&ElFrefvec[0] , n_Vrt, MPI_INT, dest, dest*2000000, comm);
            MPI_Send(&ElFovec[0]   , n_Vrt, MPI_INT, dest, dest*4000000, comm);
            MPI_Send(&n_Ele        ,     1, MPI_INT, dest, dest*5000000, comm);
            MPI_Send(&ElIDs[0]     , n_Ele, MPI_INT, dest, dest*6000000, comm);

			acull = acull + n_Vrt;
		}
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


		for(int i=0;i<expFromRank.size();i++)
		{
			int origin = expFromRank[i];
			int n_Elr;
			MPI_Recv(&n_Elr,   1, MPI_INT, origin, world_rank, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvNElVec(n_Elr);
			MPI_Recv(&recvNElVec[0], n_Elr, MPI_INT, origin, world_rank*10000, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvFElVec(n_Elr);
			MPI_Recv(&recvFElVec[0], n_Elr, MPI_INT, origin, world_rank*50000, comm, MPI_STATUS_IGNORE);
			
			std::vector<int> recvONElVec(n_Elr);
			MPI_Recv(&recvONElVec[0],n_Elr, MPI_INT, origin, world_rank*20000, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFrefElVec(n_Elr);
            MPI_Recv(&recvFrefElVec[0], n_Elr, MPI_INT, origin, world_rank*2000000, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recvFONElVec(n_Elr);
            MPI_Recv(&recvFONElVec[0], n_Elr, MPI_INT, origin, world_rank*4000000, comm, MPI_STATUS_IGNORE);
			
            int n_EleRecv;
            MPI_Recv(&n_EleRecv,   1, MPI_INT, origin, world_rank*5000000, comm, MPI_STATUS_IGNORE);

            std::vector<int> recvElIDsvec(n_EleRecv);
            MPI_Recv(&recvElIDsvec[0],   n_EleRecv, MPI_INT, origin, world_rank*6000000, comm, MPI_STATUS_IGNORE);

            collected_ElIds[origin]         = recvElIDsvec;
			collected_NIds[origin] 			= recvNElVec;
			collected_OriginalNIds[origin] 	= recvONElVec;
			collected_FIds[origin] 			= recvFElVec;
            collected_Frefs[origin]         = recvFrefElVec;
            collected_FOriginalNIds[origin] = recvFONElVec;

		}
		
		int el = 0;
        int u = 0;

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
						gl_map[gvid] = lbvid;
						ien[q]       = lbvid;
                        ien_o[q]     = gvid;
						lbvid        = lbvid + 1;
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
	}
    
    
//    for(int i=0;i<new_ien->getNrow();i++)
//    {
//        int sum = 0;
//        for(int j=0;j<new_ien->getNcol();j++)
//        {
//            sum = sum + new_ien->getVal(i,j);
//        }
//
//        if(sum == 0)
//        {
//            std::cout << world_rank << " " << i << std::endl;
//        }
//
//    }
    
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
    real_t ubvec_[]  = {1.02};
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
    
    Array<double>* xcn_ref = ReadDataSetFromFile<double>(fn_grid,"xcn");
    Array<int>* ien_ref    = ReadDataSetFromFile<int>(fn_conn,"ien");

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

    
    //std::cout << "understand bitshift " << world_rank << "  " << (world_rank & (~(1 << 3)))<< std::endl;
    //std::vector<int> LocElem     = P->getLocElem();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Domain* pDom = P->getPartitionDomain();
    std::map<int,std::vector<int> > tetras     = pDom->GTetras;
    std::set<int> ushell                       = pDom->ushell;

    
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* ien_part_map                   = P->getIENpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    Array<int>* part_global                    = P->getGlobalPartition();
    
    
    
    // Based on a partitioned hybrid mesh, we are extracting the tetra and redistribute them
    // uniformly. The hybrid mesh is partitioned without any weights which means that it might happen that the tetrahedra are initial distributed non-uniformly where a significant number of processors end up with any tetrahedra. This function distributed them equally.
    
    
    TetrahedraMesh* tmesh = ExtractTetrahedralMesh(part_global,tetras,
                                                   ief_part_map->i_map,
                                                   ifn_part_map->i_map,
                                                   ife_part_map->i_map,
                                                   if_ref_part_map->i_map,
                                                   ushell,
                                                   comm);
    
    // Once we have a uniform distribution of the tetrahedra, we determine a new OPTIMAL partitioning using the PARMETIS routine.
    PartitionInfo* partInfo = GetNewGlobalPartitioningTetrahedraMesh(tmesh,comm);

    // Determining the OPTIMAL number of elements on each rank based on the PARMETIS partitioning.
    int* newSizesOnRanks    = CommunicatePartitionLayout(partInfo->part,comm);

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

//    if(world_rank == 3)
//    {
//        for(int i=0;i<tmesh->ien_part_tetra->getNrow();i++)
//        {
//            std::cout << "Element  " << i << " -> ";
//            for(int j=0;j<tmesh->ien_part_tetra->getNcol();j++)
//            {
//                std::cout << tmesh->ien_part_tetra->getVal(i,j) << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    UpdateTetrahedraOnPartition(nTetrahedraGlob, nTetrahedra,
                                partInfo->part, tmesh,
                                us3d->xcn, xcn_pstate,
                                ife_pstate,
                                comm);

    ParallelState* ien_pstate_tet      = new ParallelState(nTetrahedraGlob,comm);

    std::map<int,int*> face2node = GetFace2EntityTetrahedraMesh(tmesh, us3d->ifn, 3, ife_pstate,nTetrahedraGlob,comm);

    //std::cout << world_rank << " number of faces = " << face2node.size() << " " << tmesh->hybF2tetF.size() << " " << tmesh->tetF2hybF.size() << std::endl;
    
//    std::map<int,int*>::iterator itst;
//
//    for(itst=face2node.begin();itst!=face2node.end();itst++)
//    {
//        std::cout << itst->first << " :: ";
//
//        for(int j=0;j<3;j++)
//        {
//            std::cout << itst->second[j] << " ";
//        }
//        std::cout << std::endl;
//    }
    
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
    std::vector<Vert*> locVs     = tmesh->LocalVerts;

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
    
    
    
    PartitionBoundary* pb = ExtractPartitionBoundary(tmesh,face2rank,face2node,comm);
    
    
    int* color_face = (int *) malloc(pb->ncomm*sizeof(int));
    int* ntifc = (int *) malloc(pb->ncomm*sizeof(int));
    int* ifc_tria_glob[pb->ncomm];
    int* ifc_tria_loc[pb->ncomm];
    
    int icomm=0;
    std::map<int,std::vector<int> >::iterator itc;

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
            
        }
        icomm++;
    }
    
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
//
//
//
//
    
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

          tensor[0] = 1.0/(0.1*0.1);
          tensor[1] = 0.0;
          tensor[2] = 0.0;
          tensor[3] = 1.0/(0.1*0.1);
          tensor[4] = 0.0;
          tensor[5] = 1.0/(0.1*0.1);

          if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],k+1)!=1)
          {
             MPI_Finalize();
             exit(EXIT_FAILURE);
          }
    }
    
    int v0,v1,v2,v3;
    int v0l,v1l,v2l,v3l;
    int teller = 0;
    
    for ( k=0; k<nTriangles; ++k )
    {
        int faceID = pb->faces4parmmg[k];

        v0 = face2node[faceID][0];
        v1 = face2node[faceID][1];
        v2 = face2node[faceID][2];
        
        v0l = tmesh->globV2locV[v0];
        v1l = tmesh->globV2locV[v1];
        v2l = tmesh->globV2locV[v2];

        if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,1.0,k+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, pb->ncomm);
    
    for(int icomm=0; icomm<pb->ncomm; icomm++ ) {

      /* Set nb. of entities on interface and rank of the outward proc */
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

      /* Set local and global index for each entity on the interface */
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

    
    /* remeshing function */
    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    
//
//
//
//

    /**/
    
    MPI_Finalize();
    
}

