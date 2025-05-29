
#include <chrono>
#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include "../../src/adapt_prismtetratrace.h"
#include "../../src/adapt_repartition.h"
#include "../../src/adapt_output_vtk.h"
#include "../../src/adapt_meshtopology_lite.h"
#include "../../src/adapt_gradreconstruct_lite.h"
#include "../../src/adapt_runparmmg.h"
#include "../../src/adapt_inputs.h"
#include "../../src/adapt_writeus3ddata.h"
#include "../../src/adapt_parmetisstate_lite.h"
#include "../../src/adapt_operations.h"
#include "../../src/adapt_partobject_lite.h"
#include "../../src/adapt.h"
// #include <cstdio>
#include <stdio.h>
//#include "/Users/dekelsch/mimic_libmesh/utilities/partitionTetrahedra/build/ThirdParty/dist/include/libmeshb7.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;



// ===================== Hex ordering =============================
//  top 
//   7------6
//   |      |
//   |      | 
//   |      |
//   4------5
//  bottom
//   3------2
//   |      |
//   |      | 
//   |      |
//   0------1
// ===================== Hex ordering =============================



void ReadSU2Mesh(MPI_Comm comm,
                const char* fm,
                std::vector<std::vector<int> > &elements,
                std::vector<std::vector<int> > &tetraElements,
                std::vector<int> &offsetsTetra,
                std::vector<int> &nlocsTetra,
                std::vector<std::vector<int> > &otherElements,
                std::vector<int> &offsetsOther,
                std::vector<int> &NlocsOther,
                std::map<int,int> &element_type,
                std::vector<int> &offsets,
                std::vector<int> &nlocs,
                std::vector<int> &nlocs_vrts,
                std::vector< std::vector<double> > &nodes,
                std::map<int,std::vector<std::vector<int> > > &bcFaces,
                std::map<int,int> &bcNFaces,
                std::map<int,std::vector<int> > &bcFaceSizes,
                std::map<int,char*> &bcTags,
                FaceSetPointer &allbcFaces)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::ifstream meshFile;
    meshFile.open(fm);
    std::string   line;
    if(!meshFile.is_open())
    {
        std::cout << "Error:: Make sure there is a mesh.su2 file in the directory where main.cpp resides. "<<std::endl;
        exit(0);
    }

    bool processHeaderElements = true;
    int lstart = 0;
    int ndim   = 0; 
    int nelem  = 0;
    int npoin  = 0;
    int c      = 0;

    //================================PROCESS ELEMENTS =====================================================
    std::map<int,std::vector<int> > elements_map;
    while (processHeaderElements) 
    {
        std::getline(meshFile, line);    // consider: while (getline(inputfile, line).good())
        std::stringstream ssline(line);

        if (line.find("NDIME=") != std::string::npos)
        {
                std::istringstream iss(line);
                std::string key;
                iss >> key >> ndim;
        }

        if (line.find("NELEM=") != std::string::npos)
        {
                std::istringstream iss(line);
                std::string key;
                iss >> key >> nelem;
                lstart = c;
                processHeaderElements = false;
        }
        c++;
    }


    int linesToSkip = c + 1;
    bool processElements    = true;
    int var                 = 0;
    int elid                = 0;
    int offset_count        = 0;
    int offset_tetra_count  = 0;
    int offset_other_count  = 0;

    while (processElements) 
    {
        std::getline(meshFile, line); 
        if (line.find("%") != std::string::npos)
        {
                processElements = false;
        }
        else
        {
            std::istringstream iss(line);
            std::vector<int> row;
            int num;
            while (iss >> num) {
                row.push_back(num);
            }
            //std::cout << "row size should be 6 " << row.size() << std::endl;
            element_type[elid] = row[0];
            offsets.push_back(offset_count);
            
            elements.push_back(row);
            elements_map[elid] = row;
            // if(row.size()!=6)
            // {
            //     std::cout << row.size() << " " << world_rank << " " << row[0] << std::endl;
            // }
            if(row[0] == 10)
            {
                tetraElements.push_back(row);
                offsetsTetra.push_back(offset_tetra_count);
                offset_tetra_count = offset_tetra_count + row.size();
                nlocsTetra.push_back(row.size());
            }
            // if(row[0] == 12)
            // {
            //     tetraElements.push_back(row);
            //     offsetsTetra.push_back(offset_tetra_count);
            //     offset_tetra_count = offset_tetra_count + row.size();
            //     nlocsTetra.push_back(row.size());
            // }
            else
            {   
                // std::cout << "row.size " << row.size() << " " << row[0] << std::endl;
                otherElements.push_back(row);
                offsetsOther.push_back(offset_other_count);
                offset_other_count = offset_other_count + row.size();
                NlocsOther.push_back(row.size());
            }

            elid = elid + 1;
            offset_count=offset_count+row.size();
        }
    }

    //================================PROCESS NODES =====================================================
    bool processHeaderNodes = true;

    while (processHeaderNodes) 
    {
        std::getline(meshFile, line);    
        std::stringstream ssline(line);

        if (line.find("NPOIN=") != std::string::npos)
        {
                std::istringstream iss(line);
                std::string key;
                iss >> key >> npoin;
                lstart = c;
                processHeaderNodes = false;
        }
        c++;
    }

    int nid  = 0;

    bool processNodes = true;
    int offset_vrts_count = 0;
    while (processNodes) 
    {
        std::getline(meshFile, line); 
        if (line.find("%") != std::string::npos)
        {
                processNodes = false;
        }
        else
        {
            std::istringstream iss(line);
            std::vector<double> row;
            double num;
            while (iss >> num) {
                row.push_back(num);
            }
            nlocs_vrts.push_back(row.size());
            offset_vrts_count=offset_vrts_count+row.size();

            nodes.push_back(row);
            nid = nid + 1;
        }
    }


    //================================PROCESS BOUNDARY FACES ===============================================

    bool processHeaderBoundaryFaces = true;
    int nmark = -1;
    while (processHeaderBoundaryFaces) 
    {
        std::getline(meshFile, line);    
        std::stringstream ssline(line);


        if (line.find("NMARK=") != std::string::npos)
        {
                std::istringstream iss(line);
                std::string key;
                iss >> key >> nmark;
                processHeaderBoundaryFaces = false;
        }
    }


    bool processBoundaryFaces = true;
    int m = 0;
    int cnt = 0;
    int mark_elem = 0;
    int nf  = 0; 
    int fid = 0;
    int ntri = 0;
    int nquad = 0;
    std::map<int,int> bcTriFaces;
    std::map<int,int> bcQuadFaces;
    std::vector<std::vector<int> > VV_bcFaces;
    while (processBoundaryFaces) 
    {
        std::getline(meshFile, line); 
        if (line.find("%") != std::string::npos)
        {
                processBoundaryFaces = false;
        }
        else
        {
            if (line.find("MARKER_TAG=") != std::string::npos)
            {
                    m++;
                    nf = 0;
                    ntri = 0;
                    nquad = 0;
                    char* tag = new char[20];
                    std::istringstream iss(line);
                    std::string key;
                    iss >> key >> tag;
                    bcTags[m+1]=tag;
            }

            else if (line.find("MARKER_ELEMS=") != std::string::npos)
            {
                    mark_elem = 0;
                    std::istringstream iss(line);
                    std::string key;
                    iss >> key >> mark_elem;
                                        
                    cnt = 0;
                    
            }
            else
            {
                std::istringstream iss(line);
                std::vector<int> row;
                int num;
                std::vector<int> face;
                int te = 0;
                while (iss >> num) {
                    row.push_back(num);
                    if(te>0)
                    {
                        face.push_back(num);
                    }
                    te++;
                }
                

                FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));
                
                FaceSetPointer::iterator testInsPointer2 = allbcFaces.find(facePointer);
                if(testInsPointer2==allbcFaces.end())
                {

                    std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer      = allbcFaces.insert(facePointer);
                    int fref            = m;
                    (*testInsPointer.first)->SetFaceRef(fref);
                    //(*testInsPointer.first)->SetFaceID(-fid);
                    bcFaces[m+1].push_back(face);
                    int nvpf            = face.size();
                    nf                  = nf + face.size();
                    bcNFaces[m+1]       = nf;
                    bcFaceSizes[m+1].push_back(nvpf); 
                    //fid++;
                }
                
                nid = nid + 1;

                if((m == (nmark)) && (cnt == mark_elem-1))
                {
                    processBoundaryFaces = false;
                }
                cnt = cnt + 1;
            }


        }
        
    }

    // FaceSetPointer m_FaceSetPointer;

    // // update faceid for the boundary faces
    // FaceSetPointer::iterator ftit;
    // for(ftit=allbcFaces.begin();ftit!=allbcFaces.end();ftit++)
    // {
    //     int ofid = (*ftit)->GetFaceID();
    //     int nfid = ofid;
    //     (*ftit)->SetFaceID(nfid);
        
    //     m_FaceSetPointer.insert((*ftit));
    // }
}





void CommunicateBoundaryFaceData(MPI_Comm comm, 
    std::map<int,std::vector<std::vector<int> > > bcFaces_onRoot, 
    std::map<int,int> bcNFaces, 
    std::map<int,std::vector<int> > bcFacesizes_onRoot, 
    FaceSetPointer &allbcFaces)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::map<int,std::vector<std::vector<int> > >::iterator itmb;

    std::vector<int> keys(bcFaces_onRoot.size(),0);
    std::vector<int> nfcs(bcFaces_onRoot.size(),0);

    int ii      = 0;
    int cntr    = 0;    

    int face_offset = 0;
    std::vector<int> face_offsets(bcFaces_onRoot.size(),0);
    int s = 0;
    for(itmb=bcFaces_onRoot.begin();itmb!=bcFaces_onRoot.end();itmb++)
    {
        int bfref       = itmb->first;
        cntr            = cntr + bcNFaces[bfref];
        face_offsets[s] = face_offset;
        face_offset     = face_offset + bcFacesizes_onRoot[bfref].size();
        s = s + 1;
    }

    std::vector<int> bc_vids(cntr,0);
    std::vector<int> bc_fids(face_offset,0);
    int offs    = 0;
    int y       = 0;
    for(itmb=bcFaces_onRoot.begin();itmb!=bcFaces_onRoot.end();itmb++)
    {
        keys[ii] = itmb->first;
        nfcs[ii] = itmb->second.size();

        for(int q=0;q<itmb->second.size();q++)
        {
            int nv      = itmb->second[q].size();
            bc_fids[y]  = nv;

            for(int p=0;p<nv;p++)
            {
                bc_vids[offs+p] = itmb->second[q][p];   
            }

            y    = y + 1;
            offs = offs + nv;
        }

        ii = ii + 1;
    }

    int bvid_before = bc_vids.size();

    int keys_size = keys.size();
    MPI_Bcast(&keys_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int bc_fids_size = bc_fids.size();
    MPI_Bcast(&bc_fids_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int bc_vids_size = bc_vids.size();
    MPI_Bcast(&bc_vids_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank != 0) 
    {
        keys.resize(keys_size);
        nfcs.resize(keys_size);
        bc_fids.resize(bc_fids_size);
        bc_vids.resize(bc_vids_size);
    }

    MPI_Bcast(keys.data(), keys_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nfcs.data(), keys_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(bc_fids.data(), bc_fids_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(bc_vids.data(), bc_vids_size, MPI_INT, 0, MPI_COMM_WORLD);
    int fid = 0;

    if(world_rank != 0)
    {
        int foffset = 0;
        int offs    = 0;
        for(int i=0;i<keys.size();i++)
        {
            int m  = keys[i];
            int nd = nfcs[i];
            for(int j=0;j<nd;j++)
            {
                int nvpf = bc_fids[foffset + j];
                std::vector<int> face(nvpf,0);
                for(int k=0;k<nvpf;k++)
                {   
                    face[k] = bc_vids[offs+k];
                }

                FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));

                FaceSetPointer::iterator testInsPointer = allbcFaces.find(facePointer);
                if(testInsPointer==allbcFaces.end())
                {
                    std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer  = allbcFaces.insert(facePointer);
                    int fref        = m;
                    (*testInsPointer.first)->SetFaceRef(fref);
                    (*testInsPointer.first)->SetFaceID(-fid);
                    fid++;
                }

                offs = offs + nvpf;
            }
            foffset = foffset + nd;
        }
    }


    // ===============================TEST OUT THE COMMUNICATION OF THE MAP ================================
    // FaceSetPointer::iterator ftit;
    // int teller = 0;
    // for(ftit=allbcFaces.begin();ftit!=allbcFaces.end();ftit++)
    // {
    //     std::vector<int> vrts = (*ftit)->GetEdgeIDs();

    //     if(world_rank == 2 && teller < 10)
    //     {
    //         std::cout << "WR = " << world_rank << " -- " << teller << " :: ";
    //         for(int i=0;i<vrts.size();i++)
    //         {
    //             std::cout << vrts[i] << " ";
    //         }
    //         std::cout << std::endl;
    //     }

    //     if(world_rank == 3 && teller < 10)
    //     {
    //         std::cout << "WR = " << world_rank << " -- " << teller << " :: ";
    //         for(int i=0;i<vrts.size();i++)
    //         {
    //             std::cout << vrts[i] << " ";
    //         }
    //         std::cout << std::endl;

    //     }

    //     teller = teller + 1;
    // }
    // ===============================TEST OUT THE COMMUNICATION OF THE MAP ================================
    
}

struct UniformMesh{
    std::map<int,std::vector<int> > e2v;
    std::map<int, std::vector<double> > vert_local;
    std::map<int,int> eltype_map;
    std::vector<int> eltype_vec;
    FaceSetPointer allbcFaces;
};


UniformMesh* ScatterElements(MPI_Comm comm, ParallelState* pstate_Elem, ParallelState* pstate_Vert, std::vector<std::vector<int> > Elements, std::vector<int> nlocs_elem, std::vector<std::vector<double> > nodes, std::vector<double> vertices_root_flatten, int nvrts, std::vector<int> nlocs_vrts)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    UniformMesh* umesh = new UniformMesh;

    std::vector<int> Elements_root_flatten;
    for (const auto& row : Elements) 
    {
            Elements_root_flatten.insert(Elements_root_flatten.end(), row.begin(), row.end());
    }
    int nElem = Elements.size();
    std::cout << "nElemnElemnElemnElem " << nElem << std::endl; 
    std::vector<int> sendcounts(world_size,0);
    std::vector<int> displs(world_size,0);

    std::vector<int> sendcountsV(world_size,0);
    std::vector<int> displsV(world_size,0);

    std::vector<int> sendcounts_faces(world_size,0);
    std::vector<int> displs_faces(world_size,0);
    std::vector<int> nlocals(world_size,0);
    std::vector<int> noffsets(world_size,0);
    std::vector<int> nlocalsv(world_size,0);
    std::vector<int> noffsetsv(world_size,0);

    if(world_rank == 0)
    {
        int end,start,index;
        int end_vrts,start_vrts,start_faces,end_faces;
        int cnt      = 0;
        int nlocal   = 0;
        int nlocalv  = 0;
        int stp      = pstate_Elem->getNlocs()[cnt];
        int stpv     = pstate_Vert->getNlocs()[cnt];

        int noffset  = 0;
        int noffsetv = 0;
        for(int i=0;i<nElem;i++)
        {
            nlocal  = nlocal  + Elements[i].size();
            // nlocalv = nlocalv + nodes[i].size();
            // std::cout << nodes[i].size() << std::endl;
            if(i==(stp-1))
            {
                nlocals[cnt]    = nlocal; 
                noffsets[cnt]   = noffset;
                cnt             = cnt + 1;
                stp             = stp + pstate_Elem->getNlocs()[cnt];
                noffset         = noffset + nlocal;
                nlocal          = 0;
            }
        }
        std::cout << "nvrts " << nvrts << "  " << nvrts*4 << " " << stpv << std::endl;
        std::cout << "nElem " << nElem << std::endl;
        cnt = 0;
        for(int i=0;i<nvrts;i++)
        {
            nlocalv = nlocalv + nodes[i].size();
            if(i==(stpv-1))
            {
                nlocalsv[cnt]  = nlocalv; 
                noffsetsv[cnt] = noffsetv;
                cnt      = cnt + 1;
                stpv     = stpv + pstate_Vert->getNlocs()[cnt];
                noffsetv = noffsetv + nlocalv;
                nlocalv  = 0;

            }
        }
    }
    
    // std::vector<int> nlocals_red(world_size,0);
    // MPI_Allreduce(nlocals.data(), 
    //             nlocals_red.data(), 
    //             world_size, MPI_INT, MPI_SUM, comm);

    // std::vector<int> nlocalsV_red(world_size,0);
    // MPI_Allreduce(nlocalsV.data(), 
    //             nlocalsV_red.data(), 
    //             world_size, MPI_INT, MPI_SUM, comm);      
    
    if(world_rank == 0)
    {
        //std::cout << "wr " << world_rank << " "  << nlocal << " " << tetraElements_root_flatten.size() << std::endl;
        int end,start,index,startV,endV;
        int end_vrts,start_vrts,start_faces,end_faces;
        for(int i=0;i<world_size;i++)
        {
            int offset          = pstate_Elem->getOffsets()[i];
            int nloc            = pstate_Elem->getNlocs()[i];
            start               = noffsets[i];
            end                 = noffsets[i]+nlocals[i];
            sendcounts[i]       = end-start;
            displs[i]           = start;


            int offsetV          = pstate_Vert->getOffsets()[i];
            int nlocV            = pstate_Vert->getNlocs()[i];
            startV               = noffsetsv[i];
            endV                 = noffsetsv[i]+nlocalsv[i];
            sendcountsV[i]       = endV-startV;
            displsV[i]           = startV;
            std::cout << endV-startV << " " << i << " " << end-start << std::endl;
            
        }
    }
  
    std::vector<int> sendcounts_red(world_size,0);
    MPI_Allreduce(sendcounts.data(), 
                  sendcounts_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);

    std::vector<int> displs_red(world_size,0);
    MPI_Allreduce(displs.data(), 
                  displs_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);
  
    std::vector<int> sendcountsV_red(world_size,0);
    MPI_Allreduce(sendcountsV.data(), 
                    sendcountsV_red.data(), 
                    world_size, MPI_INT, MPI_SUM, comm);

    std::vector<int> displsV_red(world_size,0);
    MPI_Allreduce(displsV.data(), 
                    displsV_red.data(), 
                    world_size, MPI_INT, MPI_SUM, comm);

    std::cout << "check stat " << world_rank << " :: " << sendcounts_red[world_rank] << " " << sendcountsV_red[world_rank] << std::endl;
    
    std::vector<int> elements_on_rank(sendcounts_red[world_rank],0);
    std::vector<double> verts_on_rank(sendcountsV_red[world_rank],0);

    if(world_rank == 0)
    {
        std::cout << "Scattering the mesh to other processors..." << std::endl;
    }

    
    MPI_Scatterv(Elements_root_flatten.data(), 
                 sendcounts_red.data(), displs_red.data(), MPI_INT, 
                 elements_on_rank.data(), sendcounts_red[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);
    
    MPI_Scatterv(vertices_root_flatten.data(), 
                sendcountsV_red.data(), displsV_red.data(), MPI_DOUBLE, 
                verts_on_rank.data(), sendcountsV_red[world_rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    
    std::vector<int>e2v_loc(pstate_Elem->getNlocs()[world_rank],0);

    std::vector<int>v_loc(pstate_Vert->getNlocs()[world_rank],0);
    
    MPI_Scatterv(nlocs_elem.data(), 
                 pstate_Elem->getNlocs(), pstate_Elem->getOffsets(), MPI_INT, 
                 e2v_loc.data(), pstate_Elem->getNlocs()[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(nlocs_vrts.data(), 
                 pstate_Vert->getNlocs(), pstate_Vert->getOffsets(), MPI_INT, 
                 v_loc.data(), pstate_Vert->getNlocs()[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);
    
    Elements_root_flatten.resize(0);
    std::cout << "e2v_loc " << e2v_loc.size() << std::endl;
    // std::map<int,int> eltype_map;
    // std::map<int,std::vector<int> > e2v;
    
    int off = 0;
    // std::vector<int> eltype_vec;
    for(int i = 0;i<e2v_loc.size();i++)
    {
        int nloc   = e2v_loc[i];
        //std::cout << "world rank " << world_rank << " " << nloc << " " << elements_on_rank.size() << std::endl;
        int elid   = elements_on_rank[off+nloc-1];
        int eltype = elements_on_rank[off];

        umesh->eltype_map[elid] = eltype;
        umesh->eltype_vec.push_back(2);
        std::vector<int>e2v_row(nloc-2,0);
        int c = 0;
        for(int j=1;j<nloc-1;j++)
        {
            e2v_row[c] = elements_on_rank[off+j];
            //std::cout << e2v_row[c] << " ";
            c++;
        }
        //std::cout << std::endl;
       
        umesh->e2v[elid] = e2v_row;
        off = off+nloc;

    }
    
    // std::cout << "e2v_loc " << e2v_loc.size() << std::endl;

    // std::map<int,std::vector<double> > vert_local;
    int vid = pstate_Vert->getOffsets()[world_rank];
    off = 0;
    std::map<int,int> lv2gv;
    std::map<int,int> gv2lv;
    int lvid = 0;
    int fnd = 0;
    
    for(int i = 0;i<v_loc.size();i++)
    {   
        int nloc   = v_loc[i];
        int gvid   = verts_on_rank[off+nloc-1];
        std::vector<double>v_row(3,0);
        if(nloc != 4)
        {   
            std::cout << i << " nloc " << nloc << std::endl;
        }
        

        for(int j=0;j<3;j++)
        {
            v_row[j] = verts_on_rank[off+j];
        }

        lv2gv[lvid] = gvid; 
        gv2lv[gvid] = lvid; 

        if(vid!=gvid)
        {
            std::cout << vid << " compare ids " << gvid << " " << off+nloc-1 << std::endl;
            fnd++;
        }
        

        umesh->vert_local[vid] = v_row;
        vid = vid + 1;
        lvid = lvid + 1;
        off = off + nloc;
    }
    /**/

    //renumber
    // std::map<int,std::vector<int> >::iterator re;
    // std::vector<int> new_V_offsets(world_size,0);

    // for(int i=0;i<world_size;i++)
    // {
    //     new_V_offsets[i] = pstate_Vert->getOffsets()[i]-1;
    // }

    // std::map<int,std::vector<int> > e2v2r;
    // std::map<int,std::set<int> > request_vid;
    // for(re=e2v.begin();re!=e2v.end();re++)
    // {
    //     std::vector<int> row_v2r(re->second.size(),0);


    //     for(int q=0;q<re->second.size();q++)
    //     {
    //         int vid         = re->second[q];
    //         int r           = FindRank(new_V_offsets.data(),world_size,vid);

    //         if(r!=world_rank)
    //         {
    //             request_vid[r].insert(vid);
    //         }


    //         row_v2r[q]      = r;
    //     }

    //     e2v2r[re->first]    = row_v2r;

    // }

    // int Ne = nelem;
    // int Nv = nvrts;

    return umesh;
}





int main(int argc, char** argv)
{
    clock_t start_total = clock();
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
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;


    int64_t LibIdx;
    int ver, dim, NmbVer, NmbTri, NmbTet, NmbPri;
    
    //const char* fm = "hemihyb.su2";
    const char* fm = "hemihyb_test2.su2";
    const char* fd = "hemihyb_test2.h5";

    std::vector<std::vector<int> > elements;
    std::vector<std::vector<int> > tetraElements;
    std::vector<std::vector<int> > otherElements;
    std::vector<int> offsets;
    std::vector<int> nlocs;
    std::vector<int> nlocsTetra;
    std::vector<int> nlocsOther;
    std::vector<int> nlocs_vrts;
    std::vector<std::vector<double> > nodes;
    std::map<int,std::vector<std::vector<int> > > bcFaces_onRoot;
    std::map<int,int> bcNFaces_onRoot;
    std::map<int,std::vector<int> > bcFaceSizes_onRoot;
    std::map<int,char*> bcTags;
    std::vector<int> elements_root_flatten;
    std::vector<int> tetraElements_root_flatten;
    std::vector<int> otherElements_root_flatten;
    std::vector<double> vertices_root_flatten;
    std::map<int,std::vector<double> > vertices;
    std::map<int,int> element_type;
    int nelem = 0;
    int nvrts = 0;
    double start1, end1,time_taken1;
    start1 = clock();
    // FaceSetPointer allbcFaces;
    std::vector<int> offsetsTetra;
    std::vector<int> offsetsOther;
    FaceSetPointer allbcFaces;
    int ntetraElem = 0;
    int notherElem = 0;
    if(world_rank == 0)
    {
        std::cout << "Starting to read the file on root..." << std::endl;
        
        ReadSU2Mesh(comm,fm, 
                    elements,
                    tetraElements,
                    offsetsTetra,
                    nlocsTetra,
                    otherElements,
                    offsetsOther,
                    nlocsOther,
                    element_type,
                    offsets, 
                    nlocs, 
                    nlocs_vrts, nodes, bcFaces_onRoot, bcNFaces_onRoot, bcFaceSizes_onRoot, bcTags, allbcFaces);

        for (const auto& row : nodes) 
        {
            vertices_root_flatten.insert(vertices_root_flatten.end(), row.begin(), row.end());
        }
        int nVertices = nodes.size();

        for (const auto& row : tetraElements) 
        {
                tetraElements_root_flatten.insert(tetraElements_root_flatten.end(), row.begin(), row.end());
        }
        ntetraElem = tetraElements.size();
        
        // tetraElements.clear();
        for (const auto& row : otherElements) 
        {
                // std::cout << "row " << row.size() << std::endl;
                otherElements_root_flatten.insert(otherElements_root_flatten.end(), row.begin(), row.end());
        }
        notherElem = otherElements.size();
        //otherElements.clear();

        for (const auto& row : elements) 
        {
                // std::cout << "row " << row.size() << std::endl;
                elements_root_flatten.insert(elements_root_flatten.end(), row.begin(), row.end());
        }
        nelem = elements.size();
        elements.clear();
        nvrts = nodes.size();
        std::cout << "nelem " << nelem <<  " " << nvrts << std::endl;

        

    }

    // Communicate all boundary face data to all other ranks:
    
    CommunicateBoundaryFaceData(comm, bcFaces_onRoot, bcNFaces_onRoot, bcFaceSizes_onRoot, allbcFaces);


    //std::cout << world_rank << "  allbcFaces  " << allbcFaces.size() << std::endl; 

    
    //================ End of communication boundary data to all Processors. =========================================

    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    if(world_rank == 0)
    {
        cout << setprecision(16) << "Time taken to read in the mesh.su2 on root :                       " << fixed
        << time_taken1 << setprecision(16);
        cout << " sec " << endl;
    }

    MPI_Bcast(&nelem, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nvrts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntetraElem, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&notherElem, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<std::vector<double> > t_hessian = ReadHyperSolveHessianData(fd,nvrts,comm,info);
    
    std::cout << "t_hessian " << t_hessian.size() << " " << t_hessian[0].size()<< std::endl;
    
    ParallelState* pstate_tetraElem = new ParallelState(ntetraElem,comm);
    ParallelState* pstate_otherElem = new ParallelState(notherElem,comm);
    ParallelState* pstate_Vert      = new ParallelState(nvrts,comm);

    UniformMesh* tetUniMesh = ScatterElements(comm, 
                                                pstate_tetraElem, 
                                                pstate_Vert, 
                                                tetraElements, 
                                                nlocsTetra, 
                                                nodes, 
                                                vertices_root_flatten,
                                                nvrts, 
                                                nlocs_vrts);

    UniformMesh* blUniMesh = ScatterElements(comm, 
                                                pstate_otherElem, 
                                                pstate_Vert, 
                                                otherElements, 
                                                nlocsOther, 
                                                nodes, 
                                                vertices_root_flatten,
                                                nvrts, 
                                                nlocs_vrts);


    int Ne = nelem;
    int Nv = nvrts;
    
    PartObjectLite* partitionT = new PartObjectLite(tetUniMesh->e2v, tetUniMesh->vert_local, tetUniMesh->eltype_map, tetUniMesh->eltype_vec, allbcFaces, Ne, Nv, comm);
    PartObjectLite* partitionP = new PartObjectLite(blUniMesh->e2v,   blUniMesh->vert_local,  blUniMesh->eltype_map,  blUniMesh->eltype_vec, allbcFaces, Ne, Nv, comm);

    FaceSetPointer sharedTfaces = partitionT->getAllSharedFaceMap();
    FaceSetPointer sharedPfaces = partitionP->getAllSharedFaceMap();

    tetUniMesh->e2v.clear();
    tetUniMesh->vert_local.clear();
    tetUniMesh->eltype_map.clear();
    tetUniMesh->eltype_vec.clear();

    blUniMesh->e2v.clear();
    blUniMesh->vert_local.clear();
    blUniMesh->eltype_map.clear();
    blUniMesh->eltype_vec.clear();

    FaceSetPointer InterFaceFaces;
    std::cout << "shared faces per partition T -> "<< sharedTfaces.size() << " P -> " << sharedPfaces.size() << std::endl;
    FaceSetPointer::iterator ftit;
    int shell = 0;
    if(sharedPfaces.size()>=sharedTfaces.size())
    {
        for(ftit=sharedPfaces.begin();ftit!=sharedPfaces.end();ftit++)
        {
            std::vector<int> fvid                   = (*ftit)->GetEdgeIDs();
            FaceSharedPtr facePointer               = std::shared_ptr<NekFace>(new NekFace(fvid));

            if(sharedTfaces.find(facePointer)!=sharedTfaces.end())
            {   
                std::pair<FaceSetPointer::iterator, bool> testInsPointer = InterFaceFaces.insert(facePointer);
                shell++;
            }
            
        }
    }
    else
    {
        for(ftit=sharedTfaces.begin();ftit!=sharedTfaces.end();ftit++)
        {
            std::vector<int> fvid                   = (*ftit)->GetEdgeIDs();
            FaceSharedPtr facePointer               = std::shared_ptr<NekFace>(new NekFace(fvid));

            if(sharedPfaces.find(facePointer)!=sharedPfaces.end())
            {   
                std::pair<FaceSetPointer::iterator, bool> testInsPointer = InterFaceFaces.insert(facePointer);
                shell++;
            }
        }
    }

    std::cout << "Shellfaces = " << shell << std::endl;


    MPI_Finalize();
        
}

