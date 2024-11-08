
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
#include "../../src/adapt_operations.h"
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
                std::vector<std::vector<int> > &element2face,
                std::map<int,std::vector<int> > &element2face_map,
                std::map<int,std::vector<int> > &bcfaces_map,
                std::map<int,int> &element_type,
                std::vector<int> &offsets,
                std::vector<int> &offsets_vrts,
                std::vector<int> &offsets_faces,
                std::vector<int> &nlocs,
                std::vector<int> &nlocs_vrts,
                std::vector<int> &nlocs_faces,
                std::map<int,std::vector<double> > &nodes,
                std::map<int,std::vector<std::vector<int> > > &bcFaces,
                std::map<int,char*> &bcTags)
{
    
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
    bool processElements = true;
    int var = 0;



    int elid = 0;
    int offset_count = 0;
    while (processElements) 
    {
        std::getline(meshFile, line); 
        if (line.find("%") != std::string::npos)
        {
                processElements = false;
        }
        else
        {
            // std::getline(meshFile, line);    // consider: while (getline(inputfile, line).good())
            std::istringstream iss(line);
            std::vector<int> row;
            int num;
           // std::cout << "elid " << elid << " :: ";
            while (iss >> num) {
                row.push_back(num);
                //std::cout << num << " ";
            }

            element_type[elid]  = row[0];
            //std::cout << std::endl;
            offsets.push_back(offset_count);
            // offsets_vrts.push_back(offset_count);
            nlocs.push_back(row.size());
            elements.push_back(row);
            elements_map[elid] = row;

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
            offsets_vrts.push_back(offset_vrts_count);
            nlocs_vrts.push_back(row.size());
            offset_vrts_count=offset_vrts_count+row.size();

            nodes[nid] = row;
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
    FaceSetPointer allbcFaces;
    int fid = 0;
    
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
                    char* tag = new char[20];
                    std::istringstream iss(line);
                    std::string key;
                    iss >> key >> tag;
                    bcTags[m+1]=tag;
                    std::cout << " " << m << " tag " << tag  << " end tag "<< std::endl;
            }

            else if (line.find("MARKER_ELEMS=") != std::string::npos)
            {
                    mark_elem = 0;
                    std::istringstream iss(line);
                    std::string key;
                    iss >> key >> mark_elem;
                    
                    std::cout << "bc " << m << " mark_elem " << mark_elem << std::endl;
                    
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
                    // std::cout << "tag = " << tag << std::endl;
                    std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer = allbcFaces.insert(facePointer);
                    int fref = m;
                    (*testInsPointer.first)->SetFaceRef(fref);
                    (*testInsPointer.first)->SetFaceID(-fid);
                    bcFaces[m].push_back(face);
                    bcfaces_map[fid] = face;
                    fid++;
                    
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

    std::cout << "allbcFaces " << allbcFaces.size() << std::endl;

    std::vector<std::vector<int> > e2f_map;
    FaceSetPointer m_FaceSetPointer;
    // FaceSetPointer allbcFaces;
    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    std::map<int,std::vector<int> >::iterator ite;
    int fintid          = 0;
    int offset_fcount   = 0;
    int fo = 0;
    int fint = 0;
    for(ite=elements_map.begin();ite!=elements_map.end();ite++)
    {
        int elid                = ite->first;
        int eltype              = ite->second[0];
        int nloc                = ite->second.size();
        std::vector<int> e2n_i  = ite->second;
        std::vector<int> e2v_row(nloc-2,0);
        int c = 0;
        for(int j=1;j<nloc-1;j++)
        {
            e2v_row[c] = e2n_i[j];
            c++;
        }


        switch (eltype) {
        case 10:
            e2f_map = tetra_faces;
    
            break;
        case 12:
            e2f_map = hex_faces;
       
            break;
        case 13:
            e2f_map = prism_faces;
          
            break;
        case 14:    
            e2f_map = pyramid_faces;
           
            break;
        }
        
        int nfaces = e2f_map.size();
        std:vector<int> faceids(nfaces,0);
        std::map<int,std::vector<int> > gface2bcface;
        for(int f=0;f<nfaces;f++)
        {
            int nnodes = e2f_map[f].size();

            std::vector<int> face_globvid(nnodes,0);

            for(int n=0;n<nnodes;n++)
            {
                face_globvid[n] = e2v_row[e2f_map[f][n]];
            }
            
            FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face_globvid));
            
            FaceSetPointer::iterator testInsPointer2 = m_FaceSetPointer.find(facePointer);
            if(testInsPointer2==m_FaceSetPointer.end())
            {


                // (*testInsPointer.first)->SetFaceID(fintid);
                FaceSetPointer::iterator testInsPointer3 = allbcFaces.find(facePointer);
                if(testInsPointer3!=allbcFaces.end())
                {
                    int fref    = (*testInsPointer3)->GetFaceRef();
                    int bfid    = (*testInsPointer3)->GetFaceID();
                    gface2bcface[bfid].push_back(elid);
                    faceids[f]  = bfid;
                    fo++;
                }
                else
                {
                    std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer = m_FaceSetPointer.insert(facePointer);
                    int fref   = -2;
                    (*testInsPointer.first)->SetFaceRef(fref);
                    (*testInsPointer.first)->SetFaceID(fint);
                    faceids[f] = fint;
                    fint++;
                    
                }
                fintid++;
            }
        }
        //std::cout << "elid " << elid << std::endl;
        element2face_map[elid] = faceids;
        element2face.push_back(faceids);
        nlocs_faces.push_back(faceids.size());
        offsets_faces.push_back(offset_fcount);
        offset_fcount = offset_fcount + faceids.size();
    }

    std::map<int,std::vector<int> >::iterator itr;
    std::map<int,int> bcf_o_vs_n;
    for(itr=element2face_map.begin();itr!=element2face_map.end();itr++)
    {
        for(int q=0;q<itr->second.size();q++)
        {
            int fid = itr->second[q];
            if(itr->second[q]<0)
            {
                int fnewid     = fint-fid;
                itr->second[q] = fnewid;
            }
        }
    }


    // update faceid for the boundary faces
    FaceSetPointer::iterator ftit;
    for(ftit=allbcFaces.begin();ftit!=allbcFaces.end();ftit++)
    {
        int ofid = (*ftit)->GetFaceID();
        int nfid = ofid+fint;
        (*ftit)->SetFaceID(nfid);
        
        m_FaceSetPointer.insert((*ftit));
    }


    //std::cout << "m_FaceSetPointer " << m_FaceSetPointer.size() << " " << element2face_map.size() << std::endl;
    //std::cout << "SIZE MARH " << nlocs_faces.size() << " " << elements_map.size() << " " << fo << std::endl;
    
    /**/
    // std::cout << "num_elem = " << elements.size() << " num_nodes = " << nodes.size() << " " << ndim << " " << nelem << " " << npoin << std::endl;
    // std::cout << "bcFaces " << bcFaces.size() << " " << bcTags.size() << std::endl;
    // 

    // std::map<int,std::vector<std::vector<int> > >::iterator itmb;

    // for(itmb=bcFaces.begin();itmb!=bcFaces.end();itmb++)
    // {
    //     std::cout << itmb->first << " " << itmb->second.size() << std::endl;
    //     for(int q=0;q<itmb->second.size();q++)
    //     {
    //         for(int p=0;p<itmb->second[q].size();p++)
    //         {
    //             std::cout << itmb->second[q][p] << " ";
    //         }
    //         std::cout << std::endl;
            
    //     }
    //     std::cout << std::endl;
    // }


    // std::map<int,char*>::iterator itmbb;

    // std::cout << "bcTags = " << std::endl;
    // for(itmbb=bcTags.begin();itmbb!=bcTags.end();itmbb++)
    // {
    //     std::cout << "bc tags " << itmbb->first << " " << itmbb->second << std::endl;
    // }
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
    
    const char* fm = "hemihyb.su2";
    
    std::vector<std::vector<int> > elements;
    std::vector<std::vector<int> > element2face;
    std::vector<int> offsets;
    std::vector<int> offsets_vrts;
    std::vector<int> offsets_faces;
    std::vector<int> nlocs;
    std::vector<int> nlocs_vrts;
    std::vector<int> nlocs_faces_vec;
    std::map<int,std::vector<double> > nodes;
    std::map<int,std::vector<std::vector<int> > > bcFaces;
    std::map<int,char*> bcTags;
    std::vector<int> elements_root_flatten;
    std::vector<int> element2face_root_flatten;
    std::map<int,int> element_type;
    int nelem = 0;
    int nvrts = 0;
    double start1, end1,time_taken1;
    std::map<int, std::vector<int> > element2face_map;
    start1 = clock();
    FaceSetPointer allbcFaces;
    std::map<int,std::vector<int> > bcfaces_map;
    if(world_rank == 0)
    {
        std::cout << "Starting to read the file on root..." << std::endl;

        ReadSU2Mesh(comm, fm, 
                    elements, 
                    element2face, 
                    element2face_map, 
                    bcfaces_map,
                    element_type,
                    offsets, 
                    offsets_vrts, 
                    offsets_faces, 
                    nlocs, 
                    nlocs_vrts, 
                    nlocs_faces_vec,
                    nodes, bcFaces, bcTags);

        for (const auto& row : elements) 
        {
                // std::cout << "row " << row.size() << std::endl;
                elements_root_flatten.insert(elements_root_flatten.end(), row.begin(), row.end());
        }
        nelem = elements.size();
        elements.clear();
        for (const auto& row : element2face) 
        {
                element2face_root_flatten.insert(element2face_root_flatten.end(), row.begin(), row.end());
        }
        element2face.clear();
        
        nvrts = nodes.size();

    }
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

    ParallelState* pstate_elem = new ParallelState(nelem,comm);
    ParallelState* pstate_vert = new ParallelState(nvrts,comm);

    std::vector<int> sendcounts(world_size,0);
    std::vector<int> displs(world_size,0);

    std::vector<int> sendcounts_vrts(world_size,0);
    std::vector<int> displs_vrts(world_size,0);

    std::vector<int> sendcounts_faces(world_size,0);
    std::vector<int> displs_faces(world_size,0);

    // std::cout << nlocs_faces_vec.size() << " WR " << std::endl;
    if(world_rank == 0)
    {
        int end,start,index;
        int end_vrts,start_vrts,start_faces,end_faces;
        for(int i=0;i<world_size;i++)
        {
            int offset          = pstate_elem->getOffsets()[i];
            int nloc            = pstate_elem->getNlocs()[i];
            start               = offsets[offset];
            end                 = offsets[offset+nloc-1]+nlocs[offset+nloc-1];
            sendcounts[i]       = end-start;
            displs[i]           = start;

            int offset_vrts     = pstate_elem->getOffsets()[i];
            int nloc_vrts       = pstate_elem->getNlocs()[i];
            start_vrts          = offsets_vrts[offset];
            end_vrts            = offsets_vrts[offset+nloc-1]+nlocs_vrts[offset_vrts+nloc_vrts-1];
            // std::cout << world_rank << " offsets_vrts  "<< offsets_vrts[offset] << " " << offsets_vrts[offset+nloc+1] << std::endl;
            sendcounts_vrts[i]  = end_vrts-start_vrts;
            displs_vrts[i]      = start_vrts;

            int offset_faces    = pstate_elem->getOffsets()[i];
            int nloc_faces      = pstate_elem->getNlocs()[i];
            start_faces         = offsets_faces[offset];
            end_faces           = offsets_faces[offset+nloc-1]+nlocs_faces_vec[offset_faces+nloc_faces-1];
            //std::cout << "offsets_faces  "<< offsets_faces[offset] << " " << offsets_faces[offset+nloc] << " " << start_faces << " " << end_faces << std::endl;
            sendcounts_faces[i]  = end_faces-start_faces;
            displs_faces[i]      = start_faces;
        }
    }
    std::vector<int> sendcounts_red(world_size,0);
    MPI_Allreduce(sendcounts.data(), 
                  sendcounts_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);

    std::vector<int> sendcounts_faces_red(world_size,0);
    MPI_Allreduce(sendcounts_faces.data(), 
                  sendcounts_faces_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);

    std::vector<int> displs_red(world_size,0);
    MPI_Allreduce(displs.data(), 
                  displs_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);

    std::vector<int> displs_faces_red(world_size,0);
    MPI_Allreduce(displs_faces.data(), 
                  displs_faces_red.data(), 
                  world_size, MPI_INT, MPI_SUM, comm);
    
    std::vector<int> elements_on_rank(sendcounts_red[world_rank],0);
    std::vector<int> element2faces_on_rank(sendcounts_faces_red[world_rank],0);

    if(world_rank == 0)
    {
        std::cout << "Scattering the mesh to other processors..." << std::endl;
    }
    start1 = clock();


    MPI_Scatterv(elements_root_flatten.data(), 
                sendcounts_red.data(), displs_red.data(), MPI_INT, 
                elements_on_rank.data(), sendcounts_red[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(element2face_root_flatten.data(), 
                sendcounts_faces_red.data(), displs_faces_red.data(), MPI_INT, 
                element2faces_on_rank.data(), sendcounts_faces_red[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);
    
    std::vector<int>e2v_loc(pstate_elem->getNlocs()[world_rank],0);
    std::vector<int>e2f_loc(pstate_elem->getNlocs()[world_rank],0);

    MPI_Scatterv(nlocs.data(), 
                 pstate_elem->getNlocs(), pstate_elem->getOffsets(), MPI_INT, 
                 e2v_loc.data(), pstate_elem->getNlocs()[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(nlocs_faces_vec.data(), 
                 pstate_elem->getNlocs(), pstate_elem->getOffsets(), MPI_INT, 
                 e2f_loc.data(), pstate_elem->getNlocs()[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);

    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    if(world_rank == 0)
    {
        cout << setprecision(16) << "Time to taken scatter the mesh from root to other procs :          " << fixed
        << time_taken1 << setprecision(16);
        cout << " sec " << endl;
    }

    
    elements_root_flatten.resize(0);

    std::map<int,int> eltype_map;
    std::map<int,std::vector<int> > e2v;
    std::map<int,std::vector<int> > e2f;

    // FaceSharedPtr Face2RefPointer = std::shared_ptr<NekFace>(new NekFace(face));
    
    int off = 0;
    for(int i = 0;i<e2v_loc.size();i++)
    {
        int nloc   = e2v_loc[i];
        int elid   = elements_on_rank[off+nloc-1];
        int eltype = elements_on_rank[off];
        eltype_map[elid] = eltype;

        std::vector<int>e2v_row(nloc-2,0);
        int c = 0;
        for(int j=1;j<nloc-1;j++)
        {
            e2v_row[c] = elements_on_rank[off+j];
            c++;
        }

        e2v[elid] = e2v_row;


        int nlocf = e2f_loc[i];
        std::vector<int>e2f_row(nlocf,0);
        for(int j=0;j<nlocf;j++)
        {
            e2f_row[j] = element2faces_on_rank[off+j];
        }
        e2f[elid] = e2f_row;
        off = off+nloc;

    }



    /*
    std::map<int,std::vector<int> >::iterator ite;
    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    std::vector<std::vector<int> > e2f_map;
    int fcnt = 0;
    FaceSetPointer m_FaceSetPointer;
    FaceSetPointer m_InteriorFaceSetPointer;
    int ishere = 0;
    int fintid = 0;

    for(ite=e2v.begin();ite!=e2v.end();ite++)
    {
        int elid                = ite->first;
        int eltype              = eltype_map[elid];
        std::vector<int> e2n_i  = ite->second;
        // if(world_rank == 1)
        // {
        //     std::cout << "e2v " << e2n_i.size() << std::endl;
        // }
        switch (eltype) {
        case 10:
            e2f_map = tetra_faces;
            break;
        case 12:
            e2f_map = hex_faces;
            break;
        case 13:
            e2f_map = prism_faces;
            break;
        case 14:    
            e2f_map = pyramid_faces;
            break;
        }
        
        int nfaces = e2f_map.size();
        // if(world_rank == 1)
        // {std::cout << "nfaces " << nfaces << " "  << eltype << std::endl;}
        
        for(int f=0;f<nfaces;f++)
        {
            fcnt++;
            int nnodes = e2f_map[f].size();
            std::vector<int> face_globvid(nnodes,0);

            for(int n=0;n<nnodes;n++)
            {
                face_globvid[n] = e2n_i[e2f_map[f][n]];
                // if(world_rank == 1)
                // {
                //     std::cout << "e2v " << e2n_i[e2f_map[f][n]] << std::endl;
                // }
            }
            
            FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face_globvid));
            std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_FaceSetPointer.insert(facePointer);

            if(allbcFaces.find(facePointer)==allbcFaces.end())
            {
                ishere++;
                std::pair<FaceSetPointer::iterator, bool> testInsPointer2;
                testInsPointer2 = m_InteriorFaceSetPointer.insert(facePointer);
                if(testInsPointer2.second)
                {
                    (*testInsPointer2.first)->SetFaceID(fintid);
                    (*testInsPointer2.first)->SetFaceLeftElement(elid);
                    (*testInsPointer2.first)->SetFaceRightElement(-1);
                    fintid++;
                }
                else
                {
                    (*testInsPointer2.first)->SetFaceRightElement(elid);
                }
            }
        }
    }
    
    */
    //std::cout << "wr " << world_rank << " " << e2v.size() << " " << m_FaceSetPointer.size() << std::endl;

    MPI_Finalize();
        
}

