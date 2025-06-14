
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






void WriteCompleteSU2Mesh(MPI_Comm comm, const char* fm, 
                         std::map<int,std::vector<double> > nodes_to_write, 
                         std::map<int, std::vector<int> > elements_to_write, 
                         std::map<int,std::vector<std::vector<int> > > boundaryfaces_to_write, 
                         std::map<int,char*> bctags )
{
    std::ofstream mesh_file(fm);
    std::map<int,std::vector<int> >::iterator itmi;
    std::map<int,std::vector<std::vector<int> > >::iterator itmii;
    std::map<int,std::vector<double> >::iterator itmd;
    // 3D mesh declaration


    // % Problem dimension
    // NDIME= 3    
    // %
    // % Inner element connectivity
    // %

    mesh_file << "% Problem dimension"<<std::endl;
    mesh_file << "NDIME= 3"<<std::endl;
    mesh_file << "%"<<std::endl;
    mesh_file << "% Inner element connectivity"<<std::endl;
    mesh_file << "%"<<std::endl;

     // Write elements section
    mesh_file << "NELEM= " << elements_to_write.size() << "\n";

    for(itmi=elements_to_write.begin();itmi!=elements_to_write.end();itmi++)
    {
        mesh_file << 10 << std::setw(9) << itmi->second[0]-1 << std::setw(9) << itmi->second[1]-1 << std::setw(9) << itmi->second[2]-1 << std::setw(9) << itmi->second[3]-1 << std::setw(9) << itmi->first << "\n";
    }



    mesh_file << "NPOIN= " << nodes_to_write.size() << "\n";



    for(itmd=nodes_to_write.begin();itmd!=nodes_to_write.end();itmd++)
    {
        mesh_file << std::setw(24) << setprecision(16) << itmd->second[0] << std::setw(24) << itmd->second[1] << std::setw(24) << itmd->second[2] << std::setw(9) << itmd->first-1 << std::endl;
    }

    // for (const auto& node : nodes_to_write) {
    //     mesh_file << node[0] << " " << node[1] << " " << node[2] << "\n";
    // }
    
    // Create hexahedral elements (12 = 3D hex)
    
    // std::vector<std::vector<int>> elements = {
    //     {12, 0,1,2,3,4,5,6,7}  // Single hex element for demonstration
    // };
    
   
    
    // Define boundary markers (14 = 3D quad element)


    //mesh_file << "NMARK= 3\n";
    mesh_file << "NMARK= " << bctags.size() << "\n";

    std::map<int,char*>::iterator itchar;
    for(itchar=bctags.begin();itchar!=bctags.end();itchar++)
    {
        std::cout << itchar->first << " " << itchar->second << std::endl;
    }
    for(itmii=boundaryfaces_to_write.begin();itmii!=boundaryfaces_to_write.end();itmii++)
    {
        std::cout << itmii->first << std::endl;
    }

    for(itmii=boundaryfaces_to_write.begin();itmii!=boundaryfaces_to_write.end();itmii++)
    {
        std::cout << "itmii->first " << itmii->first << " " << bctags[itmii->first]  << std::endl; 
        mesh_file << "MARKER_TAG=" << bctags[itmii->first] << "\n";
        mesh_file << "MARKER_ELEMS=" << itmii->second.size() << "\n";

        for(int q=0;q<itmii->second.size();q++)
        {
            mesh_file << 5 << std::setw(17) << itmii->second[q][0]-1 << std::setw(17) << itmii->second[q][1]-1 << std::setw(17) << itmii->second[q][2]-1 << std::endl;
        }
    }
    
    // // Inlet (front and back faces)
    // mesh_file << "MARKER_TAG= inlet\n";
    // mesh_file << "MARKER_ELEMS= 2\n";
    // mesh_file << "14 0 1 5 4\n";  // Front face
    // mesh_file << "14 3 2 6 7\n";  // Back face
    
    // // Outlet (left and right faces)
    // mesh_file << "MARKER_TAG= outlet\n";
    // mesh_file << "MARKER_ELEMS= 2\n";
    // mesh_file << "14 0 3 7 4\n";  // Left face
    // mesh_file << "14 1 2 6 5\n";  // Right face
    
    // // Walls (top and bottom faces + ends)
    // mesh_file << "MARKER_TAG= walls\n";
    // mesh_file << "MARKER_ELEMS= 3\n";
    // mesh_file << "14 4 5 6 7\n";  // Top face
    // mesh_file << "14 0 1 2 3\n";  // Bottom face
    // mesh_file << "14 2 3 7 6\n";  // End face
    
    mesh_file.close();
}








void WriteTetraSU2Mesh(MPI_Comm comm, const char* fm, std::map<int,std::vector<double> > nodes_to_write, std::map<int, std::vector<int> > elements_to_write, std::map<int,std::vector<std::vector<int> > > boundaryfaces_to_write, std::map<int,char*> bctags )
{
    std::ofstream mesh_file(fm);
    std::map<int,std::vector<int> >::iterator itmi;
    std::map<int,std::vector<std::vector<int> > >::iterator itmii;
    std::map<int,std::vector<double> >::iterator itmd;
    // 3D mesh declaration


    // % Problem dimension
    // NDIME= 3    
    // %
    // % Inner element connectivity
    // %

    mesh_file << "% Problem dimension"<<std::endl;
    mesh_file << "NDIME= 3"<<std::endl;
    mesh_file << "%"<<std::endl;
    mesh_file << "% Inner element connectivity"<<std::endl;
    mesh_file << "%"<<std::endl;

     // Write elements section
    mesh_file << "NELEM= " << elements_to_write.size() << "\n";

    for(itmi=elements_to_write.begin();itmi!=elements_to_write.end();itmi++)
    {
        mesh_file << 10 << std::setw(9) << itmi->second[0]-1 << std::setw(9) << itmi->second[1]-1 << std::setw(9) << itmi->second[2]-1 << std::setw(9) << itmi->second[3]-1 << std::setw(9) << itmi->first << "\n";
    }



    mesh_file << "NPOIN= " << nodes_to_write.size() << "\n";



    for(itmd=nodes_to_write.begin();itmd!=nodes_to_write.end();itmd++)
    {
        mesh_file << std::setw(24) << setprecision(16) << itmd->second[0] << std::setw(24) << itmd->second[1] << std::setw(24) << itmd->second[2] << std::setw(9) << itmd->first-1 << std::endl;
    }

    // for (const auto& node : nodes_to_write) {
    //     mesh_file << node[0] << " " << node[1] << " " << node[2] << "\n";
    // }
    
    // Create hexahedral elements (12 = 3D hex)
    
    // std::vector<std::vector<int>> elements = {
    //     {12, 0,1,2,3,4,5,6,7}  // Single hex element for demonstration
    // };
    
   
    
    // Define boundary markers (14 = 3D quad element)


    //mesh_file << "NMARK= 3\n";
    mesh_file << "NMARK= " << bctags.size() << "\n";

    std::map<int,char*>::iterator itchar;
    for(itchar=bctags.begin();itchar!=bctags.end();itchar++)
    {
        std::cout << itchar->first << " " << itchar->second << std::endl;
    }
    for(itmii=boundaryfaces_to_write.begin();itmii!=boundaryfaces_to_write.end();itmii++)
    {
        std::cout << itmii->first << std::endl;
    }

    for(itmii=boundaryfaces_to_write.begin();itmii!=boundaryfaces_to_write.end();itmii++)
    {
        std::cout << "itmii->first " << itmii->first << " " << bctags[itmii->first]  << std::endl; 
        mesh_file << "MARKER_TAG=" << bctags[itmii->first] << "\n";
        mesh_file << "MARKER_ELEMS=" << itmii->second.size() << "\n";

        for(int q=0;q<itmii->second.size();q++)
        {
            mesh_file << 5 << std::setw(17) << itmii->second[q][0]-1 << std::setw(17) << itmii->second[q][1]-1 << std::setw(17) << itmii->second[q][2]-1 << std::endl;
        }
    }
    
    // // Inlet (front and back faces)
    // mesh_file << "MARKER_TAG= inlet\n";
    // mesh_file << "MARKER_ELEMS= 2\n";
    // mesh_file << "14 0 1 5 4\n";  // Front face
    // mesh_file << "14 3 2 6 7\n";  // Back face
    
    // // Outlet (left and right faces)
    // mesh_file << "MARKER_TAG= outlet\n";
    // mesh_file << "MARKER_ELEMS= 2\n";
    // mesh_file << "14 0 3 7 4\n";  // Left face
    // mesh_file << "14 1 2 6 5\n";  // Right face
    
    // // Walls (top and bottom faces + ends)
    // mesh_file << "MARKER_TAG= walls\n";
    // mesh_file << "MARKER_ELEMS= 3\n";
    // mesh_file << "14 4 5 6 7\n";  // Top face
    // mesh_file << "14 0 1 2 3\n";  // Bottom face
    // mesh_file << "14 2 3 7 6\n";  // End face
    
    mesh_file.close();
}







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
                FaceSetPointer &allbcFaces,
                std::vector<std::vector<int> > &bcfaces_root,
                std::vector<int> &bcref_root,
                std::vector<int> &bcfaceNvrt_root)
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

    std::vector<int> bcref;
    std::vector<int> bcfaceNvrt;
    std::vector<std::vector<int> > bcface;
    

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
                    int fref            = m+1;
                    //std::cout << fref << " fref bc " << std::endl;
                    (*testInsPointer.first)->SetFaceRef(fref);
                    //(*testInsPointer.first)->SetFaceID(-fid);
                    bcFaces[m+1].push_back(face);
                    int nvpf            = face.size();
                    nf                  = nf + face.size();
                    bcNFaces[m+1]       = nf;
                    bcFaceSizes[m+1].push_back(nvpf); 
                    // if(fref!=1 && fref!=2 && fref!=3)
                    // {
                    //     std::cout << "Wrong fref " << fref << " " << std::endl;
                    // }
                    bcref_root.push_back(fref);
                    bcfaces_root.push_back(face);
                    bcfaceNvrt_root.push_back(face.size());
                    //std::cout << "face.size() " << face.size() << std::endl;
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


    bcTags[10]="interface";

    std::cout << "bcfaces_root " << bcfaces_root.size() << std::endl;


    // std::cout << "bcFaces " << bcFaces.size() << std::endl; 
    // std::map<int,std::vector<std::vector<int> > >::iterator bcit;

    // for(bcit=bcFaces.begin();bcit!=bcFaces.end();bcit++)
    // {
    //     std::cout << bcit->first << std::endl;
    // }
    // std::vector<int> bcface_root_flatten;
    // for (const auto& row : bcface) 
    // {
    //         bcface_root_flatten.insert(bcface_root_flatten.end(), row.begin(), row.end());
    // }

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
        // if(bfref!=1 && bfref!=2 && bfref!=3)
        // {
        // }
        std::cout << "in map " <<  bfref << " " << bcNFaces[bfref] << " " << bcFacesizes_onRoot[bfref].size() << std::endl;
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

                    // if(fref!=2 && fref!=3 && fref!=4)
                    // {
                    //     std::cout << "still wrong " << std::endl;
                    // }
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


FaceSetPointer AllGatherBoundaryFaces(MPI_Comm comm, std::vector<int> bcfaces_root, std::vector<int> bcfaceNvrt_root, std::vector<int> bcfaceref_root)
{
    FaceSetPointer allbc;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int Nbcs = 0;
    int NbcsVrt = 0;

    Nbcs         = bcfaceref_root.size();
    NbcsVrt      = bcfaces_root.size();

    
    int nBCfaces = 0;
    MPI_Allreduce(&Nbcs,  &nBCfaces, 1, MPI_INT, MPI_SUM, comm);
    int nBCfacesVrts = 0;
    MPI_Allreduce(&NbcsVrt,  &nBCfacesVrts, 1, MPI_INT, MPI_SUM, comm);

    if(world_rank!=0)
    {
        bcfaces_root.resize(nBCfacesVrts,0);
        bcfaceNvrt_root.resize(nBCfaces,0);
        bcfaceref_root.resize(nBCfaces,0);   
    }
    

    std::vector<int> bcfaces_red(nBCfacesVrts,0);
    std::vector<int> bcfaceNvrt_red(nBCfaces,0);
    std::vector<int> bcfaceref_red(nBCfaces,0);
    
    MPI_Allreduce(bcfaces_root.data(), 
                  bcfaces_red.data(), 
                  bcfaces_root.size(), MPI_INT, MPI_SUM, comm);

    MPI_Allreduce(bcfaceref_root.data(), 
                  bcfaceref_red.data(), 
                  bcfaceref_root.size(), MPI_INT, MPI_SUM, comm);

    MPI_Allreduce(bcfaceNvrt_root.data(), 
                  bcfaceNvrt_red.data(), 
                  bcfaceNvrt_root.size(), MPI_INT, MPI_SUM, comm);

    int offset = 0;
    for(int i=0;i<nBCfaces;i++)
    {
        
        int bref = bcfaceref_red[i];

        // if(bref!=2 && bref!=3 && bref!=4)
        // {
        //     std::cout << "vbref " << bref << std::endl;
        // }
        int nvrt = bcfaceNvrt_red[i];
        std::vector<int> face(nvrt,0);
        for(int j=0;j<nvrt;j++)
        {
            face[j] = bcfaces_red[offset+j];
        }

        FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));

        std::pair<FaceSetPointer::iterator, bool> testInsPointer = allbc.insert(facePointer);
        if(testInsPointer.second) // new insertion
        {
            (*testInsPointer.first)->SetFaceRef(bref);
        }
        offset=offset+nvrt;
    }
    /**/
    return allbc;
}

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
            //std::cout << endV-startV << " " << i << " " << end-start << std::endl;
            
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
    
    std::vector<int> elements_on_rank(sendcounts_red[world_rank],0);
    std::vector<double> verts_on_rank(sendcountsV_red[world_rank],0);

    // if(world_rank == 0)
    // {
    //     std::cout << "Scattering the mesh to other processors..." << std::endl;
    // }

    
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
            //std::cout << vid << " compare ids " << gvid << " " << off+nloc-1 << std::endl;
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
    const char* fm = "mesh/hemi.su2";
    //const char* fm = "testbox2.su2";
    const char* fd = "mesh/hemi.h5";

    Inputs* inputs = ReadXmlFile(comm, "metric.xml");

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
    std::vector<int> bcfaces_root_flatten;
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

    std::vector<std::vector<int> > bcfaces_root;
    std::vector<int> bcref_root;
    std::vector<int> bcfaceNvrt_root;
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
                    nlocs_vrts, nodes, bcFaces_onRoot, bcNFaces_onRoot, bcFaceSizes_onRoot, bcTags, allbcFaces, bcfaces_root, bcref_root, bcfaceNvrt_root);

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
        //std::cout << "nelem " << nelem <<  " " << nvrts << std::endl;

        for (const auto& row : bcfaces_root) 
        {
                // std::cout << "row " << row.size() << std::endl;
                bcfaces_root_flatten.insert(bcfaces_root_flatten.end(), row.begin(), row.end());
        }
    }

    // Communicate all boundary face data to all other ranks:
    FaceSetPointer allbcFacesNew = AllGatherBoundaryFaces(comm, bcfaces_root_flatten, bcfaceNvrt_root, bcref_root);

    //std::cout << "allbcFacesNew.size() " << allbcFacesNew.size() << std::endl;

    // FaceSetPointer allbcFacesNew;
    // CommunicateBoundaryFaceData(comm, bcFaces_onRoot, bcNFaces_onRoot, bcFaceSizes_onRoot, allbcFacesNew);


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

    std::map<int, std::vector<double> > t_hessian = ReadHyperSolveHessianData(fd,nvrts,comm,info);
    
    // std::cout << "read hessian " << t_hessian.size() << " " << nvrts << std::endl;
    // std::cout << "t_hessian " << t_hessian.size() << " " << t_hessian[0].size()<< std::endl;
    
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
    
    PartObjectLite* partitionT = new PartObjectLite(tetUniMesh->e2v, tetUniMesh->vert_local, tetUniMesh->eltype_map, tetUniMesh->eltype_vec, allbcFacesNew, Ne, Nv, comm);
    PartObjectLite* partitionP = new PartObjectLite(blUniMesh->e2v,   blUniMesh->vert_local,  blUniMesh->eltype_map,  blUniMesh->eltype_vec, allbcFacesNew, Ne, Nv, comm);
    std::map<int,std::vector<double> > t_hessian_new = partitionT->RedistributeVertexDataForPartition(comm, nvrts, t_hessian);

    //std::cout << "Check redistributed data " << partitionT->getLocalVertsMap().size() << " " << t_hessian_new.size() << std::endl;
    FaceSetPointer::iterator ftit;
    FaceSetPointer sharedTfaces = partitionT->getAllSharedAndInterFaceFaceMap();
    //FaceSetPointer sharedPfaces = partitionP->getAllSharedAndInterFaceFaceMap();

    tetUniMesh->e2v.clear();
    tetUniMesh->vert_local.clear();
    tetUniMesh->eltype_map.clear();
    tetUniMesh->eltype_vec.clear();

    blUniMesh->e2v.clear();
    blUniMesh->vert_local.clear();
    blUniMesh->eltype_map.clear();
    blUniMesh->eltype_vec.clear();

    std::vector<int> OwnedSharedFaceOffsets     = partitionT->getownedSharedFacesOffsets();
    std::vector<int> OwnedSharedFaceNlocs       = partitionT->getownedSharedFacesNlocs();

    std::vector<int> OwnedInteriorFaceOffsets   = partitionT->getownedInteriorFacesOffsets();
    std::vector<int> OwnedInteriorFaceNlocs     = partitionT->getownedInteriorFacesNlocs();
    
    
    FaceSetPointer OwnedSharedTfaces            = partitionT->getOwnedSharedAndInterFaceFaceMap();
    FaceSetPointer OwnedInteriorTfaces          = partitionT->getOwnedInteriorFaceFaceMap();
    FaceSetPointer ActualSharedFaceSetPointerT  = partitionT->getActualSharedFaceMap();
    FaceSetPointer ExternalInterFaceFace        = partitionT->getExternalInterFaceFaceMap();

    // FaceSetPointer OwnedInteriorPfaces          = partitionP->getOwnedInteriorFaceFaceMap();
    // FaceSetPointer ActualSharedFaceSetPointerP  = partitionP->getActualSharedFaceMap();

    int nLocalInteriorT     = OwnedInteriorTfaces.size();
    int nAllInteriorT       = 0;
    MPI_Allreduce(&nLocalInteriorT,  &nAllInteriorT,   1, MPI_INT, MPI_SUM, comm);
    int nTotalInteriorFacesT = nAllInteriorT+ActualSharedFaceSetPointerT.size();

    // int nLocalInteriorP     = OwnedInteriorPfaces.size();
    // int nAllInteriorP       = 0;
    // MPI_Allreduce(&nLocalInteriorP,  &nAllInteriorP,   1, MPI_INT, MPI_SUM, comm);
    // int nTotalInteriorFacesP = nAllInteriorP+ActualSharedFaceSetPointerP.size();

    string filename_interface = "interface_" + std::to_string(world_rank) + ".vtu";
    OutputTriMeshPartitionVTK(comm,
                             filename_interface, 
                             ExternalInterFaceFace,
                             partitionT->getLocalVertsMap());
    




    //=========================================================================================================================
    //=========================================================================================================================
    //=========================================================================================================================
    //=========================================================================================================================

    FaceSetPointer m_OwnedBoundaryFaceSetPointer        = partitionT->getOwnedBoundaryFaceFaceMap();

    PMMG_pParMesh parmesh = InitializeParMMGmeshFromHyperSolveInputs(comm, 
                                                                     partitionT, 
                                                                     ActualSharedFaceSetPointerT, 
                                                                     InterFaceFaces, 
                                                                     m_OwnedBoundaryFaceSetPointer, 
                                                                     nTotalInteriorFacesT, 
                                                                     t_hessian_new,
                                                                     inputs);

    RunParMMGAndTestPartitioningFromHyperSolveInputs(comm, parmesh, partitionT, ActualSharedFaceSetPointerT, inputs);


    //=================================================================================================================
    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;


    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                           &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }

    DistributedParallelState* pstate_AdaptedTetraElem = new DistributedParallelState(nTetrahedraOUT,comm);

    std::map<int,std::vector<double> > LocalVerts_Adapted;
    std::map<int,std::vector<double> > LocalVerts_Adapted_GlobalID;

    int nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;
    std::map<int,int> glob2locVid;

    std::map<int,int> loc2globVid_complete;
    std::map<int,int> glob2locVid_complete;
    std::vector<int> globIDs;
    
    for(int k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k;

        loc2globVid_complete[k]=nodeGloNumber;
        glob2locVid_complete[nodeGloNumber]=k;
        globIDs.push_back(nodeGloNumber);
    }

    int max_gvid = *std::max_element(globIDs.begin(), globIDs.end());
    int max_vid = 0;
    MPI_Allreduce(&max_gvid,  &max_vid, 1, MPI_INT, MPI_MAX, comm);

    //std::cout << "max vid " << max_gvid << " " << max_vid << std::endl;

    std::vector<int> required_v(nVerticesOUT,0);
    std::vector<int> corner_v(nVerticesOUT,0);
    std::vector<int> refOUT_v(nVerticesOUT,0);

    std::map<int,int> vid2ref;

    for(int i=0;i<nVerticesOUT;i++)
    {
        std::vector<double> coords(3,0.0);
        
        if ( PMMG_Get_vertex(parmesh,&(coords[0]),&(coords[1]),&(coords[2]),
                               &(refOUT_v[i]),&(corner_v[i]),&(required_v[i])) != 1 )
        {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
        }
        // if(refOUT_v[i]!=0)
        // {
        //     std::cout << "refOUT_v[k] " << refOUT_v[i] << std::endl;
        // }
        LocalVerts_Adapted[i] = coords;
        LocalVerts_Adapted_GlobalID[loc2globVid[i+1]] = coords;

    }

    // std::map<int,std::vector<int> > LocalVerts_Adapted;
    std::map<int,int> local_v_map;
    std::map<int,std::vector<int> > gE2lV_Adapted;
    std::map<int,std::vector<int> > gE2gV_Adapted;
    std::set<int> OwnedElem_Adapted;
    std::vector<int> ref(nTetrahedraOUT,0);
    std::vector<int> required(nTetrahedraOUT,0);
    int locvid = 0;
    for(int i=0;i<nTetrahedraOUT;i++)
    {
        int offset = pstate_AdaptedTetraElem->getOffsets()[world_rank];
        int elid = offset + i;
        OwnedElem_Adapted.insert(elid);

        std::vector<int> LocalVertids(4,0);
        std::vector<int> GlobalVertids(4,0);
        if ( PMMG_Get_tetrahedron(parmesh, &(LocalVertids[0]),&(LocalVertids[1]), &(LocalVertids[2]),&(LocalVertids[3]), &(ref[i]),&(required[i])) != 1 ) 
        {
          fprintf(inm,"Unable to get mesh tetra %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }

        GlobalVertids[0] = loc2globVid[LocalVertids[0]];
        GlobalVertids[1] = loc2globVid[LocalVertids[1]];
        GlobalVertids[2] = loc2globVid[LocalVertids[2]];
        GlobalVertids[3] = loc2globVid[LocalVertids[3]];

        LocalVertids[0] = LocalVertids[0]-1;
        LocalVertids[1] = LocalVertids[1]-1;
        LocalVertids[2] = LocalVertids[2]-1;
        LocalVertids[3] = LocalVertids[3]-1;

        gE2lV_Adapted[elid] = LocalVertids;
        gE2gV_Adapted[elid] = GlobalVertids;
    }
    
    std::map<int,std::vector<double> > t_hessian_on_element;
    std::map<int,std::vector<int> >::iterator itm;
    std::map<int,std::vector<int> > gE2gV_t = partitionT->getElem2VertMap();
    
    for(itm=gE2gV_t.begin();itm!=gE2gV_t.end();itm++)
    {
        int elid = itm->first;
        int nvrt = itm->second.size();

        double h00 = 0.0;double h01 = 0.0;double h02 = 0.0;
        double h10 = 0.0;double h11 = 0.0;double h12 = 0.0;
        double h20 = 0.0;double h21 = 0.0;double h22 = 0.0;
        
        int t = 0;
        for(int q=0;q<nvrt;q++)
        {
            int vid = itm->second[q];
            h00 = h00 + t_hessian_new[vid][0];
            h01 = h01 + t_hessian_new[vid][1];
            h02 = h02 + t_hessian_new[vid][2];
            h10 = h10 + t_hessian_new[vid][3];
            h11 = h11 + t_hessian_new[vid][4];
            h12 = h12 + t_hessian_new[vid][5];
            h20 = h20 + t_hessian_new[vid][6];
            h21 = h21 + t_hessian_new[vid][7];
            h22 = h22 + t_hessian_new[vid][8];

        }

        std::vector<double> hess_elem(9,0.0);
        hess_elem[0] = h00/nvrt;
        hess_elem[1] = h01/nvrt;
        hess_elem[2] = h02/nvrt;
        hess_elem[3] = h10/nvrt;
        hess_elem[4] = h11/nvrt;
        hess_elem[5] = h12/nvrt;
        hess_elem[6] = h20/nvrt;
        hess_elem[7] = h21/nvrt;
        hess_elem[8] = h22/nvrt;


        t_hessian_on_element[elid] = hess_elem;

    }

    int bc0 = 0;
    int bc1 = 0;
    int bc2 = 0;
    int bc3 = 0;
    int bc4 = 0;
    int bc10 = 0;
    int nreq2   = 0;
    int *required2  = (int*)calloc(nTrianglesOUT,sizeof(int));
    // int *faceRefs   = (int*)calloc(nTrianglesOUT,sizeof(int));
    
    std::vector<int> freq(nTrianglesOUT,0);
    std::vector<int> fref(nTrianglesOUT,0);
    std::vector<int> fsize(nTrianglesOUT,0);
    std::vector<int> fvrt;
    // std::map<int,std::vector<int> > boundaryfaces_to_write;
    // std::map<int,std::vector<int> > boundaryfaces_to_write_size;
    int reqq = 0;
    for ( k=0; k<nTrianglesOUT; k++ )
    {
        std::vector<int> face(3,0);
        if ( PMMG_Get_triangle(parmesh,&(face[0]),
                                        &(face[1]),
                                        &(face[2]),
                                        &(fref[k]),
                                        &(required2[k])) != 1 )
        {
            fprintf(inm,"Unable to get mesh triangle %d \n",k);
            ier = PMMG_STRONGFAILURE;
        }

        if(fref[k]==0)
        {
            bc0++;
        }
        if(fref[k]==1)
        {
            bc1++;
        }
        if(fref[k]==2)
        {
            bc2++;
        }
        if(fref[k]==3)
        {
            bc3++;
        }
        if(fref[k]==4)
        {
            bc4++;
        }
        if(fref[k]==10)
        {
            bc10++;
        }

        if(required2[k]==1)
        {
             for(int b=0;b<face.size();b++)
            {
                if((face[b]-1)>nVerticesOUT)
                {
                    std::cout << "(" << face[b] << " -> ref = " << refOUT_v[face[b]-1]<< ") ";
                    std::cout << std::endl;
                }
            }
            
            reqq++;
        }

        std::vector<int> face_glob(3,0);
        face_glob[0] = loc2globVid[face[0]];
        face_glob[1] = loc2globVid[face[1]];
        face_glob[2] = loc2globVid[face[2]];

        
        fsize[k] = face_glob.size();
        for(int b=0;b<face_glob.size();b++)
        {
            fvrt.push_back(face_glob[b]);
        }
        
        if ( required2 && required2[k] )  nreq2++;
    }
    std::cout << "reqq " << world_rank << " " << reqq << std::endl; 
    
    // After this point, we can start figuring out what the boundary layer mesh vertex IDs are since we have now collected to vertices that on the interface
    std::map<int,std::vector<double> > LocalVerts_BL = partitionT->getLocalVertsMap();
    std::map<int,std::vector<double> >::iterator itmb;

    
    int vert_offset = max_vid;
    for(itmb=LocalVerts_BL.begin();itmb!=LocalVerts_BL.end();itmb++)
    {
        int tagid = itmb->first;

        if(interface_vrts.find(tagid)==interface_vrts.end())
        {

        }
    }







    int nreq2_total = 0;
    MPI_Allreduce(&nreq2,  &nreq2_total, 1, MPI_INT, MPI_SUM, comm);
    int bc0_total = 0;
    MPI_Allreduce(&bc0,  &bc0_total, 1, MPI_INT, MPI_SUM, comm);
    int bc1_total = 0;
    MPI_Allreduce(&bc1,  &bc1_total, 1, MPI_INT, MPI_SUM, comm);
    int bc2_total = 0;
    MPI_Allreduce(&bc2,  &bc2_total, 1, MPI_INT, MPI_SUM, comm);
    int bc3_total = 0;
    MPI_Allreduce(&bc3,  &bc3_total, 1, MPI_INT, MPI_SUM, comm);
    int bc4_total = 0;
    MPI_Allreduce(&bc4,  &bc4_total, 1, MPI_INT, MPI_SUM, comm);
    int bc10_total = 0;
    MPI_Allreduce(&bc10,  &bc10_total, 1, MPI_INT, MPI_SUM, comm);
    std::cout << "world_rank  " << nreq2_total << " on rank " << world_rank << " bc0 " << bc0 << " " << bc0_total << " " << bc1_total << " " << bc2_total << " " << bc3_total << " " << bc4_total << " bc10 " << bc10_total << std::endl;

    std::map<int,std::vector<int> > boundaryfaces_to_write_copy;
    std::map<int,std::vector<int> > boundaryfaces_to_write_size_copy;

    std::map<int,std::vector<int> >::iterator itmm;

    DistributedParallelState* bfaces_sizes = new DistributedParallelState(fsize.size(),comm);
    DistributedParallelState* bfaces_verts = new DistributedParallelState(fvrt.size(),comm);

    std::vector<int> BfacesSizeOnRoot;
    std::vector<int> BfacesVrtsOnRoot;
    std::vector<int> BfacesRefOnRoot;

    if(world_rank == 0)
    {
        std::cout << "bfaces_sizes->getNel() " << bfaces_sizes->getNel() << std::endl;
        BfacesSizeOnRoot.resize(bfaces_sizes->getNel());
        BfacesVrtsOnRoot.resize(bfaces_verts->getNel());
        BfacesRefOnRoot.resize(bfaces_sizes->getNel());
    }

    MPI_Gatherv(&fsize.data()[0],
                fsize.size(),
                MPI_INT,
                &BfacesSizeOnRoot.data()[0],
                bfaces_sizes->getNlocs(),
                bfaces_sizes->getOffsets(),
                MPI_INT, 0, comm);

    MPI_Gatherv(&fvrt.data()[0],
                fvrt.size(),
                MPI_INT,
                &BfacesVrtsOnRoot.data()[0],
                bfaces_verts->getNlocs(),
                bfaces_verts->getOffsets(),
                MPI_INT, 0, comm);
    
    MPI_Gatherv(&fref.data()[0],
                fref.size(),
                MPI_INT,
                &BfacesRefOnRoot.data()[0],
                bfaces_sizes->getNlocs(),
                bfaces_sizes->getOffsets(),
                MPI_INT, 0, comm);



    std::map<int,std::vector<std::vector<int> > > boundaryMap_global;
    FaceSetPointer bfacesOut;

    std::map<int,std::vector<std::vector<int> > > bcFacesPerRefTetra;

    if(world_rank == 0)
    {
        int offset = 0;
        

        for(int i=0;i<BfacesSizeOnRoot.size();i++)
        {
            int fz      = BfacesSizeOnRoot[i];
            int fref    = BfacesRefOnRoot[i];

            std::vector<int> face(fz,0);

            for(int j=0;j<fz;j++)
            {
                face[j] = BfacesVrtsOnRoot[offset+j];
            }

            FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));

            FaceSetPointer::iterator testInsPointer2 = bfacesOut.find(facePointer);
            if(testInsPointer2==bfacesOut.end())
            {
                std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer      = bfacesOut.insert(facePointer);

                if(fref!=0)
                {
                    bcFacesPerRefTetra[fref].push_back(face);
                }
                
            }



            offset = offset+fz;
        }

    }

    //std::cout << "bcFacesPerRef. "<< bcFacesPerRef.size() << std::endl;
    //============================================================================================
    /**/

    // std::map<int,std::string > varnamesGrad;

    // varnamesGrad[0]     = "dU2dx2";
    // varnamesGrad[1]     = "dU2dxdy";
    // varnamesGrad[2]     = "dU2dxdz";
    // varnamesGrad[3]     = "dU2dydx";
    // varnamesGrad[4]     = "dU2dy2";
    // varnamesGrad[5]     = "dU2dydz";
    // varnamesGrad[6]     = "dU2dzdx";
    // varnamesGrad[7]     = "dU2dzdy";
    // varnamesGrad[8]     = "dU2dz2";

    // string fm_NOTadaptedTetra = "NOTadaptedTetra_" + std::to_string(world_rank) + ".vtu";
    
    // OutputTetraMeshPartitionVTK(comm,
    //                         fm_NOTadaptedTetra, 
    //                         partitionT->getLocalElemSet(), 
    //                         gE2gV_t, 
    //                         t_hessian_on_element, 
    //                         varnamesGrad, 
    //                         partitionT->getLocalVertsMap());


    // string fm_adaptedTetra = "adaptedTetra_" + std::to_string(world_rank) + ".vtu";
    // OutputTetraMeshNoSolutionPartitionVTK(comm,
    //                          fm_adaptedTetra, 
    //                          OwnedElem_Adapted,
    //                          gE2lV_Adapted,
    //                          LocalVerts_Adapted);

    // string fm_adaptedRootTetra = "adaptedTetraRoot_" + std::to_string(world_rank) + ".vtu";

    // OutputTetraMeshNoSolutionOnRootVTK(comm,
    //                          fm_adaptedRootTetra, 
    //                          OwnedElem_Adapted,
    //                          gE2gV_Adapted,
    //                          LocalVerts_Adapted_GlobalID);
    
    // string fm_NOTadaptedTetraOnRoot = "NOTadaptedTetraRoot_" + std::to_string(world_rank) + ".vtu";
    // OutputTetraMeshOnRootVTK(comm,
    //                         fm_NOTadaptedTetraOnRoot, 
    //                         partitionT->getLocalElemSet(), 
    //                         gE2gV_t, 
    //                         t_hessian_on_element, 
    //                         varnamesGrad, 
    //                         partitionT->getLocalVertsMap());

    
    // std::map<int,std::vector<int> > gE2lVOnRoot     = GatherGlobalMapOnRoot_T(gE2lV_tmp,comm);
    // std::map<int,std::vector<double> > glob_data    = GatherGlobalMapOnRoot_T(loc_data,comm);

    std::map<int,std::vector<double> > VertsOnRootTetra                 = GatherGlobalMapOnRoot_T(LocalVerts_Adapted_GlobalID,comm);
    std::map<int,std::vector<int> > TetraOnRoot                         = GatherGlobalMapOnRoot_T(gE2gV_Adapted,comm);

    std::map<int,std::vector<int> > boundaryfaces_to_write_global       = GatherGlobalMapOnRoot_T(boundaryfaces_to_write_copy,comm);
    std::map<int,std::vector<int> > boundaryfaces_to_write_size_global  = GatherGlobalMapOnRoot_T(boundaryfaces_to_write_size_copy,comm);


    std::map<int,std::vector<std::vector<int> > > BoundariesOnRoot_global;
    // std::map<int,std::vector<std::vector<int> > >::iterator itmii;
    // for(itmii=boundaryfaces_to_write.begin();itmii!=boundaryfaces_to_write.end();itmii++)
    // {
    //     std::vector<std::vector<int> > BoundariesOnRoot_i    = GatherJaggedGlobalVectorOnRoot_T(itmii->second,comm);

    //     BoundariesOnRoot_global[itmii->first] = BoundariesOnRoot_i;
    // }

    

    if(world_rank == 0)
    {
        std::cout << "VertsOnRoot " << VertsOnRootTetra.size() << " " << boundaryfaces_to_write_global.size() << " " << boundaryfaces_to_write_size_global.size() << std::endl; 

        std::map<int,std::vector<int> >::iterator tim;
        
        // for(tim=boundaryfaces_to_write_size_global.begin();tim!=boundaryfaces_to_write_size_global.end();tim++)
        // {
        //     int bref  = tim->first;
        //     int nface = tim->second.size();
        //     std::cout << "bref " << bref << " " << nface << " " << boundaryfaces_to_write_global[bref].size() << std::endl;
        //     int offset = 0;
        //     for(int q=0;q<nface;q++)
        //     {
        //         int nvrt = tim->second[q];
        //         std::vector<int> faceOnRoot(nvrt,0);
        //         for(int k=0;k<nvrt;k++)
        //         {
        //             faceOnRoot[k] = boundaryfaces_to_write_global[bref][offset+k];
        //         }
                
        //         boundaryMap_global[bref].push_back(faceOnRoot);

        //         offset = offset + nvrt;
        //     }
        // }
    }
    
    
    if(world_rank == 0)
    {
        const char* fdo = "mimic_hemi.su2";
        
        WriteTetraSU2Mesh(comm, fdo, VertsOnRootTetra, TetraOnRoot, bcFacesPerRefTetra, bcTags);
    }

    // PMMG_pParMesh InitializeParMMGmeshFromHyperSolveInputs(MPI_Comm comm, 
    //                                                     PartObjectLite* partition,
    //                                                     std::vector<int> face4parmmg,
    //                                                     FaceSetPointer allbcFaces,
    //                                                     std::vector<std::vector<double> > t_metric)


    // std::cout << "All interior faces T " << nAllInteriorT+ActualSharedFaceSetPointerT.size() << " " << nAllInteriorT << " " << ActualSharedFaceSetPointerT.size() << std::endl;
    // std::cout << "All interior faces P " << nAllInteriorP+ActualSharedFaceSetPointerP.size()<< " " << nAllInteriorP << " " << ActualSharedFaceSetPointerP.size() << std::endl;
    /*
    int m_ncomm             = m_ColorsFaces.size();
    int** m_ifc_tria_glo;
    int** m_ifc_tria_loc;
    int* m_color_face;
    int* m_ntifc;
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
            m_ifc_tria_loc[icomm][q] = m_globShF2locShF[itc->second[q]]+1;

            //std::cout << m_globShF2locShF[itc->second[q]]+1 <<  " " << itc->second[q]+1 << std::endl;
        }

        icomm++;
    }

    //std::cout << "m_ColorsFaces " << world_rank << " " << m_ColorsFaces.size() << std::endl;
    //std::cout << "fnd " << fnd << " On Rank " << world_rank << std::endl;

    int nTotalSharedFaces = OwnedSharedFaceOffsets[world_size-1]+OwnedSharedFaceNlocs[world_size-1];
    
    // for(ftit=InterFaceFaces.begin();ftit!=InterFaceFaces.end();ftit++)
    // {
    //     (*ftit)->SetFaceID(nTotalInteriorFaces+nTotalSharedFaces+fid);
    //     fid++;
    // }

*/

    MPI_Finalize();
        
}

