
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

//#include "/Users/dekelsch/mimic_libmesh/utilities/partitionTetrahedra/build/ThirdParty/dist/include/libmeshb7.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;

    // std::vector<std::vector<int> > e2f_hex(6);
    // std::vector<int> f0hex(4,0);
    // f0hex[0] = 0;
    // f0hex[1] = 3;
    // f0hex[2] = 2;
    // f0hex[3] = 1;
    // e2f_hex[0] = f0hex;
    // std::vector<int> f1hex(4,0);
    // f1hex[0] = 0;
    // f1hex[1] = 4;
    // f1hex[2] = 7;
    // f1hex[3] = 3;
    // e2f_hex[1] = f1hex;
    // std::vector<int> f2hex(4,0);
    // f2hex[0] = 0;
    // f2hex[1] = 1;
    // f2hex[2] = 5;
    // f2hex[3] = 4;
    // e2f_hex[2] = f2hex;
    // std::vector<int> f3hex(4,0);
    // f3hex[0] = 4;
    // f3hex[1] = 5;
    // f3hex[2] = 6;
    // f3hex[3] = 7;
    // e2f_hex[3] = f3hex;
    // std::vector<int> f4hex(4,0);
    // f4hex[0] = 1;
    // f4hex[1] = 2;
    // f4hex[2] = 6;
    // f4hex[3] = 5;
    // e2f_hex[4] = f4hex;
    // std::vector<int> f5hex(4,0);
    // f5hex[0] = 2;
    // f5hex[1] = 3;
    // f5hex[2] = 7;
    // f5hex[3] = 6;
    // e2f_hex[5] = f5hex;

    // std::vector<std::vector<int> > e2f_tet(4);
    // std::vector<int> f0(3,0);
    // f0[0] = 0;
    // f0[1] = 2;
    // f0[2] = 1;
    // e2f_tet[0] = f0;
    // std::vector<int> f1(3,0);
    // f1[0] = 0;
    // f1[1] = 1;
    // f1[2] = 3;
    // e2f_tet[0] = f1;
    // std::vector<int> f2(3,0);
    // f1[0] = 1;
    // f1[1] = 2;
    // f1[2] = 3;
    // e2f_tet[0] = f2;
    // std::vector<int> f3(3,0);
    // f1[0] = 0;
    // f1[1] = 3;
    // f1[2] = 2;
    // e2f_tet[0] = f3;

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
struct ElementNew {
    int id;
    std::vector<int> vertices;
};


void ReadSU2Mesh(MPI_Comm comm, const char* fm,    
                std::map<int,std::vector<int> > &elements, 
                std::map<int,int> &element_type,
                std::map<int,std::vector<double> > &nodes, 
                std::map<int,FaceSetPointer> &bcid2Face_map,
                FaceSetPointer &bcFaces,
                std::map<int,char*> &bcTags)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::ifstream meshFile;
    meshFile.open(fm);
    std::string   line;
    if(!meshFile.is_open())
    {
        std::cout << "Error:: Make sure there is a mesh.su2 file in the directory where main.cpp resides. "<<std::endl;
        exit(0);
    }

    bool processHeaderElements = true;
    int lstart  = 0;
    int ndim    = 0; 
    int nelem   = 0;
    int npoin   = 0;
    int c       = 0;

    int cols    = 6;
    //================================PROCESS ELEMENTS =====================================================

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

    if(rank == 0)
    {
        std::cout << "ndim = " << ndim << " nelem " << nelem << " lstart = " << lstart << std::endl; 
    }

    // pstate->getOffsets()

    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, fm, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    int nloc                = int(nelem/size) + ( rank < nelem%size );
    //  compute offset of rows for each proc;
    int offset              = rank*int(nelem/size) + MIN(rank, nelem%size);
    MPI_Offset start        = lstart + offset;
    MPI_Offset end          = lstart + offset+nloc;
    int local_rows          = end - start;
    int current_pos         = start;


    std::vector<std::vector<int> > local_data;
    const int MAX_LINE_LENGTH = 1024;
    char line_buffer[MAX_LINE_LENGTH];
    MPI_Status status;
    current_pos = start;
    int u = 0;
    meshFile.seekg(start);

    // std::string line;
    // std::vector<std::vector<int>> local_data;

    while (std::getline(meshFile, line)) {
        std::istringstream iss(line);
        std::vector<int> row;
        int value;
        while (iss >> value) {
            row.push_back(value);
        }
        local_data.push_back(row);
        if (meshFile.tellg() >= end) {
            break;
        }
    }

    std::cout << "local_data " << rank << " " << local_data.size() << std::endl;

    meshFile.close();
    // for (int i = 0; i < start_row; ++i) {
    //     int num_cols;
    //     file >> num_cols;
    //     file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip row data
    // }


    // while (u < nloc) {
    //     int line_length = 0;
    //     // Read a line
    //     MPI_File_read_at(file, current_pos, line_buffer, MAX_LINE_LENGTH, MPI_CHAR, &status);
    //     MPI_Get_count(&status, MPI_CHAR, &line_length);

    //     // Find the end of the line
    //     int i;
    //     for (i = 0; i < line_length; i++) {
    //         if (line_buffer[i] == '\n') {
    //             break;
    //         }
    //     }
    //     line_length = i + 1;  // Include the newline

    //     // Parse the line
    //     std::string line(line_buffer, line_length - 1);  // Exclude the newline
    //     std::vector<int> row_data;
    //     std::istringstream iss(line);
    //     int value;
    //     while (iss >> value) {
    //         row_data.push_back(value);
    //         std::cout << value << " ";
    //     }
    //     std::cout << line_length << " " << row_data.size() << std::endl;

    //     //std::cout << " row_data " << row_data.size()<< " " << current_pos << " " << end << std::endl;

    //     if (!row_data.empty()) {
    //         local_data.push_back(row_data);
    //     }
        

    //     // Move to the next line
    //     current_pos += line_length;
    //     //std::cout << i << std::endl;
    //     u = u + 1;
    // }

    //std::cout << "elements_mpi " << buffer.size() << " " << local_rows << " " << end << " " << start << " "  << " " << nloc << std::endl; 
    // std::cout << local_data.size() << " " << nloc << std::endl;
    // if(rank==0)
    // {
    //     for(int q=0;q<local_data.size();q++)
    //     {
    //         for(int p=0;p<local_data[q].size();p++)
    //         {
    //             std::cout << local_data[q][p] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
//  MPI_File_close(&file);
//     MPI_Finalize();

}




void ReadSU2MeshSerial(MPI_Comm comm, const char* fm,    
                std::map<int,std::vector<int> > &elements, 
                std::map<int,int> &element_type,
                std::map<int,std::vector<double> > &nodes, 
                std::map<int,FaceSetPointer> &bcid2Face_map,
                FaceSetPointer &bcFaces,
                std::map<int,char*> &bcTags)
{

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // ParallelState* pstate = new ParallelState(nelem,comm);
    
    std::vector<int> scatter_offsets(size,0);
    std::vector<int> scatter_offsets_red(size,0);


    if(rank == 0)
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

        std::vector<std::vector<int> > elements_root;

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
        std::cout << "ndim = " << ndim << " nelem " << nelem << " lstart = " << lstart << " " << rank << std::endl; 

        std::cout << "Reading on root..." << std::endl;
        //================================PROCESS ELEMENTS =====================================================
        int linesToSkip = c + 1;
        bool processElements = false;
        int var = 0;
        
        // std::map<int,std::vector<int> > elements;
        int elid         = 0;
        int offset_track = 0;
        std::vector<int> offset_tracker(nelem,0);

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
                while (iss >> num) {
                    row.push_back(num);
                    
                }
                //std::cout << std::endl;
                element_type[elid]  = row[0];
                std::vector<int> destVect(row.size()-1,0);
                copy(row.begin() + 1, row.end(), destVect.begin());
                elements_root.push_back(destVect);
                offset_tracker[elid] = offset_track;
                offset_track = offset_track + destVect.size();
                elid = elid + 1;
            }
        }
        std::cout << "processing header nodes..." << std::endl;
        //================================PROCESS NODES =====================================================
        bool processHeaderNodes = true;
        int y = 0;
        while (processHeaderNodes) 
        {
            std::getline(meshFile, line);    
            std::stringstream ssline(line);
            std::cout << "line processHeaderNodes " << line << std::endl; 
            if (line.find("NPOIN=") != std::string::npos)
            {
                // std::cout << "when her " << y << std::endl;
                std::istringstream iss(line);
                std::string key;
                iss >> key >> npoin;
                lstart = c;
                processHeaderNodes = false;
            }

            std::cout << "processing " << y << std::endl;
            c++;
            y++;
            
            // std::cout << processHeaderNodes << std::endl;
        }

        std::cout << "processing nodes..." << y << std::endl;
        /*
        int nid  = 0;
        // std::map<int,std::vector<double> > nodes;
        bool processNodes = true;
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
                nodes[nid] = row;
                nid = nid + 1;
                
            }
        }

        std::cout << "processing boundary faces..." << std::endl;
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

        // std::map<int,std::vector<std::vector<int> > > bcFaces;
        // std::map<int,char*> bcTags;
        bool processBoundaryFaces = true;
        int m = 0;
        int cnt = 0;
        int mark_elem = 0;
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
                        m++;
                }
                else
                {
                    std::istringstream iss(line);
                    std::vector<int> row;
                    int num;
                    
                    while (iss >> num) {
                        row.push_back(num);
                    }
                    int ftype  = row[0];
                    std::vector<int> face(row.size()-1,0);
                    copy(row.begin() + 1, row.end(), face.begin());
                    FaceSharedPtr BCFacePointer = std::shared_ptr<NekFace>(new NekFace(face));
                    std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer = bcid2Face_map[m].insert(BCFacePointer);

                    testInsPointer = bcFaces.insert(BCFacePointer);
                    nid = nid + 1;


                    if((m == (nmark)) && (cnt == mark_elem-1))
                    {
                        processBoundaryFaces = false;
                    }
                    cnt = cnt + 1;
                }
            }
        }
        */
        // std::cout << "read " << std::endl;

        // std::vector<int> scatter_offsets(size,0);
        // for(int s=0;s<size;s++)
        // {
        //     int e_nloc          = pstate->getNlocs()[s];
        //     int e_offs          = pstate->getOffsets()[s];
        //     scatter_offsets[s]  = offset_tracker[e_offs+e_nloc];
        //     std::cout << "offsets " << e_nloc << " " << scatter_offsets[s] << std::endl;
        // }
        
        
        // std::vector<int> scatter_offsets(size,0);

        // scatter_offsets[rank] = offset;
        

        // MPI_Allreduce(scatter_offsets.data(), 
        //             scatter_offsets_red.data(), 
        //             size, MPI_INT, MPI_SUM, comm);

        // for(int i=0;i<size;i++)
        // {
        //     std::cout << i << " " << scatter_offsets_red[i] << std::endl;
        // }


    }


    
    
    
    // std::vector<int> elements_root_flatten;

    // for (const auto& row : elements_root) 
    // {
    //         elements_root_flatten.insert(elements_root_flatten.end(), row.begin(), row.end());
    // }

    // MPI_Scatterv(elements_root_flatten.data(), sendcounts.data(), displacements.data(), MPI_INT,
    //              recv_buffer.data(), recv_count, MPI_INT, root, MPI_COMM_WORLD);





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
    // const char* fm = "aa.su2";
    if(world_rank == 0)
    {
        std::cout << "loading " << fm << "..." << std::endl;
    }
   
    std::map<int,std::vector<int> > elements;
    std::map<int,int> element_type;
    std::map<int,std::vector<double> > nodes;
    std::map<int,FaceSetPointer> bcid2face_map;
    std::map<int,char*> bc_tags;
    FaceSetPointer allbcFaces;
    ReadSU2MeshSerial(comm,fm,elements, element_type, nodes, bcid2face_map, allbcFaces, bc_tags);
    // std::cout << "N_elem=   "<<elements.size() << std::endl; 
    // std::cout << "N_node=   "<<nodes.size() << std::endl; 
    // std::cout << "N_bcfa=   "<<bcid2face_map.size() << std::endl; 
    /*
    std::map<int,FaceSetPointer>::iterator itbc;
    for(itbc=bcid2face_map.begin();itbc!=bcid2face_map.end();itbc++)
    {
        int bcid = itbc->first;
        FaceSetPointer bcfaces = itbc->second;
        //std::cout << "bc id " << bcid << " " << bc_tags[bcid] << " " << bcfaces.size() << std::endl;
    }


    std::map<int,std::vector<int> >::iterator ite;

    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    std::vector<std::vector<int> > e2f_map;
    FaceSetPointer m_FaceSetPointer;
    FaceSetPointer m_InteriorFaceSetPointer;
    int fcnt    = 0;
    int fintid  = 0;
    int ishere  = 0;
    for(ite=elements.begin();ite!=elements.end();ite++)
    {
        int elid                = ite->first;
        int eltype              = element_type[elid];
        std::vector<int> e2n_i  = ite->second;

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

        for(int f=0;f<nfaces;f++)
        {
            fcnt++;
            int nnodes = e2f_map[f].size();
            std::vector<int> face_globvid(nnodes,0);

            for(int n=0;n<nnodes;n++)
            {
                face_globvid[n] = e2n_i[e2f_map[f][n]];
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
                    // std::cout << (*testInsPointer2.first)->GetFaceLeftElement() << "  "<< elid << std::endl;
                    (*testInsPointer2.first)->SetFaceRightElement(elid);
                }
            }
        }
    }

    std::cout << "N_face=   " << m_FaceSetPointer.size() 
              << " " << m_InteriorFaceSetPointer.size() 
              << " " << allbcFaces.size() << " " << m_FaceSetPointer.size()-m_InteriorFaceSetPointer.size() << " " << ishere << std::endl;
    */
    // MPI_Finalize();
        
}

