
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
                std::vector<int> &offsets,
                std::vector<int> &nlocs,
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
    // std::cout << "ndim = " << ndim << " nelem " << nelem << " lstart = " << lstart << std::endl; 
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
            //std::cout << std::endl;
            offsets.push_back(offset_count);
            nlocs.push_back(row.size());
            elements.push_back(row);
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
            //std::cout << "nid " << nid << " :: ";
            while (iss >> num) {
                row.push_back(num);
                //std::cout << std::setprecision(16) << num << " ";
            }
            //std::cout << std::endl;

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
                    //std::cout << " " << m << " tag " << tag  << " end tag "<< std::endl;
                    
                    
            }

            else if (line.find("MARKER_ELEMS=") != std::string::npos)
            {
                    mark_elem = 0;
                    std::istringstream iss(line);
                    std::string key;
                    iss >> key >> mark_elem;
                    
                    // std::cout << "bc " << m << " mark_elem " << mark_elem << std::endl;
                    
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

                bcFaces[m].push_back(row);
                nid = nid + 1;


                if((m == (nmark)) && (cnt == mark_elem-1))
                {
                    processBoundaryFaces = false;
                }
                //std::cout << "cnt = " << cnt << " " << nmark << " " << mark_elem-1 << std::endl;
                cnt = cnt + 1;
            }

            // std::istringstream iss(line);
            // std::vector<double> row;
            // double num;
            // while (iss >> num) {
            //     row.push_back(num);
            // }

            // nodes[nid] = row;
            // nid = nid + 1;
        }
    }


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
    std::vector<int> offsets;
    std::vector<int> nlocs;
    std::map<int,std::vector<double> > nodes;
    std::map<int,std::vector<std::vector<int> > > bcFaces;
    std::map<int,char*> bcTags;
    std::vector<int> elements_root_flatten;
    int nelem = 0;
    double start1, end1,time_taken1;
    start1 = clock();
    if(world_rank == 0)
    {
        std::cout << "Starting to read the file on root..." << std::endl;

        ReadSU2Mesh(comm, fm, elements, offsets, nlocs, nodes, bcFaces, bcTags);

        for (const auto& row : elements) 
        {
                elements_root_flatten.insert(elements_root_flatten.end(), row.begin(), row.end());
        }
        nelem = elements.size();

        // std::cout << "Done reading the file on root..." << std::endl;
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
    
    ParallelState* pstate = new ParallelState(nelem,comm);

    std::vector<int> sendcounts(world_size,0);
    std::vector<int> displs(world_size,0);

    if(world_rank == 0)
    {
        int end,start,index;

        for(int i=0;i<world_size;i++)
        {
            int offset      = pstate->getOffsets()[i];
            int nloc        = pstate->getNlocs()[i];
            start           = offsets[offset];
            end             = offsets[offset+nloc-1]+nlocs[offset+nloc-1];
            sendcounts[i]   = end-start;
            displs[i]       = start;
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
    
    std::vector<int> elements_on_rank(sendcounts_red[world_rank],0);
    if(world_rank == 0)
    {
        std::cout << "Scattering the mesh to other processors..." << std::endl;
    }
    start1 = clock();
    MPI_Scatterv(elements_root_flatten.data(), 
                sendcounts_red.data(), displs_red.data(), MPI_INT, 
                elements_on_rank.data(), sendcounts_red[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);

    std::vector<int>e2v_loc(pstate->getNlocs()[world_rank],0);
    MPI_Scatterv(nlocs.data(), 
                pstate->getNlocs(), pstate->getOffsets(), MPI_INT, 
                e2v_loc.data(), pstate->getNlocs()[world_rank], MPI_INT,
                 0, MPI_COMM_WORLD);
    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    if(world_rank == 0)
    {
        cout << setprecision(16) << "Time to taken scatter the mesh from root to other procs :          " << fixed
        << time_taken1 << setprecision(16);
        cout << " sec " << endl;
    }
    // if(world_rank == 0)
    // {
    //     std::cout << "Done Scattering the mesh to other processors..." << std::endl;
    // }
    elements_root_flatten.resize(0);

    // if(world_rank == 1)
    // {
    //     int off = 0;
    //     for(int i = 0;i<e2v_loc.size();i++)
    //     {
    //         int nloc = e2v_loc[i];

    //         for(int j=0;j<nloc;j++)
    //         {
    //             std::cout << elements_on_rank[off+j] <<" ";
    //         }
    //         std::cout << std::endl;

    //         off = off+nloc;
    //     }
    // }


    MPI_Finalize();
        
}

