
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
//#include "/Users/dekelsch/mimic_libmesh/utilities/partitionTetrahedra/build/ThirdParty/dist/include/libmeshb7.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;




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



    // Open the mesh file for reading
    LibIdx = GmfOpenMesh( "hemi.meshb", GmfRead, &ver, &dim );


    // Get the number of vertices and triangles
    NmbVer = GmfStatKwd( LibIdx, GmfVertices  );
    NmbTri = GmfStatKwd( LibIdx, GmfTriangles );
    NmbTet = GmfStatKwd( LibIdx, GmfTetrahedra);
    NmbPri = GmfStatKwd( LibIdx, GmfPrisms);

    std::vector<std::vector<double> > Coords(NmbVer);
    std::vector<std::vector<int> > Nodes(NmbTri);
    std::vector<std::vector<int> > Nodes_Tet(NmbTet);
    std::vector<int> Domains(NmbVer,0.0);


    GmfGotoKwd(LibIdx, GmfVertices);

    for(i=0;i<NmbVer;i++)
    {
        std::vector<double> Coords_row(3,0.0);
        GmfGetLin( LibIdx, GmfVertices, &Coords_row[0], &Coords_row[1], &Coords_row[2], &Domains[i] );
        std::cout << Coords_row[0] << " " << Coords_row[1] << " " << Coords_row[2] << " " << Domains[i] << std::endl;
        Coords[i] = Coords_row;
    }

    // // Move the file pointer to the triangles keyword
    GmfGotoKwd( LibIdx, GmfTriangles );

    // // Read each line of triangle data into your own data structures
    for(i=0;i<NmbTri;i++)
    {
        std::vector<int> Nodes_row(4,0);
        GmfGetLin( LibIdx, GmfTriangles, &Nodes_row[0], &Nodes_row[1], &Nodes_row[2], &Nodes_row[3] );
        std::cout << " tri " << Nodes_row[0] << " " << Nodes_row[1] << " " << Nodes_row[2] << " " << Nodes_row[3] << std::endl;
        Nodes[i] = Nodes_row;
    }


    GmfGotoKwd( LibIdx, GmfTetrahedra );

     // // Read each line of triangle data into your own data structures
    for(i=0;i<NmbTet;i++)
    {
        std::vector<int> Nodes_Tet_row(4,0);
        GmfGetLin( LibIdx, GmfTetrahedra, &Nodes_Tet_row[0], &Nodes_Tet_row[1], &Nodes_Tet_row[2], &Nodes_Tet_row[3] );
        std::cout << " tet " << Nodes_Tet_row[0] << " " << Nodes_Tet_row[1] << " " << Nodes_Tet_row[2] << " " << Nodes_Tet_row[3] << std::endl;
        Nodes_Tet[i] = Nodes_Tet_row;
    }



    // // Close the mesh file !
    GmfCloseMesh( LibIdx );


    std::cout << "stats = " << "nverts " << NmbVer << " ntri " << NmbTri << " ntet " << NmbTet << " npri " << NmbPri << std::endl;

    MPI_Finalize();
        
}

