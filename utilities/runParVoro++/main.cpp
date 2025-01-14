
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
#include <random>

#include <voronoi/VoronoiKD.hpp>
// #include <cstdio>
#include <stdio.h>
//#include "/Users/dekelsch/mimic_libmesh/utilities/partitionTetrahedra/build/ThirdParty/dist/include/libmeshb7.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

typedef ParticleBlock               Block;



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
    
    const char* fn_grid = "inputs/grid101010.h5";
    const char* fn_conn = "inputs/grid101010.h5";

    Inputs* inputs = ReadXmlFile(comm, "inputs/metric.xml");

    mesh* meshRead = ReadUS3DMesh(fn_conn,fn_grid,
                                        inputs->ReadFromStats,
                                        inputs->StateVar,
                                        comm,info);

    double xmin = -1.0/* minimum x coordinate */;
    double xmax =  1.0/* maximum x coordinate */;
    double ymin = -1.0/* minimum y coordinate */;
    double ymax =  1.0/* maximum y coordinate */;
    double zmin = -1.0/* minimum z coordinate */;
    double zmax =  1.0/* maximum z coordinate */;
    int nx = 10/* number of grid cells in x direction */;
    int ny = 10/* number of grid cells in y direction */;
    int nz = 10/* number of grid cells in z direction */;

    //container con(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, false, false, false, 8);

    diy::mpi::environment     env(argc, argv);
    diy::mpi::communicator    world;

    int rank = world.rank();      // MPI usual
    int size = world.size();      // MPI usual
    int tot_blocks = size;        // total number of blocks in the domain
    int num_threads = 1;          // number of threads diy can use
    int mem_blocks = -1;          // number of blocks to keep in memory
    int num_files = -1;           // number of output files for each timestep
    
    bool kdtree, wrap, debug, help;
    // Create VoronoiKD object
    //int nx=10, ny=10, nz=10;
    diy::Master  master(world, num_threads, mem_blocks, &Block::create_block, &Block::destroy_block);

    Voronoi_kd(master, wrap, nx, ny, nz);
    // Compute Voronoi tessellation in parallel
    

    
    MPI_Finalize();
        
}

