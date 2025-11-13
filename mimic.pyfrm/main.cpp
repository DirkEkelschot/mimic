#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include "src/adapt.h"
#include "src/adapt_io.h"

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

    // const char* us3d_grid_file = "mesh/grid.h5";

    std::vector<std::vector<double> > pyfr_xcn;
    std::map<int, std::vector<int> > pyfr_ien_tet;
    std::map<int, std::vector<int> > pyfr_ien_pri;
    std::vector<std::vector<double> > us3d_xcn;
    std::vector<std::vector<std::vector<float> > > pyfr_data;
    // us3d_xcn = ReadDataSetFromFileInParallel_Lite<double>(us3d_grid_file,"xcn",comm,info);

    const char* pyfr_grid_file = "inputs/mesh.pyfrm";
    const char* pyfr_data_file = "inputs/data.pyfrs";

    pyfr_xcn        = ReadVerticesFromPyFRMeshFileInParallel_Lite<double>(pyfr_grid_file,"nodes",comm,info);
    // pyfr_ien_pri = ReadElementsFromPyFRMeshFileInParallel_Lite(pyfr_grid_file,"eles","pri",comm,info);
    pyfr_ien_tet    = ReadElementsFromPyFRMeshFileInParallel_Lite(pyfr_grid_file,"eles","tet",comm,info);
    pyfr_data       = ReadSolutionFromPyFRFileInParallel_Lite(pyfr_data_file, "tet", 1, comm, info);

    // std::cout << "pyfr_ien_tet " << pyfr_ien_pri.size() << "pyfr_ien_tet " << pyfr_ien_tet.size()  << std::endl;
    PartObjectLite* partitionP = new PartObjectLite(pyfr_ien_pri, pyfr_xcn, tetUniMesh->eltype_map, tetUniMesh->eltype_vec, allbcFacesNew, Ne, Nv, comm);
    // // PartObjectLite* partitionT = new PartObjectLite(pyfr_ien_tet, pyfr_xcn, tetUniMesh->eltype_map, tetUniMesh->eltype_vec, allbcFacesNew, Ne, Nv, comm);


    // std::cout << "rank : " << world_rank << " reads in " << us3d_xcn.size() << " " << pyfr_xcn.size() << " vertices " << std::endl;
    // std::cout << "rank : " << world_rank << " reads in " << pyfr_ien_pri.size() << " " << pyfr_ien_tet.size() << " vertices " << std::endl;

    MPI_Finalize();

}