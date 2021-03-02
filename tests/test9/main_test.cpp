#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include <iomanip>



int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    double* Ptri = new double[3*3];
    Ptri[0*3+0] = 0.0;Ptri[0*3+1] = 0.0;Ptri[0*3+2] = 0.0;
    Ptri[1*3+0] = 1.0;Ptri[1*3+1] = 0.0;Ptri[1*3+2] = 0.0;
    Ptri[2*3+0] = 0.0;Ptri[2*3+1] = 1.0;Ptri[2*3+2] = 0.0;
    
    double dst = ComputeTriSurfaceArea(Ptri);

    if(world_rank == 0)
    {
        if(fabs(dst-0.5)<1.0e-12)
        {
            std::cout << "Test for computing the area of a triangles has PASSED."<<std::endl;
        }
        else{
            std::cout << dst << " Test for computing the area of a triangles has FAILED."<<std::endl;

        }
    }
    
    double* Pquad = new double[4*3];
    Pquad[0*3+0] = 0.0;Pquad[0*3+1] = 0.0;Pquad[0*3+2] = 0.0;
    Pquad[1*3+0] = 1.0;Pquad[1*3+1] = 0.0;Pquad[1*3+2] = 0.0;
    Pquad[2*3+0] = 1.0;Pquad[2*3+1] = 1.0;Pquad[2*3+2] = 0.0;
    Pquad[3*3+0] = 0.0;Pquad[3*3+1] = 1.0;Pquad[3*3+2] = 0.0;

    double dsq = ComputeQuadSurfaceArea(Pquad);

    if(world_rank == 0)
    {
        if(fabs(dsq-1.0)<1.0e-12)
        {
            std::cout << "Test for computing the area of a quadrilateral has PASSED."<<std::endl;
        }
        else
        {
            std::cout << "Test for computing the area of a quadrilateral has FAILED."<<std::endl;

        }
    }
    
    
    double* Ppr = new double[6*3];
    Ppr[0*3+0] = 0.0;Ppr[0*3+1] = 0.0;Ppr[0*3+2] = 0.0;
    Ppr[1*3+0] = 1.0;Ppr[1*3+1] = 0.0;Ppr[1*3+2] = 0.0;
    Ppr[2*3+0] = 0.0;Ppr[2*3+1] = 1.0;Ppr[2*3+2] = 0.0;
    Ppr[3*3+0] = 0.0;Ppr[3*3+1] = 0.0;Ppr[3*3+2] = 1.0;
    Ppr[4*3+0] = 1.0;Ppr[4*3+1] = 0.0;Ppr[4*3+2] = 1.0;
    Ppr[5*3+0] = 0.0;Ppr[5*3+1] = 1.0;Ppr[5*3+2] = 1.0;
    double Vpr = ComputeVolumePrismCell(Ppr);
    if(world_rank == 0)
    {
        if(fabs(Vpr-0.5)<1.0e-12)
        {
            std::cout << "Test for computing the volume of a prism has PASSED."<<std::endl;
        }
        else{
            std::cout << Vpr << std::endl;
            std::cout << "Test for computing the volume of a prism has FAILED."<<std::endl;

        }
    }
    
    double* Pte = new double[4*3];
    Pte[0*3+0] = 0.0;Pte[0*3+1] = 0.0;Pte[0*3+2] = 0.0;
    Pte[1*3+0] = 1.0;Pte[1*3+1] = 0.0;Pte[1*3+2] = 0.0;
    Pte[2*3+0] = 0.0;Pte[2*3+1] = 1.0;Pte[2*3+2] = 0.0;
    Pte[3*3+0] = 0.0;Pte[3*3+1] = 0.0;Pte[3*3+2] = 1.0;
    
    double Vte = ComputeVolumeTetCell(Ppr);
    if(world_rank == 0)
    {
        if(fabs(Vte-(1.0/6.0))<1.0e-12)
        {
            std::cout << "Test for computing the volume of a tetrahedra has PASSED."<<std::endl;
        }
        else{
            std::cout << "Test for computing the volume of a tetrahedra has FAILED."<<std::endl;

        }
    }
  
    
//    delete v0;
//    delete v1;
    MPI_Finalize();
    
}

