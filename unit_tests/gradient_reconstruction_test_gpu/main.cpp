
#include <chrono>
#include <iostream>
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
//#include "Kokkos_Core.hpp"
#include <Kokkos_Core.hpp>
//#include <KokkosBatched_QR_Decl.hpp>
//#include <KokkosBatched_QR_Impl.hpp>
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;




int main() {
    // Initialize Kokkos

  //  MPI_Init(NULL, NULL);
    FILE            *inm;
  //  MPI_Comm comm = MPI_COMM_WORLD;
  //  MPI_Info info = MPI_INFO_NULL;
    int world_size;
   // MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
   // MPI_Comm_rank(comm, &world_rank);


   Kokkos::initialize();
   {
        // Define a simple 3x3 matrix A
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> A("Matrix A", 3, 3);
   
        // Initialize matrix A with some values
        A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
        A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
        A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;
        
        // Perform QR decomposition using KokkosBatched
        // This step involves using the functions from KokkosBatched_QR_Decl.hpp
        // For demonstration, assume we have a function performQR that takes A and returns Q and R
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> Q("Matrix Q", 3, 3);
        Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> R("Matrix R", 3, 3);
        
        // Example call to perform QR decomposition
        // performQR(A, Q, R); // Implement this function using KokkosBatched_QR_Decl.hpp
        
        // Print Q and R to verify the decomposition
        std::cout << "Matrix Q:\n";
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << Q(i,j) << " ";
            }
            std::cout << "\n";
        }
        
        std::cout << "Matrix R:\n";
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << R(i,j) << " ";
            }
            std::cout << "\n";
        }
   }  
    Kokkos::finalize();

//    MPI_Finalize();
    return 0;
}

