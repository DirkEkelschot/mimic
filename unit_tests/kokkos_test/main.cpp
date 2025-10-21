#include <Kokkos_Core.hpp>
#include <KokkosBatched_QR_Decl.hpp>

int main() {
    Kokkos::initialize();
    {
        constexpr int numMatrices = 3;
        constexpr int rows = 4;
        constexpr int cols = 3;

        using execution_space = Kokkos::DefaultExecutionSpace;
        using memory_space = execution_space::memory_space;

        // Use LayoutRight for column-major storage
        using matrix_view_type = Kokkos::View<double***, Kokkos::LayoutRight, 
            Kokkos::Device<execution_space, memory_space>>;

        matrix_view_type A("A", numMatrices, rows, cols);
        matrix_view_type Q("Q", numMatrices, rows, rows);
        matrix_view_type R("R", numMatrices, rows, cols);

        // Initialize matrices
        auto h_A = Kokkos::create_mirror_view(A);
        for (int m = 0; m < numMatrices; ++m) {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    h_A(m, i, j) = (i + j + m + 1.0);
                }
            }
        }
        Kokkos::deep_copy(A, h_A);

        // Team policy for batched execution
        Kokkos::TeamPolicy<execution_space> policy(numMatrices, Kokkos::AUTO);
        Kokkos::parallel_for(
            "BatchedQR", policy,
            KOKKOS_LAMBDA(const Kokkos::TeamPolicy<execution_space>::member_type& member) {
                const int m = member.league_rank();
                
                auto A_sub = Kokkos::subview(A, m, Kokkos::ALL(), Kokkos::ALL());
                auto Q_sub = Kokkos::subview(Q, m, Kokkos::ALL(), Kokkos::ALL());
                auto R_sub = Kokkos::subview(R, m, Kokkos::ALL(), Kokkos::ALL());

                KokkosBatched::QR<
                    Kokkos::TeamPolicy<execution_space>::member_type,  // MemberType
                    KokkosBatched::Mode::Team,                         // Execution mode
                    KokkosBatched::Algo::QR::Unblocked                 // Algorithm
                >::invoke(
                    member,  // Team member handle
                    A_sub,   // Input matrix (overwritten with Q)
                    R_sub,   // R output (upper triangular)
                    Q_sub    // Q output (orthogonal matrix)
                );
            }
        );

        // Verify results
        auto h_Q = Kokkos::create_mirror_view(Q);
        auto h_R = Kokkos::create_mirror_view(R);
        Kokkos::deep_copy(h_Q, Q);
        Kokkos::deep_copy(h_R, R);

        for (int m = 0; m < numMatrices; ++m) {
            printf("\nMatrix %d:\n", m);
            
            printf("Q:\n");
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < rows; ++j) {
                    printf("%8.4f ", h_Q(m, i, j));
                }
                printf("\n");
            }
            
            printf("R:\n");
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    printf("%8.4f ", h_R(m, i, j));
                }
                printf("\n");
            }
        }
    }
    Kokkos::finalize();
    return 0;
}
