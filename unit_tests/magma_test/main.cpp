#include <magma_v2.h>
#include <iostream>
#include <vector>

const int BATCH_SIZE = 1000000;
const int M = 3;
const int N = 3;

int main() {
    // Initialize MAGMA
    magma_init();
    magma_queue_t queue;
    magma_queue_create(0, &queue);

    // Host template matrix (column-major)
    std::vector<double> h_A = {1.0, 4.0, 7.0,
                               2.0, 5.0, 8.0,
                               3.0, 6.0, 9.0};

    // Device memory allocation
    magmaDouble_ptr d_A[BATCH_SIZE];
    magmaDouble_ptr d_tau[BATCH_SIZE];
    magma_int_t info[BATCH_SIZE];
    
    // Create batch array pointers
    double** d_A_array;
    double** d_tau_array;
    magma_malloc((void**)&d_A_array, BATCH_SIZE * sizeof(double*));
    magma_malloc((void**)&d_tau_array, BATCH_SIZE * sizeof(double*));

    // Initialize batch matrices
    for(int i = 0; i < BATCH_SIZE; i++) {
        magma_dmalloc(&d_A[i], M*N);
        magma_dmalloc(&d_tau[i], N);
        magma_dsetmatrix(M, N, h_A.data(), M, d_A[i], M, queue);
        d_A_array[i] = d_A[i];
        d_tau_array[i] = d_tau[i];
    }

    // Perform batched QR factorization
    magma_dgeqrf_batched(
        M, N,
        d_A_array, M,
        d_tau_array,
        info,
        BATCH_SIZE,
        queue
    );

    // Verify first result
    std::vector<double> h_R(M*N);
    magma_dgetmatrix(M, N, d_A[0], M, h_R.data(), M, queue);

    std::cout << "First R matrix:\n";
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < N; j++) {
            std::cout << h_R[i + j*M] << "\t";
        }
        std::cout << "\n";
    }

    // Cleanup
    for(int i = 0; i < BATCH_SIZE; i++) {
        magma_free(d_A[i]);
        magma_free(d_tau[i]);
    }
    magma_free(d_A_array);
    magma_free(d_tau_array);
    magma_queue_destroy(queue);
    magma_finalize();

    return 0;
}
