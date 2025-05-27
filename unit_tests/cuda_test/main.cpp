#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <iostream>
#include <vector>

const int NUM_STREAMS = 32;  // Optimal for most GPUs
const int m = 3;
const int n = 3;
const int BATCH_SIZE = 1000000;

int main() {
    // Original matrix (column-major)
    std::vector<double> h_template = {1.0, 4.0, 7.0,  // Column 1
                                      2.0, 5.0, 8.0,  // Column 2
                                      3.0, 6.0, 9.0}; // Column 3

    // Device memory allocation
    double* d_A;       // All matrices (BATCH_SIZE * m*n elements)
    double* d_tau;     // Tau values (BATCH_SIZE * n elements)
    int* d_info;       // Info codes (BATCH_SIZE elements)
    
    cudaMalloc(&d_A, BATCH_SIZE * m*n * sizeof(double));
    cudaMalloc(&d_tau, BATCH_SIZE * n * sizeof(double));
    cudaMalloc(&d_info, BATCH_SIZE * sizeof(int));

    // Initialize matrices (copy template to all positions)
    for(int i = 0; i < BATCH_SIZE; i++) {
        cudaMemcpy(d_A + i*m*n, h_template.data(), 
                  m*n*sizeof(double), cudaMemcpyHostToDevice);
    }

    // Create streams and handles
    cudaStream_t streams[NUM_STREAMS];
    cusolverDnHandle_t handles[NUM_STREAMS];
    for(int i = 0; i < NUM_STREAMS; i++) {
        cudaStreamCreate(&streams[i]);
        cusolverDnCreate(&handles[i]);
        cusolverDnSetStream(handles[i], streams[i]);
    }

    // Workspace setup
    int lwork;
    double* d_work[NUM_STREAMS];
    cusolverDnDgeqrf_bufferSize(handles[0], m, n, d_A, m, &lwork);
    for(int i = 0; i < NUM_STREAMS; i++) {
        cudaMalloc(&d_work[i], lwork * sizeof(double));
    }

    // Process batch using streams
    for(int i = 0; i < BATCH_SIZE; i++) {
        int stream_id = i % NUM_STREAMS;
        double* matrix_ptr = d_A + i*m*n;
        double* tau_ptr = d_tau + i*n;
        int* info_ptr = d_info + i;
        // std::cout << stream_id << " " << handles[stream_id] << std::endl;
        cusolverDnDgeqrf(handles[stream_id], m, n,
                        matrix_ptr, m, tau_ptr,
                        d_work[stream_id], lwork, info_ptr);
    }

    // Synchronize and check errors
    std::vector<int> h_info(BATCH_SIZE);
    for(int i = 0; i < NUM_STREAMS; i++) {
        cudaStreamSynchronize(streams[i]);
    }
    cudaMemcpy(h_info.data(), d_info, BATCH_SIZE*sizeof(int), cudaMemcpyDeviceToHost);

    // Verify first result
    std::vector<double> h_first_result(m*n);
    cudaMemcpy(h_first_result.data(), d_A, m*n*sizeof(double), cudaMemcpyDeviceToHost);

    std::cout << "First R matrix:" << std::endl;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << h_first_result[i + j*m] << "\t";
        }
        std::cout << std::endl;
    }

    // Cleanup
    for(int i = 0; i < NUM_STREAMS; i++) {
        cudaFree(d_work[i]);
        cusolverDnDestroy(handles[i]);
        cudaStreamDestroy(streams[i]);
    }
    cudaFree(d_A);
    cudaFree(d_tau);
    cudaFree(d_info);

    return 0;
}
