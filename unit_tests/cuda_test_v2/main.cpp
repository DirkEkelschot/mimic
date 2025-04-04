#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>

int main() {
    const int m = 2;        // Rows per matrix
    const int n = 2;        // Columns per matrix
    const int batchSize = 2; // Number of matrices
    const int ltau = (m < n) ? m : n; // Tau vector length

    // Create CUBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);

    // Host memory allocation
    double h_A[m*n*batchSize] = {
        // First matrix (column-major)
        1.0, 3.0,   // Column 1
        2.0, 4.0,   // Column 2
        
        // Second matrix (column-major)
        5.0, 7.0,   // Column 1
        6.0, 8.0    // Column 2
    };
    
    double h_Tau[ltau*batchSize] = {0};

    // Device memory allocation
    double *d_A, *d_Tau;
    cudaMalloc(&d_A, m*n*batchSize*sizeof(double));
    cudaMalloc(&d_Tau, ltau*batchSize*sizeof(double));

    // Copy matrices to device
    cudaMemcpy(d_A, h_A, m*n*batchSize*sizeof(double), cudaMemcpyHostToDevice);

    // Create device pointer arrays
    double **d_Aarray, **d_TauArray;
    cudaMalloc(&d_Aarray, batchSize*sizeof(double*));
    cudaMalloc(&d_TauArray, batchSize*sizeof(double*));

    // Host-side pointer arrays
    double *h_Aarray[batchSize];
    double *h_TauArray[batchSize];
    for(int i=0; i<batchSize; i++){
        h_Aarray[i] = d_A + i*m*n;
        h_TauArray[i] = d_Tau + i*ltau;
    }

    // Copy pointer arrays to device
    cudaMemcpy(d_Aarray, h_Aarray, batchSize*sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(d_TauArray, h_TauArray, batchSize*sizeof(double*), cudaMemcpyHostToDevice);

    // Perform batched QR
    int info;
    cublasStatus_t status = cublasDgeqrfBatched(
        handle, 
        m, 
        n, 
        d_Aarray, 
        m,    // Leading dimension (lda)
        d_TauArray, 
        &info, 
        batchSize
    );

    if (status != CUBLAS_STATUS_SUCCESS || info != 0) {
        std::cerr << "QR decomposition failed! Error: " << status << std::endl;
        return 1;
    }

    // Copy results back to host
    cudaMemcpy(h_A, d_A, m*n*batchSize*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Tau, d_Tau, ltau*batchSize*sizeof(double), cudaMemcpyDeviceToHost);

    // Print results
    std::cout.precision(4);
    for(int b=0; b<batchSize; b++){
        std::cout << "\nMatrix " << b+1 << " results:\nR matrix (upper triangular):\n";
        const double* mat = h_A + b*m*n;
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                std::cout << (j >= i ? mat[i + j*m] : 0.0) << "\t";
            }
            std::cout << "\n";
        }

        std::cout << "Householder vectors (below diagonal):\n";
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                std::cout << (j < i ? mat[i + j*m] : 0.0) << "\t";
            }
            std::cout << "\n";
        }

        std::cout << "Tau values: ";
        const double* tau = h_Tau + b*ltau;
        for(int i=0; i<ltau; i++){
            std::cout << tau[i] << " ";
        }
        std::cout << "\n";
    }

    // Cleanup
    cudaFree(d_A);
    cudaFree(d_Tau);
    cudaFree(d_Aarray);
    cudaFree(d_TauArray);
    cublasDestroy(handle);

    return 0;
}
