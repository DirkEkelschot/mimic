#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <map>
#include <mpi.h>

std::map<int, std::map<int, std::vector<int> > > readDataInt(const char* filename)
{

    std::ifstream dataFile;
    dataFile.open(filename);
    if(!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file." << std::endl;
    }

    std::map<int, std::map<int, std::vector<int> > > data; // Container to store all rows

    std::string line;
    int lineentry = 0;
    std::map<int,int> id_2_rowsize;
    while (std::getline(dataFile, line)) { // Read each line from the file
        std::istringstream iss(line);  // Create a string stream for the line
        std::vector<int> row;      // Container for the current row
        double number;
        while (iss >> number) {       // Extract numbers from the line
            row.push_back(number);
        }
        int rowsize = row.size(); // Get the number of columns in the current row
        // data[rowsize].push_back(row);
        data[rowsize][lineentry] = row;
        lineentry++;
    }

    std::cout << "Read in data size is " << data.size() << std::endl;

    dataFile.close();

    return data;

}


std::map<int, std::map<int, std::vector<double> > > readDataMat(const char* filename,
                                                            std::map<int, std::map<int, std::vector<int> > > &matsizes)
{

    std::ifstream dataFile;
    dataFile.open(filename);
    if(!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file." << std::endl;
    }

    std::map<int, std::map<int, std::vector<double> > > data; // Container to store all rows

    std::string line;
    int lineentry = 0;
    std::map<int,int> id_2_rowsize;
    while (std::getline(dataFile, line)) { // Read each line from the file
        std::istringstream iss(line);  // Create a string stream for the line
        std::vector<double> row;      // Container for the current row
        std::vector<int> rowsizes(2);
        double number;
        while (iss >> number) {       // Extract numbers from the line
            row.push_back(number);
        }
        int rowsize = row.size(); // Get the number of columns in the current row
        // data[rowsize].push_back(row);
        auto start = row.begin();
        auto end   = row.end() - 2;  
        std::vector<double> subset(start, end);

        data[rowsize][lineentry]     = subset;

        rowsizes[0] = (int) row[rowsize-2];
        rowsizes[1] = (int) row[rowsize-1];
        matsizes[rowsize][lineentry] = rowsizes;
        lineentry++;
    }

    std::cout << "Read in data size is " << data.size() << std::endl;

    dataFile.close();

    return data;

}


std::map<int, std::map<int, std::vector<double> > > readDataDouble(const char* filename)
{

    std::ifstream dataFile;
    dataFile.open(filename);
    if(!dataFile.is_open()) {
        std::cerr << "Error: Unable to open file." << std::endl;
    }

    std::map<int, std::map<int, std::vector<double> > > data; // Container to store all rows

    std::string line;
    int lineentry = 0;
    std::map<int,int> id_2_rowsize;
    while (std::getline(dataFile, line)) { // Read each line from the file
        std::istringstream iss(line);  // Create a string stream for the line
        std::vector<double> row;      // Container for the current row
        double number;
        while (iss >> number) {       // Extract numbers from the line
            row.push_back(number);
        }
        int rowsize = row.size(); // Get the number of columns in the current row
        // data[rowsize].push_back(row);
        data[rowsize][lineentry] = row;
        lineentry++;
    }

    std::cout << "Read in data size is " << data.size() << std::endl;

    dataFile.close();

    return data;

}

int main() {

    clock_t start, end;
    double dur_max,time_taken;

    std::map<int, std::map<int, std::vector<int> > > matsizes_map;
    std::map<int, std::map<int, std::vector<double> > > amat = readDataMat("Amats.dat", matsizes_map);
    std::map<int, std::map<int, std::vector<double> > > bvec = readDataDouble("bvecs.dat");

    // std::map<int, std::map<int, std::vector<int> > > amm = readDataInt("Am.dat");
    // std::map<int, std::map<int, std::vector<int> > > amn = readDataInt("An.dat");

    auto itA = amat.begin();
    auto itb = bvec.begin();

    std::map<int, std::map<int, std::vector<double> > >::iterator dit;
    std::vector<int> sizes;
    std::vector<int> matsizes;
    for(dit=amat.begin();dit!=amat.end();dit++)
    {   
        //std::cout << dit->first << " " << dit->second.size() << std::endl;
        sizes.push_back(dit->second.size());
        matsizes.push_back(dit->first);
        //matrows.push_back(amm[]);
    }

    auto maxElementIter = std::max_element(sizes.begin(), sizes.end());
    int maxElemIndex    = std::distance(sizes.begin(), maxElementIter);

    std::cout << sizes[maxElemIndex] << " " << matsizes[maxElemIndex] << " " << sizes.size() << std::endl;
    int batchSize = sizes[maxElemIndex];
    int matsize = matsizes[maxElemIndex];
    int n = 9;
    int m = (int) matsize/n;
    
    std::cout << "==========sizes==========" << std::endl;
    std::cout << "m = " << m << " " << n << std::endl;
    std::cout << "=========================" << std::endl;
    // const int m = 2;        // Rows per matrix
    // const int n = 2;        // Columns per matrix
    // const int batchSize = 2; // Number of matrices
    // const int ltau = (m < n) ? m : n; // Tau vector length

    // Create CUBLAS handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    // Host memory allocation
    
    for(int i=0;i<sizes.size();i++)
    {
        matsize     = matsizes[i];
        batchSize   =  sizes[i];


        // double* h_A = new double[matsize*batchSize];

        // size_t index = 0;
        // std::map<int, std::map<int, std::vector<double> > >::iterator its;

        // std::cout << amat.size() << std::endl;
        // std::map<int, std::vector<double> >::iterator itmat;

        // for(itmat=amat[matsize].begin();itmat!=amat[matsize].end();itmat++)
        // {
        //     int elid = itmat->first;
        //     int nq   = itmat->second.size();

        //     for(int q=0;q<nq;q++)
        //     {
        //         h_A[index] = itmat->second[q];
        //         index++;
        //     }

        // }
        // start = clock();
        // int ltau = m;
        // double h_Tau[ltau*batchSize] = {0};

        // // Device memory allocation
        // double *d_A, *d_Tau;
        // cudaMalloc(&d_A, matsize*batchSize*sizeof(double));
        // cudaMalloc(&d_Tau, ltau*batchSize*sizeof(double));

        // // Copy matrices to device
        // cudaMemcpy(d_A, h_A, matsize*batchSize*sizeof(double), cudaMemcpyHostToDevice);

        // // Create device pointer arrays
        // double **d_Aarray, **d_TauArray;
        // cudaMalloc(&d_Aarray, batchSize*sizeof(double*));
        // cudaMalloc(&d_TauArray, batchSize*sizeof(double*));

        // // Host-side pointer arrays
        // double *h_Aarray[batchSize];
        // double *h_TauArray[batchSize];
        // for(int i=0; i<batchSize; i++){
        //     h_Aarray[i]   = d_A + i*matsize;
        //     h_TauArray[i] = d_Tau + i*ltau;
        // }

        // // Copy pointer arrays to device
        // cudaMemcpy(d_Aarray, h_Aarray, batchSize*sizeof(double*), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_TauArray, h_TauArray, batchSize*sizeof(double*), cudaMemcpyHostToDevice);
        // end = clock();

        // time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
        // std::cout << "time_taken loading data " << time_taken << std::endl;
        // // Perform batched QR
        // int info;
        // start = clock();
        // std::cout << "Before Batched" << std::endl;
        // cublasStatus_t status = cublasDgeqrfBatched(
        //     handle, 
        //     m, 
        //     9, 
        //     d_Aarray, 
        //     m,    // Leading dimension (lda)
        //     d_TauArray, 
        //     &info, 
        //     batchSize
        // );

        // std::cout << "After Batched" << std::endl;
        // end = clock();

        // time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
        // std::cout << "time_taken " << time_taken << std::endl;
        // // MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        // // if(world_rank == 0)
        // // {
        // //     cout << setprecision(3) << "Time taken to compute " << batchSize << " QR compositions on GPU :          " << dur_max << setprecision(3);
        // //     cout << " sec " << endl;
        // // }


        // // if (status != CUBLAS_STATUS_SUCCESS || info != 0) {
        // //     std::cerr << "QR decomposition failed! Error: " << status << std::endl;
        // //     return 1;
        // // }

        // // Copy results back to host
        // cudaMemcpy(h_A, d_A, m*n*batchSize*sizeof(double), cudaMemcpyDeviceToHost);
        // cudaMemcpy(h_Tau, d_Tau, ltau*batchSize*sizeof(double), cudaMemcpyDeviceToHost);
        
        // // Print results
        
        // // std::cout.precision(4);
        // // for(int b=0; b<batchSize; b++){
        // //     std::cout << "\nMatrix " << b+1 << " results:\nR matrix (upper triangular):\n";
        // //     const double* mat = h_A + b*m*n;
        // //     for(int i=0; i<m; i++){
        // //         for(int j=0; j<n; j++){
        // //             std::cout << (j >= i ? mat[i + j*m] : 0.0) << "\t";
        // //         }
        // //         std::cout << "\n";
        // //     }

        // //     std::cout << "Householder vectors (below diagonal):\n";
        // //     for(int i=0; i<m; i++){
        // //         for(int j=0; j<n; j++){
        // //             std::cout << (j < i ? mat[i + j*m] : 0.0) << "\t";
        // //         }
        // //         std::cout << "\n";
        // //     }

        // //     std::cout << "Tau values: ";
        // //     const double* tau = h_Tau + b*ltau;
        // //     for(int i=0; i<ltau; i++){
        // //         std::cout << tau[i] << " ";
        // //     }
        // //     std::cout << "\n";
        // // }
        
        // //Cleanup
        // cudaFree(d_A);
        // cudaFree(d_Tau);
        // cudaFree(d_Aarray);
        // cudaFree(d_TauArray);
        // cublasDestroy(handle);


        // delete[] h_A;
    }
    

    

    return 0;
}
