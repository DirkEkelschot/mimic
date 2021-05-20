#include "adapt.h"
#include "adapt_array.h"
#include "adapt_datatype.h"
#ifndef ADAPT_IO_H
#define ADAPT_IO_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline hid_t hid_from_type(const int &)
{
    return H5T_NATIVE_INT;
}

inline hid_t hid_from_type(const double &)
{
    return H5T_NATIVE_DOUBLE;
}

inline hid_t hid_from_type(const char &)
{
    return H5T_STRING;
}

template <typename T>
inline hid_t hid_from_type() {
    return hid_from_type(T());
}

inline hid_t h5tools_get_native_type(hid_t type)
{
    hid_t p_type;
    H5T_class_t type_class;

    type_class = H5Tget_class(type);
    if (type_class == H5T_BITFIELD)
        p_type = H5Tcopy(type);
    else
        p_type = H5Tget_native_type(type, H5T_DIR_DEFAULT);

return(p_type);
}

double* ReadDataSetDoubleFromFile(const char* file_name, const char* dataset_name);

template<typename T>
Array<T>* ReadDataSetFromFile(const char* file_name, const char* dataset_name)
{
    
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
    hid_t dspace          = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    
    int nrow            = dims[0];
    int ncol            = dims[1];
    
    hid_t memspace_id   = H5Screate_simple( 2, dims, NULL );
    hid_t file_space_id = H5Dget_space(dset_id);
    
    Array<T>* A_t = new Array<T>(nrow,ncol);
    
    hid_t status = H5Dread(dset_id, hid_from_type<T>(), memspace_id, file_space_id, H5P_DEFAULT, A_t->data);
    
    status = H5Dclose(dset_id);
    
    return A_t;
}


template<typename T>
Array<T>* ReadDataSetFromGroupFromFile(const char* file_name, const char* group_name, const char* dataset_name)
{
    hid_t status;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id       = H5Gopen(file_id,group_name,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(group_id,dataset_name,H5P_DEFAULT);
    
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    
    int nrow            = dims[0];
    int ncol            = dims[1];
    
    hid_t memspace_id   = H5Screate_simple( 2, dims, NULL );
    hid_t file_space_id = H5Dget_space(dset_id);

    Array<T>* A_t;
    
    if(hid_from_type<T>()==H5T_STRING)
    {
        hid_t type = H5Dget_type(dset_id);
        hid_t native_type = h5tools_get_native_type(type);
        int n_element = dims[0];
        size_t type_size = std::max(H5Tget_size(type), H5Tget_size(native_type));
        
        A_t = new Array<T>(n_element,type_size);
        
        status = H5Dread(dset_id, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, A_t->data);
        
        status = H5Tclose(native_type);
        status = H5Tclose(type);
    }
    else{
        
        A_t = new Array<T>(nrow,ncol);
        
        status = H5Dread(dset_id, hid_from_type<T>(), memspace_id, file_space_id, H5P_DEFAULT, A_t->data);
    }
     
    status = H5Dclose(dset_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
    return A_t;
}


template<typename T>
Array<T>* ReadDataSetFromRunInFile(const char* file_name, const char* run_name,const char* dataset_name)
{

    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id       = H5Gopen(file_id,"solution",H5P_DEFAULT);
    hid_t run_id         = H5Gopen(group_id,run_name,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(run_id,dataset_name,H5P_DEFAULT);
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    
    int nrow             = dims[0];
    int ncol             = dims[1];
    
    hid_t memspace_id    = H5Screate_simple( 2, dims, NULL );
    hid_t file_space_id  = H5Dget_space(dset_id);
    
    Array<T>* A_t = new Array<T>(nrow,ncol);
    std::clock_t start;
    double duration;
    start = std::clock();
    hid_t status = H5Dread(dset_id, hid_from_type<T>(), memspace_id, file_space_id, H5P_DEFAULT, A_t->data);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "timer_serial = " << duration << std::endl;
    
    status = H5Dclose(dset_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
    
    return A_t;
}

template<typename T>
Array<T>* ReadUS3DGhostCellsFromRun(const char* file_name, const char* run_name,const char* dataset_name, int Nel)
{
    
    herr_t ret;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id       = H5Gopen(file_id,"solution",H5P_DEFAULT);
    hid_t run_id         = H5Gopen(group_id,run_name,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(run_id,dataset_name,H5P_DEFAULT);
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow             = dims[0];
    int ncol             = dims[1];
    int N = nrow-Nel;
    int g_offset = Nel;
    
    hsize_t offset[2];   // hyperslab offset in the file
    hsize_t count[2];    // size of the hyperslab in the file
    offset[0] = g_offset;
    offset[1] = 0;
    count[0]  = N;
    count[1]  = ncol;
    
    ret = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    /*
    * Define the memory dataspace.
    */
    hsize_t     dimsm[2];              /* memory space dimensions */
    dimsm[0] = N;
    dimsm[1] = ncol;
    hid_t memspace = H5Screate_simple (2, dimsm, NULL);
    /*
    * Define memory hyperslab.
    */
    hsize_t      offset_out[2];   // hyperslab offset in memory
    hsize_t      count_out[2];    // size of the hyperslab in memory
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0]  = N;
    count_out[1]  = ncol;
    
    ret = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL,count_out, NULL);
    
    Array<T>* A_t = new Array<T>(N,ncol);
    
    ret = H5Dread (dset_id, hid_from_type<T>(), memspace, dspace, H5P_DEFAULT, A_t->data);
    
    return A_t;
}

template<typename T>
ParArray<T>* ReadDataSetFromRunInFileInParallel(const char* file_name, const char* run_name,const char* dataset_name, int g, int Nel, MPI_Comm comm, int rank, int size, MPI_Info info)
{
    herr_t ret;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id       = H5Gopen(file_id,"solution",H5P_DEFAULT);
    hid_t run_id         = H5Gopen(group_id,run_name,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(run_id,dataset_name,H5P_DEFAULT);
    
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow             = dims[0];
    int ncol             = dims[1];
    int N = nrow;
    int g_offset = 0;
    if (strcmp(dataset_name, "interior") == 0)
    {
        if(g == 1)
        {
            g_offset = Nel;
            N = nrow-Nel;
        }
        else{
            g_offset = 0;
            N = Nel;
        }
        
    }
    
    if (strcmp(dataset_name, "stats-mean") == 0)
    {
        if(g == 1)
        {
            g_offset = Nel;
            N = nrow-Nel;
        }
        else{
            g_offset = 0;
            N = Nel;
        }
        
    }
    //std::cout << "Nel = " << N << std::endl;
    ParArray<T>* A_ptmp = new ParArray<T>(N,ncol,comm,rank,size);
    
    hsize_t              offsets[2];
    hsize_t              counts[2];
    offsets[0]           = A_ptmp->getOffset(rank);
    offsets[1]           = 0;
    counts[0]            = A_ptmp->getNloc(rank);
    counts[1]            = ncol;
    
    ret=H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    
    hsize_t     dimsm[2];
    dimsm[0] = A_ptmp->getNloc(rank);
    dimsm[1] = ncol;
    
    hid_t memspace = H5Screate_simple (2, dimsm, NULL);
    
    hsize_t     offsets_out[2];
    hsize_t     counts_out[2];
    offsets_out[0] = 0;
    offsets_out[1] = 0;
    counts_out[0]  = A_ptmp->getNloc(rank);
    counts_out[1]  = ncol;
    
    ret = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets_out, NULL,counts_out, NULL);
        
    ret = H5Dread (dset_id, hid_from_type<T>(), memspace, dspace, H5P_DEFAULT, A_ptmp->data);
    
    H5Dclose(dset_id);
    H5Fclose(file_id);
    
    return A_ptmp;
}

template<typename T>
ParArray<T>* ReadDataSetFromFileInParallel(const char* file_name, const char* dataset_name, MPI_Comm comm, int rank, int size, MPI_Info info)
{
    
    hid_t acc_tpl1          = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    acc_tpl1                = H5P_DEFAULT;
    // Open file and data set to get dimensions of array;
    
    hid_t file_id           = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t dset_id           = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
    hid_t dspace            = H5Dget_space(dset_id);
    int ndims               = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow                = dims[0];
    int ncol                = dims[1];
    int N                   = nrow;
    int nloc                = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset              = rank*int(N/size) + MIN(rank, N%size);
    ParArray<T>* PA         = new ParArray<T>(N,ncol,comm,rank,size);
    
    
    hsize_t              offsets[2];
    hsize_t              counts[2];
    offsets[0]           = offset;
    offsets[1]           = 0;
    counts[0]            = nloc;
    counts[1]            = ncol;
    
    ret=H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    
    hsize_t     dimsm[2];
    dimsm[0] = nloc;
    dimsm[1] = ncol;
    
    hid_t memspace = H5Screate_simple (2, dimsm, NULL);
    
    hsize_t     offsets_out[2];
    hsize_t     counts_out[2];
    offsets_out[0] = 0;
    offsets_out[1] = 0;
    counts_out[0]  = nloc;
    counts_out[1]  = ncol;
     
    ret = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets_out, NULL,counts_out, NULL);
    
    ret = H5Dread (dset_id, hid_from_type<T>(), memspace, dspace, H5P_DEFAULT, PA->data);

    H5Sclose(dspace);
    H5Sclose(memspace);

    H5Dclose(dset_id);
    H5Fclose(file_id);

    
    return PA;
}

template<typename T>
Array<T>* ReadDataSetFromFileInParallelToRoot(const char* file_name, const char* dataset_name, MPI_Comm comm, int rank, int size, MPI_Info info)
{
    
    hid_t acc_tpl1       = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    
    acc_tpl1             = H5P_DEFAULT;
    // Open file and data set to get dimensions of array;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow             = dims[0];
    int ncol             = dims[1];
    int N                = nrow;
    
    ParArray<T>* parA     = new ParArray<T>(N,ncol,comm,rank,size);
    
    hsize_t              offsets[2];
    hsize_t              counts[2];
    offsets[0]           = parA->getOffset(rank);
    offsets[1]           = 0;
    counts[0]            = parA->getNloc(rank);
    counts[1]            = ncol;
    
    ret=H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    
    hsize_t     dimsm[2];
    dimsm[0] = parA->getNloc(rank);
    dimsm[1] = ncol;
    
    hid_t memspace = H5Screate_simple (2, dimsm, NULL);
    
    hsize_t     offsets_out[2];
    hsize_t     counts_out[2];
    offsets_out[0] = 0;
    offsets_out[1] = 0;
    counts_out[0]  = parA->getNloc(rank);
    counts_out[1]  = ncol;
    
    ret = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets_out, NULL,counts_out, NULL);
    
    ret = H5Dread (dset_id, hid_from_type<T>(), memspace, dspace, H5P_DEFAULT, parA->data);
    
    Array<T>* A_ptot = new Array<T>(N,ncol);
    
    int* nlocs_tmp   = new int[size];
    int* offsets_tmp = new int[size];
    
    for(int i=0;i<size;i++)
    {
        nlocs_tmp[i] = parA->getNloc(i)*ncol;
        offsets_tmp[i] = parA->getOffset(i)*ncol;
    }
    
    if(hid_from_type<T>()==H5T_NATIVE_INT)
    {
        MPI_Gatherv(parA->data, parA->getNloc(rank)*ncol, MPI_INT, A_ptot->data, nlocs_tmp, offsets_tmp, MPI_INT, 0, comm);
    }
    if(hid_from_type<T>()==H5T_NATIVE_DOUBLE)
    {
        MPI_Gatherv(parA->data, parA->getNloc(rank)*ncol, MPI_DOUBLE, A_ptot->data, nlocs_tmp, offsets_tmp, MPI_DOUBLE, 0, comm);
    }
    

    H5Dclose(dset_id);
    H5Fclose(file_id);
    
    delete parA;
    delete[] nlocs_tmp;
    delete[] offsets_tmp;
    
    return A_ptot;
}

template<typename T>
Array<T>* ReadDataSetFromFileInParallelToAll(const char* file_name, const char* dataset_name, MPI_Comm comm, int rank, int size, MPI_Info info)
{
    
    hid_t acc_tpl1       = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    
    acc_tpl1             = H5P_DEFAULT;
    // Open file and data set to get dimensions of array;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t dset_id        = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
    hid_t dspace         = H5Dget_space(dset_id);
    int ndims            = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow             = dims[0];
    int ncol             = dims[1];
    int N                = nrow;
    
    ParArray<T>* parA     = new ParArray<T>(N,ncol,comm,rank,size);
    
    
    hsize_t              offsets[2];
    hsize_t              counts[2];
    offsets[0]           = parA->getOffset(rank);
    offsets[1]           = 0;
    counts[0]            = parA->getNloc(rank);
    counts[1]            = ncol;
    
    ret=H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    
    hsize_t     dimsm[2];
    dimsm[0] = parA->getNloc(rank);
    dimsm[1] = ncol;
    
    hid_t memspace = H5Screate_simple (2, dimsm, NULL);
    
    hsize_t     offsets_out[2];
    hsize_t     counts_out[2];
    offsets_out[0] = 0;
    offsets_out[1] = 0;
    counts_out[0]  = parA->getNloc(rank);
    counts_out[1]  = ncol;
    
    ret = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets_out, NULL,counts_out, NULL);
    
    ret = H5Dread (dset_id, hid_from_type<T>(), memspace, dspace, H5P_DEFAULT, parA->data);
    
    Array<T>* A_ptot = new Array<T>(N,ncol);
    
    int* nlocs_tmp   = new int[size];
    int* offsets_tmp = new int[size];
    
    for(int i=0;i<size;i++)
    {
        nlocs_tmp[i] = parA->getNloc(i)*ncol;
        offsets_tmp[i] = parA->getOffset(i)*ncol;
    }
    
    if(hid_from_type<T>()==H5T_NATIVE_INT)
    {
        MPI_Allgatherv(parA->data, parA->getNloc(rank)*ncol, MPI_INT, A_ptot->data, nlocs_tmp, offsets_tmp, MPI_INT, comm);
    }
    if(hid_from_type<T>()==H5T_NATIVE_DOUBLE)
    {
        MPI_Allgatherv(parA->data, parA->getNloc(rank)*ncol, MPI_DOUBLE, A_ptot->data, nlocs_tmp, offsets_tmp, MPI_DOUBLE, comm);
    }
    

    H5Dclose(dset_id);
    H5Fclose(file_id);
    
    delete parA;
    return A_ptot;
}

std::vector<double> ReadMetricInputs(const char* fn_metric);

void WriteUS3DGridFromMMG_it0(MMG5_pMesh mmgMesh,MMG5_pSol mmgSol, US3D* us3d);

void WriteUS3DGridFromMMG_itN(MMG5_pMesh mmgMesh,MMG5_pSol mmgSol, US3D* us3d);

US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, int readFromStats, MPI_Comm comm,  int size, int rank, MPI_Info info);

#endif
