#include "adapt.h"
//#include "adapt_part_func.h"
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
ParArray<T>* ReadDataSetFromRunInFileInParallel(const char* file_name, const char* run_name,const char* dataset_name, int g, int Nel, MPI_Comm comm, MPI_Info info)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
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
    //std::cout << "Nel = " << N << std::endl;
    ParArray<T>* A_ptmp = new ParArray<T>(N,ncol,comm);
    
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
    /*
    std::cout << "============filled?===" << rank << "=======" << std::endl;
    for(int i=offset;i<offset+1;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            std::cout << A_ptmp.getVal(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "=============================" << std::endl;
    
    Array<T> A_ptot(nrow,ncol);
    */
    
    //MPI_Reduce(A_ptmp.data, A_ptot.data, nrow*ncol, MPI_DOUBLE, MPI_SUM, 0,comm);
    /*
    if (rank != 0)
    {
        MPI_Send(A_ptmp.data, nloc*ncol, MPI_DOUBLE, 0, 1234, MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(A_ptmp.data, nloc*ncol, MPI_DOUBLE, rank, 1234, MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    }
    */
    
    
    
    //ret=H5Pclose(acc_tpl1);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    
    return A_ptmp;
}











template<typename T>
ParArray<T>* ReadDataSetFromFileInParallel(const char* file_name, const char* dataset_name, MPI_Comm comm, MPI_Info info)
{
    
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::cout << rank << " " << size << std::endl;
    
    hid_t acc_tpl1          = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    //herr_t ret            = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    //herr_t ret            = H5Pset_dxpl_mpio(,comm,info);
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
    ParArray<T>* PA         = new ParArray<T>(N,ncol,comm);
    //ParArray<T>* parA     = new ParArray<T>(N,ncol,comm);
    //ParallelState* pstate = parA->getParallelState();
    
    
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
Array<T>* ReadDataSetFromFileInParallelToRoot(const char* file_name, const char* dataset_name, MPI_Comm comm, MPI_Info info)
{
    
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::cout << rank << " " << size << std::endl;
    
    hid_t acc_tpl1       = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    //herr_t ret           = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    
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
    
    //Partition the array;
    // compute local number of rows per proc;
//    ParVar* pv = CreateParallelData(N, comm);
//
//    Array<T>* A_pt = new Array<T>(pv->nlocs[rank],ncol);
//
//    A_pt->offset = pv->offsets[rank];
//    A_pt->nglob  = nrow;
    ParArray<T>* parA     = new ParArray<T>(N,ncol,comm);
    
    
    
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
Array<T>* ReadDataSetFromFileInParallelToAll(const char* file_name, const char* dataset_name, MPI_Comm comm, MPI_Info info)
{
    
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::cout << rank << " " << size << std::endl;
    
    hid_t acc_tpl1       = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    //herr_t ret           = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    
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
    
    //Partition the array;
    // compute local number of rows per proc;
//    ParVar* pv = CreateParallelData(N, comm);
//
//    Array<T>* A_pt = new Array<T>(pv->nlocs[rank],ncol);
//
//    A_pt->offset = pv->offsets[rank];
//    A_pt->nglob  = nrow;
    ParArray<T>* parA     = new ParArray<T>(N,ncol,comm);
    
    
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


void WriteUS3DGridFromMMG(MMG5_pMesh mmgMesh, US3D* us3d,
                          std::map<int,std::vector<int> > bnd_face_map,
                          std::set<std::set<int> > u_fset);

//US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, MPI_Comm comm, MPI_Info info);

US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, MPI_Comm comm, MPI_Info info);









//US3dData ReadUS3Ddata(const char* filename, const char* run){
    //
    
    //hid_t status;
    
    //US3dData V;
    
   // Array<double>* C          = ReadDataSetFromFile<double>("../../grids/grid_kl4.h5","xcn");
    //Array<double>* boundaries = ReadDataSetFromRunInFile<double>("../../grids/data_r.h5","run_1","boundaries");
    //Array<double>* interior   = ReadDataSetFromRunInFile<double>("../../grids/data_r.h5","run_1","interior");
    
    //int* dim = boundaries->getDim();
   
    //std::cout << dim[0] << " " << dim[1] << std::endl;
    
    /*
    for(int i=0;i<dim[0];i++)
    {
      for(int j=0;j<dim[1];j++)
      {
          std::cout << boundaries.getVal(i,j) << " ";
      }
      std::cout << std::endl;
    }
    */
    //for(int i;i<)
    
    /*
    Array< double > Array_d(2,3);
    Array_d.SetVal(1,1,1.12453);
    std::cout << Array_d.GetVal(1,1) << std::endl;
    */
    
    //Array<double>A;
    //A.data=XCN;
    //A.nrow=265850;
    //A.ncol=3;
    //PrintMat(A);
    //std::cout << ncol << " " << nrow << std::endl;
    //status = H5Dclose(xcn_id);
    /*
    
    DataSet dataset_grd           = file_grd.openDataSet("xcn");
    DataSpace dataspace_grd       = dataset_grd.getSpace();
    
    DataSet dataset_grd_ifn           = file_grd.openDataSet("ifn");
    DataSpace dataspace_grd_ifn       = dataset_grd_ifn.getSpace();
    
    hsize_t dims_grd_ifn_out[2];
    int ndims_grd_ifn = dataspace_grd_ifn.getSimpleExtentDims( dims_grd_ifn_out, NULL);
    double *  ifn     = new double[dims_grd_ifn_out[0]*dims_grd_ifn_out[1]];
    
    H5File file_conn("conn.h5", H5F_ACC_RDWR);
    DataSet dataset_conn           = file_conn.openDataSet("iee");
    DataSpace dataspace_conn       = dataset_conn.getSpace();
    
    DataSet dataset_conn_ief       = file_conn.openDataSet("ief");
    DataSpace dataspace_conn_ief   = dataset_conn_ief.getSpace();
    
    DataSet dataset_conn_ife       = file_conn.openDataSet("ife");
    DataSpace dataspace_conn_ife   = dataset_conn_ife.getSpace();
    
    DataSet dataset_conn_ien       = file_conn.openDataSet("ien");
    DataSpace dataspace_conn_ien   = dataset_conn_ien.getSpace();
    
    Group grp = file_conn.openGroup ("zones");
    DataSet dataset_zdefs = grp.openDataSet("zdefs");
    DataSpace dataspace_zdefs   = dataset_zdefs.getSpace();
    
    hsize_t dims_zdefs_out[2];
    int ndims_zdefs = dataspace_zdefs.getSimpleExtentDims( dims_zdefs_out, NULL);
    std::cout << dims_zdefs_out[0] << " " << dims_zdefs_out[1] << std::endl;
    */
     
     
     
     
    /*
    H5File file(filename, H5F_ACC_RDWR);
    Group grp = file.openGroup ("solution");
    Group grp2 = grp.openGroup(run);
    DataSet grp_bound_dataset           = grp2.openDataSet("boundaries");
    DataSpace grp_bound_dataspace       = grp_bound_dataset.getSpace();
    
    DataSet grp_interior_dataset        = grp2.openDataSet("interior");
    DataSpace grp_interior_dataspace    = grp_interior_dataset.getSpace();
    */
    // Read the boundary data from the us3d data.h5 file.

    
    
    /*
    hsize_t dims_grd_out[2];
    int ndims_grd = dataspace_grd.getSimpleExtentDims( dims_grd_out, NULL);
    
    hsize_t dims_conn_out[2];
    int ndims_conn = dataspace_conn.getSimpleExtentDims( dims_conn_out, NULL);
    
    hsize_t dims_conn_ief_out[2];
    int ndims_conn_ief = dataspace_conn_ief.getSimpleExtentDims( dims_conn_ief_out, NULL);
    
    hsize_t dims_conn_ife_out[2];
    int ndims_conn_ife = dataspace_conn_ife.getSimpleExtentDims( dims_conn_ife_out, NULL);
    
    hsize_t dims_conn_ien_out[2];
    int ndims_conn_ien = dataspace_conn_ien.getSimpleExtentDims( dims_conn_ien_out, NULL);
    */
     
     
     
     
    /*hsize_t dims_bound_out[2];
    int ndims_bound = grp_bound_dataspace.getSimpleExtentDims( dims_bound_out, NULL);
    
    hsize_t dims_interior_out[2];
    int ndims_interior_out = grp_interior_dataspace.getSimpleExtentDims( dims_interior_out, NULL);
    */
    // Construct US3DData structure.
    
    
    /*
    V.rows_grd=dims_grd_out[0];
    V.cols_grd=dims_grd_out[1];
    
    V.rows_conn=dims_conn_out[0];
    V.cols_conn=dims_conn_out[1];
    
    //V.rows_bound=dims_bound_out[0];
    //V.cols_bound=dims_bound_out[1];
    
    //V.rows_interior=dims_interior_out[0];
    //V.cols_interior=dims_interior_out[1];
    
    //std::cout << "Connection " << dims_conn_out[0]<<" " <<dims_conn_out[1] << std::endl;
    
    V.Coordinates  = new double[dims_grd_out[0]*dims_grd_out[1]];
    V.Connection   = new int[dims_conn_out[0]*dims_conn_out[1]];
    
    int *  ief     = new int[dims_conn_ief_out[0]*dims_conn_ief_out[1]];
    int *  ife     = new int[dims_conn_ife_out[0]*dims_conn_ife_out[1]];
    int *  ien     = new int[dims_conn_ien_out[0]*dims_conn_ien_out[1]];
    int *  zdefs_i     = new int[dims_zdefs_out[0]*dims_zdefs_out[1]];
    //V.Boundaries = new double[dims_bound_out[0]*dims_bound_out[1]];
    //V.Interior   = new double[dims_interior_out[0]*dims_interior_out[1]];

    const int RANK_OUT = 2;
    
    hsize_t dimsm_zdefs[2];
    dimsm_zdefs[0] = dims_zdefs_out[0];
    dimsm_zdefs[1] = dims_zdefs_out[1];
    DataSpace memspace_zdefs( RANK_OUT, dimsm_zdefs );
    
    hsize_t dimsm_grd[2];
    dimsm_grd[0] = dims_grd_out[0];
    dimsm_grd[1] = dims_grd_out[1];
    DataSpace memspace_grd( RANK_OUT, dimsm_grd );
    
    hsize_t dimsm_grd_ifn[2];
    dimsm_grd_ifn[0] = dims_grd_ifn_out[0];
    dimsm_grd_ifn[1] = dims_grd_ifn_out[1];
    DataSpace memspace_grd_ifn( RANK_OUT, dimsm_grd_ifn );
    
    hsize_t dimsm_conn[2];
    dimsm_conn[0] = dims_conn_out[0];
    dimsm_conn[1] = dims_conn_out[1];
    DataSpace memspace_conn( RANK_OUT, dimsm_conn );
    
    hsize_t dimsm_conn_ief[2];
    dimsm_conn_ief[0] = dims_conn_ief_out[0];
    dimsm_conn_ief[1] = dims_conn_ief_out[1];
    DataSpace memspace_conn_ief( RANK_OUT, dimsm_conn_ief );
    
    hsize_t dimsm_conn_ife[2];
    dimsm_conn_ife[0] = dims_conn_ife_out[0];
    dimsm_conn_ife[1] = dims_conn_ife_out[1];
    DataSpace memspace_conn_ife( RANK_OUT, dimsm_conn_ife );
    
    hsize_t dimsm_conn_ien[2];
    dimsm_conn_ien[0] = dims_conn_ien_out[0];
    dimsm_conn_ien[1] = dims_conn_ien_out[1];
    DataSpace memspace_conn_ien( RANK_OUT, dimsm_conn_ien );
     */
    /*
    hsize_t dimsm_bound[2];
    dimsm_bound[0] = dims_bound_out[0];
    dimsm_bound[1] = dims_bound_out[1];
    DataSpace memspace_bound( RANK_OUT, dimsm_bound );
    
    hsize_t dimsm_interior[2];
    dimsm_interior[0] = dims_interior_out[0];
    dimsm_interior[1] = dims_interior_out[1];
    
    DataSpace memspace_interior( RANK_OUT, dimsm_interior );
    */
    //DataType dt_org = grp_bound_dataset.getDataType();
    
    
    
    /*
    DataType dt_org = dt_double;
    if(dt_org==dt_int)
    {
        dt_org = dt_int;
    }
    else if(dt_org==dt_double)
    {
        dt_org = dt_double;
    }
    
    dataset_grd.read( V.Coordinates, dt_org, memspace_grd, dataspace_grd );
    dataset_grd_ifn.read( ifn, dt_org, memspace_grd_ifn, dataspace_grd_ifn );
    dataset_conn.read( V.Connection, dt_int, memspace_conn, dataspace_conn );
    dataset_conn_ief.read( ief, dt_int, memspace_conn_ief, dataspace_conn_ief );
    dataset_conn_ife.read( ife, dt_int, memspace_conn_ife, dataspace_conn_ife );
    dataset_conn_ien.read( ien, dt_int, memspace_conn_ien, dataspace_conn_ien );
    dataset_zdefs.read( zdefs_i, dt_int, memspace_zdefs, dataspace_zdefs );

    for(int i=0;i<dims_zdefs_out[0];i++)
    {
        for(int j=0;j<dims_zdefs_out[1];j++)
        {
            std::cout << zdefs_i[i*dims_zdefs_out[1]+j] << " ";
        }
        std::cout << " " << std::endl;
        
    }
    //grp_bound_dataset.read( V.Boundaries, dt_org, memspace_bound, grp_bound_dataspace );
    //grp_interior_dataset.read( V.Interior, dt_org, memspace_interior, grp_interior_dataspace );
    
    // generate the appropriate face2element and element2face map
    int Nfaces    = dims_conn_ife_out[0];
    int Nelement  = dims_conn_ief_out[0];
    int Nnodes    = dims_grd_out[0];
        
    std::vector<std::vector<int>> F2E(Nfaces);
    std::vector<std::vector<int>> E2F(Nelement);
    
    std::vector<std::vector<int>> F2N(Nfaces);
    std::vector<std::vector<int>> N2F(Nnodes);

    std::vector<std::unordered_set<int>> N2N(Nnodes);
    
    for(int i=0;i<Nfaces;i++)
    {
        for(int j=1;j<dims_grd_ifn_out[1]-3;j++)
        {
            F2N[i].push_back(ifn[i*dims_grd_ifn_out[1]+j]);
            N2F[ifn[i*dims_grd_ifn_out[1]+j]].push_back(i);
        }
    }
    
    for(int i=0;i<F2N.size();i++)
    {
        //std::cout << "new = map ";
        for(int j=0;j<F2N[i].size();j++)
        {
            
            //0-->1 3
            //1-->2 0
            //2-->3 1
            //3-->0 2
            
            
            N2N[F2N[i][0]].insert(F2N[i][1]);
            N2N[F2N[i][0]].insert(F2N[i][3]);
            
            N2N[F2N[i][1]].insert(F2N[i][2]);
            N2N[F2N[i][1]].insert(F2N[i][0]);
            
            N2N[F2N[i][2]].insert(F2N[i][3]);
            N2N[F2N[i][2]].insert(F2N[i][1]);
            
            N2N[F2N[i][3]].insert(F2N[i][0]);
            N2N[F2N[i][3]].insert(F2N[i][2]);
        }
    }
    */
    /*
    for(int i=0;i<N2N.size();i++)
    {
        std::unordered_set<int>::iterator it = N2N[i].begin();
        
        while (it != N2N[i].end())
        {
            // Print the element
            std::cout << (*it) << " ";
            //Increment the iterator
            it++;
        }
    }
    */
    
    /*
    for(int i=0;i<Nelement;i++)
    {
        for(int j=1;j<dims_conn_ief_out[1];j++)
        {
            E2F[i].push_back(fabs(ief[i*dims_conn_ief_out[1]+j]));
            F2E[fabs(ief[i*dims_conn_ief_out[1]+j])].push_back(i);
        }
    }
    
    V.face2element = F2E;
    V.element2face = E2F;
    
    V.face2node = F2N;
    V.node2face = N2F;

    std::vector<std::vector<int>> N2E(Nnodes);
    std::vector<std::vector<int>> E2N(dims_conn_out[0]);
    std::vector<int> ElType(Nelement);
    int tel =0;
    for(int i=0;i<Nelement;i++)
    {
        ElType[i]=ien[i*9+0];
        
        for(int j=0;j<8;j++)
        {
            E2N[i].push_back(ien[i*9+j+1]);
            
            if(ien[i*9+j+1] < Nnodes)
            {
                N2E[ien[i*9+j+1]].push_back(i);
            }
        }
    }
    
    V.node2element = N2E;
    V.element2node = E2N;
    V.ElType = ElType;
    
    file_grd.close();
    file_conn.close();
    //file.close();
    */
//    return V;
//}

#endif
