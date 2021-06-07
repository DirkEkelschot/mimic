#include "adapt.h"
#include "adapt_partition.h"
#include "adapt_datatype.h"
#include "adapt_io.h"

#ifndef ADAPT_PAROPS_H
#define ADAPT_PAROPS_H

using namespace std;




inline MPI_Datatype get_mpi_datatype(const int &)
{
    return MPI_INT;
}

inline MPI_Datatype get_mpi_datatype(const double &)
{
    return MPI_DOUBLE;
}

template <typename T>
inline MPI_Datatype get_mpi_datatype() {
    return get_mpi_datatype(T());
}

template<typename T>
Array<T>* GatherArrayOnRoot(Array<T>* A,MPI_Comm comm, MPI_Info info)
{
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    MPI_Datatype mpi_dtype = get_mpi_datatype<T>();
    
    int nglob = 0;
    int ncol  = A->getNcol();
    Array<T>* gA;
    
    int* G_nlocs          = new int[size];
    int* red_G_nlocs      = new int[size];
    int* G_offsets        = new int[size];

    for(int i=0;i<size;i++)
    {
        G_nlocs[i] = 0;
        
        if(i==rank)
        {
            G_nlocs[i] = A->getNrow()*A->getNcol();
        }
        else
        {
            G_nlocs[i] = 0;
        }
    }
    MPI_Allreduce(G_nlocs, red_G_nlocs, size, mpi_dtype, MPI_SUM, comm);

    int offset = 0;
    for(int i=0;i<size;i++)
    {
        G_offsets[i] = offset;
        offset       = offset+red_G_nlocs[i];
        nglob        = nglob+red_G_nlocs[i];
    }
    nglob = nglob/A->getNcol();
    
    if(rank == 0)
    {
        gA = new Array<T>(nglob,ncol);
    }
    else
    {
        gA = new Array<T>(1,1);
    }
    
    MPI_Gatherv(&A->data[0],
                A->getNrow()*ncol,
                MPI_INT,
                &gA->data[0],
                red_G_nlocs,
                G_offsets,
                MPI_INT, 0, comm);
    
    return gA;
}


Array<double>* GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, std::map<int,Array<double>*> mv_map, MPI_Comm comm);


Mesh* ReduceMeshToRoot(ParArray<int>* ien,
                       ParArray<int>* ief,
                       ParArray<double>* xcn,
                       ParArray<int>* ifn,
                       ParArray<int>* ife,
                       ParArray<int>* if_ref,
                       MPI_Comm comm, MPI_Info info);

MeshTransfer* GetMeshTransfer(std::map<int,std::vector<int> > Ate, std::vector<int> ufaces_vec, std::vector<int> uverts_vec, MPI_Comm comm, MPI_Info info);

Array<int>* GatherTetrahedraOnRoot(std::map<int,std::vector<int> > Ate, MPI_Comm comm, MPI_Info info);

std::map<int,std::vector<int> > GatherElementsOnRoot(std::map<int,std::vector<int> >Apr, std::map<int,std::vector<int> > Ate, MPI_Comm comm, MPI_Info info);
#endif
