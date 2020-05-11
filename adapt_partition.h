#ifndef ADAPT_PARTITION_H
#define ADAPT_PARTITION_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include "adapt.h"
#include "adapt_datastruct.h"
#include "adapt_operations.h"
#include "adapt_parstate.h"
#include "adapt_parmetisstate.h"
#include "adapt_array.h"

struct GathervObject
{
    int size;
    int* nlocs;
    int* offsets;
    int* data;
    int length;
};

std::vector<int> GetAdjacencyForUS3D_V4(ParArray<int>* ief, MPI_Comm comm);

void Example3DPartitioning(MPI_Comm comm);


ParVar* CreateParallelData(int N, MPI_Comm comm);

// This function computes the eptr and eind array which are required for most ParMetis APIs.
// This is done based on the global element2node (e2n) array, the number of elements, the communicator
// and the type of elements we are currently using. For now the type of elements is fixed to an integer
// assuming that we use one type of element throughout the whole mesh. However this needs to become
// an int* in order to allow for hybrid meshes.
// e2n has the Nvert per element stored consecutively for each element. Hence this array is Nel*NvertPerElement long.
Array<int> GatherArrayToAll(Array<int> locarr, MPI_Comm comm);

GathervObject* GetGathervObject(int nloc, MPI_Comm comm);

ParVar_ParMetis* CreateParallelDataParmetis(ParArray<int>* e2n, MPI_Comm comm, int type);

int* GetPartitionInfo(ParArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm);

Array<int>* DeterminePartitionLayout(ParArray<int>* ien, Array<int>* ien_root, MPI_Comm comm);

Partition* CollectVerticesPerRank(ParArray<int>* ien, Array<double>* xcn_r, MPI_Comm comm);

Partition* CollectElementsPerRank(ParArray<int>* ien, Array<int>* ien_root, MPI_Comm comm);

void DivideElements(Array<int>* part_on_root, Array<int>* ien_on_root, MPI_Comm comm);

#endif
