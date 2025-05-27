#include "adapt_partition.h"



DistributedParallelState::DistributedParallelState(int nloc, MPI_Comm comm)
{
    int i;
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process;
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int* nlocs_tmp      = new int[world_size];
    nlocs               = new int[world_size];
    offsets             = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        nlocs_tmp[i] = 0;
        nlocs[i]     = 0;

        if(i==world_rank)
        {
            nlocs_tmp[i] = nloc;
        }
        else
        {
            nlocs_tmp[i] = 0;
        }
    }
    
    MPI_Allreduce(nlocs_tmp, nlocs, world_size, MPI_INT, MPI_SUM, comm);

    int o = 0;

    for(i=0;i<world_size;i++)
    {
        offsets[i] = o;
        o          = o+nlocs[i];
    }
    
    Nel = o;
     
}// This is the constructor

DistributedParallelState::~DistributedParallelState()
{
    delete[] nlocs;
    delete[] offsets;
}

int* DistributedParallelState::getOffsets( void )
{
    return offsets;
}

int* DistributedParallelState::getNlocs( void )
{
    return nlocs;
}

int DistributedParallelState::getOffset( int rank )
{
    return offsets[rank];
}

int DistributedParallelState::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

int DistributedParallelState::getNel( void )
{
  return Nel;
}

