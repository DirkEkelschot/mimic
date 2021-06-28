#include "adapt.h"
#ifndef ADAPT_DISTRI_PARSTATE_H
#define ADAPT_DISTRI_PARSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class DistributedParallelState {
   public:
    DistributedParallelState(int nloc, MPI_Comm comm);
    int* getOffsets( void );
    int* getNlocs( void );
    int getNloc( int rank );
    int getOffset (int rank );
    int getNel( void );
    
      
   private:
      int Nel;
      MPI_Comm comm;
      int* offsets;
      int* nlocs;
};

inline DistributedParallelState::DistributedParallelState(int nloc, MPI_Comm comm)
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
        o             = o+nlocs[i];
    }
    
    Nel = o;
     
}// This is the constructor

inline int* DistributedParallelState::getOffsets( void )
{
    return offsets;
}

inline int* DistributedParallelState::getNlocs( void )
{
    return nlocs;
}

inline int DistributedParallelState::getOffset( int rank )
{
    return offsets[rank];
}

inline int DistributedParallelState::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

inline int DistributedParallelState::getNel( void )
{
  return Nel;
}

#endif
