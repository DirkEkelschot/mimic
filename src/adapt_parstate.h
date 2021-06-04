#include "adapt.h"
#ifndef ADAPT_PARSTATE_H
#define ADAPT_PARSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class ParallelState {
   public:
    ParallelState(int N, MPI_Comm c, int world_rank, int world_size);
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

inline ParallelState::ParallelState(int N, MPI_Comm c, int rank, int size)
{
    std::cout << "here 1 in ParallelState::ParallelState" << std::endl;
    Nel  = N;
    comm = c;

    int nloc             = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(N/size) + MIN(rank, N%size);
    
    std::cout << "here 2 in ParallelState::ParallelState" << std::endl;
    int* proc_nlocs                 = new int[size];
    int* proc_offset                = new int[size];
    nlocs                           = new int[size];
    offsets                         = new int[size];
    for(int i=0;i<size;i++)
    {
        nlocs[i]   = 0;
        offsets[i] = 0;
         
        if(i==rank)
        {
            proc_nlocs[i]  = nloc;
            proc_offset[i] = offset;
        }
        else
        {
            proc_nlocs[i]  = 0;
            proc_offset[i] = 0;
        }
    }
    std::cout << "here 3 in ParallelState::ParallelState";
    std::cout << " size=" << size << " rank=" << rank << std::endl;
    std::cout << "proc_nlocs=" <<proc_nlocs<< " proc_offset" << proc_offset << std::endl;
    MPI_Allreduce(proc_nlocs,  nlocs,   size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(proc_offset, offsets, size, MPI_INT, MPI_SUM, comm);
    std::cout << "here 4 in ParallelState::ParallelState" << std::endl;
     
}// This is the constructor

inline int* ParallelState::getOffsets( void )
{
    return offsets;
}

inline int* ParallelState::getNlocs( void )
{
    return nlocs;
}

inline int ParallelState::getOffset( int rank )
{
    return offsets[rank];
}

inline int ParallelState::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

inline int ParallelState::getNel( void )
{
  return Nel;
}

#endif
