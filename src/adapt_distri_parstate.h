#include "adapt.h"
#ifndef ADAPT_DISTRI_PARSTATE_H
#define ADAPT_DISTRI_PARSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class DistributedParallelState {
   public:
        DistributedParallelState(){};
        DistributedParallelState(int nloc, MPI_Comm comm);
        ~DistributedParallelState();
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



#endif
