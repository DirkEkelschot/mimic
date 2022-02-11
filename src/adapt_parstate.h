#include "adapt.h"
#ifndef ADAPT_PARSTATE_H
#define ADAPT_PARSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

class ParallelState {
   public:
    ParallelState(){};
    ParallelState(int N, MPI_Comm c);
    ~ParallelState();
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
