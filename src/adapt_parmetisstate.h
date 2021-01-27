#include "adapt.h"
#ifndef ADAPT_PARMETISSTATE_H
#define ADAPT_PARMETISSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class ParallelState_Parmetis {
   public:
    ParallelState_Parmetis(ParArray<int>* e2n, MPI_Comm comm, int type);
    int* getNlocs( void );
    int* getElmdist( void );
    int getNloc( int rank );
    int getElmdistAtRank (int rank );
    int getNpolocAtRank (int rank );
    int getNel( void );
    int* getNpolocs( void );
    int* getEptr( void );
    int* getEind( void );
      
   private:
      int  Nel;
      int* elmdist;
      int* nlocs;
      int* npo_locs;
      int* eptr;
      int* eind;
};

inline ParallelState_Parmetis::ParallelState_Parmetis(ParArray<int>* e2n, MPI_Comm comm, int type)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int Nel = e2n->getNglob();
    //std::cout << "number of elements = " << Nel;
    int nloc             = int(Nel/size) + ( rank < Nel%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(Nel/size) + MIN(rank, Nel%size);
    
    int npo_loc = 0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += type;
    }
    
    int* nlocs_tmp         = new int[size];
    nlocs             = new int[size];
    int* npo_locs_tmp      = new int[size];
    npo_locs          = new int[size];

    for(int i=0;i<size;i++)
    {
        nlocs[i]        = 0;
        npo_locs[i]     = 0;
        
        if(i==rank)
        {
            nlocs_tmp[i]        = nloc;
            npo_locs_tmp[i]     = npo_loc;
        }
        else
        {
            nlocs_tmp[i]        = 0;
            npo_locs_tmp[i]     = 0;
        }
    }
    
    int* elm_dist_tmp          = new int[size+1];
    int* npo_offset_tmp        = new int[size+1];
    elmdist                    = new int[size+1];
    int* npo_offset            = new int[size+1];
    
    for(int i=0;i<size+1;i++)
    {
        elmdist[i]   = 0;
        npo_offset[i] = 0;
        if(i==rank)
        {
            elm_dist_tmp[i]   = offset;
            npo_offset_tmp[i] = offset*type;
        }
        else
        {
            elm_dist_tmp[i]  = 0;
            npo_offset_tmp[i] = 0;
        }
    }
    
    
    MPI_Allreduce(nlocs_tmp,        nlocs,      size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs_tmp,     npo_locs,   size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(elm_dist_tmp,     elmdist,    size+1,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_offset_tmp,   npo_offset, size+1,   MPI_INT, MPI_SUM, comm);

    elmdist[size] = Nel;
    npo_offset[size] = Nel*type;

    eptr = new int[nloc+1];
    eind = new int[npo_loc];
    eptr[0]  = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+type;
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = e2n->data[j];
        }
    }
//    eind = e2n->data;
    
    // The constructor builds the following arrays:
    // elmdist
    // nlocs
    // npo_locs
    // npo_offset
    // eptr
    // eind
    
}// This is the constructor

inline int* ParallelState_Parmetis::getElmdist( void )
{
    return elmdist;
}

inline int* ParallelState_Parmetis::getNlocs( void )
{
    return nlocs;
}

inline int* ParallelState_Parmetis::getNpolocs( void )
{
    return npo_locs;
}

inline int* ParallelState_Parmetis::getEind( void )
{
    return eind;
}

inline int* ParallelState_Parmetis::getEptr( void )
{
    return eptr;
}

inline int ParallelState_Parmetis::getElmdistAtRank( int rank )
{
    return elmdist[rank];
}

inline int ParallelState_Parmetis::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

inline int ParallelState_Parmetis::getNpolocAtRank( int rank )
{
    return npo_locs[rank];
}


inline int ParallelState_Parmetis::getNel( void )
{
  return Nel;
}

#endif
