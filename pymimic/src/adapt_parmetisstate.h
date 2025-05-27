#include "adapt.h"
#ifndef ADAPT_PARMETISSTATE_H
#define ADAPT_PARMETISSTATE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class ParallelState_Parmetis {
   public:
    //ParallelState_Parmetis(ParArray<int>* e2n, MPI_Comm comm, int type);
    ParallelState_Parmetis(ParArray<int>* e2n, Array<int>* elTypes, ParArray<int>* iet, MPI_Comm comm);
    int* getNlocs( void );
    int* getElmdist( void );
    int* getElmWgt( void );
    int getNloc( int rank );
    int getElmdistAtRank (int rank );
    int getNpolocAtRank (int rank );
    int getNel( void );
    int* getNpolocs( void );
    int* getEptr( void );
    int* getEind( void );
    int getNcommonNodes(void);
      
   private:
      int  Nel;
      int* elmdist;
      int* nlocs;
      int* npo_locs;
      int* eptr;
      int* eind;
      int ncommonnodes;
      int* elmwgt;
};



inline ParallelState_Parmetis::ParallelState_Parmetis(ParArray<int>* e2n, Array<int>* elTypes, ParArray<int>* ie_Nv, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int Nel = e2n->getNglob();

    int nloc             = int(Nel/size) + ( rank < Nel%size );
    //  compute offset of rows for each proc;

    int npo_loc     = 0;
    int npo_loc_tot = 0;

    if(elTypes->getVal(0,0) == 1 || elTypes->getVal(1,0) == 1)
    {
        ncommonnodes = 3;
    }
    else if(elTypes->getVal(2,0) == 1)
    {
        ncommonnodes = 4;
    }
    
    elmwgt = new int[nloc];
    
    for(int i=0;i<nloc;i++)
    {
        npo_loc += ie_Nv->getVal(i,0);
        if(ie_Nv->getVal(i,0)==4)
        {
            elmwgt[i] = 1;
        }
        if(ie_Nv->getVal(i,0)==6)
        {
            elmwgt[i] = 1;
        }
        if(ie_Nv->getVal(i,0)==8)
        {
            elmwgt[i] = 1;
        }
        
        
    }
    
    MPI_Allreduce(&npo_loc, &npo_loc_tot, 1, MPI_INT, MPI_SUM, comm);

    int* nlocs_tmp          = new int[size];
    nlocs                   = new int[size];
    int* npo_locs_tmp       = new int[size];
    npo_locs                = new int[size];

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

    MPI_Allreduce(nlocs_tmp,        nlocs,      size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs_tmp,     npo_locs,   size,     MPI_INT, MPI_SUM, comm);

    //int* npo_offset            = new int[size+1];
    elmdist                    = new int[size+1];

    int nelOffset = 0;
    int npoOffset = 0;
    for(int i=0;i<size;i++)
    {
        elmdist[i]      = nelOffset;
        npoOffset       = npoOffset+npo_locs[i];
        nelOffset       = nelOffset+nlocs[i];
    }

    elmdist[size]       = Nel;

//    for(int i=0;i<size+1;i++)
//    {
//        std::cout << elmdist[i] << " " << npo_offset[i] << std::endl;
//    }

    eptr     = new int[nloc+1];
    eind     = new int[npo_loc];
    eptr[0]  = 0;
    int k    = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+ie_Nv->getVal(i,0);
        k = 0;
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            //eind[j]  = e2n->data[j];
            eind[j]    = e2n->getVal(i,k);
            k++;
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

inline int* ParallelState_Parmetis::getElmWgt( void )
{
    return elmwgt;
}


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

inline int ParallelState_Parmetis::getNcommonNodes( void )
{
    return ncommonnodes;
}

#endif
