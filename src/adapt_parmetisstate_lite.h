#include "adapt.h"
#ifndef ADAPT_PARMETISSTATE_LITE_H
#define ADAPT_PARMETISSTATE_LITE_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class ParallelState_Parmetis_Lite {
   public:
    //ParallelState_Parmetis_Lite(ParArray<int>* e2n, MPI_Comm comm, int type);
    ParallelState_Parmetis_Lite(std::map<int,std::vector<int> > e2n, 
                           std::vector<int> elTypes, 
                           MPI_Comm comm);
    std::vector<int> getNlocs( void );
    std::vector<int> getElmdist( void );
    std::vector<int> getElmWgt( void );
    int getNloc( int rank );
    int getElmdistAtRank (int rank );
    int getNpolocAtRank (int rank );
    int getNtotalElem( void );
    std::vector<int> getNpolocs( void );
    std::vector<int> getEptr( void );
    std::vector<int> getEind( void );
    int getNcommonNodes(void);
      
   private:
      int  nElemTotal;
      std::vector<int> elmdist;
      std::vector<int> nlocs;
      std::vector<int> npo_locs;
      std::vector<int> eptr;
      std::vector<int> eind;
      int ncommonnodes;
      std::vector<int> elmwgt;
};



inline ParallelState_Parmetis_Lite::ParallelState_Parmetis_Lite(std::map<int,std::vector<int> > e2n, 
                           std::vector<int> elTypes, 
                           MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);

    int nloc             = e2n.size();
    nElemTotal       = 0;
    MPI_Allreduce(&nloc, &nElemTotal, 1, MPI_INT, MPI_SUM, comm);

    int npo_loc     = 0;
    int npo_loc_tot = 0;

    if(elTypes[0] == 1 || elTypes[1] == 1)
    {
        ncommonnodes = 3;
    }
    else if(elTypes[2] == 1)
    {
        ncommonnodes = 4;
    }
    
    
    elmwgt = std::vector<int>(nloc,0);
    std::map<int,std::vector<int> >::iterator itmiv;
    int i = 0;
    for(itmiv=e2n.begin();itmiv!=e2n.end();itmiv++)
    {
        int el_gid = itmiv->first;
        npo_loc   += e2n[el_gid].size();
        
        if(e2n[el_gid].size()==4)
        {
            elmwgt[i] = 1;
        }
        if(e2n[el_gid].size()==6)
        {
            elmwgt[i] = 1;
        }
        if(e2n[el_gid].size()==8)
        {
            elmwgt[i] = 1;
        }
        
        i++;
    }

    
    MPI_Allreduce(&npo_loc, &npo_loc_tot, 1, MPI_INT, MPI_SUM, comm);

    std::vector<int> nlocs_tmp(size,0);
    nlocs = std::vector<int>(size,0);
    std::vector<int> npo_locs_tmp(size,0);
    npo_locs = std::vector<int>(size,0);

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

    MPI_Allreduce(nlocs_tmp.data(),        nlocs.data(),      size,     MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs_tmp.data(),     npo_locs.data(),   size,     MPI_INT, MPI_SUM, comm);

    elmdist = std::vector<int>(size+1,0);

    int nelOffset = 0;
    int npoOffset = 0;
    for(int i=0;i<size;i++)
    {
        elmdist[i]      = nelOffset;
        npoOffset       = npoOffset+npo_locs[i];
        nelOffset       = nelOffset+nlocs[i];
    }

    elmdist[size]       = nElemTotal;

//    for(int i=0;i<size+1;i++)
//    {
//        std::cout << elmdist[i] << " " << npo_offset[i] << std::endl;
//    }

    eptr = std::vector<int>(nloc+1,0);
    eind = std::vector<int>(npo_loc,0);

    eptr[0]  = 0;
    int k    = 0;
    i        = 0;

    for(itmiv=e2n.begin();itmiv!=e2n.end();itmiv++)
    {
        int el_gid = itmiv->first;
        eptr[i+1]  = eptr[i]+e2n[el_gid].size();
        k = 0;

        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j]    = itmiv->second[k];
            k++;
        }

        i++;
    }
}// This is the constructor

inline std::vector<int> ParallelState_Parmetis_Lite::getElmWgt( void )
{
    return elmwgt;
}


inline std::vector<int> ParallelState_Parmetis_Lite::getElmdist( void )
{
    return elmdist;
}

inline std::vector<int> ParallelState_Parmetis_Lite::getNlocs( void )
{
    return nlocs;
}

inline std::vector<int> ParallelState_Parmetis_Lite::getNpolocs( void )
{
    return npo_locs;
}

inline std::vector<int> ParallelState_Parmetis_Lite::getEind( void )
{
    return eind;
}

inline std::vector<int> ParallelState_Parmetis_Lite::getEptr( void )
{
    return eptr;
}

inline int ParallelState_Parmetis_Lite::getElmdistAtRank( int rank )
{
    return elmdist[rank];
}

inline int ParallelState_Parmetis_Lite::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

inline int ParallelState_Parmetis_Lite::getNpolocAtRank( int rank )
{
    return npo_locs[rank];
}


inline int ParallelState_Parmetis_Lite::getNtotalElem( void )
{
  return nElemTotal;
}

inline int ParallelState_Parmetis_Lite::getNcommonNodes( void )
{
    return ncommonnodes;
}

#endif
