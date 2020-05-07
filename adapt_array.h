#include "adapt.h"
#include "adapt_parstate.h"

#ifndef ADAPT_ARRAY_H
#define ADAPT_ARRAY_H

template <typename T> class Array {
    public:
        T *data;
        Array(){}
    
        Array(int r, int c)
        {
            spanArray(r,c);
        }
    
    
        void setVal(int i, int j, T val)
        {
            data[i*ncol+j] = val;
        }
        T getVal(int i, int j)
        {
            return  data[i*ncol+j];
        }
        int getNrow( void )
        {
            return nrow;
        }
        int getNcol( void )
        {
            return ncol;
        }
        void spanArray( int r, int c)
        {
            nrow = r;
            ncol = c;
            int length = nrow*ncol;
            data = new T[length];
        }
        int* getDim()
        {
            int* dim = new int[2];
            dim[0] = nrow;
            dim[1] = ncol;
            return  dim;
        }
    
    private:
        int nrow;
        int ncol;
};
    

template <typename T> class ParArray : public Array<T>
{
    public:
        ParArray(int N, int c, MPI_Comm comm): Array<T>()
        {
            pstate = new ParallelState(N,comm);
            int size;
            MPI_Comm_size(comm, &size);
            int rank;
            MPI_Comm_rank(comm, &rank);
            
            int r = pstate->getNloc(rank);
            
            this->spanArray(r,c);
            
            nglob = N;
        }
        ParallelState* getParallelState( void )
        {
            return pstate;
        }
        int getNglob( void )
        {
            return nglob;
        }
    private:
        ParallelState* pstate;
        int nglob;
};


#endif
