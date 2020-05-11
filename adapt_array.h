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


template <typename T> class JagArray {
    public:
        T *data;
        JagArray(){}
    
        JagArray(int r, int* c)
        {
            spanJagArray(r,c);
        }
    
    
        void setVal(int i, int j, T val)
        {
            data[offset[i]+j] = val;
        }
        T getVal(int i, int j)
        {
            return  data[offset[i]+j];
        }
        int getNrow( void )
        {
            return nrow;
        }
        int getNcol( int i )
        {
            return ncol[i];
        }
        void spanJagArray( int r, int* c)
        {
            nrow = r;
            ncol = c;
            int length=0;
            int* offset = new int[r];
            offset[0] = 0;
            
            for(int i=0;i<r;i++)
            {
                length = length+c[i];
                if(i>0)
                {
                    offset[i] = offset[i-1]+c[i];
                }
            }
        std::cout << length << std::endl;
            data = new T[length];
            /*
            data = new T[length];
             */
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
        int* ncol;
    int* offset;
};

#endif
