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
            //pstate = new ParallelState(N,comm);
            int size;
            MPI_Comm_size(comm, &size);
            int rank;
            MPI_Comm_rank(comm, &rank);
            
            nloc             = int(N/size) + ( rank < N%size );
            //  compute offset of rows for each proc;
            offset           = rank*int(N/size) + MIN(rank, N%size); 
            //pstate = new ParallelState(N,comm); 
            this->spanArray(nloc,c);
            
            nglob = N;
        }
        int getNglob( void )
        {
            return nglob;
        }
        int getOffset(int rank)
	{
	    return offset;
	}
        int getNloc(int rank)
        {
 	    return nloc;
        }
    private:
        ParallelState* pstate;
        int nloc;
        int offset;
        int nglob;
};


template <typename T> class JagArray {
    public:
        T *data;
        JagArray(){}
    
        JagArray(int r, int* c)
        {
            this->spanJagArray(r,c);
        }
        void setVal(int i, int j, T val)
        {
            data[oset[i]+j] = val;
        }
        T getVal(int i, int j)
        {
            return  data[oset[i]+j];
        }
        int getNrow( void )
        {
            return nrow;
        }
        int getNcol( int i )
        {
            return ncol[i];
        }
        int getRowOffset( int i )
        {
            return oset[i];
        }
        int* getOffsets( void )
        {
            return oset;
        }
        int getLength()
        {
            return length;
        }
        void spanJagArray( int r, int* c)
        {
            
            nrow        =   r;
            ncol        =   c;
            length      =   0;
            oset        = new int[nrow];
            oset[0]     = 0;

            for(int i=0;i<r;i++)
            {
                length = length+c[i];
                if(i<r-1)
                {
                    oset[i+1]=oset[i]+c[i];
                }

            }
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
        int* ncol;
        int* oset;
        int length;
        
};

#endif
