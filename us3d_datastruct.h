//#include <string>
//#include <vector>
//#include <map>
//#include <iostream>

#include "adapt.h"

struct ParVar
{
    int size;
    int* nlocs;
    int* offsets;
};

struct ParVar_ParMetis
{
    int size;
    int* nlocs;
    int* elmdist;
    int* npo_locs;
    int* npo_offset;
    int* eptr;
    int* eind;
};



template <typename T> class JaggedArray {
private:

public:
    
    int nloc;
    int* ncols;
    int* offset;
    T *data;
    JaggedArray(){}
    
    JaggedArray(int r, int* nc)
    {
        
        nloc = r;
        ncols = nc;
        
        int size = 0;
        for(int i=0;i<r;i++)
        {
            size = size+ncols[i];
        }
        
        data = new T[size];
    
    }
};


template <typename T> class Array {
    private:
    
    public:
    
    int nglob;
    int nloc;
    int offset;
    int ncol;
    
    T *data;
    Array(){}
    
    Array(int r, int c)
    {
        
        nloc = r;
        ncol = c;
        int size = nloc*ncol;
        
        data = new T[size];
        
        /*
        for(int i=0;i<size;i++)
        {
            data[i] = 0;
        }
        */
    }
    
    Array(int r, int c, int o)
    {
        nloc   = r;
        ncol   = c;
        offset = o;
        int size = nloc*ncol;
        
        data = new T[size];
        
        /*
        for(int i=0;i<size;i++)
        {
            data[i] = 0;
        }
        */
        
    }
    
    void setVal(int i, int j, T val)
    {
        data[i*ncol+j] = val;
    }
    T getVal(int i, int j)
    {
        return  data[i*ncol+j];
    }
    int getOffset()
    {
        return  offset;
    }
    int* getDim()
    {
        int* dim = new int[2];
        dim[0] = nloc;
        dim[1] = ncol;
        return  dim;
    }
};
    

template <typename T> class ParallelArray : public Array<T>
{
    public:
        ParVar* pv;
    
    ParallelArray(int r, int c, ParVar* pv_): Array<T>(r,c)
    {
        pv = pv_;
    }
};

    
template<typename T>
struct Array2
{
    T* data;
    int nrow;
    int ncol;
};


struct Vert
{
    double x;
    double y;
    double z;
};




struct TmpStruct
{
    int* data;
    int* offsets;
    int* nlocs;
    
    int* sizing;
    int* offsets_sizing;
    int* nlocs_sizing;
};

struct LocalPartitionData
{
    std::map<int, int> loc2glob_el;
    std::map<int, int> glob2loc_el;
    
    std::map<int, int> loc2glob_vrt;
    std::map<int, int> glob2loc_vrt;
    
    Array<int>* ien_loc;
};

//template<typename T>
struct US3dData
{
    int rows_grd;
    int cols_grd;
    double* Coordinates;
    
    int rows_conn;
    int cols_conn;
    int* Connection;
    
    int rows_bound;
    int cols_bound;
    double* Boundaries;
    
    int rows_interior;
    int cols_interior;
    double* Interior;

    std::vector<std::vector<int> > element2face;
    std::vector<std::vector<int> > face2element;
    
    std::vector<std::vector<int> > element2node;
    std::vector<std::vector<int> > node2element;
    
    std::vector<std::vector<int> > face2node;
    std::vector<std::vector<int> > node2face;
    
    std::vector<int> ElType;
    
};


template <class T>
void printArray(Array<T> A)
{
    int m = A.nrow;
    int n = A.ncol;
    std::cout << " " << std::endl;
    std::cout << "[";
    for(int i=0;i<m;i++)
    {
        std::cout << "[";
        for(int j=0;j<n;j++)
        {
            std::cout << A.data[i*n+j] << ", ";
        }
        
        if (i == m-1)
        {
            std::cout << "]";
        }
        else
        {
            std::cout << "]," << std::endl;
        }
    }
    std::cout << "]" << std::endl;
    std::cout << " " << std::endl;
};




template <class T>
void PrintUS3dData(US3dData Vec, char tag)
{
    if (tag == 'g')
    {
        std::cout << " ===============Grid Coordinate Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_grd;i++)
        {
            for(int j=0;j<Vec.cols_grd;j++)
            {
                std::cout << Vec.Coordinates[i*Vec.cols_grd+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'c')
    {
        std::cout << " ===============Connection Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_conn;i++)
        {
            for(int j=0;j<Vec.cols_conn;j++)
            {
                std::cout << Vec.Connection[i*Vec.cols_conn+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'b')
    {
        std::cout << " ===============Boundary Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_bound;i++)
        {
            for(int j=0;j<Vec.cols_bound;j++)
            {
                std::cout << Vec.Boundaries[i*Vec.cols_bound+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'i')
    {
        std::cout << " ===============Interior Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_interior;i++)
        {
            for(int j=0;j<Vec.cols_interior;j++)
            {
                std::cout << Vec.Interior[i*Vec.cols_interior+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }


};

