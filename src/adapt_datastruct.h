#include "adapt.h"
#include "adapt_array.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#ifndef ADAPT_DATASTRUCT_H
#define ADAPT_DATASTRUCT_H

struct Domain
{
    Array<int>* LocElem2LocNode;
    std::vector<int> loc_part_verts;
    std::vector<int> glob_part_verts;
    std::map<int,int> gv2lpv;
    std::map<int,int> lv2gpv;
    
    std::map<int,std::vector<int> > vert2elem;
    
};



struct i_part_map
{
    std::map<int,std::vector<int> > i_map;
    std::map<int,std::vector<int> > i_inv_map;
};

struct ParVar
{
    int size;
    int* nlocs;
    int* offsets;
};

struct ParArrayOnRoot
{
    int size;
    int* nlocs;
    int* offsets;
    int* data;
    int length;
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


struct Vert
{
    double x=0.0;
    double y=0.0;
    double z=0.0;
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




struct Partition_old
{
    
    int ndim;
    
    std::map<int, int> loc2glob_Vmap;
    std::map<int, int> glob2loc_Vmap;
    
    Array<int>* loc2glob_Varr;
    Array<int>* glob2loc_Varr;
    
    Array<double>* Verts;
    Array<int>* ien;
    
    int* xadj;
    int* adjncy;
};



struct PartitionStruct
{
    int ndim;
    std::vector<Vert> Verts;
    Array<int>* loc_elem2verts_glob;
    Array<int>* loc_elem2verts_loc;
    std::map<int,int> v_loc2glob;
    std::map<int,int> v_glob2loc;
    int* xadj;
    int* adjncy;
    Array<double>* rho_elem;
    Array<double>* rho_vert;
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


/*
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
*/

#endif
