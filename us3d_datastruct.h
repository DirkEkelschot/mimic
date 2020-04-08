#include <string>
#include <vector>




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

struct ParVar
{
    int size;
    int* nlocs;
    int* offsets;
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
