#include "adapt.h"
#include "adapt_array.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#ifndef ADAPT_DATASTRUCT_H
#define ADAPT_DATASTRUCT_H



struct Inputs{
    double hgrad;
    double hmin;
    double hmax;
    double MetScale;
    double hausd;
    int ReadFromStats;
    int RunWakRefinement;
    double hwake;
    int niter;
    int recursive;
    int extended;
    int StateVar;
};


struct Element{
    
    int globID;
    std::vector<int> GlobalNodes;
    
    std::map<int,std::vector<int> > GlobalFace2GlobalNode;
    std::map<int,std::vector<int> > GlobalFace2LocalNode;
    
    std::map<int,std::vector<int> > LocalFace2GlobalNode;
    std::map<int,std::vector<int> > LocalFace2LocalNode;
};

struct Vert
{
    //double x=0.0;
    //double y=0.0;
    //double z=0.0;
    double x;
    double y;
    double z;
};


// struct Mesh{
//     Array<double>* xcn;
//     Array<int>* ief;
//     Array<int>* ien;
//     Array<int>* iet;
//     Array<int>* if_ref;
//     Array<int>* ifn;
//     Array<int>* ife;
//     Array<int>* if_Nv;
// };



struct US3D{
    
    ParArray<double>* xcn;
    Array<int>* elTypes;
    ParArray<int>* ien;
    ParArray<int>* ief;
    ParArray<int>* iee;
    ParArray<int>* iet;
    ParArray<int>* ie_Nv;
    ParArray<int>* ie_Nf;
    ParArray<int>* if_Nv;
    
    ParArray<int>* ifn;
    ParArray<int>* ife;
    ParArray<int>* if_ref;
    Array<int>* ie_tetCnt;
    std::map<std::set<int>,int> tria_ref_map;
    std::map<std::set<int>,int> quad_ref_map;
    std::map<int,int> vert_ref_map;
    
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,int> zone2bcref;
    std::map<int,char*> zone2name;
    ParArray<double>* interior;
    Array<double>* ghost;
    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;
    std::map<int,std::vector<int> > bref2zone;
    Array<char>* znames;
    Array<int>* zdefs;
    std::vector<int> bnd_m;
    int* bnd_map;
    std::map<int,std::vector<int> > face_map_gen;
    std::map<int,std::vector<int> > bnd_face_map;
    int nBnd;
};


struct mesh{
    
    int nElem;
    int nFace;
    int nVert;
    std::map<int,std::vector<double> > xcn;
    std::vector<int> elTypes;
    std::map<int,std::vector<int> > ien;
    std::map<int,std::vector<int> > ief;
    std::map<int,std::vector<int> > iee;
    std::map<int,std::vector<int> > ifn;
    std::map<int,std::vector<int> > ife;

    std::map<int,int> ie_Nv;
    std::map<int,int> ie_Nf;
    std::map<int,std::vector<int> > if_Nv;
    std::map<int,std::vector<int> > if_ref;
    std::map<int,int> iet;
    std::map<int,std::vector<double> > interior;
    std::map<int,std::vector<double> > ghost;

    std::vector<int> element2rank;
    std::vector<int> ie_tetCnt;

    std::map<std::set<int>,int> tria_ref_map;
    std::map<std::set<int>,int> quad_ref_map;
    std::map<int,int> vert_ref_map;
    
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,std::vector<int> > ranges_ref;
    std::map<int,int> zone2bcref;
    std::map<int,char*> zone2name;

    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;

    std::map<int,std::vector<int> > bref2zone;
    std::vector<std::vector<char> > znames;
    std::vector<std::vector<int> > zdefs;
    std::vector<int> bnd_m;
    int* bnd_map;
    std::map<int,std::vector<int> > face_map_gen;
    std::map<int,std::vector<int> > bnd_face_map;
    int nBnd;
    int ntetra;
    int nhexahedral;
    int nprism;
    int npyramid;
};



struct Vec3D
{
    double c0;
    double c1;
    double c2;
};


struct Domain
{
    std::map<int,std::vector<int> > ushell;
    std::map<int,std::vector<double> > ushell_centroid;
    int ncomm;
    std::vector<int> faces_ref;
    std::vector<std::vector<int> > faces_part;
//    int* ntifc;
//    int* color_face;
//    std::vector<int* > ifc_tria_loc;
//    std::vector<int* > ifc_tria_glob;
    std::map<int,std::vector<int> > Elements;
    std::map<int,std::vector<int> > Hexes;
    std::map<int,std::vector<int> > Prisms;
    std::map<int,std::vector<int> > Tetras;
    std::map<int,std::vector<int> > Pyramids;
    std::map<int,std::vector<int> > ref2bcface;
    
    std::map<int,std::vector<int> > GHexes;
    std::map<int,std::vector<int> > GPrisms;
    std::map<int,std::vector<int> > GTetras;
    std::map<int,std::vector<int> > GPyramids;
    Array<int>* LocElem2LocNode;
    std::vector<int> loc_part_verts;
    std::vector<int> glob_part_verts;
    std::map<int,int> gv2lpv;
    std::map<int,int> lv2gpv;
    std::map<int,std::vector<int> > vert2elem;
    std::map<int,int> gv2lpartv;
    std::map<int,int> lpartv2gv;
};

struct i_part_map
{
    std::map<int,std::vector<int> > i_map;
    std::map<int,std::vector<int> > i_inv_map;
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
