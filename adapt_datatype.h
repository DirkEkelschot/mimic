#include "adapt_array.h"

#ifndef ADAPT_DATATYPE_H
#define ADAPT_DATATYPE_H

struct US3D{
    
    ParArray<double>* xcn;
    
    ParArray<int>* ien;
    ParArray<int>* ief;
    ParArray<int>* iee;
    
    Array<int>* ifn;
    Array<int>* ifn_ref;
    std::map<std::set<int>,int> tria_ref_map;
    std::vector<int*> tria_ref;
    Array<int>* ife;
    
    ParArray<double>* interior;
    Array<double>* ghost;
    
    Array<char>* znames;
    Array<int>* zdefs;
    std::vector<int> bnd_m;
    int* bnd_map;
    int nBnd;
};

struct Vec3D
{
    double c0;
    double c1;
    double c2;
};
#endif
