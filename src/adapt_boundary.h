#include "adapt_array.h"

#ifndef ADAPT_BOUNDARY_H
#define ADAPT_BOUNDARY_H

using namespace std;

class BoundaryMap {
    public:
        BoundaryMap(){};
        BoundaryMap(Array<int>* ifn, Array<int>* if_ref);
        std::map<int,std::vector<int> > getBfaceMap();
        std::map<std::set<int>,int> getTriaRefMap();
        std::map<std::set<int>,int> getQuadRefMap();
        std::map<int,int> getNodeRefMap();
    private:
        std::map<int,std::vector<int> > bnd_face_map;
        std::map<std::set<int>,int> tria_ref_map;
        std::map<std::set<int>,int> quad_ref_map;
        std::map<int,int> node_ref_map;
};

#endif
