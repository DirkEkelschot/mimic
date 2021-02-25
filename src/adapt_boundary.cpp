#include "adapt_boundary.h"

BoundaryMap::BoundaryMap(Array<int>* ifn, Array<int>* if_ref)
{

    std::set<int> node_ref_set;
    
    std::set<int> tria0;
    std::set<int> tria00;
    std::set<int> tria1;
    std::set<int> tria11;
    
    std::set<int> quad;
    int faceid;
    int nodeid;
    int nrow_ifn = ifn->getNrow();
    //Array<int>* ifn_ref  = new Array<int>(nrow_ifn,1);
    int ref;
    
    std::map<int,std::vector<int> >::iterator bmit;
    
    int r2=0;int r10=0;int r36=0;int r3=0;
    for(int i=0;i<nrow_ifn;i++)
    {
        ref = if_ref->getVal(i,0);
        //ifn_ref->setVal(i,0,ref);
        
        faceid = i;
        if(ref != 2)
        {
            bnd_face_map[ref].push_back(faceid);
        }
        
        for(int j=0;j<4;j++)
        {
            nodeid = ifn->getVal(i,j); // This is actually node ID!!!!
            //ifn_copy->setVal(i,j,ifn_g->getVal(i,j+1)-1);
            
            if(ref!=2)
            {
                if(node_ref_set.find(nodeid)==node_ref_set.end())
                {
                    node_ref_set.insert(nodeid);
                    node_ref_map[nodeid] = if_ref->getVal(i,0);
                }
            }
        }
        
        tria0.insert(ifn->getVal(i,0));
        tria0.insert(ifn->getVal(i,1));
        tria0.insert(ifn->getVal(i,2));
        
        tria00.insert(ifn->getVal(i,0));
        tria00.insert(ifn->getVal(i,2));
        tria00.insert(ifn->getVal(i,3));
        
        tria1.insert(ifn->getVal(i,0));
        tria1.insert(ifn->getVal(i,1));
        tria1.insert(ifn->getVal(i,3));
        
        tria11.insert(ifn->getVal(i,1));
        tria11.insert(ifn->getVal(i,2));
        tria11.insert(ifn->getVal(i,3));
        
        quad.insert(ifn->getVal(i,0));
        quad.insert(ifn->getVal(i,1));
        quad.insert(ifn->getVal(i,2));
        quad.insert(ifn->getVal(i,3));
        
        if(tria_ref_map.find(tria0)==tria_ref_map.end() && ref!=2)
        {
            tria_ref_map[tria0]  = ref;
            tria_ref_map[tria00] = ref;
        }
        if(tria_ref_map.find(tria1)==tria_ref_map.end() && ref!=2)
        {
            tria_ref_map[tria1]  = ref;
            tria_ref_map[tria11] = ref;
        }
        
        if(quad_ref_map.find(quad)==quad_ref_map.end() && ref!=2)
        {
            quad_ref_map[quad] = ref;
        }
        
        tria0.clear();
        tria00.clear();
        tria1.clear();
        tria11.clear();
        quad.clear();
    }
}


std::map<int,std::vector<int> > BoundaryMap::getBfaceMap()
{
    return bnd_face_map;
}

std::map<std::set<int>,int> BoundaryMap::getTriaRefMap()
{
    return tria_ref_map;
}

std::map<std::set<int>,int> BoundaryMap::getQuadRefMap()
{
    return quad_ref_map;
}

std::map<int,int> BoundaryMap::getNodeRefMap()
{
    return node_ref_map;
}


