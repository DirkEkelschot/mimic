#include "adapt.h"
#include "adapt_array.h"
#include "adapt_datastruct.h"
#include "adapt_distri_parstate.h"
#include "adapt_schedule.h"
#include "adapt_io.h"
#ifndef ADAPT_DISTRI_PRISMATICLAYER_H
#define ADAPT_DISTRI_PRISMATICLAYER_H
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))



class PrismaticLayer {
   public:
        PrismaticLayer(){};
        PrismaticLayer(std::map<int,std::vector<int> > elements,
                       std::map<int,std::vector<int> > ief_part_map,
                       std::map<int,std::vector<int> > ifn_part_map,
                       std::map<int,std::vector<int> > ife_part_map,
                       std::map<int,int> if_ref_part_map,
                       std::map<int,int> if_Nv_part_map,
					   std::map<int,std::vector<int> > if_Erank_part_map,
                       std::map<int,std::vector<int> > ushell,
                       std::map<int,int> tag2locV,
                       std::vector<Vert*> locVerts,
                       std::map<int,int> shellvert2ref_glob,
					   std::map<int,std::vector<int> > range_id, 
                       MPI_Comm comm);
        
        ~PrismaticLayer();
    
        std::map<int,int> getLeftElementGlobalIDForFaces();
        std::map<int,int> getRightElementGlobalIDForFaces();
    
        std::map<int,std::vector<int> > getBoundaryFace2NodeMap();
        std::map<int,std::vector<int> > getInternalFace2NodeMap();
        std::map<int,std::vector<int> > getOwnedSharedFace2NodeMap();

        std::map<int,int> getGlobal2TagElementMap();
        std::map<int,int> getTag2GlobalElementMap();
    
        Array<double>* getInternalCoordinates();
        Array<double>* getSharedCoordinates();
        
        std::map<int,int> getNotOwnedSharedVerticesMap();
        std::map<int,int> getVertexTag2GlobalMap();
        std::map<int,int> getTag2Element4TetPrismInterface();
        std::map<int,std::vector<int> > getBoundaryCondition2FaceID();
        std::map<int,std::vector<int> > getBoundaryConditionZone2FaceID();
        Array<int>* getElementType();
        std::map<int,int> getSharedVertexMap();
   private:
        std::map<int,std::vector<int> > SharedFace2Node;
        std::map<int,std::vector<int> > InternalFace2Node;
        std::map<int,std::vector<int> > BoundaryFace2Node;
        std::map<int,int> tagE2gE;
        std::map<int,int> gE2tagE;
        std::map<int,int> lhp;
        std::map<int,int> rhp;
        std::map<int,std::vector<int> > ref2bcface;
        std::map<int,std::vector<int> > zone2bcface;
        std::map<int,int> tag2element_TetPrismInterface;
        std::map<int,int> sharedVmap;
        Array<double>* xcn_int;
        Array<double>* xcn_shared;
        std::map<int,int> SharedVertsNotOwned;
        std::map<int,int> tagV2globV;
        std::map<int,int> globV2tagV;
        
        Array<int>* iet;
};



#endif
