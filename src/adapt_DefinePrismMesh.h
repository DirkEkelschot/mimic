#include "adapt_recongrad.h"
#include "adapt_io.h"
#include "adapt_parops.h"
#include "adapt_output.h"
#include "adapt_boundary.h"
#include "adapt_distri_parstate.h"
#include "adapt_redistribute.h"
#include "adapt_geometry.h"

#ifndef ADAPT_DEFINEPRISMMESH_H
#define ADAPT_DEFINEPRISMMESH_H

struct newNumberingNodesFaces
{
    std::map<int,int> internalVertMap;
    std::map<int,std::vector<int> > sharedFace2Node;
    std::map<int,std::vector<int> > shellFace2Node;
    std::map<int,std::vector<int> > shFace2Node;
    std::map<int,std::vector<int> > intFace2Node;
    std::map<int,std::vector<int> > bcFace2Node;
    std::map<int,int> glob2locF;
    std::map<int,int> loc2globF;
    std::map<int,int> tagE2gE;
    std::map<int,int> gE2tagE;
    std::map<int,int> sharedVmap;
    std::map<int,int> glob2locF_sha;
    std::map<int,int> loc2globF_sha;
    std::map<int,std::vector<int> > face2node;
    std::map<int,int > ifref;
    std::map<int,std::vector<int> > ien;
    std::map<int,std::vector<int> > ief;
    DistributedParallelState* dist;
    std::map<int,int> lh;
    std::map<int,int> rh;
    std::map<int,int> lh_sha;
    std::map<int,int> rh_sha;
    std::map<int,std::vector<int> > pbcmap;
    std::map<int,int> tag2ElementID;
    std::map<int,int> localV2tagV;
    std::map<int,int> tagV2localV;
    Array<double>* xcn;
    Array<double>* xcn_int;
    Array<double>* xcn_shared;
    int nTotUniqueNonShellVerts;
    std::map<int,int> sharedvert2rank;
    std::map<int,Vert*> tag2coords;
    std::map<int,int> local2globalVertMap;
    std::map<int,int> local2globalVertMapReal;
    Array<int>* iet;
    std::map<int,int> SharedVertsOwned;
    std::map<int,int> NonSharedVertsOwned;
    std::map<int,int> SharedVertsNotOwned;
    std::map<int,int> tag2glob_prism;
    std::map<int,int> glob2tag_prism;
};

double CheckFaceOrientation(Vert* VcF, std::vector<Vert*> Vface, Vert* Vijk);

newNumberingNodesFaces* DetermineNewNumberingOfElementSubset_Test2(Array<int>* part_global,
                                          std::map<int,std::vector<int> > elements,
                                          std::map<int,std::vector<int> > ief_part_map,
                                          std::map<int,std::vector<int> > ifn_part_map,
                                          std::map<int,std::vector<int> > ife_part_map,
                                          std::map<int,std::vector<int> > if_ref_part_map,
                                          std::map<int,std::vector<int> > if_Nv_part_map,
                                          std::map<int,std::vector<int> > ushell,
                                          std::map<int,int> tag2locV,
                                          std::vector<Vert*> locVerts,
                                          std::map<int,int> shellvert2ref_glob,
                                                                   MPI_Comm comm );

newNumberingNodesFaces* DetermineNewNumberingOfElementSubset_Test(Array<int>* part_global,
                                          std::map<int,std::vector<int> > elements,
                                          std::map<int,std::vector<int> > ief_part_map,
                                          std::map<int,std::vector<int> > ifn_part_map,
                                          std::map<int,std::vector<int> > ife_part_map,
                                          std::map<int,std::vector<int> > if_ref_part_map,
                                          std::map<int,std::vector<int> > if_Nv_part_map,
                                          std::map<int,std::vector<int> > ushell,
                                          std::map<int,int> tag2locV,
                                          std::vector<Vert*> locVerts,
                                          std::map<int,int> shellvert2ref_glob,
                                                                  MPI_Comm comm );

newNumberingNodesFaces* DetermineNewNumberingOfElementSubset(Array<int>* part_global,
                                          std::map<int,std::vector<int> > elements,
                                          std::map<int,std::vector<int> > ief_part_map,
                                          std::map<int,std::vector<int> > ifn_part_map,
                                          std::map<int,std::vector<int> > ife_part_map,
                                          std::map<int,std::vector<int> > if_ref_part_map,
                                          std::map<int,std::vector<int> > if_Nv_part_map,
                                          std::map<int,std::vector<int> > ushell,
                                          std::map<int,int> tag2locV,
                                          std::vector<Vert*> locVerts,
                                          std::map<int,int> shellvert2ref_glob,
                                                             MPI_Comm comm );

#endif

