#include "adapt_geometry.h"
#include "adapt_partition.h"
#include "adapt_geometry.h"
#include "adapt_compute.h"

#ifndef ADAPT_TOPOLOGY_H
#define ADAPT_TOPOLOGY_H

using namespace std;

struct Element{
    
    int globID;
    std::vector<int> GlobalNodes;
    
    std::map<int,std::vector<int> > GlobalFace2GlobalNode;
    std::map<int,std::vector<int> > GlobalFace2LocalNode;
    
    std::map<int,std::vector<int> > LocalFace2GlobalNode;
    std::map<int,std::vector<int> > LocalFace2LocalNode;
};


struct Mesh_Topology_BL{
    std::map<int,std::vector<int> > BLlayers;
    std::vector<std::vector<int> > BndFaces;
    int Nprisms;
    std::map<int,std::vector<std::vector<int> > > BLlayersPrisms;
    std::map<int,std::vector<Element*> > BLlayersElements;
    std::map<int,std::vector<int> > GlobalElement2GlobalNode;
    std::map<int,std::vector<int> > LocalElement2GlobalNode;
    std::vector<int> exteriorElIDs;
    std::map<int,vector<int> > exteriorVertIDs;
    std::map<int,vector<double> > exteriorVerts;
    std::map<int,int> verts_g2l_ex;
    std::map<int,std::vector<double> > local_ex_verts;
    std::vector<int> outer_shell_faces;
    std::vector<int> elements;
    std::map<int,std::vector<std::vector<int> > > bcQuad;
    std::map<int,std::vector<std::vector<int> > > bcTria;
};


class Mesh_Topology {
    public:
        Mesh_Topology(){};
        Mesh_Topology(Partition* Pa, MPI_Comm comm);
        void DetermineBoundaryLayerElements(Partition* Pa, Array<int>* ife_in, int nLayer, int bID, MPI_Comm comm);
        std::map<int,vector<Vec3D*> > getNormals();
        std::map<int,vector<Vec3D*> > getRvectors();
        std::map<int,vector<Vec3D*> > getdXfXc();
        std::map<int,vector<double> > getdr();
        std::map<int,vector<double> > getdS();
        Array<int>* getIFN();
        std::map<int,double> getVol();
        std::vector<double> ReduceUToVertices(Domain* dom, std::map<int,double> Uelem);
        std::map<int,int> getFace2Ref();
        std::map<int,std::vector<int> > getRef2Face();
        std::map<int,int> getVert2Ref();
        std::map<int,std::vector<int> > getRef2Vert();
        Mesh_Topology_BL* getBLMeshTopology();
    private:
        Array<double>* cc;
        std::map<int,vector<Vec3D*> > normals;
        std::map<int,vector<Vec3D*> > rvector;
        std::map<int,vector<Vec3D*> > dxfxc;
        std::map<int,vector<double> > dr;
        std::map<int,vector<double> > dS;
        std::map<int,double> Vol;
        Partition* P;
        Array<int>* ifn;
        std::map<int,int> face2ref;
        std::map<int,std::vector<int> > ref2face;
        std::map<int,int> vert2ref;
        std::map<int,std::vector<int> > ref2vert;
        std::map<int,int> Bface2Element;
        std::map<int,int> Bface2LocID;
        std::map<int,Vec3D*> Bface2Normal;
        std::map<int,std::vector<int> > BLlayers; 
        Mesh_Topology_BL* mesh_topo_bl;
        
    
        MPI_Comm c;
};


#endif
