#include "adapt_geometry.h"
#include "adapt_partition.h"
#include "adapt_geometry.h"
#include "adapt_compute.h"

#ifndef ADAPT_TOPOLOGY_H
#define ADAPT_TOPOLOGY_H

using namespace std;

class Mesh_Topology {
    public:
        Mesh_Topology(){};
        Mesh_Topology(Partition* Pa, Array<int>* ifn_in, Array<double>* ghost, std::map<int,double> U, MPI_Comm comm);
        std::map<int,vector<Vec3D*> > getNormals();
        std::map<int,vector<Vec3D*> > getRvectors();
        std::map<int,vector<Vec3D*> > getdXfXc();
        std::map<int,vector<double> > getdr();
        std::map<int,vector<double> > getdS();
        Array<int>* getIFN();
        std::map<int,double> getVol();
        std::vector<double> ReduceUToVertices(Array<double>* Uelem);
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
        MPI_Comm c;
};


#endif
