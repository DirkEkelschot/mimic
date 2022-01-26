#include "adapt_topology.h"
#include "adapt_output.h"

Mesh_Topology::Mesh_Topology(Partition* Pa, MPI_Comm comm)
{
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int nlocElem, start, end, offset, nloc, np, loc_vid, size, rank, lid;
    int vf0, vf1, vf2, vf3, vf4, vf5, vf6, vf7, fid;
    double wi, ds0, ds1 ,ds2, ds3, ds4, ds5, u_po,orient0,orient1,orient2,orient3,orient4,orient5,L0,L1,L2,L3,L4,L5;
    
    //ifn = ifn_in;
    c = comm;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
    
    int Nel = Pa->getGlobalPartition()->getNrow();
    
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> > gE2lV 	= Pa->getGlobElem2LocVerts();
    std::vector<Vert*> locVerts  			= Pa->getLocalVerts();
    std::map<int,std::vector<int> > gE2gV 	= Pa->getGlobElem2GlobVerts();

    std::vector<int> Loc_Elem             	= Pa->getLocElem();
    int nLocElem                          	= Loc_Elem.size();
    std::map<int,int> gV2lV               	= Pa->getGlobalVert2LocalVert();
        
    i_part_map* ifn_part_map    = Pa->getIFNpartmap();
    i_part_map* ief_part_map    = Pa->getIEFpartmap();
    i_part_map* iee_part_map    = Pa->getIEEpartmap();
    i_part_map* if_Nv_part_map  = Pa->getIF_Nvpartmap();
    
    std::vector<int> vijkIDs;
    std::map<int,int> LocElem2Nv = Pa->getLocElem2Nv();
    int tel     = 0;
    std::vector<Vert*> face;
    std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
    double volume = 0.0;
    for(int i=0;i<nLocElem;i++)
    {
        int gEl = Loc_Elem[i];
        int NvEl = LocElem2Nv[gEl];
        vijkIDs = gE2lV[gEl];
        
        double* Pijk = new double[NvEl*3];

        for(int k=0;k<vijkIDs.size();k++)
        {
           loc_vid     = vijkIDs[k];
           Pijk[k*3+0] = locVerts[loc_vid]->x;
           Pijk[k*3+1] = locVerts[loc_vid]->y;
           Pijk[k*3+2] = locVerts[loc_vid]->z;
        }

        Vert* Vijk     = ComputeCentroidCoord(Pijk, vijkIDs.size());
        
        int NfPEl      = LocElem2Nf[gEl];

        if(NfPEl == 6)
        {
            volume  = ComputeVolumeHexCell(Pijk);
        }
        if(NfPEl == 4)
        {
            volume  = ComputeVolumeTetCell(Pijk);
        }
        if(NfPEl == 5)
        {
            volume  = ComputeVolumePrismCell(Pijk);
        }
        
        Vol[gEl] = volume;
        std::set<int> vs;
        std::vector<int> vrts;
        
        for(int s=0;s<NfPEl;s++)
        {
            int adjID = iee_part_map->i_map[gEl][s];
            int Nvadj = LocElem2Nv[adjID];
            
            if(adjID<Nel)
            {
                double* Po  = new double[gE2lV[adjID].size()*3];

                for(int k=0;k<gE2lV[adjID].size();k++)
                {
                    loc_vid     = gE2lV[adjID][k];
                    Po[k*3+0] = locVerts[loc_vid]->x;
                    Po[k*3+1] = locVerts[loc_vid]->y;
                    Po[k*3+2] = locVerts[loc_vid]->z;
                }

                Vert* Vpo = ComputeCentroidCoord(Po,gE2lV[adjID].size());

                double d = sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                                (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                                (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));

                Vec3D* rf = new Vec3D;
                rf->c0    = (Vpo->x-Vijk->x)/d;
                rf->c1    = (Vpo->y-Vijk->y)/d;
                rf->c2    = (Vpo->z-Vijk->z)/d;

                rvector[gEl].push_back(rf);
                dr[gEl].push_back(d);
                delete Vpo;
                delete[] Po;
                
                for(int k=0;k<Nvadj;k++)
                {
                   int gV = gE2gV[adjID][k];
                   if(vs.find(gV)==vs.end())
                   {
                     vs.insert(gV);
                     vrts.push_back(gV);
                   }
                }
                
                int faceid = ief_part_map->i_map[gEl][s];
                Vert* Vface = new Vert;
                int NvPerF = if_Nv_part_map->i_map[faceid][0];
                double* F = new double[NvPerF*3];
                
                for(int r=0;r<NvPerF;r++)
                {
                    int gvid = ifn_part_map->i_map[faceid][r];
                    int lvid = gV2lV[gvid];
                    
                    //vert2ref[gvid] = ref;
                    //ref2vert[ref].push_back(gvid);

                    Vert* V = new Vert;
                    V->x    = locVerts[lvid]->x;
                    V->y    = locVerts[lvid]->y;
                    V->z    = locVerts[lvid]->z;
                    
                    F[r*3+0] = V->x;
                    F[r*3+1] = V->y;
                    F[r*3+2] = V->z;
                    
                    Vface->x = Vface->x+locVerts[lvid]->x;
                    Vface->y = Vface->y+locVerts[lvid]->y;
                    Vface->z = Vface->z+locVerts[lvid]->z;

                    face.push_back(V);
                }
                
                if(NvPerF==3) // triangle
                {
                    Vface->x = Vface->x/NvPerF;
                    Vface->y = Vface->y/NvPerF;
                    Vface->z = Vface->z/NvPerF;
                    
                    Vec3D* r0 = new Vec3D;
                    //double Lr = ComputeEdgeLength(Vface,Vijk);

                    r0->c0 = (Vface->x-Vijk->x);///Lr;
                    r0->c1 = (Vface->y-Vijk->y);///Lr;
                    r0->c2 = (Vface->z-Vijk->z);///Lr;
                    
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[2]->x-face[0]->x;
                    v1->c1 = face[2]->y-face[0]->y;
                    v1->c2 = face[2]->z-face[0]->z;
                    
                    Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    vfacevector[gEl].push_back(Vface);
                    ds0 = ComputeTriSurfaceArea(F);
                    dS[gEl].push_back(ds0);
                    normals[gEl].push_back(n0);
                    dxfxc[gEl].push_back(r0);
                    
                }
                if(NvPerF==4) // quad
                {
                    Vface->x = Vface->x/NvPerF;
                    Vface->y = Vface->y/NvPerF;
                    Vface->z = Vface->z/NvPerF;
                    
                    Vec3D* r0 = new Vec3D;
                    r0->c0 = (Vface->x-Vijk->x);///Lr
                    r0->c1 = (Vface->y-Vijk->y);///Lr
                    r0->c2 = (Vface->z-Vijk->z);///Lr
                    
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[3]->x-face[0]->x;
                    v1->c1 = face[3]->y-face[0]->y;
                    v1->c2 = face[3]->z-face[0]->z;
                    
                    Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    vfacevector[gEl].push_back(Vface);
                    ds0 = ComputeQuadSurfaceArea(F);
                    dS[gEl].push_back(ds0);
                    normals[gEl].push_back(n0);
                    dxfxc[gEl].push_back(r0);
                    
                }
                
                //rvector[gEl].push_back(rf);

                //delete Vface;
                delete[] F;
                face.clear();
            }
            else // If boundary face then search data in the correct ghost cell;
            {
                fid = ief_part_map->i_map[gEl][s];
                Vert* Vface = new Vert;
                Vface->x = 0.0;Vface->y = 0.0;Vface->z = 0.0;

                int NvPerF = if_Nv_part_map->i_map[fid][0];
                double* F = new double[NvPerF*3];
                for(int s=0;s<NvPerF;s++)
                {
                    //int gvid = ifn->getVal(fid,s);
                    int gvid = ifn_part_map->i_map[fid][s];

                    int lvid = gV2lV[gvid];

                    Vface->x = Vface->x+locVerts[lvid]->x;
                    Vface->y = Vface->y+locVerts[lvid]->y;
                    Vface->z = Vface->z+locVerts[lvid]->z;
                    
                    Vert* V = new Vert;
                    V->x    = locVerts[lvid]->x;
                    V->y    = locVerts[lvid]->y;
                    V->z    = locVerts[lvid]->z;
                    
                    F[s*3+0] = V->x;
                    F[s*3+1] = V->y;
                    F[s*3+2] = V->z;
                    
                    face.push_back(V);
                    
                }


                if(NvPerF==3)
                {
                    Vface->x = Vface->x/3.0;
                    Vface->y = Vface->y/3.0;
                    Vface->z = Vface->z/3.0;

                    double d = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)+
                                 (Vface->y-Vijk->y)*(Vface->y-Vijk->y)+
                                 (Vface->z-Vijk->z)*(Vface->z-Vijk->z));

                    Vec3D* rff = new Vec3D;
                    rff->c0    = (Vface->x-Vijk->x)/d;
                    rff->c1    = (Vface->y-Vijk->y)/d;
                    rff->c2    = (Vface->z-Vijk->z)/d;
                    
                    Vec3D* r0 = new Vec3D;
                    //double Lr = ComputeEdgeLength(Vface,Vijk);

                    r0->c0 = (Vface->x-Vijk->x);///Lr;
                    r0->c1 = (Vface->y-Vijk->y);///Lr;
                    r0->c2 = (Vface->z-Vijk->z);///Lr;
                    
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[2]->x-face[0]->x;
                    v1->c1 = face[2]->y-face[0]->y;
                    v1->c2 = face[2]->z-face[0]->z;
                    
                    Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    Vert* npos = new Vert;
                    
                    npos->x = Vijk->x+2.0*(r0->c0);
                    npos->y = Vijk->y+2.0*(r0->c1);
                    npos->z = Vijk->z+2.0*(r0->c2);
                    
                    vfacevector[gEl].push_back(npos);
                    rvector[gEl].push_back(rff);
                    dr[gEl].push_back(d*2.0);
                    normals[gEl].push_back(n0);
                    dxfxc[gEl].push_back(r0);
                    //delete rf;
                }

                if(NvPerF==4)
                {
                    Vface->x = Vface->x/4.0;
                    Vface->y = Vface->y/4.0;
                    Vface->z = Vface->z/4.0;

                    double d = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)+
                                 (Vface->y-Vijk->y)*(Vface->y-Vijk->y)+
                                 (Vface->z-Vijk->z)*(Vface->z-Vijk->z));

                    Vec3D* rff = new Vec3D;
                    rff->c0    = (Vface->x-Vijk->x)/d;
                    rff->c1    = (Vface->y-Vijk->y)/d;
                    rff->c2    = (Vface->z-Vijk->z)/d;
                    
                    Vec3D* r0 = new Vec3D;
                    r0->c0 = (Vface->x-Vijk->x);///Lr
                    r0->c1 = (Vface->y-Vijk->y);///Lr
                    r0->c2 = (Vface->z-Vijk->z);///Lr
                    
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[3]->x-face[0]->x;
                    v1->c1 = face[3]->y-face[0]->y;
                    v1->c2 = face[3]->z-face[0]->z;
                    
                    Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    Vert* npos = new Vert;
                    
                    npos->x = Vijk->x+2.0*(r0->c0);
                    npos->y = Vijk->y+2.0*(r0->c1);
                    npos->z = Vijk->z+2.0*(r0->c2);
                    
                    vfacevector[gEl].push_back(npos);
                    rvector[gEl].push_back(rff);
                    dr[gEl].push_back(d*2.0);
                    normals[gEl].push_back(n0);
                    dxfxc[gEl].push_back(r0);
                    //delete rf;

                }
            }
        }
        vs.clear();
        E2V_scheme[gEl] = vrts;
        delete[] Pijk;
        vijkIDs.clear();
    }
    
//    delete ifn_part_map;
//    delete ief_part_map;
//    delete iee_part_map;
//    delete if_Nv_part_map;
//    delete if_ref_part_map;
//    
//    delete v0,v1,fc0,fc1,fc2,fc3,fc4,fc5;
//    gE2lV.clear();
//    locVerts.clear();
//    gE2gV.clear();
//    gV2lV.clear();
//    Loc_Elem.clear();
//    LocElem2Nv.clear();
//    LocElem2Nf.clear();
    
    delete v0;
    delete v1;
    
}


Mesh_Topology::~Mesh_Topology()
{
    
    //E2V_scheme.clear();
	
    std::map<int,vector<Vec3D*> >::iterator itdes;
    for(itdes=normals.begin();itdes!=normals.end();itdes++)
    {
        for(int q=0;q<itdes->second.size();q++)
        {
            delete itdes->second[q];
        }
    }
//    
    for(itdes=rvector.begin();itdes!=rvector.end();itdes++)
    {
        for(int q=0;q<itdes->second.size();q++)
        {
            delete itdes->second[q];
        }
    }
//    
    for(itdes=dxfxc.begin();itdes!=dxfxc.end();itdes++)
    {
        for(int q=0;q<itdes->second.size();q++)
        {
            delete itdes->second[q];
        }
    }
    
    std::map<int,std::vector<double> >::iterator itdesVecDouble;
    for(itdesVecDouble=dr.begin();itdesVecDouble!=dr.end();itdesVecDouble++)
    {
        itdesVecDouble->second.clear();
    }

    for(itdesVecDouble=dS.begin();itdesVecDouble!=dS.end();itdesVecDouble++)
    {
        itdesVecDouble->second.clear();
    }
//
    face2ref.clear();
    
    std::map<int,std::vector<int> >::iterator itdesVecInt;
    for(itdesVecInt=ref2face.begin();itdesVecInt!=ref2face.end();itdesVecInt++)
    {
        itdesVecInt->second.clear();
    }
    for(itdesVecInt=ref2vert.begin();itdesVecInt!=ref2vert.end();itdesVecInt++)
    {
        itdesVecInt->second.clear();
    }
}



std::map<int,std::vector<Vert*> > Mesh_Topology::getVfacevector()
{
    return vfacevector;
}
std::map<int,std::vector<int> > Mesh_Topology::getScheme_E2V()
{
    return E2V_scheme;
}

std::map<int,vector<Vec3D*> > Mesh_Topology::getNormals()
{
    return normals;
}
std::map<int,vector<Vec3D*> > Mesh_Topology::getRvectors()
{
    return rvector;
}
std::map<int,vector<Vec3D*> > Mesh_Topology::getdXfXc()
{
    return dxfxc;
}
std::map<int,vector<double> > Mesh_Topology::getdr()
{
    return dr;
}
std::map<int,vector<double> > Mesh_Topology::getdS()
{
    return dS;
}
std::map<int,double > Mesh_Topology::getVol()
{
    return Vol;
}
Array<int>* Mesh_Topology::getIFN()
{
    return ifn;
}
std::map<int,int> Mesh_Topology::getFace2Ref()
{
    return face2ref;
}
std::map<int,std::vector<int> > Mesh_Topology::getRef2Face()
{
    return ref2face;
}
std::map<int,int> Mesh_Topology::getVert2Ref()
{
    return vert2ref;
}
std::map<int,std::vector<int> > Mesh_Topology::getRef2Vert()
{
    return ref2vert;
}
Mesh_Topology_BL* Mesh_Topology::getBLMeshTopology()
{
    return mesh_topo_bl;
}
