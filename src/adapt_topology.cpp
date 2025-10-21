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
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    int Nel = Pa->getGlobalPartition()->getNrow();
    
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> > gE2lV 	         = Pa->getGlobElem2LocVerts();
    std::vector<std::vector<double> > locVerts       = Pa->getLocalVerts();
    std::map<int,std::vector<int> > gE2gV 	         = Pa->getGlobElem2GlobVerts();

    std::vector<int> Loc_Elem             	         = Pa->getLocElem();
    int nLocElem                          	         = Loc_Elem.size();
    std::map<int,int> gV2lV               	         = Pa->getGlobalVert2LocalVert();
        
    i_part_map* ifn_part_map    = Pa->getIFNpartmap();
    i_part_map* ief_part_map    = Pa->getIEFpartmap();
    i_part_map* iee_part_map    = Pa->getIEEpartmap();
    i_part_map* if_Nv_part_map  = Pa->getIF_Nvpartmap();
    i_part_map* ien_part_map    = Pa->getIENpartmap();

    std::vector<int> vijkIDs;
    std::map<int,int> LocElem2Nv = Pa->getLocElem2Nv();
    int tel     = 0;
    std::vector<std::vector<double> > face;
    std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
    double volume = 0.0;
    std::map<int,std::vector<int> >::iterator iee_it;
    //for(iee_it=iee_part_map->i_map.begin();
//        iee_it!=iee_part_map->i_map.end();
//        iee_it++)
    int nLoc_Elem = Loc_Elem.size();
    
    
    for(int q=0;q<nLoc_Elem;q++)
    {
        //int gEl  = iee_it->first;
        int gEl  = Loc_Elem[q];

        int NvEl  = ien_part_map->i_map[gEl].size();
        int NfPEl = iee_part_map->i_map[gEl].size();

        int nadj_stored  = LocElem2Nf[gEl];
        
        if(NfPEl!=nadj_stored)
        {
            std::cout << "Big error " << std::endl;
        }
        
        std::vector<double> Pijk(NvEl*3);
        for(int k=0;k<NvEl;k++)
        {
           int global_vid = ien_part_map->i_map[gEl][k];
           loc_vid     = gV2lV[global_vid];
           Pijk[k*3+0] = locVerts[loc_vid][0];
           Pijk[k*3+1] = locVerts[loc_vid][1];
           Pijk[k*3+2] = locVerts[loc_vid][2];
        }
        
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk, NvEl);
        
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
        
        for(int o=0;o<NfPEl;o++)
        {
            int adjID = iee_part_map->i_map[gEl][o];
            
            if(adjID<Nel)
            {
                int Nvadj = ien_part_map->i_map[adjID].size();

                //double* Po  = new double[Nvadj*3];
                std::vector<double> Po(Nvadj*3);
                for(int k=0;k<Nvadj;k++)
                {
                    int global_vid  = ien_part_map->i_map[adjID][k];
                    loc_vid         = gV2lV[global_vid];
                    Po[k*3+0]       = locVerts[loc_vid][0];
                    Po[k*3+1]       = locVerts[loc_vid][1];
                    Po[k*3+2]       = locVerts[loc_vid][2];
                }

                std::vector<double> Vpo = ComputeCentroidCoord(Po,Nvadj);

                double d            = sqrt((Vpo[0]-Vijk[0])*(Vpo[0]-Vijk[0])+
                                           (Vpo[1]-Vijk[1])*(Vpo[1]-Vijk[1])+
                                           (Vpo[2]-Vijk[2])*(Vpo[2]-Vijk[2]));
                
                std::vector<double> rf(3);
                
                rf[0]              = (Vpo[0]-Vijk[0])/d;
                rf[1]              = (Vpo[1]-Vijk[1])/d;
                rf[2]              = (Vpo[2]-Vijk[2])/d;

                rvector[gEl].push_back(rf);
                dr[gEl].push_back(d);
            
                //delete[] Po;
                
                for(int k=0;k<Nvadj;k++)
                {
                   int gV = ien_part_map->i_map[adjID][k];
                   if(vs.find(gV)==vs.end())
                   {
                     vs.insert(gV);
                     vrts.push_back(gV);
                   }
                }
                
                
                int faceid = ief_part_map->i_map[gEl][o];
                std::vector<double> Vface(3);
                int NvPerF = if_Nv_part_map->i_map[faceid][0];
                double* F = new double[NvPerF*3];
                
                for(int r=0;r<NvPerF;r++)
                {
                    int gvid = ifn_part_map->i_map[faceid][r];
                    int lvid = gV2lV[gvid];
                    
                    //vert2ref[gvid] = ref;
                    //ref2vert[ref].push_back(gvid);

                    std::vector<double> V(3);
                    V[0]    = locVerts[lvid][0];
                    V[1]    = locVerts[lvid][1];
                    V[2]    = locVerts[lvid][2];
                    
                    F[r*3+0] = V[0];
                    F[r*3+1] = V[1];
                    F[r*3+2] = V[2];
                    
                    Vface[0] = Vface[0]+locVerts[lvid][0];
                    Vface[1] = Vface[1]+locVerts[lvid][1];
                    Vface[2] = Vface[2]+locVerts[lvid][2];

                    face.push_back(V);
                }
                
                if(NvPerF==3) // triangle
                {
                    Vface[0] = Vface[0]/NvPerF;
                    Vface[1] = Vface[1]/NvPerF;
                    Vface[2] = Vface[2]/NvPerF;
                    
                    std::vector<double> r0(3);

                    r0[0] = (Vface[0]-Vijk[0]);///Lr;
                    r0[1] = (Vface[1]-Vijk[1]);///Lr;
                    r0[2] = (Vface[2]-Vijk[2]);///Lr;
                    
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                    
                    std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
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
                    Vface[0] = Vface[0]/NvPerF;
                    Vface[1] = Vface[1]/NvPerF;
                    Vface[2] = Vface[2]/NvPerF;
                
                    std::vector<double> r0(3);
                    r0[0] = (Vface[0]-Vijk[0]);///Lr
                    r0[1] = (Vface[1]-Vijk[1]);///Lr
                    r0[2] = (Vface[2]-Vijk[2]);///Lr
                    
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[3][0]-face[0][0];
                    v1[1] = face[3][1]-face[0][1];
                    v1[2] = face[3][2]-face[0][2];
                    
                    std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
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
                fid = ief_part_map->i_map[gEl][o];
                int NvPerF = if_Nv_part_map->i_map[fid][0];
                
                std::vector<double> Vface(3);
                Vface[0] = 0.0;
                Vface[1] = 0.0;
                Vface[2] = 0.0;
                
                for(int s=0;s<NvPerF;s++)
                {
                    //int gvid = ifn->getVal(fid,s);
                    int gvid = ifn_part_map->i_map[fid][s];
                    int lvid = gV2lV[gvid];

                    Vface[0] = Vface[0]+locVerts[lvid][0];
                    Vface[1] = Vface[1]+locVerts[lvid][1];
                    Vface[2] = Vface[2]+locVerts[lvid][2];
                    
                    std::vector<double> V(3);
                    V[0]    = locVerts[lvid][0];
                    V[1]    = locVerts[lvid][1];
                    V[2]    = locVerts[lvid][2];
                    
                    face.push_back(V);
                }

                Vface[0] = Vface[0]/NvPerF;
                Vface[1] = Vface[1]/NvPerF;
                Vface[2] = Vface[2]/NvPerF;
                
                std::vector<double> r0(3);
                r0[0] = (Vface[0]-Vijk[0]);
                r0[1] = (Vface[1]-Vijk[1]);
                r0[2] = (Vface[2]-Vijk[2]);
                
                double d = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])+
                             (Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])+
                             (Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));

 
                std::vector<double> rff(3);
                rff[0]    = (Vface[0]-Vijk[0])/d;
                rff[1]    = (Vface[1]-Vijk[1])/d;
                rff[2]    = (Vface[2]-Vijk[2])/d;
            
                if(NvPerF==3)
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                }

                if(NvPerF==4)
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[3][0]-face[0][0];
                    v1[1] = face[3][1]-face[0][1];
                    v1[2] = face[3][2]-face[0][2];
                }
                
                std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
                orient0   = DotVec3D(r0,n0);
                
                if(orient0<0.0)
                {
                    NegateVec3D(n0);
                }
                
                double rdotn = DotVec3D(r0,n0);
                
                std::vector<double> reflect(3);
                reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                reflect[2] = r0[2]-2.0*(rdotn)*n0[2];

                std::vector<double> Vghost(3);
                Vghost[0] = Vface[0] - reflect[0];
                Vghost[1] = Vface[1] - reflect[1];
                Vghost[2] = Vface[2] - reflect[2];
                
                ghostVerts[adjID] = Vghost;
                
                std::vector<double> npos(3);
                
                npos[0] = Vface[0];
                npos[1] = Vface[1];
                npos[2] = Vface[2];
                
                vfacevector[gEl].push_back(npos);
                rvector[gEl].push_back(rff);
                dr[gEl].push_back(d);
                normals[gEl].push_back(n0);
                dxfxc[gEl].push_back(r0);
                
                face.clear();
            }
        }
        vs.clear();
        E2V_scheme[gEl] = vrts;
        //delete[] Pijk;
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

}


Mesh_Topology::~Mesh_Topology()
{
    
    //E2V_scheme.clear();
	
//    std::map<int,vector<std::vector<double> > >::iterator itdes;
//    for(itdes=normals.begin();itdes!=normals.end();itdes++)
//    {
//        for(int q=0;q<itdes->second.size();q++)
//        {
//            delete itdes->second[q];
//        }
//    }
////
//    for(itdes=rvector.begin();itdes!=rvector.end();itdes++)
//    {
//        for(int q=0;q<itdes->second.size();q++)
//        {
//            delete itdes->second[q];
//        }
//    }
////
//    for(itdes=dxfxc.begin();itdes!=dxfxc.end();itdes++)
//    {
//        for(int q=0;q<itdes->second.size();q++)
//        {
//            delete itdes->second[q];
//        }
//    }
    
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



std::map<int,std::vector<std::vector<double> > > Mesh_Topology::getVfacevector()
{
    return vfacevector;
}
std::map<int,std::vector<int> > Mesh_Topology::getScheme_E2V()
{
    return E2V_scheme;
}


std::map<int,vector<std::vector<double> > > Mesh_Topology::getNormals()
{
    return normals;
}
std::map<int,vector<std::vector<double> > > Mesh_Topology::getRvectors()
{
    return rvector;
}
std::map<int,vector<std::vector<double> > > Mesh_Topology::getdXfXc()
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
std::map<int,std::vector<double> > Mesh_Topology::getGhostVerts()
{
    return ghostVerts;
}
