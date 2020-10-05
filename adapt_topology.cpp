#include "adapt_topology.h"
#include "adapt_output.h"
Mesh_Topology::Mesh_Topology(Partition* Pa, Array<int>* ifn_in, std::map<int,double> U, int* bnd_map, std::map<int,std::vector<int> > bnd_face_map, int nBnd, MPI_Comm comm)
{
    int nlocElem, start, end, offset, nloc, np, loc_vid, size, rank, lid;
    int vf0, vf1, vf2, vf3, vf4, vf5, vf6, vf7, fid;
    double wi, ds0, ds1 ,ds2, ds3, ds4, ds5, u_po,orient0,orient1,orient2,orient3,orient4,orient5,L0,L1,L2,L3,L4,L5;
    
    ifn = ifn_in;
    P = Pa;
    c = comm;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
    int Nel = Pa->getGlobalPartition()->getNrow();
    
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    np                                    = 8;
    std::map<int,std::vector<int> > gE2lV = Pa->getGlobElem2LocVerts();
    std::vector<Vert> locVerts            = Pa->getLocalVerts();
    std::map<int,int> gE2lE               = Pa->getGlobalElement2LocalElement();
    std::map<int,int> lE2gE               = Pa->getLocalElement2GlobalElement();
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    
    double* Pijk                          = new double[np*3];
    std::vector<int> Loc_Elem             = Pa->getLocElem();
    int nLocElem                          = Loc_Elem.size();
    std::map<int,int> gV2lV               = Pa->getGlobalVert2LocalVert();
    
    std::vector<int> ElemPart             = Pa->getLocElem();
    i_part_map* ief_part_map = Pa->getIEFpartmap();
    i_part_map* iee_part_map = Pa->getIEEpartmap();
    std::vector<int> vijkIDs;
    
    cc                     = new Array<double>(nLocElem,3);
    double* Po  = new double[np*3];
    int tel     = 0;
    Vert* fc0 = new Vert;
    Vert* fc1 = new Vert;
    Vert* fc2 = new Vert;
    Vert* fc3 = new Vert;
    Vert* fc4 = new Vert;
    Vert* fc5 = new Vert;
    std::vector<Vert*> face;
    
    int ref   = 0;
    
    //bnd_m.push_back(zdefs->getVal(zdefs->getNrow()-1,4));
    int fint = bnd_map[0];
    for(int i=0;i<nLocElem;i++)
    {
        int gEl = Loc_Elem[i];
        
        vijkIDs = gE2lV[gEl];

        for(int k=0;k<vijkIDs.size();k++)
        {
           loc_vid     = vijkIDs[k];
           Pijk[k*3+0] = locVerts[loc_vid].x;
           Pijk[k*3+1] = locVerts[loc_vid].y;
           Pijk[k*3+2] = locVerts[loc_vid].z;
        }

        Vert* Vijk     = ComputeCenterCoord(Pijk, np);
        double volume  = ComputeVolumeHexCell(Pijk);
        Vol[gEl]       = volume;
        
        for(int s=0;s<6;s++)
        {
            int faceid = ief_part_map->i_map[gEl][s];
            
            Vert* Vface = new Vert;
            
            for(int r=0;r<4;r++)
            {
                int gvid = ifn->getVal(faceid,r);
                int lvid = gV2lV[gvid];
                
                vert2ref[gvid] = ref;
                ref2vert[ref].push_back(gvid);
                
                Vert* V = new Vert;
                V->x    = locVerts[lvid].x;
                V->y    = locVerts[lvid].y;
                V->z    = locVerts[lvid].z;
                
                Vface->x = Vface->x+locVerts[lvid].x;
                Vface->y = Vface->y+locVerts[lvid].y;
                Vface->z = Vface->z+locVerts[lvid].z;
                
                face.push_back(V);
            }

            Vface->x = Vface->x/4.0;
            Vface->y = Vface->y/4.0;
            Vface->z = Vface->z/4.0;
            
//            L0 = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)
//                     +(Vface->y-Vijk->y)*(Vface->y-Vijk->y)
//                     +(Vface->z-Vijk->z)*(Vface->z-Vijk->z));
            
            Vec3D* r0 = new Vec3D;
            
            r0->c0 = (Vface->x-Vijk->x);
            r0->c1 = (Vface->y-Vijk->y);
            r0->c2 = (Vface->z-Vijk->z);
            
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
            
            //double orientaft = DotVec3D(r0,n0);
            if(faceid<fint) // identify the internal face;
            {
                ref = 0;
                face2ref[faceid] = ref;
                ref2face[ref].push_back(faceid);
            }
            else // identify the boundary interface and based on bnd_map, determine the reference value.
            {
                ref = FindBoundaryID(bnd_map,nBnd,faceid)+1;
                face2ref[faceid]        = ref;
                ref2face[ref].push_back(faceid);
                Bface2Element[faceid]   = gEl;
                Bface2Normal[faceid]    = n0;
                Bface2LocID[faceid]     = s;
            }
            
            ds0 = ComputeSurfaceArea(v0,v1);
            dS[gEl].push_back(ds0);
            normals[gEl].push_back(n0);
            dxfxc[gEl].push_back(r0);
            
            face.clear();
        }
        
        tel = 0;
       
        for(int j=0;j<6;j++)
        {
           int adjID = iee_part_map->i_map[gEl][j];

           if(adjID<Nel)// If internal element;
           {
//               lid  = gE2lE[adjID];
//               u_po = U[adjID];
               
               for(int k=0;k<gE2lV[adjID].size();k++)
               {
                   loc_vid     = gE2lV[adjID][k];
                   Po[k*3+0] = locVerts[loc_vid].x;
                   Po[k*3+1] = locVerts[loc_vid].y;
                   Po[k*3+2] = locVerts[loc_vid].z;
               }
               
               Vert* Vpo = ComputeCenterCoord(Po,8);

               double d = sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                               (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                               (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));
            
               Vec3D* rf = new Vec3D;
               rf->c0    = (Vpo->x-Vijk->x)/d;
               rf->c1    = (Vpo->y-Vijk->y)/d;
               rf->c2    = (Vpo->z-Vijk->z)/d;
               
               rvector[gEl].push_back(rf);
               dr[gEl].push_back(d);
           }
           else // If boundary face then search data in the correct ghost cell;
           {
               fid = ief_part_map->i_map[gEl][j];
               Vert* Vpo = new Vert;
               Vpo->x = 0.0;Vpo->y = 0.0;Vpo->z = 0.0;
               for(int s=0;s<4;s++)
               {
                   int gvid = ifn->getVal(fid,s);
                   int lvid = gV2lV[gvid];

                   Vpo->x = Vpo->x+locVerts[lvid].x;
                   Vpo->y = Vpo->y+locVerts[lvid].y;
                   Vpo->z = Vpo->z+locVerts[lvid].z;
               }

               Vpo->x = Vpo->x/4.0;
               Vpo->y = Vpo->y/4.0;
               Vpo->z = Vpo->z/4.0;

               double d = 2.0*sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                            (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                            (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));
                              
               Vec3D* rf = new Vec3D;
               rf->c0    = (Vpo->x-Vijk->x)/d;
               rf->c1    = (Vpo->y-Vijk->y)/d;
               rf->c2    = (Vpo->z-Vijk->z)/d;

               rvector[gEl].push_back(rf);
               dr[gEl].push_back(d);
               delete Vpo;
           }
           tel++;
        }
    }
    
    
    std::map<int,std::vector<int> > layers = DetermineBoundaryLayerElements(Pa,10,4,c);
    
    if(layers.size()!=0)
    {
        std::vector<int> elements;
        std::set<int> un_elements;

        std::map<int,std::vector<int> >::iterator itt;

        for(itt=layers.begin();itt!=layers.end();itt++)
        {
            for(int q=0;q<itt->second.size();q++)
            {
                if(un_elements.find(itt->second[q])==un_elements.end())
                {
                    un_elements.insert(itt->second[q]);
                    elements.push_back(itt->second[q]);
                }
            }
        }

        OutputBLElements(Pa, elements, c);
    }
    
    
    
//
//    if(rank == 0)
//    {
//        std::map<int,std::vector<int> >::iterator itt;
//
//        for(itt=layers.begin();itt!=layers.end();itt++)
//        {
//            std::cout << itt->first << " -> ";
//            for(int q=0;q<itt->second.size();q++)
//            {
//                std::cout << itt->second[q] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
//
    
    delete[] Pijk;
    delete[] Po;
    
    delete v0;
    delete v1;
    
    delete fc0;
    delete fc1;
    delete fc2;
    delete fc3;
    delete fc4;
    delete fc5;
    
}

std::map<int,std::vector<int> > Mesh_Topology::DetermineBoundaryLayerElements(Partition* Pa, int nLayer, int bID, MPI_Comm comm)
{
    int size;
    int rank;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    int fid_start,gEl,t,gElnew,min_index,lid,fid_new,lEl;
    double newf_id;
    std::vector<double> dp(5);
    std::map<int,std::vector<int> > BLelements;
    std::vector<int> gElvec;
    std::map<int,std::vector<int> > gF2gE = Pa->getglobFace2GlobalElements();
    i_part_map* gE2gF = Pa->getIEFpartmap();
    std::map<int,int> gE2lE               = Pa->getGlobalElement2LocalElement();
    std::vector<std::vector<int> > loc_elem2verts_loc = Pa->getLocalElem2LocalVert();
    std::vector<Vert> LVerts                          =  Pa->getLocalVerts();

    //std::map<int,std::vector<int> >::iterator it;
    Array<int>* gPart = Pa->getGlobalPartition();
    int tel = 0;
    int gElvec0,gElvec1;
    std::vector<int> ElLayer;
    for(int i=0;i<ref2face[bID].size();i++)
    {
        fid_start = ref2face[bID][i];
        
        ;
        Vec3D* nb = Bface2Normal[fid_start];
        //NegateVec3D(nb);
        lid = Bface2LocID[fid_start];
        gEl = Bface2Element[fid_start];
        ElLayer.push_back(gEl);
        for(int k=0;k<nLayer;k++)
        {
            t = 0;
            for(int j=0;j<6;j++)
            {
                Vec3D* nt = normals[gEl][j];
                if(j!=lid)
                {
                    dp[t] = DotVec3D(nb,nt);
                    
                    //std::cout << t << " "  << dp[t]  << std::endl;
                                      
                    t++;
                }
            }
            
            min_index = std::min_element(dp.begin(),dp.end())-dp.begin();
            double min_val = *std::min_element(dp.begin(),dp.end());
            fid_new   = gE2gF->i_map[gEl][min_index];
            nb        = normals[gEl][min_index];
            NegateVec3D(nb);

            if(gE2gF->i_inv_map[fid_new].size()==2)
            {
                gElvec0    = gE2gF->i_inv_map[fid_new][0];
                gElvec1    = gE2gF->i_inv_map[fid_new][1];
                if(gE2gF->i_map.find(gEl)==gE2gF->i_map.end())
                {
                    std::cout << "not in map 1" << std::endl;
                }
                if(gE2gF->i_inv_map.find(fid_new)==gE2gF->i_inv_map.end())
                {
                    std::cout << "not in map 2" << std::endl;
                }
                
                if(gElvec0==gEl)
                {
                    gElnew = gElvec1;
                }
                else if(gElvec1==gEl)
                {
                    gElnew = gElvec0;
                }
                else{
                    std::cout << rank << " "  << gEl << " :: " << gElvec0 << " ->  " << gPart->getVal(gElvec0,0) << " " << gElvec1 << " -> " << gPart->getVal(gElvec1,0) << " " << gE2gF->i_inv_map[fid_new].size() << " not on partition" << std::endl;
                }

                ElLayer.push_back(gElnew);
                gEl=gElnew;
            }
            else
            {
                gElvec0    = gE2gF->i_inv_map[fid_new][0];
                std::cout << rank << " "  << gEl << " :: " << gElvec0 << std::endl;
            
            }
        }
        //std::cout << std::endl;
        BLelements[fid_start] = ElLayer;
        ElLayer.clear();
    }
    
    return BLelements;
}


std::vector<double> Mesh_Topology::ReduceUToVertices(Array<double>* Uelem)
{
    std::vector<double> Uelem_all = P->PartitionAuxilaryData(Uelem, c);
    std::map<int,std::vector<double> > collect_Ui;
    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
    
    for(int i=0;i<Uelem_all.size();i++)
    {
        double uinew = Uelem_all[i];
        int loc_v;
        
        for(int j=0;j<8;j++)
        {
            loc_v = loc_elem2verts_loc[i][j];
            collect_Ui[loc_v].push_back(uinew);
            
            //collect_Vi[loc_v].push_back(Vol[gEl]);
        }
    }
    std::map<int,std::vector<double> >::iterator it_rhos;
    std::vector<double> uivert;
    
    for(it_rhos=collect_Ui.begin();it_rhos!=collect_Ui.end();it_rhos++)
    {
        double sum_u = 0;
        
        for(int q = 0;q<it_rhos->second.size();q++)
        {
            sum_u    = sum_u + it_rhos->second[q];
        }
        uivert.push_back(sum_u/it_rhos->second.size());
        
    }
        
    return uivert;
    
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
