#include "adapt_output.h"

using namespace std;


void OutputBoundaryLayerPrisms(Array<double>* xcn_g, Mesh_Topology_BL* BLmesh, MPI_Comm comm,string fname)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::map<int,std::vector<Element* > >::iterator iter;

    std::set<int> unique_prism_verts;
    std::vector<int> u_prism_v;

    Array<int>* local_prisms = new Array<int>(BLmesh->Nprisms,6);
    std::map<int,int> gv2lv_prisms;
    int npr =0;
    int lvid=0;
    std::set<std::set<int> > unique_faces;
    std::vector<std::vector<int> > ufaces;
    int fc = 0;
    std::map<int,int> lhface;
    std::map<int,int> rhface;
    int TotNumFaceNodes = 0;
    int numit = 0;
    std::map<std::set<int>,int> face2ID;
    for(iter=BLmesh->BLlayersElements.begin();iter!=BLmesh->BLlayersElements.end();iter++)
    {
        numit=iter->second.size();
        
        for(int p=0;p<numit;p++)
        {
            
            std::vector<int> prism = iter->second[p]->GlobalNodes;
            
            int gElID = iter->second[p]->globID;
            
            for(int q=0;q<prism.size();q++)
            {
                if(unique_prism_verts.find(prism[q])==unique_prism_verts.end())
                {
                    unique_prism_verts.insert(prism[q]);
                    u_prism_v.push_back(prism[q]);
                    local_prisms->setVal(npr,q,lvid);
                    gv2lv_prisms[prism[q]] = lvid;
                    lvid++;
                }
                else
                {
                    local_prisms->setVal(npr,q,gv2lv_prisms[prism[q]]);
                }
            }
           
            std::map<int,std::vector<int> > LF2GN = iter->second[p]->LocalFace2GlobalNode;
//
            std::map<int,std::vector<int> >::iterator itm;
            for(itm=LF2GN.begin();itm!=LF2GN.end();itm++)
            {
                std::set<int> face;
                for(int q=0;q<itm->second.size();q++)
                {
                    face.insert(itm->second[q]);
                }
                if(unique_faces.find(face)==unique_faces.end())
                {
                    unique_faces.insert(face);
                    std::set<int>::iterator its;
                    std::vector<int> tmp;
                    face2ID[face] = fc;
                    for(int q=0;q<itm->second.size();q++)
                    {
                        int local_id = gv2lv_prisms[itm->second[q]];
                        tmp.push_back(local_id+1); //maintaining face ordering.
                    }
                    lhface[fc]=gElID+1;
                    ufaces.push_back(tmp);
                    TotNumFaceNodes=TotNumFaceNodes+tmp.size();
                    tmp.clear();
                    fc++;
                }
                else
                {
                    int fcid = face2ID[face];
                    rhface[fcid]=gElID+1;
                }
                face.clear();
            }
            npr++;
        }
    }

    std::map<int,int>::iterator itmap;
    int tel = 0;
    for(itmap=lhface.begin();itmap!=lhface.end();itmap++)
    {
        if(rhface.find(itmap->first)==rhface.end())
        {
            rhface[itmap->first] = 0;
        }
        tel++;
    }
//
    string filename = fname + std::to_string(world_rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile << "Zone" << std::endl;
    myfile << "ZoneType=FEPolyhedron" << std::endl;
    myfile <<"Nodes="<< u_prism_v.size() << std::endl;
    myfile <<"Faces="<< ufaces.size() << std::endl;
    myfile <<"Elements="<< npr << std::endl;
    myfile <<"TotalNumFaceNodes="<<TotNumFaceNodes << std::endl;
    myfile << "NumConnectedBoundaryFaces="<< 0 << std::endl;
    myfile << "TotalNumBoundaryConnections="<< 0 << std::endl;
    myfile << "Datapacking=BLOCK" << std::endl;
    myfile << "Varlocation=(4=CellCentered)" << std::endl;
    myfile << "DT=(DOUBLE DOUBLE DOUBLE)" << std::endl;
    int Nlim = 6;
    for(int i=0;i<u_prism_v.size();i++)
    {
        myfile <<setprecision(6)<< xcn_g->getVal(u_prism_v[i],0) << " ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
    }
    myfile << std::endl;
    Nlim = 6;
    for(int i=0;i<u_prism_v.size();i++)
    {
        
        myfile <<setprecision(6)<< xcn_g->getVal(u_prism_v[i],1) << " ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
        
    }
    myfile << std::endl;
    Nlim = 6;
    for(int i=0;i<u_prism_v.size();i++)
    {
        myfile <<setprecision(6)<< xcn_g->getVal(u_prism_v[i],2) << " ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
    }
    myfile << std::endl;
    myfile << "#node count per face"<<std::endl;
    Nlim = 6;
    for(int i=0;i<ufaces.size();i++)
    {
        myfile << ufaces[i].size()<<" ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
    }
    myfile << std::endl;
    myfile << "#face nodes"<<std::endl;
    for(int i=0;i<ufaces.size();i++)
    {
        for(int q=0;q<ufaces[i].size();q++)
        {
            myfile << ufaces[i][q] << " ";
        }
        myfile<<std::endl;
    }

    myfile << "#left elements (negative indicates boundary connnection)" << std::endl;
    Nlim = 6;
    int i = 0;
    for(itmap=lhface.begin();itmap!=lhface.end();itmap++)
    {
        
        myfile <<  itmap->second << " ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
        i++;
    }
    myfile << std::endl;
    myfile << "#right elements" << std::endl;
    Nlim = 6;
    i = 0;
    for(itmap=rhface.begin();itmap!=rhface.end();itmap++)
    {
        myfile <<  itmap->second << " ";
        if(i>Nlim)
        {
            myfile << std::endl;
            Nlim = Nlim+6;
        }
        i++;
    }
    myfile << std::endl;
    myfile << "#boundary connection counts" <<std::endl;
    myfile << "#boundary connection elements" << std::endl;
    myfile << "#boundary connection zones" << std::endl;
    myfile.close();
    
    
    int bndv=0;
    Array<int>* BndFace_Output = new Array<int>(BLmesh->BndFaces.size(),3);
    std::map<int,int> gbnd2lbnd;
    std::set<int> U_BndNodes_Set;
    std::vector<int> U_BndNodes_Vec;
    for(int i=0;i<BLmesh->BndFaces.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            int gbndvert = BLmesh->BndFaces[i][j];

            if(U_BndNodes_Set.find(gbndvert)==U_BndNodes_Set.end())
            {
                U_BndNodes_Set.insert(gbndvert);
                U_BndNodes_Vec.push_back(gbndvert);
                BndFace_Output->setVal(i,j,bndv);
                gbnd2lbnd[gbndvert]=bndv;
                bndv++;
            }
            else
            {
                BndFace_Output->setVal(i,j,gbnd2lbnd[gbndvert]);
            }
        }
    }
    
    string filename2 = "boundary_BL.dat";
    ofstream myfile2;
    myfile2.open(filename2);
    myfile2 << "TITLE=\"boundary.tec\"" << std::endl;
    myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
    myfile2 <<"ZONE N = " << U_BndNodes_Vec.size() << ", E = " << BLmesh->BndFaces.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
    for(int i=0;i<U_BndNodes_Vec.size();i++)
    {
      myfile2 << xcn_g->getVal(U_BndNodes_Vec[i],0)
        << " " << xcn_g->getVal(U_BndNodes_Vec[i],1)
        << " " << xcn_g->getVal(U_BndNodes_Vec[i],2) << std::endl;
    }

    for(int i=0;i<BLmesh->BndFaces.size();i++)
    {
      myfile2 << BndFace_Output->getVal(i,0)+1
        << " " << BndFace_Output->getVal(i,1)+1
        << " " << BndFace_Output->getVal(i,2)+1 << std::endl;
    }
    myfile2.close();
}

void OutputMesh_MMG(MMG5_pMesh mmgMesh)
{
    std::ofstream myfile;
    myfile.open("mmgMesh_latest.dat");
    myfile << "TITLE=\"new_volume.tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << mmgMesh->np << ", E = " << mmgMesh->ne << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<mmgMesh->np;i++)
    {
        myfile << mmgMesh->point[i+1].c[0] << " " <<mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] <<  std::endl;
    }
    for(int i=1;i<=mmgMesh->ne;i++)
    {
        myfile << mmgMesh->tetra[i].v[0] << " " << mmgMesh->tetra[i].v[1] << " " << mmgMesh->tetra[i].v[2] << " " << mmgMesh->tetra[i].v[3] << std::endl;
    }
    myfile.close();
}






void OutputBoundaryID_MMG(MMG5_pMesh mmgMesh, std::map<int,std::vector<int> > ref2bface, int bndID)
{
    std::cout << "Writing boundary from mmgMesh "<< bndID << " to boundary_" << bndID<<  ".dat" << std::endl;
    map< int, int > Loc2GlobBound;
    map< int, Vert> BC_verts;

    std::vector<int> bfaceIDs = ref2bface[bndID];
    
    int n_bc_faces  =  bfaceIDs.size();
    int* Loc        = new int[n_bc_faces*3];
    int cnt = 0;
    Vert V;
    int tel = 0;
    std::cout << bndID << " " << bfaceIDs.size() << std::endl;
    for(int j=0;j<bfaceIDs.size();j++)
    {
        int fid = bfaceIDs[j];
        for(int k=0;k<3;k++)
        {
            int val = mmgMesh->tria[fid].v[k];
            if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
            {
                Loc[tel*3+k]=Loc2GlobBound[val];
            }
            else
            {
                Loc2GlobBound[val] = cnt;
                Loc[tel*3+k]=cnt;
                V.x = mmgMesh->point[val].c[0];
                V.y = mmgMesh->point[val].c[1];
                V.z = mmgMesh->point[val].c[2];
                BC_verts[cnt] = V;
                cnt++;
            }
        }
        tel++;
    }
//
    ofstream myfile;

    string filename = "boundary_mmg_" + std::to_string(bndID) + ".dat";

    myfile.open(filename);
    myfile << "TITLE=\"boundary.tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
    myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
    for(int i=0;i<BC_verts.size();i++)
    {
      myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << std::endl;
    }

    for(int i=0;i<n_bc_faces;i++)
    {
      myfile << Loc[i*3+0]+1 << "    " << Loc[i*3+1]+1 << "   " << Loc[i*3+2]+1 << std::endl;
    }
    myfile.close();

    delete[] Loc;
    
}


void OutputBoundaryID(Partition* Pa, Mesh_Topology* meshTopo, US3D* us3d, int bndID)
{
    std::cout << "Writing boundary "<< bndID << " to boundary_" << bndID<<  ".dat" << std::endl;
    map< int, int > Loc2GlobBound;
    map< int, Vert> BC_verts;

    std::map<int,int> gV2lV               = Pa->getGlobalVert2LocalVert();
    std::vector<Vert> locVerts            = Pa->getLocalVerts();
    std::map<int,int> face2ref            = meshTopo->getFace2Ref();
    std::map<int,std::vector<int> > ref2face            = meshTopo->getRef2Face();
    std::map<int,std::vector<int> > ref2vert            = meshTopo->getRef2Vert();
    int* bnd_map = us3d->bnd_map;

    std::vector<int> bfaceIDs = ref2face[bndID];
    
    int n_bc_faces  =  bfaceIDs.size();
    int* Loc        = new int[n_bc_faces*4];
    int cnt = 0;
    Vert V;
    int tel = 0;
    std::cout << bfaceIDs[0] << " " << bfaceIDs[bfaceIDs.size()-1] << " " << us3d->nBnd << " " << n_bc_faces << std::endl;
    
    cout << "\nMin Element = "
    << *min_element(bfaceIDs.begin(), bfaceIDs.end());
    cout << "\nMax Element = "
    << *max_element(bfaceIDs.begin(), bfaceIDs.end());
    
    for(int j=0;j<bfaceIDs.size();j++)
    {
        int fid = bfaceIDs[j];
        for(int k=0;k<4;k++)
        {
            int val = us3d->ifn->getVal(fid,k);
            int lvid = gV2lV[val];
            if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
            {
                Loc[tel*4+k]=Loc2GlobBound[val];
            }
            else
            {
                Loc2GlobBound[val] = cnt;
                Loc[tel*4+k]=cnt;
                V.x = locVerts[lvid].x;
                V.y = locVerts[lvid].y;
                V.z = locVerts[lvid].z;
                BC_verts[cnt] = V;
                cnt++;
            }
        }
        tel++;
    }
//
    ofstream myfile;

    string filename = "boundary_" + std::to_string(bndID) + ".dat";

    myfile.open(filename);
    myfile << "TITLE=\"boundary.tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
    myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
    for(int i=0;i<BC_verts.size();i++)
    {
      myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << std::endl;
    }

    for(int i=0;i<n_bc_faces;i++)
    {
      myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
    }
    myfile.close();

    delete[] Loc;
    
}


void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs)
{
    int nrow = zdefs->getNrow();
    int ncol = znames->getNcol();

    std::cout << "printing boundary data..." << nrow << " " << zdefs->getNcol() << std::endl;
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            //std::cout << " (" << i << "," << j << ") ";
            std::cout << znames->getVal(i,j) << "";
        }
        std::cout << " :: ";
        for(int j=0;j<zdefs->getNcol();j++)
        {
            std::cout << zdefs->getVal(i,j) << " ";
        }
        std::cout << std::endl;
    }
    
}




void OutputPartition(Partition* part, ParArray<int>* ien, Array<double>* H,  MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<std::vector<int> > loc_elem2verts_loc = part->getLocalElem2LocalVert();
    std::map<int,int> globV2locV = part->getGlobalVert2LocalVert();
    int nloc = ien->getNrow();
    int ncol = ien->getNcol();
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\",  \"drhodx\",  \"drhody\",  \"drhodz\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    std::vector<Vert> LVerts =  part->getLocalVerts();



        
    Array<int>* ien_local= new Array<int>(nloc,ncol); 
    std::set<int> v_used;
    std::vector<int> LocalVerticesID;
    for(int i=0;i<nloc;i++)
    {
        for(int j=0;j<ncol;j++)
        {
	    int g_v_id = ien->getVal(i,j);
	    int lv_id = globV2locV[g_v_id];
            ien_local->setVal(i,j,lv_id);
            if(v_used.find(lv_id)==v_used.end())
            {
                v_used.insert(lv_id);
                LocalVerticesID.push_back(lv_id);
            }
        }
    }
    
    
    delete ien_local;



    
    int nvert = LocalVerticesID.size();
    myfile <<"ZONE N = " << nvert << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    //std::cout << rank << " number of nodes -> " << nvert << " " << H->getNrow() << std::endl;
    /*for(int i=0;i<nvert;i++)
    {
       myfile << LVerts[LocalVerticesID[i]].x << "   " << LVerts[LocalVerticesID[i]].y << "   " << LVerts[LocalVerticesID[i]].z << "	" << H->getVal(LocalVerticesID[i],0) << "	"<< H->getVal(LocalVerticesID[i],1) << "	" << H->getVal(LocalVerticesID[i],2) <<  std::endl;
    }*/
    for(int i=0;i<nvert;i++)
    {   
       myfile << LVerts[LocalVerticesID[i]].x << "   " << LVerts[LocalVerticesID[i]].y << "   " << LVerts[LocalVerticesID[i]].z <<  std::endl;
    }
    for(int i=0;i<nloc;i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    
    
    myfile.close();
}


void OutputBLElements(Partition* part, Mesh_Topology_BL* mesh_topology_bl,  MPI_Comm comm, string fname)
{
    std::vector<int> bl_elem;
    std::set<int> un_bl_elem;

    std::map<int,std::vector<int> >::iterator itt;

    for(itt=mesh_topology_bl->BLlayers.begin();itt!= mesh_topology_bl->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            if(un_bl_elem.find(itt->second[q])==un_bl_elem.end())
            {
                un_bl_elem.insert(itt->second[q]);
                bl_elem.push_back(itt->second[q]);
            }
        }
    }
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::cout << "checkcheck " << rank << " " << mesh_topology_bl->exteriorVertIDs.size() <<  " " << mesh_topology_bl->exteriorElIDs.size() << std::endl;
    std::vector<Vert> LVerts                          =  part->getLocalVerts();
    std::vector<std::vector<int> > loc_elem2verts_loc =  part->getLocalElem2LocalVert();
    std::map<int,int> gE2lE                           =  part->getGlobalElement2LocalElement();
    i_part_map* ien_part_map                          =  part->getIENpartmap();
    std::set<int> v_used;
    std::vector<int> LocalVerticesID;
    int vid,gEl,lEl;
    std::set<int> gvert_used;
    std::map<int,int> gvert2lvert;
    std::set<int> vert_used;
    std::vector<int> vert_plot;
    std::map<int,int> gv2lv;
    Array<int>* ien_bl    = new Array<int>(bl_elem.size(),8);
    Array<int>* ien_bl_ex = new Array<int>(mesh_topology_bl->exteriorElIDs.size(),8);
    int lvid = 0;
    int gvid;
    for(int i=0;i<bl_elem.size();i++)
    {
        gEl = bl_elem[i];
        lEl = gE2lE[gEl];
        
        for(int j=0;j<8;j++)
        {
            vid = loc_elem2verts_loc[lEl][j];
            gvid = ien_part_map->i_map[gEl][j];
            if(gvert_used.find(gvid)==gvert_used.end())
            {
                gvert_used.insert(gvid);
                gvert2lvert[gvid]=vid;
            }
            if(vert_used.find(vid)==vert_used.end())
            {
                vert_used.insert(vid);
                vert_plot.push_back(vid);
                gv2lv[vid] = lvid;
                ien_bl->setVal(i,j,lvid);
                lvid++;
            }
            else
            {
                ien_bl->setVal(i,j,gv2lv[vid]);
            }
        }
    }
    std::vector<int> vert_plot_ex;
    std::set<int> gvert_used2;
    std::map<int,int> gv2lv_ex;
    int lvid2=0;
    for(int i=0;i<mesh_topology_bl->exteriorElIDs.size();i++)
    {
	int el_id = mesh_topology_bl->exteriorElIDs[i];
        for(int j=0;j<8;j++)
        {
            gvid = mesh_topology_bl->exteriorVertIDs[el_id][j];
	    std::cout << gvid << " ";
            if(gvert_used2.find(gvid)==gvert_used2.end())
            {
		//std::cout << lvid2 << std::endl;
                gvert_used2.insert(gvid);
                vert_plot_ex.push_back(mesh_topology_bl->verts_g2l_ex[gvid]);
                gv2lv_ex[gvid] = lvid2;
                ien_bl_ex->setVal(i,j,lvid2);
                lvid2++;
            }
            else
            {
		//std::cout << "gv2lv_ex[gvid] " << gv2lv_ex[gvid] << std::endl;
                ien_bl_ex->setVal(i,j,gv2lv_ex[gvid]);
            }
        }
std::cout << std::endl;
    }
     
    string filename = fname + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"BL_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << vert_plot.size() << ", E = " << bl_elem.size() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl; 
 

    for(int i=0;i<vert_plot.size();i++)
    {    
        myfile << LVerts[vert_plot[i]].x << " "
               << LVerts[vert_plot[i]].y << " "
               << LVerts[vert_plot[i]].z << std::endl;
    }   
    for(int i=0;i<bl_elem.size();i++)
    {    
        myfile << ien_bl->getVal(i,0)+1 << "  " <<
                 ien_bl->getVal(i,1)+1 << "  " <<
                 ien_bl->getVal(i,2)+1 << "  " <<
                 ien_bl->getVal(i,3)+1 << "  " <<
                 ien_bl->getVal(i,4)+1 << "  " <<
                 ien_bl->getVal(i,5)+1 << "  " <<
                 ien_bl->getVal(i,6)+1 << "  " <<
                 ien_bl->getVal(i,7)+1 << std::endl;
    }   
   myfile.close();
    if(vert_plot_ex.size()!=0)
    {
    	string filename2 = fname + std::to_string(rank) + "_exterior.dat";
    	ofstream myfile2;
    	myfile2.open(filename2);
   	 myfile2 << "TITLE=\"BL_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    	myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    	myfile2 <<"ZONE N = " << vert_plot_ex.size() << ", E = " << ien_bl_ex->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    //myfile <<"ZONE N = " << vert_plot.size() << ", E = " << bl_elem.size() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl; 
    	for(int i=0;i<vert_plot_ex.size();i++)
    	{
        	myfile2 << mesh_topology_bl->local_ex_verts[vert_plot_ex[i]][0] << " "
               	<< mesh_topology_bl->local_ex_verts[vert_plot_ex[i]][1] << " "
              	 << mesh_topology_bl->local_ex_verts[vert_plot_ex[i]][2] << std::endl;
    	}
    
    	for(int i=0;i<ien_bl_ex->getNrow();i++)
   	 {
       		myfile2 << ien_bl_ex->getVal(i,0)+1 << "  " <<
                 ien_bl_ex->getVal(i,1)+1 << "  " <<
                 ien_bl_ex->getVal(i,2)+1 << "  " <<
                 ien_bl_ex->getVal(i,3)+1 << "  " <<
                 ien_bl_ex->getVal(i,4)+1 << "  " <<
                 ien_bl_ex->getVal(i,5)+1 << "  " <<
                 ien_bl_ex->getVal(i,6)+1 << "  " <<
                 ien_bl_ex->getVal(i,7)+1 << std::endl;
    	}
   	 myfile2.close();
    }    
}


void OutputBLElementsOnRoot(Array<double>* xcn_root, Array<int>* ien_root, std::vector<int> elements,  MPI_Comm comm, string fname)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::set<int> v_used;
    std::vector<int> LocalVerticesID;
    int vid,gEl,lEl;
    std::set<int> vert_used;
    std::vector<int> vert_plot;
    std::map<int,int> gv2lv;
    Array<int>* ien_bl = new Array<int>(elements.size(),8);
    int lvid = 0;
    for(int i=0;i<elements.size();i++)
    {
        gEl = elements[i];
        
        for(int j=0;j<8;j++)
        {
            vid = ien_root->getVal(gEl,j);
            if(vert_used.find(vid)==vert_used.end())
            {
                vert_used.insert(vid);
                vert_plot.push_back(vid);
                gv2lv[vid] = lvid;
                ien_bl->setVal(i,j,lvid);
                lvid++;
            }
            else
            {
                ien_bl->setVal(i,j,gv2lv[vid]);
            }
        }
    }
    
    string filename = fname + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"BL_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << vert_plot.size() << ", E = " << elements.size() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    
    for(int i=0;i<vert_plot.size();i++)
    {
        //myfile << LVerts[vert_plot[i]].x << "   " << LVerts[vert_plot[i]].y << "   " << LVerts[vert_plot[i]].z << std::endl;
        myfile << xcn_root->getVal(vert_plot[i],0) << " " << xcn_root->getVal(vert_plot[i],1) << " " << xcn_root->getVal(vert_plot[i],2) << std::endl;
    }
    
    for(int i=0;i<elements.size();i++)
    {
       myfile << ien_bl->getVal(i,0)+1 << "  " <<
                 ien_bl->getVal(i,1)+1 << "  " <<
                 ien_bl->getVal(i,2)+1 << "  " <<
                 ien_bl->getVal(i,3)+1 << "  " <<
                 ien_bl->getVal(i,4)+1 << "  " <<
                 ien_bl->getVal(i,5)+1 << "  " <<
                 ien_bl->getVal(i,6)+1 << "  " <<
                 ien_bl->getVal(i,7)+1 << std::endl;
    }
    myfile.close();
}


void OutputCompletePartition(Partition* part, ParArray<int>* ien, Array<double>*H, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<std::vector<int> > loc_elem2verts_loc = part->getLocalElem2LocalVert();
    int nloc = ien->getNrow();
    int ncol = 8;
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\",  \"drhodx\",  \"drhody\",  \"drhodz\"" << std::endl;
    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    std::vector<Vert> LVerts =  part->getLocalVerts();
    int nvert = LVerts.size();
    myfile <<"ZONE N = " << nvert << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    //std::cout << rank << " number of nodes -> " << nvert << " " << H->getNrow() << std::endl;
    for(int i=0;i<nvert;i++)
    {
       myfile << LVerts[i].x << "   " << LVerts[i].y << "   " << LVerts[i].z << " " << H->getVal(i,0) << " " << H->getVal(i,1) << " " << H->getVal(i,2) << std::endl;
    }
    
    
    for(int i=0;i<nloc;i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    
    
    myfile.close();
}

//void OutputGradient(Partition* parttn, Array<double>* H, ParallelState* pstate, MPI_Comm comm)
//{
//
//    int size;
//    MPI_Comm_size(comm, &size);
//    // Get the rank of the process
//    int rank;
//    MPI_Comm_rank(comm, &rank);
//
//    int nloc = parttn->getPart()->getNrow();
//    int ncol = 8;
//
//    int loc_v_id, glob_v_id;
//    std::vector<std::vector<int> > loc_elem2verts_loc = new Array<int>(nloc,ncol);
//    int offset = pstate->getOffset(rank);
//    std::map<int,int> GlobalElement2LocalElement = parttn->getGlobalElement2LocalElement();
//    std::map<int,std::vector<int> > GlobElem2LocVerts = parttn->getGlobElem2LocVerts();
//    std::vector<int> loc_vert;
//    std::set<int> vert_used;
//    std::vector<int> vert_plot;
//    for(int i=0;i<nloc;i++)
//    {
//        int g_el_id = i+offset;
//        loc_vert = GlobElem2LocVerts[g_el_id];
//        for(int j=0;j<8;j++)
//        {
//            loc_v_id = loc_elem2verts_loc[i][j];
//            if(vert_used.find(loc_v_id)==vert_used.end())
//            {
//                vert_used.insert(loc_v_id);
//                vert_plot.push_back(loc_v_id);
//            }
//
//        }
//        loc_vert.erase(loc_vert.begin(),loc_vert.end());
//    }
//
//
//    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
//    ofstream myfile;
//    myfile.open(filename);
//    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
//    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"rho\", \"drhox\", \"drhoy\", \"drhoz\"" << std::endl;
//    std::vector<Vert> LVerts =  parttn->getLocalVerts();
//    int nvert = vert_plot.size();
//    myfile <<"ZONE N = " << nvert << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
//    Array<double>* U0 = parttn->getUvert();
//    //std::cout << rank << " number of nodes -> " << nvert << " " << H->getNrow() << std::endl;
//    for(int i=0;i<vert_plot.size();i++)
//    {
//       myfile << LVerts[vert_plot[i]].x << "   " << LVerts[vert_plot[i]].y << "   " << LVerts[vert_plot[i]].z << "   " << U0->getVal(vert_plot[i],0) << " " << H->getVal(vert_plot[i],0) << " " << H->getVal(vert_plot[i],1) << " " << H->getVal(vert_plot[i],2) << std::endl;
//    }
//
//
//    for(int i=0;i<nloc;i++)
//    {
//       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
//                 loc_elem2verts_loc[i][1]+1 << "  " <<
//                 loc_elem2verts_loc[i][2]+1 << "  " <<
//                 loc_elem2verts_loc[i][3]+1 << "  " <<
//                 loc_elem2verts_loc[i][4]+1 << "  " <<
//                 loc_elem2verts_loc[i][5]+1 << "  " <<
//                 loc_elem2verts_loc[i][6]+1 << "  " <<
//                 loc_elem2verts_loc[i][7]+1 << std::endl;
//    }
//
//
//    myfile.close();
//}



void OutputZone(Partition* part, Array<double>* H, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<std::vector<int> > loc_elem2verts_loc = part->getLocalElem2LocalVert();
    int nloc = loc_elem2verts_loc.size();
    int ncol = 8;
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"rho\", \"drhox\", \"drhoy\", \"drhoz\"" << std::endl;
    std::vector<Vert> LVerts =  part->getLocalVerts();
    int nvert = LVerts.size();
    myfile <<"ZONE N = " << nvert << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    Array<double>* U0 = part->getUvert();
    //std::cout << rank << " number of nodes -> " << nvert << " " << H->getNrow() << std::endl;
    for(int i=0;i<nvert;i++)
    {
       myfile << LVerts[i].x << "   " << LVerts[i].y << "   " << LVerts[i].z << "   " << U0->getVal(i,0) << " " << H->getVal(i,0) << " " << H->getVal(i,1) << " " << H->getVal(i,2) << std::endl;
    }
    
    
    for(int i=0;i<nloc;i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    
    
    myfile.close();
}

void OutputQuantityPartition(Partition_old* pa, Array<double>* Quan, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = pa->ien->getNrow();
    int ncol = pa->ien->getNcol();
    int nloc = nrow;
    
    int gid;
    int lid;
    Vert V;
    
    //set<int> gid_set;
    std::map<int,Vert> vert_out;
    std::map<int,double> quan_out;
    int v=0;int u=0;int el=0;
    int* l_vert_id = new int[nrow*(ncol-1)];
    map< int, int > gid_set;

    double Q;
    for(int i=0;i<nloc;i++)
    {
        Q = Quan->getVal(i,0);
        for(int j=0;j<ncol-1;j++)
        {
            gid = pa->ien->getVal(i,j+1)-1;
            lid = pa->glob2loc_Vmap[gid];
            
            if ( gid_set.find( gid ) != gid_set.end() )
            {
                l_vert_id[el*(ncol-1)+j]=gid_set[gid];
            }
            else
            {
                l_vert_id[el*(ncol-1)+j]=v;
                
                V.x = pa->Verts->getVal(lid,0);
                V.y = pa->Verts->getVal(lid,1);
                V.z = pa->Verts->getVal(lid,2);
                
                vert_out[u] = V;
                quan_out[u] = Q;
                u++;
            }
            v++;
        }
        el++;
    }
    
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dJ\"" << std::endl;
    myfile <<"ZONE N = " << vert_out.size() << ", E = " << nrow << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    
    std::cout << "number of nodes -> " << vert_out.size() << std::endl;
    
    for(int i=0;i<vert_out.size();i++)
    {
       myfile << vert_out[(i)].x << "   " << vert_out[(i)].y << "   " << vert_out[(i)].z << "   " << quan_out[(i)] << std::endl;
    }

    for(int i=0;i<nrow;i++)
    {
       myfile << l_vert_id[i*8+0]+1 << "  " <<
                 l_vert_id[i*8+1]+1 << "  " <<
                 l_vert_id[i*8+2]+1 << "  " <<
                 l_vert_id[i*8+3]+1 << "  " <<
                 l_vert_id[i*8+4]+1 << "  " <<
                 l_vert_id[i*8+5]+1 << "  " <<
                 l_vert_id[i*8+6]+1 << "  " <<
                 l_vert_id[i*8+7]+1 << std::endl;
    }
    
    
    myfile.close();
    delete[] l_vert_id;

}

/*
void OutputPartionVolumes(ParArray<int>* ien, Array<double>* xcn_on_root, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    Partition_old* pv = CollectVerticesPerRank(ien,xcn_on_root,comm);
    
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    
    int gid;
    int lid;
    Vert V;
    
    //set<int> gid_set;
    std::map<int,Vert> vert_out;
    int v=0;int u=0;int el=0;
    int* l_vert_id = new int[nrow*(ncol-1)];
    map< int, int > gid_set;
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            gid = ien->getVal(i,j+1)-1;
            lid = pv->glob2loc_Vmap[gid];
            
            if ( gid_set.find( gid ) != gid_set.end() )
            {
                l_vert_id[el*(ncol-1)+j]=gid_set[gid];
            }
            else
            {
                l_vert_id[el*(ncol-1)+j]=v;
                
                V.x = pv->Verts->getVal(lid,0);
                V.y = pv->Verts->getVal(lid,1);
                V.z = pv->Verts->getVal(lid,2);
                
                vert_out[u]=V;
                
                u++;
            }
            v++;
        }
        el++;
    }
    
    string filename = "volume_per_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << vert_out.size() << ", E = " << nrow << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    
    std::cout << "number of nodes -> " << vert_out.size() << std::endl;
    
    for(int i=0;i<vert_out.size();i++)
    {
       myfile << vert_out[(i)].x << "   " << vert_out[(i)].y << "   " << vert_out[(i)].z << std::endl;
    }

    for(int i=0;i<nrow;i++)
    {
       myfile << l_vert_id[i*8+0]+1 << "  " <<
                 l_vert_id[i*8+1]+1 << "  " <<
                 l_vert_id[i*8+2]+1 << "  " <<
                 l_vert_id[i*8+3]+1 << "  " <<
                 l_vert_id[i*8+4]+1 << "  " <<
                 l_vert_id[i*8+5]+1 << "  " <<
                 l_vert_id[i*8+6]+1 << "  " <<
                 l_vert_id[i*8+7]+1 << std::endl;
    }
    myfile.close();
    delete[] l_vert_id;
}
*/

/*
void OutputPartitionFaces()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    const char* fn_conn="grids/adept/conn.h5";
    const char* fn_grid="grids/adept/grid.h5";
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);

    std::vector<int> partfaces = GetAdjacencyForUS3D_V4(ief, comm);
    
    if(rank == 0)
    {
        Array<double>* xcn = ReadDataSetFromFile<double>(fn_grid,"xcn");
        Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
        int n_bc_faces = partfaces.size();
        Vert V;
        int* Loc = new int[n_bc_faces*4];
        map< int, int > Loc2GlobBound;
        map< int, Vert> P_verts;
        int cnt = 0;
        int tel = 0;
        int teller = 0;

        for(int j=0;j<partfaces.size();j++)
        {
            int face = partfaces[j];
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(face,k+1);

                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    
                    P_verts[cnt] = V;
                    
                    cnt++;
                }
                
                teller=teller+1;
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "partfaces_" + std::to_string(rank) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"partitionfaces.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << P_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        std::cout << "number of boundary nodes -> " << P_verts.size() << std::endl;
        
        for(int i=0;i<P_verts.size();i++)
        {
           myfile << P_verts[(i)].x << "   " << P_verts[(i)].y << "   " << P_verts[(i)].z << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "  "
                  << Loc[i*4+1]+1 << "  "
                  << Loc[i*4+2]+1 << "  "
                  << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
        delete ifn;
        delete xcn;
     }
}
*/


void WriteBoundaryDataInSerial3(Array<double>* xcn)
{
    string filename = "boundary_nodes.dat";
    ofstream myfile;
    myfile.open(filename);
    
    
    for(int i=0;i<3097156;i++)
    {
        myfile << xcn->getVal(i,0) << " " << xcn->getVal(i,1) << " " << xcn->getVal(i,2) << std::endl;
    }
    myfile.close();
}


void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts)
{
    for(int bc=3;bc<zdefs->getNrow();bc++)
    {
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    J_verts[cnt]  = detJ_verts[val-1];
                    V_verts[cnt]  = vol_verts[val-1];
                    Jno_verts[cnt]  = Jnorm_verts[val-1];
                    cnt++;
                }
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dJ\", \"Vol\", \"dJ/Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << J_verts[(i)] << "   " << V_verts[(i)] << "   " << Jno_verts[(i)] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
    }
}



void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts)
{
    for(int bc=3;bc<zdefs->getNrow();bc++)
    {
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        //int* pltJ = new int[n_bc_faces*4];
        int teller = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                //pltJ[teller] = val;
                //std::cout << val << " ";
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    V_verts[cnt]   = vol_verts[val-1];
                    cnt++;
                }
                
                teller=teller+1;
            }
            //std::cout << std::endl;
            tel++;
        }
        //cout << "\nlargest id = " << largest(pltJ,n_bc_faces*4) << std::endl;
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        std::cout << "number of boundary nodes -> " << BC_verts.size() << std::endl;
        
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << V_verts[i] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
    }
}
