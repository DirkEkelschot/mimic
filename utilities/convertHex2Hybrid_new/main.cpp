#include "../../src/adapt_io.h"
#include "../../src/adapt_boundary.h"
#include "../../src/hex2tet.h"
#include "../../src/adapt_bltopology.h"
#include "../../src/NekFace.h"

#include "cgnslib.h"




std::vector<double> ComputeCentroid(std::vector<std::vector<double> > elemcoords)
{
    std::vector<double> cent(3);
    cent[0] = 0.0;cent[1] = 0.0;cent[2] = 0.0;
    for(int i = 0;i<elemcoords.size();i++)
    {
        cent[0] = cent[0] + elemcoords[i][0];
        cent[1] = cent[1] + elemcoords[i][1];
        cent[2] = cent[2] + elemcoords[i][2];
    }
    
    cent[0] = cent[0] / elemcoords.size();
    cent[1] = cent[1] / elemcoords.size();
    cent[2] = cent[2] / elemcoords.size();
    
    return cent;
}



struct Tetra{
    
    int globID;
    int* VertIDs;
    
    std::map<int,std::vector<int> > GlobalFace2GlobalVert;
    std::map<int,std::vector<int> > GlobalFace2LocalVert;
    
    std::map<int,std::vector<int> > LocalFace2GlobalVert;
    std::map<int,std::vector<int> > LocalFace2LocalVert;
};



struct Prism{
    
    int globID;
    int* VertIDs;
    
    std::map<int,std::vector<int> > GlobalFace2GlobalVert;
    std::map<int,std::vector<int> > GlobalFace2LocalVert;
    
    std::map<int,std::vector<int> > LocalFace2GlobalVert;
    std::map<int,std::vector<int> > LocalFace2LocalVert;
};

struct Prisms
{
	std::vector<Prism* > prisms;
	FaceSetPointer m_FaceSet;
	std::set<int> handledHexes;	
	std::set<int> handledQuads;
	FaceSetPointer m_QuadFaceHandled;	
};





struct MeshPoint
{
    //double x=0.0;
    //double y=0.0;
    //double z=0.0;
	int id;
    double x;
    double y;
    double z;
    int ref;
};




struct CutHex
{
	int cutType; // cutType=0 is the Rising cut, cutType=1 is Falling cut.
	std::vector<int> nodes;
	std::vector<std::vector<int> > prisms;
    std::vector<std::vector<int> > tetras;
};





void OutputBoundaryLayerPrisms(Array<double>* xcn_g, Prisms* bp, MPI_Comm comm,string fname)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::vector<Prism* >::iterator iter;

    std::set<int> unique_prism_verts;
    std::vector<int> u_prism_v;

    //Array<int>* local_prisms = new Array<int>(BLmesh->Nprisms,6);
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
    
    for(int i=0;i<bp->prisms.size();i++)
    {
       
		int* prismVrts = bp->prisms[i]->VertIDs;
		
		
		
		int gElID = i;
		
//		std::vector<std::vector<double> > elemCoords;
		
		for(int q=0;q<6;q++)
		{
//			std::vector<double> coords(3);
//			coords[0] = xcn_g->getVal(prismVrts[q],0);
//			coords[1] = xcn_g->getVal(prismVrts[q],1);
//			coords[2] = xcn_g->getVal(prismVrts[q],2);
			if(unique_prism_verts.find(prismVrts[q])==unique_prism_verts.end())
			{
				unique_prism_verts.insert(prismVrts[q]);
				u_prism_v.push_back(prismVrts[q]);
				gv2lv_prisms[prismVrts[q]] = lvid;
				lvid++;
			}
			
//			elemCoords.push_back(coords);

		}
		
		
//		std::vector<double> cent = ComputeCentroid(elemCoords);

		std::map<int,std::vector<int> > LF2GN = bp->prisms[i]->LocalFace2GlobalVert;
//
		std::map<int,std::vector<int> >::iterator itm;
		for(itm=LF2GN.begin();itm!=LF2GN.end();itm++)
		{
			std::set<int> face;
			
			std::vector<std::vector<double> > fcoords;
			for(int q=0;q<itm->second.size();q++)
			{
				face.insert(itm->second[q]);
//				std::vector<double> coords(3);
//				coords[0] = xcn_g->getVal(prismVrts[q],0);
//				coords[1] = xcn_g->getVal(prismVrts[q],1);
//				coords[2] = xcn_g->getVal(prismVrts[q],2);
//				fcoords.push_back(coords);
			}
			
//			std::vector<double> Fcent = ComputeCentroid(fcoords);
//			
//			if(fcoords.size()==3)
//			{
//				Vec3D* vecf0_0 = new Vec3D;
//				vecf0_0->c0 = fcoords[1][0]-fcoords[0][0];
//				vecf0_0->c1 = fcoords[1][1]-fcoords[0][1];
//				vecf0_0->c2 = fcoords[1][2]-fcoords[0][2];
//				Vec3D* vecf0_1 = new Vec3D;
//				vecf0_1->c0 = fcoords[2][0]-fcoords[0][0];
//				vecf0_1->c1 = fcoords[2][1]-fcoords[0][1];
//				vecf0_1->c2 = fcoords[2][2]-fcoords[0][2];
//				
//				Vec3D* nf0 = ComputeSurfaceNormal(vecf0_0,vecf0_1);
//				
//		        Vec3D* r0 = new Vec3D;
//		        r0->c0 = (Fcent[0]-cent[0]);
//		        r0->c1 = (Fcent[1]-cent[1]);
//		        r0->c2 = (Fcent[2]-cent[2]);
//		        
//		        double orient0 = DotVec3D(r0,nf0);
//		        std::cout << fcoords.size() << " " << orient0 << std::endl;
//
//		        
//			}
//			if(fcoords.size()==4)
//			{
//				Vec3D* vecf0_0 = new Vec3D;
//				vecf0_0->c0 = fcoords[1][0]-fcoords[0][0];
//				vecf0_0->c1 = fcoords[1][1]-fcoords[0][1];
//				vecf0_0->c2 = fcoords[1][2]-fcoords[0][2];
//				Vec3D* vecf0_1 = new Vec3D;
//				vecf0_1->c0 = fcoords[3][0]-fcoords[0][0];
//				vecf0_1->c1 = fcoords[3][1]-fcoords[0][1];
//				vecf0_1->c2 = fcoords[3][2]-fcoords[0][2];
//				
//				Vec3D* nf0 = ComputeSurfaceNormal(vecf0_0,vecf0_1);
//				
//		        Vec3D* r0 = new Vec3D;
//		        r0->c0 = (Fcent[0]-cent[0]);
//		        r0->c1 = (Fcent[1]-cent[1]);
//		        r0->c2 = (Fcent[2]-cent[2]);
//		        
//		        double orient0 = DotVec3D(r0,nf0);
//		        
//		        std::cout << fcoords.size() << " " << orient0 << std::endl;
//			}
//			
			
			
			
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
//		std::cout<< "LF2GN.soze  " << LF2GN.size() << std::endl;
		npr++;
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
    
    /*
    int bndv=0;
    Array<int>* BndFace_Output = new Array<int>(bp->BndFaces.size(),3);
    std::map<int,int> gbnd2lbnd;
    std::set<int> U_BndNodes_Set;
    std::vector<int> U_BndNodes_Vec;
    for(int i=0;i<bp->BndFaces.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            int gbndvert = bp->BndFaces[i][j];

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
    myfile2 <<"ZONE N = " << U_BndNodes_Vec.size() << ", E = " << bp->BndFaces.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
    for(int i=0;i<U_BndNodes_Vec.size();i++)
    {
      myfile2 << xcn_g->getVal(U_BndNodes_Vec[i],0)
        << " " << xcn_g->getVal(U_BndNodes_Vec[i],1)
        << " " << xcn_g->getVal(U_BndNodes_Vec[i],2) << std::endl;
    }

    for(int i=0;i<bp->BndFaces.size();i++)
    {
      myfile2 << BndFace_Output->getVal(i,0)+1
        << " " << BndFace_Output->getVal(i,1)+1
        << " " << BndFace_Output->getVal(i,2)+1 << std::endl;
    }
    myfile2.close();
    
    delete BndFace_Output;
    gbnd2lbnd.clear();
    
    U_BndNodes_Set.clear();
    U_BndNodes_Vec.clear();
    */
    
}

//=========================================================================
//=========================================================================
//=========================================================================















int FindOppositeLocalFaceHex(int fref)
{
	int O_fref = -1;
	if(fref == 0)
	{
		O_fref=1;
	}
	if(fref == 1)
	{
		O_fref=0;
	}
	if(fref == 2)
	{
		O_fref=3;
	}
	if(fref == 3)
	{
		O_fref=2;
	}
	if(fref == 4)
	{
		O_fref=5;
	}
	if(fref == 5)
	{
		O_fref=4;
	}
	
	return O_fref;
}





std::map<int,FaceSharedPtr> FindLocalFaceInHexV2(FaceSharedPtr faceTestPtr, 
												int* hnodes, 
												int* hfaces, int face, 
												int &locFid, int &fref)
{
	
	int m_rhorient[6][4] = {{0,3,2,1},
						    {4,5,6,7},
		                    {0,4,7,3},
		                    {1,2,6,5},
		                    {0,1,5,4},
					        {3,7,6,2}};
		
		int m_lhorient[6][4] = {{0,1,2,3},
						        {4,7,6,5},
		                        {0,3,7,4},
		                        {1,5,6,2},
		                        {0,4,5,1},
								{3,2,6,5}};
		
		
	std::map<int,FaceSharedPtr> facemap;
			
	for(int s=0;s<6;s++)
	{
		if(hfaces[s]==face)
		{
			locFid=s;
		}
	
		std::vector<int> v_face(4);
		
		for(int t=0;t<4;t++)
		{
			int lid = m_rhorient[s][t];
			
			v_face[t] = hnodes[lid];
		}
		
		FaceSharedPtr facePtr = std::shared_ptr<NekFace>(new NekFace(v_face));
		
		NekFace* fa = new NekFace(v_face);
		NekFace* fo = new NekFace(faceTestPtr->GetEdgeIDs());
		
		if(facePtr==faceTestPtr)
		{
			fref = s;
			faceTestPtr->SetFaceRef(s);
		}
		
		facemap[s] = facePtr;

	}
	if(locFid == -1)
	{
		std::cout << "ERROR: bfaceid not in element. " << face << std::endl; 
	}
	
	
	return facemap;
}





std::map<int,FaceSharedPtr> GetLocalFaceMap_Hex(FaceSharedPtr faceTestPtr, Array<int>* ief, Array<int>* ien, int elid)
{
	
	int m_rhorient[6][4] = {{0,3,2,1},
						    {4,5,6,7},
		                    {0,4,7,3},
		                    {1,2,6,5},
		                    {0,1,5,4},
					        {3,7,6,2}};
		
		int m_lhorient[6][4] = {{0,1,2,3},
						        {4,7,6,5},
		                        {0,3,7,4},
		                        {1,5,6,2},
		                        {0,4,5,1},
								{3,2,6,5}};
		
		
	std::map<int,FaceSharedPtr> id2face;
			
	for(int s=0;s<6;s++)
	{
		std::vector<int> v_face(4);
		
		for(int t=0;t<4;t++)
		{
			int lid = m_rhorient[s][t];
			
			v_face[t] = ien->getVal(elid,lid);
		}
		
		FaceSharedPtr facePtr = std::shared_ptr<NekFace>(new NekFace(v_face));
		
		id2face[s] = facePtr;
		
	}
	
	return id2face;
}





std::map<int,FaceSharedPtr> GetLocalFaceMapHex(std::vector<int> hnodes_ids)
{
	
	int m_rhorient[6][4] = {{0,3,2,1},
						    {4,5,6,7},
		                    {0,4,7,3},
		                    {1,2,6,5},
		                    {0,1,5,4},
					        {3,7,6,2}};
		
		int m_lhorient[6][4] = {{0,1,2,3},
						        {4,7,6,5},
		                        {0,3,7,4},
		                        {1,5,6,2},
		                        {0,4,5,1},
								{3,2,6,5}};
		
		
	std::map<int,FaceSharedPtr> facemap;
			
	for(int s=0;s<6;s++)
	{
		std::vector<int> v_face(4);
		
		for(int t=0;t<4;t++)
		{
			int lid = m_rhorient[s][t];
			
			v_face[t] = hnodes_ids[lid];
		}
		
		FaceSharedPtr facePtr = std::shared_ptr<NekFace>(new NekFace(v_face));
		
		facemap[s] = facePtr;
		
	}
	
	return facemap;
}


std::vector<FaceSharedPtr> CutQuadIntoTris(FaceSharedPtr QuadfacePtr,int cutType)
{
	
	std::vector<FaceSharedPtr> tris;
	if(cutType == 0)
	{
		std::vector<int> cutF0(3);
		cutF0[0] = QuadfacePtr->GetEdgeIDs()[0];
		cutF0[1] = QuadfacePtr->GetEdgeIDs()[1];
		cutF0[2] = QuadfacePtr->GetEdgeIDs()[2];
		FaceSharedPtr Tri0FacePtr = std::shared_ptr<NekFace>(new NekFace(cutF0));

		std::vector<int> cutF1(3);
		cutF1[0] = QuadfacePtr->GetEdgeIDs()[0];
		cutF1[1] = QuadfacePtr->GetEdgeIDs()[2];
		cutF1[2] = QuadfacePtr->GetEdgeIDs()[3];	
		FaceSharedPtr Tri1FacePtr = std::shared_ptr<NekFace>(new NekFace(cutF1));
		
		tris.push_back(Tri0FacePtr);
		tris.push_back(Tri1FacePtr);
	}
	else if(cutType == 1)
	{
		std::vector<int> cutF0(3);
		cutF0[0] = QuadfacePtr->GetEdgeIDs()[0];
		cutF0[1] = QuadfacePtr->GetEdgeIDs()[1];
		cutF0[2] = QuadfacePtr->GetEdgeIDs()[3];
		FaceSharedPtr Tri0FacePtr = std::shared_ptr<NekFace>(new NekFace(cutF0));

		std::vector<int> cutF1(3);
		cutF1[0] = QuadfacePtr->GetEdgeIDs()[1];
		cutF1[1] = QuadfacePtr->GetEdgeIDs()[2];
		cutF1[2] = QuadfacePtr->GetEdgeIDs()[3];		
		FaceSharedPtr Tri1FacePtr = std::shared_ptr<NekFace>(new NekFace(cutF1));
	}
	
	return tris;
}

std::vector<Tetra*> CutPrismInTets(Prism* p, FaceSetPointer &cutFaces, FaceSetPointer &tetFaces, FaceSetPointer &m_CutFaceSet, FaceSetPointer HexFaces)
{

	int N = cutFaces.size();
	std::vector<Tetra*> tets;
	
	int m_faceOrient[4][3] = {{1,3,2},
							  {0,1,3},
			                  {0,2,3},
			                  {0,2,1}};
	
	std::vector<int> m_DefaultCuts(3);
	
	if(N==0)
	{
		int map[3][4] = {{0,1,2,4},
						 {1,3,4,0},
						 {1,4,3,5}};
		
		int cutTypes[3] = {0,1,0};
		
		for(int i=0;i<3;i++)
		{
			// Make sure that the quad faces are marked as cut.
			std::vector<int> Q    = p->LocalFace2GlobalVert[i+2];
			FaceSharedPtr facePtr = std::shared_ptr<NekFace>(new NekFace(Q));
			facePtr->SetFaceRef(cutTypes[0]);
			Tetra* t = new Tetra;
			t->VertIDs = new int[4];
			for(int j=0;j<4;j++)
			{
				int lid = map[i][j];
				t->VertIDs[j] = p->VertIDs[lid];	
			}
			
			if(HexFaces.find(facePtr)!=HexFaces.end())
			{
				m_CutFaceSet.insert(facePtr);
			}
			
			for(int k=0;k<4;k++)
			{
				for(int l=0;l<3;l++)
				{
					int gid = t->VertIDs[m_faceOrient[k][l]];
					t->LocalFace2GlobalVert[k].push_back(gid);
				}
//				
				FaceSharedPtr TriFacePtr = std::shared_ptr<NekFace>(new NekFace(t->LocalFace2GlobalVert[k]));
				
				tetFaces.insert(TriFacePtr);
			}			
			tets.push_back(t);
		}
	}
	if(N==1)
	{		
		int cutType;
		int locFid = -1;
		for(int i=0;i<3;i++)
		{
			std::vector<int> Q    = p->LocalFace2GlobalVert[i+2];
			FaceSharedPtr facePtr = std::shared_ptr<NekFace>(new NekFace(Q));
			FaceSetPointer::iterator FPointer 	=  cutFaces.find(facePtr);
			locFid = 2+i;
			if(FPointer!=cutFaces.end())
			{
				cutType = (*FPointer)->GetFaceRef();
			}
			else
			{
				cutType = 0;
				facePtr->SetFaceRef(0);
				cutFaces.insert(facePtr);
			}
		}
		
		
		// Different configurations:
		//[0 1 0]
		int tetConfig010[3][4] = {{0,1,2,4},
							      {1,3,4,0},
						          {1,4,3,5}};
		//[0 0 1]
		int tetConfig001[3][4] = {{0,1,2,5},
								  {0,5,4,3},
				                  {0,4,5,2}};
		//[0 1 1]
		int tetConfig011[3][4] = {{0,1,2,4},
							  	  {0,1,4,5},
								  {3,4,5,1}};
		//[1 1 0]
		int tetConfig110[3][4] = {{0,1,2,3},
							      {1,4,3,5},
						          {1,3,4,2}};
		//[1 0 1]
		int tetConfig101[3][4] = {{0,1,2,5},
								  {2,3,5,4},
						          {2,5,3,0}};
		//[1 0 0]
		int tetConfig100[3][4] = {{0,1,2,3},
							      {3,5,2,4},
						          {3,2,5,1}};

		std::vector<std::vector<int> > map;
		
		if(locFid == 2 && cutType == 0)// [0 0 1] [0 1 0] [0 1 1]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig001[i][j];
				}
				map.push_back(entry);
			}
		}
		if(locFid == 2 && cutType == 1)// [1 0 0] [1 1 0] [1 0 1]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig100[i][j];
				}
				map.push_back(entry);
			}
		}
		if(locFid == 3 && cutType == 0)// [0 0 1] [1 0 0] [1 0 1]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig001[i][j];
				}
				map.push_back(entry);
			}	
		}
		if(locFid == 3 && cutType == 1)// [0 1 1] [1 1 0] [0 1 0]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig011[i][j];
				}
				map.push_back(entry);
			}	
		}
		if(locFid == 4 && cutType == 0)// [1 0 0] [1 1 0] [0 1 0]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig100[i][j];
				}
				map.push_back(entry);
			}
		}
		if(locFid == 4 && cutType == 1)// [0 0 1] [1 0 1] [0 1 1]
		{
			for(int i=0;i<3;i++)
			{
				std::vector<int> entry(4);
				for(int j=0;j<4;j++)
				{
					entry[j] = tetConfig001[i][j];
				}
				map.push_back(entry);
			}	
		}
		
		for(int i=0;i<3;i++)
		{
			Tetra* t = new Tetra;
			t->VertIDs = new int[4];
			for(int j=0;j<4;j++)
			{
				int lid = map[i][j];
				t->VertIDs[j] = p->VertIDs[lid];	
			}
			
			for(int k=0;k<4;k++)
			{
				for(int l=0;l<3;l++)
				{
					int gid = t->VertIDs[m_faceOrient[k][l]];
					t->LocalFace2GlobalVert[k].push_back(gid);
				}
				
				FaceSharedPtr TriFacePtr = std::shared_ptr<NekFace>(new NekFace(t->LocalFace2GlobalVert[k]));
				
				tetFaces.insert(TriFacePtr);
			}
//			
			tets.push_back(t);
		}
		/**/
	}
	
	return tets;

}



std::vector<Prism*> CutHexInPrism(int* hnodes, 
								  int* hfaces, 
								  FaceSharedPtr QuadFacePtr, 
								  int face,
								  int cutType, int &locFid, FaceSetPointer &m_CutFaceSet)
{

	int locFid_reorient = -1;
	std::map<int,FaceSharedPtr> facemap = FindLocalFaceInHexV2(QuadFacePtr, hnodes, hfaces, face, locFid, locFid_reorient);
	
	Prism* bprsm0 = new Prism;
	Prism* bprsm1 = new Prism;
	
	bprsm0->VertIDs = new int[6];
	bprsm1->VertIDs = new int[6];
	
	std::vector<int>cutF0(3);
	std::vector<int>cutF1(3);	
	QuadFacePtr->SetFaceRef(cutType);
	m_CutFaceSet.insert(QuadFacePtr);
	
	if(cutType == 0)
	{
		cutF0[0] = facemap[locFid_reorient]->GetEdgeIDs()[0];
		cutF0[1] = facemap[locFid_reorient]->GetEdgeIDs()[1];
		cutF0[2] = facemap[locFid_reorient]->GetEdgeIDs()[2];

		bprsm0->VertIDs[0] = cutF0[0];
		bprsm0->VertIDs[1] = cutF0[1];
		bprsm0->VertIDs[2] = cutF0[2];
		bprsm0->LocalFace2GlobalVert[0]=cutF0;

		cutF1[0] = facemap[locFid_reorient]->GetEdgeIDs()[0];
		cutF1[1] = facemap[locFid_reorient]->GetEdgeIDs()[2];
		cutF1[2] = facemap[locFid_reorient]->GetEdgeIDs()[3];

		bprsm1->VertIDs[0] = cutF1[0];
		bprsm1->VertIDs[1] = cutF1[1];
		bprsm1->VertIDs[2] = cutF1[2];
		bprsm1->LocalFace2GlobalVert[0]=cutF1;
		
	}
	else if(cutType == 1)
	{
		cutF0[0] = facemap[locFid_reorient]->GetEdgeIDs()[0];
		cutF0[1] = facemap[locFid_reorient]->GetEdgeIDs()[1];
		cutF0[2] = facemap[locFid_reorient]->GetEdgeIDs()[3];

		bprsm0->VertIDs[0] = cutF0[0];
		bprsm0->VertIDs[1] = cutF0[1];
		bprsm0->VertIDs[2] = cutF0[2];
		bprsm0->LocalFace2GlobalVert[0]=cutF0;

		
		cutF1[0] = facemap[locFid_reorient]->GetEdgeIDs()[1];
		cutF1[1] = facemap[locFid_reorient]->GetEdgeIDs()[2];
		cutF1[2] = facemap[locFid_reorient]->GetEdgeIDs()[3];		
		
		bprsm1->VertIDs[0] = cutF1[0];
		bprsm1->VertIDs[1] = cutF1[1];
		bprsm1->VertIDs[2] = cutF1[2];
		bprsm1->LocalFace2GlobalVert[0]=cutF1;
	}
	
	int O_locFid_reorient = FindOppositeLocalFaceHex(locFid_reorient);
	
	std::vector<int>O_cutF0(3);
	std::vector<int>O_cutF1(3);		
	
	std::vector<int> v_faceTest(4);	
	v_faceTest[0] = facemap[O_locFid_reorient]->GetEdgeIDs()[0];
	v_faceTest[1] = facemap[O_locFid_reorient]->GetEdgeIDs()[1];
	v_faceTest[2] = facemap[O_locFid_reorient]->GetEdgeIDs()[2];
	v_faceTest[3] = facemap[O_locFid_reorient]->GetEdgeIDs()[3];
	
	FaceSharedPtr OppositeQuadFacePtr = std::shared_ptr<NekFace>(new NekFace(v_faceTest));
	OppositeQuadFacePtr->SetFaceRef(cutType);
	m_CutFaceSet.insert(OppositeQuadFacePtr);
	
	//std::cout << "OppositeQuadFacePtr " << OppositeQuadFacePtr->GetFaceRef()<< " " << cutType << std::endl;

	
	
	if(cutType == 0)
	{	
		O_cutF0[0]  = facemap[O_locFid_reorient]->GetEdgeIDs()[0];
		O_cutF0[1]  = facemap[O_locFid_reorient]->GetEdgeIDs()[2];
		O_cutF0[2]  = facemap[O_locFid_reorient]->GetEdgeIDs()[3];

		bprsm0->VertIDs[3] = O_cutF0[0];
		bprsm0->VertIDs[4] = O_cutF0[1];
		bprsm0->VertIDs[5] = O_cutF0[2];
		bprsm0->LocalFace2GlobalVert[1]=O_cutF0;
		
		O_cutF1[0]  = facemap[O_locFid_reorient]->GetEdgeIDs()[0];
		O_cutF1[1]  = facemap[O_locFid_reorient]->GetEdgeIDs()[1];
		O_cutF1[2]  = facemap[O_locFid_reorient]->GetEdgeIDs()[2];

		bprsm1->VertIDs[3] = O_cutF1[0];
		bprsm1->VertIDs[4] = O_cutF1[1];
		bprsm1->VertIDs[5] = O_cutF1[2];
		bprsm1->LocalFace2GlobalVert[1]=O_cutF1;
		
	}
	else if(cutType == 1)
	{
		O_cutF0[0]  = facemap[O_locFid_reorient]->GetEdgeIDs()[0];
		O_cutF0[1]  = facemap[O_locFid_reorient]->GetEdgeIDs()[1];
		O_cutF0[2]  = facemap[O_locFid_reorient]->GetEdgeIDs()[3];

		bprsm0->VertIDs[3] = O_cutF0[0];
		bprsm0->VertIDs[4] = O_cutF0[1];
		bprsm0->VertIDs[5] = O_cutF0[2];
		bprsm0->LocalFace2GlobalVert[1]=O_cutF0;
		
		O_cutF1[0]  = facemap[O_locFid_reorient]->GetEdgeIDs()[3];
		O_cutF1[1]  = facemap[O_locFid_reorient]->GetEdgeIDs()[1];
		O_cutF1[2]  = facemap[O_locFid_reorient]->GetEdgeIDs()[2];	

		bprsm1->VertIDs[3] = O_cutF1[0];
		bprsm1->VertIDs[4] = O_cutF1[1];
		bprsm1->VertIDs[5] = O_cutF1[2];
		bprsm1->LocalFace2GlobalVert[1]=O_cutF1;
	}
	
	
	std::vector<int> p0quad0(4);
	p0quad0[0] = bprsm0->VertIDs[0];
	p0quad0[1] = bprsm0->VertIDs[2];
	p0quad0[2] = bprsm0->VertIDs[4];
	p0quad0[3] = bprsm0->VertIDs[3];
	std::vector<int> p1quad0(4);
	p1quad0[0] = bprsm1->VertIDs[0];
	p1quad0[1] = bprsm1->VertIDs[2];
	p1quad0[2] = bprsm1->VertIDs[4];
	p1quad0[3] = bprsm1->VertIDs[3];

	bprsm0->LocalFace2GlobalVert[2] = p0quad0;
	bprsm1->LocalFace2GlobalVert[2] = p1quad0;
			
	std::vector<int> p0quad1(4);
	p0quad1[0] = bprsm0->VertIDs[1];
	p0quad1[1] = bprsm0->VertIDs[5];
	p0quad1[2] = bprsm0->VertIDs[4];
	p0quad1[3] = bprsm0->VertIDs[2];
	std::vector<int> p1quad1(4);
	p1quad1[0] = bprsm1->VertIDs[1];
	p1quad1[1] = bprsm1->VertIDs[5];
	p1quad1[2] = bprsm1->VertIDs[4];
	p1quad1[3] = bprsm1->VertIDs[2];
	
	bprsm0->LocalFace2GlobalVert[3] = p0quad1;
	bprsm1->LocalFace2GlobalVert[3] = p1quad1;
	
	std::vector<int> p0quad2(4);
	p0quad2[0] = bprsm0->VertIDs[0];
	p0quad2[1] = bprsm0->VertIDs[3];
	p0quad2[2] = bprsm0->VertIDs[5];
	p0quad2[3] = bprsm0->VertIDs[1];
	std::vector<int> p1quad2(4);
	p1quad2[0] = bprsm1->VertIDs[0];
	p1quad2[1] = bprsm1->VertIDs[3];
	p1quad2[2] = bprsm1->VertIDs[5];
	p1quad2[3] = bprsm1->VertIDs[1];
	
	bprsm0->LocalFace2GlobalVert[4] = p0quad2;
	bprsm1->LocalFace2GlobalVert[4] = p1quad2;
	
	std::vector<Prism*> ps;
	
	ps.push_back(bprsm0);
	ps.push_back(bprsm1);
	
	return ps;
}






Prisms* CreateBoundaryLayerPrisms(std::map<int,std::vector<int> > bnd_face_map,
		Array<int>* ien,
		Array<int>* ief,
		Array<int>* ife,
		Array<int>* ifn,
		std::map<std::set<int>,int> tria_ref,
		std::map<std::set<int>,int> quad_ref,
		int wall_id,
		int cutType,
		int nLayer,
		FaceSetPointer &m_CutFaceSet)
{
	
	Prisms* bp = new Prisms;
	
	int Nel = ien->getNrow();
	int elid;
	
	int fnd = 0;
	std::set<int> handledHexes;
	
	
	int m_rhorient[6][4] = {{0,3,2,1},
					        {4,5,6,7},
	                        {0,4,7,3},
	                        {1,2,6,5},
	                        {0,1,5,4},
							{3,7,6,2}};
	
	int m_lhorient[6][4] = {{0,1,2,3},
					        {4,7,6,5},
	                        {0,3,7,4},
	                        {1,5,6,2},
	                        {0,4,5,1},
							{3,2,6,5}};
	
	FaceSetPointer element2face;
	FaceSet fset;
	int tel = 0;
	int tt=0;
	int nfo=0;
	for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
	{
		int bfaceid = bnd_face_map[wall_id][bf];
		int face 	= bfaceid;
		int elid0   = ife->getVal(face,0);
		int elid1   = ife->getVal(face,1);
					
		if(elid0<Nel)
		{
			elid = elid0;
		}
		else
		{
			elid = elid1;
		}

		for(int c=0;c<nLayer;c++)
		{

			int locFid = -1;
			int fref   = -1;
			
			std::vector<int> v_faceTest(4);	
			v_faceTest[0] = ifn->getVal(face,0);
			v_faceTest[1] = ifn->getVal(face,1);
			v_faceTest[2] = ifn->getVal(face,2);
			v_faceTest[3] = ifn->getVal(face,3);
			
			FaceSharedPtr faceTestPtr = std::shared_ptr<NekFace>(new NekFace(v_faceTest));
			faceTestPtr->SetFaceRef(cutType);
			int offset_ien = elid*8;
			int* tmp_ien   = ien->data + offset_ien;
			int offset_ief = elid*6;
			int* tmp_ief   = ief->data + offset_ief;
			
			std::vector<Prism*> prsms = CutHexInPrism(tmp_ien, tmp_ief,faceTestPtr,face,cutType,locFid, m_CutFaceSet);
			
			if(c == (nLayer-1))
			{ 				
				FaceSharedPtr facePtr0 = std::shared_ptr<NekFace>(new NekFace(prsms[0]->LocalFace2GlobalVert[1]));	
				bp->m_FaceSet.insert(facePtr0);
				FaceSharedPtr facePtr1 = std::shared_ptr<NekFace>(new NekFace(prsms[1]->LocalFace2GlobalVert[1]));	
				bp->m_FaceSet.insert(facePtr1);
				
				int face_loc = FindOppositeLocalFaceHex(locFid);
				face 		 = ief->getVal(elid,face_loc);
				
				std::vector<int> v_faceTestFinal(4);	
				v_faceTestFinal[0] = ifn->getVal(face,0);
				v_faceTestFinal[1] = ifn->getVal(face,1);
				v_faceTestFinal[2] = ifn->getVal(face,2);
				v_faceTestFinal[3] = ifn->getVal(face,3);
				
				FaceSharedPtr faceTestFinalPtr = std::shared_ptr<NekFace>(new NekFace(v_faceTestFinal));
				//bp->m_QuadFaceHandled.insert(faceTestFinalPtr);
			}		
						
			bp->prisms.push_back(prsms[0]);
			bp->prisms.push_back(prsms[1]);
			
			bp->handledHexes.insert(elid);
			bp->handledQuads.insert(face);
			//bp->m_QuadFaceHandled.insert(faceTestPtr);
			
			int face_loc = FindOppositeLocalFaceHex(locFid);
			
			face = ief->getVal(elid,face_loc);
			
			int elid_prev 	= elid;
			elid0   		= ife->getVal(face,0);
			elid1   		= ife->getVal(face,1);
			
			if(elid0!=elid_prev)
			{
				elid = elid0;
			}
			else
			{
				elid = elid1;
			}
		}
	}
		
	return bp;
	
}







int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;
    
    if(world_size > 1)
    {
        if(world_rank == 0)
        {
            std::cout << "This application only runs in serial so try to rerun using: mpirun -np 1 ./ConvertHex2HybridMesh " << std::endl;
        }
        
        MPI_Finalize();
        
        return -1;
    }
    
    const char* fn_grid   = "inputs/grid.h5";
    const char* fn_conn   = "inputs/conn.h5";
    int ReadFromStats = 0;
    US3D* us3d    = ReadUS3DGrid(fn_conn,fn_grid,ReadFromStats,comm,info);

    //US3D* us3d  = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve       = us3d->xcn->getNglob();
    int Nel_part  = us3d->ien->getNrow();
    int Nel_glob  = us3d->ien->getNglob();
    
    std::cout << "Starting to extract the boundary map... " << std::endl;
    BoundaryMap* bmap = new BoundaryMap(us3d->ifn, us3d->if_ref);
    std::cout << "Finished extracting the boundary map... " << std::endl;
    
    
    int nPrmsNormal = 0;
    
    int wall_id     = 3;
    
    if (argc < 2)
    {
        std::cout << "ERROR: Please provide the desired number of prisms in normal direction with respect to the wall." << std::endl;
        return 1;
    }
    else
    {
        nPrmsNormal = atoi(argv[1]);
    }

	std::map<int,std::vector<int> > bnd_face_map = bmap->getBfaceMap();
	std::map<std::set<int>,int> tria_ref_map     = bmap->getTriaRefMap();
	std::map<std::set<int>,int> quad_ref_map     = bmap->getQuadRefMap();
	std::cout << "bnd_face_map " << bnd_face_map.size() << std::endl;
	std::map<int,int> vert_ref_map = bmap->getNodeRefMap();
       
	FaceSetPointer m_CutFaceSet;
	
	Prisms* bp0 = CreateBoundaryLayerPrisms(bnd_face_map, us3d->ien, 
    		us3d->ief, us3d->ife, us3d->ifn,
			tria_ref_map, quad_ref_map, wall_id,0,nPrmsNormal,m_CutFaceSet);   
   
	
	
    OutputBoundaryLayerPrisms(us3d->xcn, bp0, comm, "cut0_");
    
    //Extract all Hexes that need to be converted into tets;
    std::map<int,std::vector<int> > ToBeTets;
    
    for(int i=0;i<us3d->ien->getNrow();i++)
    {
    	if(bp0->handledHexes.find(i)==bp0->handledHexes.end())
    	{
    		std::vector<int> en(8);
    		for(int j=0;j<8;j++)
    		{
    			en[j] = us3d->ien->getVal(i,j);
    		}
    		ToBeTets[i] = en;
    	}
    }
    
	int m_rhorient[6][4] = {{0,3,2,1},
								{4,5,6,7},
								{0,4,7,3},
								{1,2,6,5},
								{0,1,5,4},
								{3,7,6,2}};
			
	int m_lhorient[6][4] = {{0,1,2,3},
							{4,7,6,5},
							{0,3,7,4},
							{1,5,6,2},
							{0,4,5,1},
							{3,2,6,5}};
    
    
    std::map<int,std::vector<int> >::iterator ittets;
    
    int flag2 = 0;
    
    

    Prisms* np = new Prisms;
    FaceSetPointer tetFaces;
    std::vector<Tetra*> Tts;
    int howmany = 0;
    std::cout << "ToBeTets  " << ToBeTets.size() << std::endl;
    for(ittets=ToBeTets.begin();ittets!=ToBeTets.end();ittets++)
    {
    	int hexID 			   = ittets->first;
    	std::vector<int> nodes = ittets->second;
    	int offset_ien 				= hexID*8;
		int* tmp_ien   				= us3d->ien->data + offset_ien;
		int offset_ief 				= hexID*6;
		int* tmp_ief   				= us3d->ief->data + offset_ief;
		
		
		std::map<int,FaceSharedPtr> facemap;
		FaceSetPointer cutFaces;
		std::map<int,FaceSharedPtr> cutFaces2;
		std::map<int,FaceSharedPtr> cutFacesLoc;
		std::map<int,int> cutLoc2FaceID;
		int flag = 0;
		std::map<int,FaceSharedPtr> locFaces;
		int N = 0;
		std::vector<int> fcut;
		
		FaceSetPointer hexFaces;
		
		
		// Analyze the hex:
		std::vector<int> orientations(6);
		for(int s=0;s<6;s++)
		{
			orientations[s] = -1;
			std::vector<int> v_face(4);
			int face = us3d->ief->getVal(hexID,s);
			
			for(int t=0;t<4;t++)
			{
				int lid 	= m_rhorient[s][t];
				v_face[t] 	= us3d->ien->getVal(hexID,lid);
			}
			
			FaceSharedPtr facePtr 				= std::shared_ptr<NekFace>(new NekFace(v_face));
			facemap[s] 							= facePtr;
			FaceSetPointer::iterator FPointer 	= m_CutFaceSet.find(facePtr);			
			
			locFaces[face] = facePtr;
			hexFaces.insert(facePtr);
			if(FPointer!=m_CutFaceSet.end()) // If face is found here, it means that the orieintation is reversed with respect to current element.
			{
				fcut.push_back(s);
				
				cutFaces.insert((*FPointer));
				cutFaces2[face]=(*FPointer);
				cutFacesLoc[s]=(*FPointer);
				cutLoc2FaceID[s] = face;
				orientations[s] = (*FPointer)->GetFaceRef();
				//std::cout << "orientations[s] From Face " << orientations[s] << std::endl;
				N++;
			}
			//std::cout << "orientations[s] " << orientations[s] << std::endl;
		}
		
		//std::cout << std::endl;
		// Find opposite face orientations:
		std::vector<std::vector<int> > pairOrientations(3);
		std::vector<int> pairOrientTag(3);
		for(int s=0;s<3;s++)
		{
			int f0 = orientations[s];
			int f1 = orientations[s+3];
			std::vector<int> pa(2);
			pa[0] = f0;
			pa[1] = f1;
			pairOrientations[s] = pa;
			pairOrientTag[s] = -1;
			
			
			if(pairOrientations[s][0]==0 && pairOrientations[s][1]==0)	
			{
				pairOrientTag[s] = 0;
			}
			if(pairOrientations[s][0]==1 && pairOrientations[s][1]==1)	
			{
				pairOrientTag[s] = 1;
			}
			if(pairOrientations[s][0]==0 && pairOrientations[s][1]==1)	
			{
				pairOrientTag[s] = 2;
			}
			if(pairOrientations[s][0]==1 && pairOrientations[s][1]==0)	
			{
				pairOrientTag[s] = 3;
			}
			// if cut is opposite and orientation is same
			
			// if cut is solo
			if(pairOrientations[s][0]==-1 && pairOrientations[s][1]==0 ||
			   pairOrientations[s][0]==0 && pairOrientations[s][1]==-1)
			{
				pairOrientTag[s] = 4;
			}
			if(pairOrientations[s][0]==-1 && pairOrientations[s][1]==1 ||
			   pairOrientations[s][0]==1 && pairOrientations[s][1]==-1)
			{
				pairOrientTag[s] = 5;
			}
			// no face in this direction is cut at all
			if(pairOrientations[s][0]==-1 && pairOrientations[s][1]==-1)
			{
				pairOrientTag[s] = 6;
			}
			
			if(pairOrientTag[s]==-1)
			{
				std::cout <<  "huh " << pairOrientations[s][0] << " " << pairOrientations[s][1] <<  " " << s+3 <<  std::endl;
				std::cout << "f 0 " << f0 << " " << f1 << std::endl;
			}
		}
		
		
		//std::cout << "N = " << N << std::endl;
		//std::cout << "N = " << N << " fcut " << fcut.size() << " " << m_CutFaceSet.size() << std::endl;
		int cutFaces_before = m_CutFaceSet.size();
		
		if(N==0)
		{
			int face			 		= us3d->ief->getVal(hexID,0);
			FaceSharedPtr facePtr   	= locFaces[face];
			int cutType 				= 0;
			int locFid					= -1;
							
			std::vector<Prism*> prsms 	= CutHexInPrism(tmp_ien, tmp_ief,facePtr,face,cutType,locFid,m_CutFaceSet);
			std::vector<Tetra*> tets    = CutPrismInTets(prsms[0],cutFaces,tetFaces,m_CutFaceSet,hexFaces);
			
		}
		if(N==1)
		{
//			CutPrismInTets(Prism* p, FaceSetPointer &cutFaces, FaceSetPointer &tetFaces)
			
			FaceSetPointer::iterator itc2;
			std::map<int,FaceSharedPtr>::iterator itc;
			for(itc=cutFaces2.begin();itc!=cutFaces2.end();itc++)
			{	
				int face			 		= itc->first;
				FaceSharedPtr facePtr   	= itc->second;
				int cutType 				= facePtr->GetFaceRef();
				int locFid					= -1;
				std::vector<Prism*> prsms 	= CutHexInPrism(tmp_ien, tmp_ief,facePtr,face,cutType,locFid,m_CutFaceSet);
				std::vector<Tetra*> tets0    = CutPrismInTets(prsms[0],cutFaces,tetFaces,m_CutFaceSet,hexFaces);
				Tts.push_back(tets0[0]);
				Tts.push_back(tets0[1]);
				std::vector<Tetra*> tets1    = CutPrismInTets(prsms[1],cutFaces,tetFaces,m_CutFaceSet,hexFaces);
				Tts.push_back(tets1[0]);
				Tts.push_back(tets1[1]);
				howmany++;
			}			
		}
		if(N==2)
		{
			int oppo_present = -1;
			int which_facePair = -1;
			for(int s=0;s<3;s++)
			{
				int tag = pairOrientTag[s];

				if(tag == 0 || tag == 1)
				{
					oppo_present = tag;
					which_facePair = s;
				}
				
				if(which_facePair != -1)
				{
					int cutType = tag;
					int locFid					= -1;
					int face 					= cutLoc2FaceID[which_facePair];
					FaceSharedPtr facePtr       = cutFacesLoc[which_facePair];
					std::vector<Prism*> prsms 	= CutHexInPrism(tmp_ien, tmp_ief,facePtr,face,cutType,locFid,m_CutFaceSet);
					std::vector<Tetra*> tets    = CutPrismInTets(prsms[0],cutFaces,tetFaces,m_CutFaceSet,hexFaces);
				}
				
			}
			
			
		}
		if(N==3)
		{

		}
		int cutFaces_after = m_CutFaceSet.size();
		int diffcutFaces   = cutFaces_after-cutFaces_before;
		//if(diffcutFaces > 1)
		//{
			//std::cout << "diff cut faces " << diffcutFaces << " " << N << std::endl;

		//}
		//std::cout << "tetFaces == " << tetFaces.size() << std::endl;
	
    }
    //std::cout << "Tts.size " << Tts.size()<< " " << howmany << std::endl;
    //std::cout << "np new " << np->prisms.size() << std::endl;
    MPI_Finalize();
    
}
