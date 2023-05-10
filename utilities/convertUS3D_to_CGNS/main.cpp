#include "../../src/adapt_io.h"
#include "../../src/adapt_boundary.h"
#include "../../src/hex2tet.h"
#include "../../src/adapt_bltopology.h"
#include "../../src/NekFace.h"

#include "cgnslib.h"






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
	
    const char* fn_grid   = "inputs_cgns/grid.h5";
    const char* fn_conn   = "inputs_cgns/conn.h5";
    int ReadFromStats = 0;
    US3D* us3d       = ReadUS3DGrid(fn_conn,fn_grid,ReadFromStats,comm,info);
    int nvrt_og      = us3d->xcn->getNrow();
    int nel 	     = us3d->ien->getNrow();
    int nf 	         = us3d->ife->getNrow();
    std::map<int, std::vector<int> > tetra;
    std::map<int, std::vector<int> > prism;
    
    for(int i=0;i<nel;i++)
    {
    	int type = us3d->iet->getVal(i,0);
    	
    	if(type == 6)// prism
    	{
    		std::vector<int> pri(6);
    		
    		for(int j=0;j<6;j++)
    		{
    			pri[j] = us3d->ien->getVal(i,j);
    		}
    		
    		prism[i] = pri;
    		
    	}
    	if(type == 2)// tet
		{
    		std::vector<int> tet(4);
    		for(int j=0;j<4;j++)
    		{
    			tet[j] = us3d->ien->getVal(i,j);
    		}
    		tetra[i] = tet;
		}
    }
    
    FaceSetPointer m_T_P_Interface;
    std::map<int,int> interface_nodes;
    int loc_vint = 0;
    std::map<int,int> bndMap;
    std::map<int,std::set<int> > bndMapRef;
    int intRef = 0;
    int intfid = 0;
    std::map<int,std::vector<int> > interface_Fcs;
    std::map<int,int> interface_map;
    for(int i=0;i<nf;i++)
	{
    	int el_l = us3d->ife->getVal(i,0);
    	int el_r = us3d->ife->getVal(i,1);
    	
    	if(el_l < nel && el_r < nel)
    	{
    		int el_l_type = us3d->iet->getVal(el_l,0);
			int el_r_type = us3d->iet->getVal(el_r,0);
			
			if((el_l_type==2 && el_r_type==6) ||
			   (el_l_type==6 && el_r_type==2))
			{
				std::vector<int> Fcs(3);
				for(int j=0;j<3;j++)
				{
					int vint = us3d->ifn->getVal(i,j);
					
					if(interface_nodes.find(vint)==interface_nodes.end())
					{
						interface_nodes[vint] = loc_vint;
						
						loc_vint++;
					}
					else
					{
						int loc_vint_n = interface_nodes[vint];
					}
					
					Fcs[j] = vint;
				}
				
				interface_map[i] = intfid;
				interface_Fcs[i] = Fcs;
				FaceSharedPtr RefFacePointer = std::shared_ptr<NekFace>(new NekFace(Fcs));
				RefFacePointer->SetFaceRef(intRef);
				pair<FaceSetPointer::iterator, bool> testInsPointer;
				testInsPointer = m_T_P_Interface.insert(RefFacePointer);
				intfid++;
				intRef++;
			}
    	}

    	int ref = us3d->if_ref->getVal(i,0);
    	
    	if(ref != 2)
    	{
    		int Nv = us3d->if_Nv->getVal(i,0);
    		for(int j=0;j<Nv;j++)
			{
				int vI = us3d->ifn->getVal(i,j);
				
				if(bndMapRef[ref].find(vI)==bndMapRef[ref].end())
				{
					bndMapRef[ref].insert(vI);
				}
				
//				if(bndMap.find(vI)==bndMap.end())
//				{
//					bndMap[vI]=ref;
//				}	
			}
    	}
	}

    delete us3d->ife;
    delete us3d->iet;
    delete us3d->ghost;
    delete us3d->ien;
    delete us3d->iee;
    delete us3d->ie_Nv;


//  delete us3d->xcn;
    
    std::map<int,std::vector<int> >::iterator itm;
//    std::map<int,int>::iterator itr;
//    for(itr=interface_nodes.begin();itr!=interface_nodes.end();itr++)
//    {
//    	std::cout << "interface_nodes "<<itr->first << " " << itr->second << std::endl;
//    }
    
	//=============================================================
	//=============================================================
	//=============================================================
    

    cgsize_t n_prism_start = 1;
    cgsize_t * prism_write = new cgsize_t[prism.size()*6];
    //cgsize_t prism_write[prism.size()][6];
    std::map<int,int> uvrt_g2l_prism;
    int loc_vrt_id_prism = 0;
    int loc_prism = 0;
    
    std::vector<int> coordsPrisms;
    //std::vector<std::vector<double> > coordsPrisms;
    int* interface_prism = new int[interface_nodes.size()];
    std::map<int,std::vector<int> > bnd_prism;
    std::set<int> bnd_set_p;
    std::map<int,std::vector<int> > IntFc_Tri_Prism;
    std::map<int,std::vector<std::vector<int> > > BndFc_Tri_Prism;
    std::map<int,std::vector<std::vector<int> > > BndFc_Quad_Prism;
    
	for(itm = prism.begin();itm!=prism.end();itm++)
	{
		for(int j=0;j<itm->second.size();j++)
		{
			if(uvrt_g2l_prism.find(itm->second[j])==uvrt_g2l_prism.end())
			{
				uvrt_g2l_prism[itm->second[j]] = loc_vrt_id_prism+1;
				prism_write[loc_prism*6+j] = loc_vrt_id_prism+1;
				
//				std::vector<double> coords(3);
//				coords[0] = us3d->xcn->getVal(itm->second[j],0);
//				coords[1] = us3d->xcn->getVal(itm->second[j],1);
//				coords[2] = us3d->xcn->getVal(itm->second[j],2);
				coordsPrisms.push_back(itm->second[j]);
				
				if(interface_nodes.find(itm->second[j])!=interface_nodes.end())
				{
					int loc_i = interface_nodes[itm->second[j]];
					interface_prism[loc_i] = loc_vrt_id_prism+1;
				}
				
				loc_vrt_id_prism++;
			}
			else
			{
				int loc_v = uvrt_g2l_prism[itm->second[j]];
				prism_write[loc_prism*6+j] = loc_v;
			}
		}
		
		loc_prism++;
	}
	
	
	
	
	std::vector<int> IntFace_Prism(m_T_P_Interface.size()*3);
	
	int tell = 0;
	for(itm = prism.begin();itm!=prism.end();itm++)
	{
		int Nf = us3d->ie_Nf->getVal(itm->first,0);
		if(Nf != 5)
		{
			std::cout << "ERROR: rhis is not a prism." << std::endl;
		}	
		for(int k=0;k<Nf;k++)
		{
			int fid        = us3d->ief->getVal(itm->first,k);
			int fref       = us3d->if_ref->getVal(fid,0);
			
			std::vector<int> IntFc_Prism(3);
			
			if(interface_Fcs.find(fid)!=interface_Fcs.end())
			{
				int intloc_fid = interface_map[fid];
				
				for(int p=0;p<3;p++)
				{
					int vid           			  = interface_Fcs[fid][p];
					int local_v       			  = uvrt_g2l_prism[vid];
					IntFace_Prism[intloc_fid*3+p] = local_v;
					IntFc_Prism[p]    			  = local_v;
				}
				
				IntFc_Tri_Prism[fid] = IntFc_Prism;
			}
			
			
			if(fref != 2)
			{
				int f_Nv = us3d->if_Nv->getVal(fid,0);
				std::vector<int> BndFc(f_Nv);
				
				for(int p=0;p<f_Nv;p++)
				{
					int vid = us3d->ifn->getVal(fid,p);
					int local_v = uvrt_g2l_prism[vid];
					BndFc[p] = local_v;
					
					if(bndMapRef[fref].find(vid)!=bndMapRef[fref].end())
					{
						bnd_prism[fref].push_back(local_v);
					}
				}
				if(f_Nv == 3)
				{
					BndFc_Tri_Prism[fref].push_back(BndFc);
					tell++;
				}
				if(f_Nv == 4)
				{
					BndFc_Quad_Prism[fref].push_back(BndFc);
				}
			}
		}
	}
	
	
	std::map<int,std::vector<std::vector<int> > >::iterator itbt;

	std::map<int,std::vector<int> >::iterator itbb;
		

	//=============================================================
	//=============================================================
	//=============================================================
	
	
	cgsize_t n_tetra_start = 0;
	//cgsize_t tetra_write[tetra.size()][4];
	cgsize_t* tetra_write = new cgsize_t[tetra.size()*4];
	
	std::map<int,int> uvrt_g2l_tetra;
	int loc_vrt_id_tetra = 0;
	int loc_tetra = 0;
	std::vector<int> coordstetras;
	//std::vector<std::vector<double> > coordstetras;
    int* interface_tetra = new int[interface_nodes.size()];
    std::map<int,std::vector<int> > bnd_tetra;
    std::map<int,std::vector<int> > IntFc_Tri_Tetra;

    std::set<int> bnd_set_t;
    std::map<int,std::vector<std::vector<int> > > BndFc_Tri_Tetra;
    
	for(itm = tetra.begin();itm!=tetra.end();itm++)
	{
		for(int j=0;j<itm->second.size();j++)
		{
			if(uvrt_g2l_tetra.find(itm->second[j])==uvrt_g2l_tetra.end())
			{
				uvrt_g2l_tetra[itm->second[j]] = loc_vrt_id_tetra+1;
				tetra_write[loc_tetra*4+j] = loc_vrt_id_tetra+1;
				
//				std::vector<double> coords(3);
//				coords[0] = us3d->xcn->getVal(itm->second[j],0);
//				coords[1] = us3d->xcn->getVal(itm->second[j],1);
//				coords[2] = us3d->xcn->getVal(itm->second[j],2);
				coordstetras.push_back(itm->second[j]);
				if(interface_nodes.find(itm->second[j])!=interface_nodes.end())
				{
					int loc_i = interface_nodes[itm->second[j]];
					interface_tetra[loc_i] = loc_vrt_id_tetra+1;
				}
				
				loc_vrt_id_tetra++;
			}
			else
			{
				int loc_v = uvrt_g2l_tetra[itm->second[j]];
				tetra_write[loc_tetra*4+j] = loc_v;
			}
			
			if(bndMap.find(itm->second[j])!=bndMap.end() &&
					bnd_set_t.find(itm->second[j])==bnd_set_t.end())
			{
				bnd_set_t.insert(itm->second[j]);
				int loc_v = uvrt_g2l_tetra[itm->second[j]];
				int ref   = bndMap[itm->second[j]];
				bnd_tetra[ref].push_back(loc_v);
			}
		}
		
				
		loc_tetra++;
	}
	
	
	std::vector<int> IntFace_Tetra(m_T_P_Interface.size()*3);
	
	for(itm = tetra.begin();itm!=tetra.end();itm++)
	{
		int Nf = us3d->ie_Nf->getVal(itm->first,0);
		if(Nf != 4)
		{
			std::cout << "ERROR: rhis is not a tetra." << std::endl;
		}	
		
		
		
		for(int k=0;k<Nf;k++)
		{
			int fid  = us3d->ief->getVal(itm->first,k);
			int fref = us3d->if_ref->getVal(fid,0);
			
			if(interface_Fcs.find(fid)!=interface_Fcs.end())
			{
				int intloc_fid = interface_map[fid];
				std::vector<int> IntFc_Tetra(3);
				for(int p=0;p<3;p++)
				{
					int vid           			  = interface_Fcs[fid][p];
					int local_v       			  = uvrt_g2l_tetra[vid];
					IntFace_Tetra[intloc_fid*3+p] = local_v;
					IntFc_Tetra[p]    			  = local_v;
				}
				
				IntFc_Tri_Tetra[fid] = IntFc_Tetra;
			}
			
			
			if(fref != 2)
			{
				int f_Nv = us3d->if_Nv->getVal(fid,0);
				std::vector<int> BndFc(f_Nv);
				
				for(int p=0;p<f_Nv;p++)
				{
					int vid = us3d->ifn->getVal(fid,p);
					int local_v = uvrt_g2l_tetra[vid];
					BndFc[p] = local_v;
				}
				
				BndFc_Tri_Tetra[fref].push_back(BndFc);
				
			}			
		}
	}
	
	
	delete us3d->ifn;
    delete us3d->if_Nv;
	delete us3d->ie_Nf;
	delete us3d->if_ref;
	
	cgsize_t n_prism_end = loc_prism;
		
	int nvrt_prism = uvrt_g2l_prism.size();

//	double coordXPrism[nvrt_prism];
//	double coordYPrism[nvrt_prism];
//	double coordZPrism[nvrt_prism];
	
	double* coordXPrism = new double[nvrt_prism];
	double* coordYPrism = new double[nvrt_prism];
	double* coordZPrism = new double[nvrt_prism];

	for(int i=0;i<nvrt_prism;i++)
	{
		coordXPrism[i] = us3d->xcn->getVal(coordsPrisms[i],0);
		coordYPrism[i] = us3d->xcn->getVal(coordsPrisms[i],1);
		coordZPrism[i] = us3d->xcn->getVal(coordsPrisms[i],2);
	}
		
		
		
	
	cgsize_t n_tetra_end = loc_tetra;
		
	int nvrt_tetra = uvrt_g2l_tetra.size();
	double* coordXtetra = new double[nvrt_tetra];
	double* coordYtetra = new double[nvrt_tetra];
	double* coordZtetra = new double[nvrt_tetra];

	for(int i=0;i<nvrt_tetra;i++)
	{
		
		coordXtetra[i] = us3d->xcn->getVal(coordstetras[i],0);
		coordYtetra[i] = us3d->xcn->getVal(coordstetras[i],1);
		coordZtetra[i] = us3d->xcn->getVal(coordstetras[i],2);
//		
//		coordXtetra[i] = coordstetras[i][0];
//		coordYtetra[i] = coordstetras[i][1];
//		coordZtetra[i] = coordstetras[i][2];
	}
	
	std::map<int,int>::iterator it;

	

	
	

	int n_tria_tet_tot = 0;
	
	for(itbt=BndFc_Tri_Tetra.begin();itbt!=BndFc_Tri_Tetra.end();itbt++)
	{
		n_tria_tet_tot = n_tria_tet_tot+itbt->second.size();
	}
	
	
	int n_tria_prism_tot = 0;
	
	std::map<int,std::vector<int> >::iterator itbi;

	n_tria_prism_tot = n_tria_prism_tot + IntFc_Tri_Prism.size();
	for(itbt=BndFc_Tri_Prism.begin();itbt!=BndFc_Tri_Prism.end();itbt++)
	{
		n_tria_prism_tot = n_tria_prism_tot+itbt->second.size();
	}
	int n_quad_prism_tot = 0;
	for(itbt=BndFc_Quad_Prism.begin();itbt!=BndFc_Quad_Prism.end();itbt++)
	{
		n_quad_prism_tot = n_quad_prism_tot+itbt->second.size();
	}
	
	std::map<int,std::vector<int> > BndFc_Tri_Tetra_Conv;
	for(itbt=BndFc_Tri_Tetra.begin();itbt!=BndFc_Tri_Tetra.end();itbt++)
	{
		int ref = itbt->first;
		std::vector<int> tria_per_ref(itbt->second.size()*3);
		
		for(int q=0;q<itbt->second.size();q++)
		{
			
			for(int p=0;p<3;p++)
			{
				tria_per_ref[q*3+p] = itbt->second[q][p];
				//std::cout << itbt->second[q][p] << " ";
			}
			//std::cout << std::endl;
		}
		
		BndFc_Tri_Tetra_Conv[ref] = tria_per_ref;
	}
	
	
	
	
	
	std::map<int,std::vector<int> > BndFc_Tri_Prism_Conv;
	for(itbt=BndFc_Tri_Prism.begin();itbt!=BndFc_Tri_Prism.end();itbt++)
	{
		int ref = itbt->first;
		std::vector<int> tria_per_ref(itbt->second.size()*3);
		
		for(int q=0;q<itbt->second.size();q++)
		{
			for(int p=0;p<3;p++)
			{
				tria_per_ref[q*3+p] = itbt->second[q][p];
				//std::cout << itbt->second[q][p] << " ";
			}
			//std::cout << std::endl;
		}
		
		BndFc_Tri_Prism_Conv[ref] = tria_per_ref;
	}
	
	std::map<int,std::vector<int> > BndFc_Quad_Prism_Conv;
	for(itbt=BndFc_Quad_Prism.begin();itbt!=BndFc_Quad_Prism.end();itbt++)
	{
		int ref = itbt->first;
		std::vector<int> quad_per_ref(itbt->second.size()*4);
		
		for(int q=0;q<itbt->second.size();q++)
		{
			for(int p=0;p<4;p++)
			{
				quad_per_ref[q*4+p] = itbt->second[q][p];
			}
		}
		BndFc_Quad_Prism_Conv[ref] = quad_per_ref;
	}
	

	//=============================================================
	//=============================================================

	  	  
	   int size_prism[3] = {nvrt_prism,loc_prism, n_tria_prism_tot+n_quad_prism_tot};
	   
	   //std::cout <<  "size_prism " <<   size_prism[0] << " " << size_prism[1] << " " << size_prism[2] << std::endl;
//	   size_prism[0] = nvrt_prism;
//	   size_prism[1] = loc_prism-1;
//       size_prism[2] = 0;

//	   std::cout << "Zone dimensions " << isize[0] << " " << isize[1] << " " << isize[2] << std::endl;
	   int index_file,icelldim,iphysdim,index_base;
	   int index_zone_p,index_coordx_p,index_coordy_p,index_coordz_p;
	   int index_zone_t,index_coordx_t,index_coordy_t,index_coordz_t;
	   char basename[33],zonename[33];
	   int indexe;
	   

	   if (cg_open("MIMIC_mesh.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
	   	   
	   strcpy(basename,"Base");
	   icelldim=3;
	   iphysdim=3;
	   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
	   int ier;
	   
	   cg_goto(index_file,index_base,NULL);
//	   std::string version = "Dirk wrote this";
	   cg_descriptor_write("Dirk wrote this", "Information");
	   std::map<int,std::vector<int> >::iterator itbc;
	   float exp[5];
	   for (int n = 0; n < 5; n++)
	   {
		   exp[n] = (float)0.0;
	   }
	   exp[1] = (float)1.0;
	   int n_tria_offset = 0;//n_prism_end;
	   int nbnd=0;
	   int indextr3;
	   int indextr7;
	   int indextr10;
	   int indextr36;
	   
	   int index_fam3;
	   int bc_index3;
	   
	   int index_fam7;
	   int bc_index7;
	   
	   int index_fam10;
	   int bc_index10;
	   
	   int index_fam36;
	   int bc_index36;
	   std::map<int,std::vector<int> >::iterator itb;
	   int indextrint;
	   
	   int index_conn = 0;
	   int* transform = new int[3];
	   transform[0] = 1;
	   transform[1] = 2;
	   transform[2] = 3;
	   int icount = interface_nodes.size();
	   
	   strcpy(zonename,"blk-1");
	   cg_zone_write(index_file, index_base, zonename, size_prism, Unstructured, &index_zone_p);	
//	   
	   cg_coord_write(index_file,index_base,index_zone_p,RealDouble,"CoordinateX", coordXPrism,&index_coordx_p);
	   cg_goto(index_file, index_base, "blk-1", 0, "GridCoordinates", 0, "CoordinateX", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   cg_coord_write(index_file,index_base,index_zone_p,RealDouble,"CoordinateY", coordYPrism,&index_coordy_p);
	   cg_goto(index_file, index_base, "blk-1", 0, "GridCoordinates", 0, "CoordinateY", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   cg_coord_write(index_file,index_base,index_zone_p,RealDouble,"CoordinateZ", coordZPrism,&index_coordz_p);
	   cg_goto(index_file, index_base, "blk-1", 0, "GridCoordinates", 0, "CoordinateZ", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   
	   cg_goto(index_file, index_base, "blk-1", 0, "GridCoordinates_t", 1, NULL);
	   cg_dataclass_write(Dimensional);
	   cg_units_write(MassUnitsUserDefined, LengthUnitsUserDefined,
						  TimeUnitsUserDefined, TemperatureUnitsUserDefined,
					  AngleUnitsUserDefined);
	   	   
	   	   
	   cg_section_write(index_file, index_base, index_zone_p, 
			   "PrismElements", PENTA_6, n_prism_start, n_prism_end, 0,
			   prism_write, &indexe);
	   
	   //std::cout << "range for PrismElements :: " <<  n_prism_start << " - " << n_prism_end << std::endl;
	   

	   
	   int ntria_int_prism = IntFc_Tri_Prism.size();
	   std::vector<int> int_prism_vrts(ntria_int_prism*3);
	   int p=0;
	   for(itb=IntFc_Tri_Prism.begin();itb!=IntFc_Tri_Prism.end();itb++)
	   {
		   int ref = itb->first;

		   for(int j=0;j<itb->second.size();j++)
		   {
			   int_prism_vrts[p*3+j] = itb->second[j];
		   }  
		   p++;
	   }
	   
	

	   

	   //std::cout << "ntria_int_prism " << ntria_int_prism << std::endl;
	   

	   cg_section_write(index_file, index_base, index_zone_p, "P-T_Interface", 
	   					TRI_3, n_tria_offset+1, n_tria_offset+ntria_int_prism, 0, IntFace_Prism.data(), &indextrint);

	   n_tria_offset = n_tria_offset+ntria_int_prism;
	   

	   
	   for(itb=BndFc_Tri_Prism_Conv.begin();itb!=BndFc_Tri_Prism_Conv.end();itb++)
	   {

		   int ref = itb->first;
		   int n_tria_bc = BndFc_Tri_Prism[ref].size();
		   const char* bcname;
		   std::vector<int> bc_interface = bnd_prism[ref];
		   int bc_intface_size = bc_interface.size();
		   
//		   std::cout << "ref - " << ref << " -> ";
//		   for(int q=0;q<itb->second.size();q++)
//		   {
//			   std::cout << itb->second[q] << " ";
//		   }
//		   
//		   std::cout << std::endl;
		   
		   std::vector<int>ElList(n_tria_bc);
		   int tid = n_tria_offset+1;
		   for(int q=0;q<n_tria_bc;q++)
		   {
			   ElList[q] = tid;
			   tid++;
		   }
		   
		   if(ref == 3)
		   {
			   int indextr3;
			   bcname="BCWallViscous_PT";
				
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCWallViscous,
					   ElementList, n_tria_bc, ElList.data(), &indextr3);
						   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr3, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr3, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam3);
			   cg_fambc_write(index_file, index_base, index_fam3, "FamBC", BCWallViscous, &bc_index3);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc, 
								0, itb->second.data(), &indextr3);
   //			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << " bc_intface_size " << bc_intface_size << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   

		   }
		   if(ref == 7)
		   {
			   bcname="BCSymmetryPlane_PT";
			   
			   int indextr7;
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCSymmetryPlane,
					   ElementList, n_tria_bc, ElList.data(), &indextr7);
			   			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr7, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr7, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam7);
			   cg_fambc_write(index_file, index_base, index_fam7, "FamBC", BCSymmetryPlane, &bc_index7);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
								TRI_3,n_tria_offset+1, n_tria_offset+n_tria_bc,
			   	   	   	   	    0, itb->second.data(), &indextr7);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   

			   
		   }
		   if(ref == 10)
		   {
			   int indextr10;
			   bcname="BCInflow_PT";
				
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCInflow,
					   ElementList, n_tria_bc, ElList.data(), &indextr10);
			   			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr10, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr10, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam10);
			   cg_fambc_write(index_file, index_base, index_fam10, "FamBC", BCInflow, &bc_index10);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc, 
			   					0, itb->second.data(), &indextr10);
//			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << " bc_intface_size " << bc_intface_size << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   

		   }
		   if(ref == 36)
		   {
			   int indextr36;
			   bcname="BCOutflow_PT";
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCOutflow,
					   	   	   ElementList, n_tria_bc, ElList.data(), &indextr36);
			   			   			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr36, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr36, "end");
			   
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam36);
			   cg_fambc_write(index_file, index_base, index_fam36, "FamBC", BCOutflow, &bc_index36);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc,
			   	   	   	   	    0, itb->second.data(), &indextr36);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   

		   }


		   
		   
	   }
	   
	   
	   
	   for(itb=BndFc_Quad_Prism_Conv.begin();itb!=BndFc_Quad_Prism_Conv.end();itb++)
	   {

		   int ref = itb->first;
		   int n_tria_bc = BndFc_Quad_Prism[ref].size();
		   std::vector<int> bc_interface = bnd_prism[ref];
		   int bc_intface_size = bc_interface.size();
		   const char* bcname;
		   
		   std::vector<int>ElList(n_tria_bc);
		   int tid = n_tria_offset+1;
		   for(int q=0;q<n_tria_bc;q++)
		   {
			   ElList[q] = tid;
			   tid++;
		   }
		   			   
		   			   
		   if(ref == 3)
		   {
			   bcname="BCWallViscous_PQ";
			   int indextr3;
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCWallViscous,
			   			      ElementList, n_tria_bc, ElList.data(), &indextr3);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr3, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr3, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam3);
			   cg_fambc_write(index_file, index_base, index_fam3, "FamBC", BCWallViscous, &bc_index3);
			   
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
					   QUAD_4, n_tria_offset+1, n_tria_offset+n_tria_bc, 0, itb->second.data(), &indextr3);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   			   
//			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   
		   }
		   if(ref == 7)
		   {
			   bcname="BCSymmetryPlane_PQ";
			   
			   int indextr7;
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCSymmetryPlane,
					   	   	   ElementList, n_tria_bc, ElList.data(), &indextr7);
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr7, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr7, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam7);
			   cg_fambc_write(index_file, index_base, index_fam7, "FamBC", BCSymmetryPlane, &bc_index7);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
					   QUAD_4,n_tria_offset+1, n_tria_offset+n_tria_bc, 0, itb->second.data(), &indextr7);
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;

		   }
		   if(ref == 10)
		   {
			   
			   bcname="BCInflow_PQ";
			   
			   int indextr10;
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCInflow,
					   	   	   ElementList, n_tria_bc, ElList.data(), &indextr10);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr10, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr10, "end");
			   
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam10);
			   cg_fambc_write(index_file, index_base, index_fam10, "FamBC", BCInflow, &bc_index10);
				
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
					   QUAD_4, n_tria_offset+1, n_tria_offset+n_tria_bc, 0, itb->second.data(), &indextr10);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   n_tria_offset = n_tria_offset+n_tria_bc;
		   }
		   if(ref == 36)
		   {
			   bcname="BCOutflow_PQ";
	  //		
			   int indextr36;
			   
			   
			   
			   
			   cg_boco_write (index_file, index_base, index_zone_p, bcname, BCOutflow,
					   	   	   ElementList, n_tria_bc, ElList.data(), &indextr36);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_p, indextr36, CGNS_ENUMV(Vertex));
			   
			   cg_goto(index_file, 1, "blk-1", 0, "ZoneBC_t", 1, "BC_t", indextr36, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam36);
			   cg_fambc_write(index_file, index_base, index_fam10, "FamBC", BCOutflow, &bc_index36);
			   
			   cg_section_write(index_file, index_base, index_zone_p, bcname, 
					   QUAD_4, n_tria_offset+1, n_tria_offset+n_tria_bc, 0, itb->second.data(), &indextr36);
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << " bc_intface_size " << bc_intface_size << std::endl;

			   n_tria_offset = n_tria_offset+n_tria_bc;
			   
			   nbnd++;
		   }
	   
	   }
	   
	   //=================================================================
	   //=================================================================
	   //=================================================================
	   //=================================================================
	   
	   
	   int size_tetra[3] = {nvrt_tetra, n_tetra_end, n_tria_tet_tot};
	   //std::cout <<  "size_tetra " <<   size_tetra[0] << " " << size_tetra[1] << " " << size_tetra[2] << std::endl;

	   //std::cout << " nvrt_tetra " <<  nvrt_tetra << " " << n_tetra_end << " " << n_tria_tet_tot  << std::endl;
	   strcpy(zonename,"blk-2");
	  
	   cg_zone_write(index_file, index_base, zonename, size_tetra, Unstructured, &index_zone_t);
	   
	   cg_coord_write(index_file,index_base,index_zone_t,RealDouble,"CoordinateX", coordXtetra,&index_coordx_t);
	   cg_goto(index_file, index_base, "blk-2", 0, "GridCoordinates", 0, "CoordinateX", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   cg_coord_write(index_file,index_base,index_zone_t,RealDouble,"CoordinateY", coordYtetra,&index_coordy_t);
	   cg_goto(index_file, index_base, "blk-2", 0, "GridCoordinates", 0, "CoordinateY", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   cg_coord_write(index_file,index_base,index_zone_t,RealDouble,"CoordinateZ", coordZtetra,&index_coordz_t);
	   cg_goto(index_file, index_base, "blk-2", 0, "GridCoordinates", 0, "CoordinateZ", 0, NULL);
	   cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
	   
	   cg_goto(index_file, index_base, "blk-2", 0, "GridCoordinates_t", 1, NULL);
	   cg_dataclass_write(Dimensional);
	   cg_units_write(MassUnitsUserDefined, LengthUnitsUserDefined,
						  TimeUnitsUserDefined, TemperatureUnitsUserDefined,
					  AngleUnitsUserDefined);
		   
	   cg_section_write(index_file, index_base, index_zone_t, 
			   "TetraElements", TETRA_4, 1, n_tetra_end, 0,
			   tetra_write, &indexe);
	   
	   n_tria_offset = n_tetra_end;
	   
	   int ntria_int_tetra = IntFc_Tri_Tetra.size();
	   //std::cout << "ntria_int_tetra " << ntria_int_tetra << std::endl;
	   cg_section_write(index_file, index_base, index_zone_t, "T-P_Interface", 
	   	   				TRI_3, n_tria_offset+1, n_tria_offset+ntria_int_tetra, 0, IntFace_Tetra.data(), &indextrint);
	   
	   n_tria_offset = n_tria_offset+ntria_int_tetra;
	   
	   
	   for(itb=BndFc_Tri_Tetra_Conv.begin();itb!=BndFc_Tri_Tetra_Conv.end();itb++)
	   {
		   int ref = itb->first;
		   int n_tria_bc = BndFc_Tri_Tetra[ref].size();
		   const char* bcname;
		   std::vector<int> bc_interface = bnd_prism[ref];
		   int bc_intface_size = bc_interface.size();

		   std::vector<int>ElList(n_tria_bc);
		   int tid = n_tria_offset+1;
		   for(int q=0;q<n_tria_bc;q++)
		   {
			   ElList[q] = tid;
			   tid++;
		   }
		   
		   if(ref == 3)
		   {
			   int indextr3;
			   bcname="BCWallViscous_TT";
				
			   cg_boco_write (index_file, index_base, index_zone_t, 
					   	   	  bcname, BCWallViscous,
							  ElementList, n_tria_bc, ElList.data(), 
							  &indextr3);
						   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, indextr3, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", indextr3, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam3);
			   cg_fambc_write(index_file, index_base, index_fam3, "FamBC", BCWallViscous, &bc_index3);
			   
			   cg_section_write(index_file, index_base, index_zone_t, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc, 
								0, itb->second.data(), &indextr3);
	  //			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << " bc_intface_size " << bc_intface_size << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;

		   }
		   if(ref == 7)
		   {
			   bcname="BCSymmetryPlane_TT";
			   
			   int indextr7;
			   
			   cg_boco_write (index_file, index_base, index_zone_t, bcname, BCSymmetryPlane,
					   ElementList, n_tria_bc, ElList.data(), &indextr7);
						   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, indextr7, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", indextr7, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam7);
			   cg_fambc_write(index_file, index_base, index_fam7, "FamBC", BCSymmetryPlane, &bc_index7);
			   
			   cg_section_write(index_file, index_base, index_zone_t, bcname, 
								TRI_3,n_tria_offset+1, n_tria_offset+n_tria_bc,
								0, itb->second.data(), &indextr7);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;

		   }
		   if(ref == 10)
		   {
			   int indextr10;
			   bcname="BCInflow_TT";
				
			   cg_boco_write (index_file, index_base, index_zone_t, bcname, BCInflow,
					   ElementList, n_tria_bc, ElList.data(), &indextr10);
						   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, indextr10, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", indextr10, "end");
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam10);
			   cg_fambc_write(index_file, index_base, index_fam10, "FamBC", BCInflow, &bc_index10);
			   
			   cg_section_write(index_file, index_base, index_zone_t, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc, 
								0, itb->second.data(), &indextr10);
   //			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << " bc_intface_size " << bc_intface_size << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;

		   }
		   if(ref == 36)
		   {
			   int indextr36;
			   bcname="BCOutflow_TT";
			   
			   cg_boco_write (index_file, index_base, index_zone_t, bcname, BCOutflow,
							   ElementList, n_tria_bc, ElList.data(), &indextr36);
									   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, indextr36, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", indextr36, "end");
			   
			   
			   cg_family_write(index_file, index_base,bcname,&index_fam36);
			   cg_fambc_write(index_file, index_base, index_fam36, "FamBC", BCOutflow, &bc_index36);
			   
			   cg_section_write(index_file, index_base, index_zone_t, bcname, 
								TRI_3, n_tria_offset+1, n_tria_offset+n_tria_bc,
								0, itb->second.data(), &indextr36);
			   
			   std::cout << "range for " << bcname << " :: " << n_tria_offset+1 << " - " << n_tria_offset+n_tria_bc << std::endl;
			   
			   n_tria_offset = n_tria_offset+n_tria_bc;
			   

		   }
	   }
//	   int index_conn = 0;
//	   int* transform = new int[3];
//	   transform[0] = 1;
//	   transform[1] = 2;
//	   transform[2] = 3;
	   
	   icount = interface_nodes.size();
	   //cg_1to1_write(index_file,index_base,index_zone_t,"Interior","blk-1",interface_tetra,interface_prism,transform,&index_conn);
	   ier = cg_conn_write(index_file,index_base,index_zone_t,"Interface",Vertex,Abutting1to1,
	                 PointList,icount,interface_tetra,"blk-1",Unstructured,
	                 PointListDonor,Integer,icount,interface_prism,&index_conn);
	   
	   /*
	   std::cout << "indextr10 " << indextr10 << " indextr36 " << indextr36 << std::endl;  
	   
	   int index_conn = 0;
	   int* transform = new int[3];
	   transform[0] = 1;
	   transform[1] = 2;
	   transform[2] = 3;
	   int icount = interface_nodes.size();
	   //cg_1to1_write(index_file,index_base,index_zone_t,"Interior","blk-1",interface_tetra,interface_prism,transform,&index_conn);
	   ier = cg_conn_write(index_file,index_base,index_zone_t,"Interface",Vertex,Abutting1to1,
	                 PointList,icount,interface_tetra,"blk-1",Unstructured,
	                 PointListDonor,Integer,icount,interface_prism,&index_conn);
	   
	  
	   std::cout << "indexe " << indexe << std::endl;
	   
	   for(itbc=bnd_tetra.begin();itbc!=bnd_tetra.end();itbc++)
	   {
		   int cgbc_wall;
		   
		   if(itbc->first == 3)
		   {
			   cg_boco_write (index_file, index_base, index_zone_t, "BCWallViscous", BCWallViscous,
								   PointList, itbc->second.size(), itbc->second.data(), &cgbc_wall);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, cgbc_wall, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", cgbc_wall, "end");
			   //cg_user_data_write("FamBC_UserId", CGNS_ENUMV(Character), &index_dim, "FamBC_UserId_WallViscous");
			   cg_descriptor_write("FamBC_UserId", "FamBC_UserId_WallViscous");

			   
		   }
		   if(itbc->first == 7)
		   {
			   cg_boco_write (index_file, index_base, index_zone_t, "BCSymmetryPlane", BCSymmetryPlane,
								   PointList, itbc->second.size(), itbc->second.data(), &cgbc_wall);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, cgbc_wall, CGNS_ENUMV(Vertex));
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", cgbc_wall, "end");
			   //cg_user_data_write("FamBC_UserId", CGNS_ENUMV(Character), &index_dim, "FamBC_UserId_SymmetryPlane");
			   cg_descriptor_write("FamBC_UserId", "FamBC_UserId_SymmetryPlane");
		   }
		   if(itbc->first == 10)
		   {
			   cg_boco_write (index_file, index_base, index_zone_t, "BCInflow", BCInflow,
								   PointList, itbc->second.size(), itbc->second.data(), &cgbc_wall);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, cgbc_wall, CGNS_ENUMV(Vertex));
			   int index_dim = 1;
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", cgbc_wall, "end");
			   //cg_user_data_write("FamBC_UserId", CGNS_ENUMV(Character), &index_dim, "FamBC_UserId_Inflow");
			   cg_descriptor_write("FamBC_UserId", "FamBC_UserId_Inflow");

		   }
		   if(itbc->first == 36)
		   {
			   cg_boco_write (index_file, index_base, index_zone_t, "BCOutflow", BCOutflow,
								   PointList, itbc->second.size(), itbc->second.data(), &cgbc_wall);
			   
			   cg_boco_gridlocation_write(index_file, index_base, index_zone_t, cgbc_wall, CGNS_ENUMV(Vertex));
			   int index_dim = 1;
			   cg_goto(index_file, 1, "blk-2", 0, "ZoneBC_t", 1, "BC_t", cgbc_wall, "end");
			   //cg_user_data_write("FamBC_UserId", CGNS_ENUMV(Character), &index_dim, "FamBC_UserId_Outflow");
			   cg_descriptor_write("FamBC_UserId", "FamBC_UserId_Outflow");

		   }   	   
		   
		   std::cout << "cgbc_wall " << cgbc_wall << std::endl;
	   }
	   
	   int index_conn = 0;
	   int* transform = new int[3];
	   transform[0] = 1;
	   transform[1] = 2;
	   transform[2] = 3;
	   int icount = interface_nodes.size();
	   
	   //cg_1to1_write(index_file,index_base,index_zone_t,"Interior","blk-1",interface_tetra,interface_prism,transform,&index_conn);
	   ier = cg_conn_write(index_file,index_base,index_zone_t,"Interface",Vertex,Abutting1to1,
	                 PointList,icount,interface_tetra,"blk-1",Unstructured,
	                 PointListDonor,Integer,icount,interface_prism,&index_conn);
	   
	   for(int q=0;q<icount;q++)
	   {
		   std::cout << "interface = " << interface_tetra[q] << " " << interface_prism[q] << std::endl;
	   }
	   /* */
	   //std::cout << "ier = " << ier << std::endl;
	   cg_close(index_file);
	   MPI_Finalize();  
	   std::cout << "\nSuccessfully wrote grid to file mesh2.cgns\n" << std::endl;
	   /**/
	   
	   return 0;
	   
}
