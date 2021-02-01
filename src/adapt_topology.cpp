#include "adapt_topology.h"
#include "adapt_output.h"

Mesh_Topology::Mesh_Topology(Partition* Pa, MPI_Comm comm)
{
    int nlocElem, start, end, offset, nloc, np, loc_vid, size, rank, lid;
    int vf0, vf1, vf2, vf3, vf4, vf5, vf6, vf7, fid;
    double wi, ds0, ds1 ,ds2, ds3, ds4, ds5, u_po,orient0,orient1,orient2,orient3,orient4,orient5,L0,L1,L2,L3,L4,L5;
    
    //ifn = ifn_in;
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
    i_part_map* ifn_part_map = Pa->getIFNpartmap();
    i_part_map* ife_part_map = Pa->getIFEpartmap();
    i_part_map* ief_part_map = Pa->getIEFpartmap();
    i_part_map* iee_part_map = Pa->getIEEpartmap();
    i_part_map* if_ref_part_map = Pa->getIFREFpartmap();
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
    //int fint = bnd_map[0];
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
                int gvid = ifn_part_map->i_map[faceid][r];
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
            
//          L0 = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)
//                   +(Vface->y-Vijk->y)*(Vface->y-Vijk->y)
//                   +(Vface->z-Vijk->z)*(Vface->z-Vijk->z));
            
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
            
            face2ref[faceid]        = ref;
            ref2face[ref].push_back(faceid);
            int ref = if_ref_part_map->i_map[faceid][0];
            if(ref!=2)
            {
                Bface2Element[faceid]   = gEl;
                Bface2Normal[faceid]    = n0;
                Bface2LocID[faceid]     = s;
            }
            
            if(ref!=2)
            {
                face2ref[faceid]        = ref;
                ref2face[ref].push_back(faceid);
            }
//            //double orientaft = DotVec3D(r0,n0);
//            if(faceid<fint) // identify the internal face;
//            {
//                ref = 0;
//                face2ref[faceid] = ref;
//                ref2face[ref].push_back(faceid);
//            }
//            else // identify the boundary interface and based on bnd_map, determine the reference value.
//            {
//                ref = FindBoundaryID(bnd_map,nBnd,faceid)+1;
//                face2ref[faceid]        = ref;
//                ref2face[ref].push_back(faceid);
//                //Bface2Element[faceid]   = gEl;
//                //Bface2Normal[faceid]    = n0;
//                //Bface2LocID[faceid]     = s;
//            }
            
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
                   //int gvid = ifn->getVal(fid,s);
                   int gvid = ifn_part_map->i_map[fid][s];

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
    
    //DetermineBoundaryLayerElements(Pa,ife_in,15,3,c);
    

//    delete[] Pijk;
//    delete[] Po;
//
//    delete v0;
//    delete v1;
//
//    delete fc0;
//    delete fc1;
//    delete fc2;
//    delete fc3;
//    delete fc4;
//    delete fc5;
    
}


/*
void Mesh_Topology::DetermineBoundaryLayerElements(Partition* Pa, Array<int>* ife_in, int nLayer, int bID, MPI_Comm comm)
{
        int size;
    int rank;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    int fid_start,gEl,t,gElnew,min_index,lid,fid_new,lEl;
    double newf_id;
    std::vector<double> dp(6);
    std::map<int,std::vector<int> > BLelements;
    std::vector<int> gElvec;
    std::map<int,std::vector<int> > gF2gE               = Pa->getglobFace2GlobalElements();
    i_part_map* gE2gF                                   = Pa->getIEFpartmap();
    std::map<int,int> gE2lE                             = Pa->getGlobalElement2LocalElement();
    std::vector<std::vector<int> > loc_elem2verts_loc   = Pa->getLocalElem2LocalVert();
    std::vector<Vert> LVerts                            = Pa->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV               = Pa->getGlobElem2LocVerts();
    i_part_map* ien_part_map = Pa->getIENpartmap();
    //std::map<int,std::vector<int> >::iterator it;
    Array<int>* gPart = Pa->getGlobalPartition();
    int tel = 0;
    int gElvec0,gElvec1;
    std::vector<int> ElLayer;
    int tel2 = 0;
    double min_val;
    i_part_map* ife_part_map = Pa->getIFEpartmap();
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_elem;
    int ui=0;
    mesh_topo_bl = new Mesh_Topology_BL;
    for(int i=0;i<ref2face[bID].size();i++)
    {
        fid_start = ref2face[bID][i];
    
        Vec3D* nb = Bface2Normal[fid_start];
        //NegateVec3D(nb);
        gEl = Bface2Element[fid_start];
        ElLayer.push_back(gEl);
        for(int k=0;k<nLayer;k++)
        {
            t = 0;
            for(int j=0;j<6;j++)
            {
                Vec3D* nt = normals[gEl][j];
                dp[j] = DotVec3D(nb,nt);
            }
            
            min_index       = std::min_element(dp.begin(),dp.end())-dp.begin();
            min_val  = *std::min_element(dp.begin(),dp.end());
            
            fid_new         = gE2gF->i_map[gEl][min_index];
            nb              = normals[gEl][min_index];
            NegateVec3D(nb);

            if(gE2gF->i_inv_map[fid_new].size()==2)
            {
                gElvec0    = gE2gF->i_inv_map[fid_new][0];
                gElvec1    = gE2gF->i_inv_map[fid_new][1];
                
                if(gElvec0==gEl)
                {
                    gElnew = gElvec1;
                }
                if(gElvec1==gEl)
                {
                    gElnew = gElvec0;
                }

                ElLayer.push_back(gElnew);
                gEl=gElnew;
            }
            else
            {

//                gElvec0 = ife_in->getVal(fid_new,0);
//                gElvec1 = ife_in->getVal(fid_new,1);
                gElvec0 = ife_part_map[fid_new][0];
                gElvec1 = ife_part_map[fid_new][1];
                int rank_id0 = gPart->getVal(gElvec0,0);
                int rank_id1 = gPart->getVal(gElvec1,0);
                if(gElvec0==gEl)
                {
                    rank2req_elem[rank_id1].push_back(gElvec1);
                }
                if(gElvec1==gEl)
                {
                    rank2req_elem[rank_id0].push_back(gElvec0);
                }
                ui++;
            }
        }
        mesh_topo_bl->BLlayers[fid_start] = ElLayer;
        ElLayer.clear();
    }
    
    
    ScheduleObj* iee_schedule = DoScheduling(rank2req_elem,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;

    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_elem.begin(); it != rank2req_elem.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                //for(int n=0;n<n_req;n++)
		//{
		//	std::cout << rank << " requested " << it->second[n] << std::endl;
		//}
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+20*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 9876*2*7654+dest*40, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+20*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*40, comm, MPI_STATUS_IGNORE);
            
            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }


    std::cout << rank << " " << ui << " " << rank2req_elem.size() << std::endl;
    //std::cout << rank << " " << ui << " " << rank2req_elem.size() << std::endl;
    
    
    //std::map<int,std::vector<int> >::iterator iter;

    //for(iter=rank2req_elem.begin();iter!=rank2req_elem.end();iter++)
    //{
	//std::cout << "Rank" << rank << " requests from " << iter->first << " for " << iter->second.size() << " IDs"<< std::endl;  
    //}
    //for(iter=reqstd_E_IDs_per_rank.begin();iter!=reqstd_E_IDs_per_rank.end();iter++)
    //{
//	std::cout << " Rank "<< rank << " received request from " << iter->first << " for " << iter->second.size() << std::endl;
 //   }
    
    
    std::map<int,int > recv_back_Niee;
    std::map<int,int* > recv_back_el_ids;
    std::map<int,int* > recv_back_vrt_ids;
    std::map<int,double* > recv_back_vcrds;
    
    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_E_IDs_per_rank.begin(); it != reqstd_E_IDs_per_rank.end(); it++)
            {
                int ne_send_el  = it->second.size();
                int ne_send_vrt = it->second.size()*8;
                int ne_send_vrt_coord = it->second.size()*8*3;
                int* iee_send_v = new int[ne_send_vrt];
                double* iee_send_v_crd = new double[ne_send_vrt_coord];
                for(int u=0;u<it->second.size();u++)
                {
                    for(int w=0;w<8;w++)
                    {
                        int lV = gE2lV[it->second[u]][w];
		        iee_send_v[u*8+w] = ien_part_map->i_map[it->second[u]][w];
                        iee_send_v_crd[u*8+w+0]=LVerts[lV].x;
                        iee_send_v_crd[u*8+w+1]=LVerts[lV].y;
                        iee_send_v_crd[u*8+w+2]=LVerts[lV].z;
			std::cout << "rank = " << rank << " local ID " << lV << " " << LVerts[lV].x << " " << LVerts[lV].y << " " << LVerts[lV].z << std::endl;
                    }
                }

                int dest = it->first;
                MPI_Send(&ne_send_el, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second[0], ne_send_el, MPI_INT, dest, 9876*7777+dest*888, comm);
                
                MPI_Send(&iee_send_v[0], ne_send_vrt, MPI_INT, dest, 9876*8888+dest*888, comm);
                
                MPI_Send(&iee_send_v_crd[0], ne_send_vrt_coord, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                //delete[] iee_send_e;
                //delete[] iee_send_v;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
	    int n_recv_back;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            double* recv_back_vcrds_arr   = new double[n_recv_back*8*3];
            int*    recv_back_vids_arr  = new int[n_recv_back*8];
            int*    recv_back_ids_arr   = new int[n_recv_back];
            MPI_Recv(&recv_back_ids_arr[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_vids_arr[0], n_recv_back*8, MPI_INT, q, 9876*8888+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_vcrds_arr[0], n_recv_back*8*3, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Niee[q]     = n_recv_back;
            recv_back_el_ids[q]   = recv_back_ids_arr;
            recv_back_vrt_ids[q]   = recv_back_vids_arr;
            recv_back_vcrds[q]      = recv_back_vcrds_arr;
         }
    }

    
    std::map<int,int >::iterator iter;
    int ntotal=0;
    int lid_vert = 0;
    std::set<int> unique_vert_ex;
    for(iter=recv_back_Niee.begin();iter!=recv_back_Niee.end();iter++)
    {
        int L = iter->second;
        for(int s=0;s<L;s++)
        {
            int el_id = recv_back_el_ids[iter->first][s];
            mesh_topo_bl->exteriorElIDs.push_back(el_id);
            
            for(int p=0;p<8;p++)
            {
 //		std::cout << "rank " << rank << " " << recv_back_vrt_ids[iter->first][s*8+p] << std::endl;
                mesh_topo_bl->exteriorVertIDs[el_id].push_back(recv_back_vrt_ids[iter->first][s*8+p]);
		
		if(unique_vert_ex.find(recv_back_vrt_ids[iter->first][s*8+p])==unique_vert_ex.end())
		{
			std::vector<double> vert(3);
			vert[0]=recv_back_vcrds[iter->first][s*8+p+0];
			vert[1]=recv_back_vcrds[iter->first][s*8+p+1];
			vert[2]=recv_back_vcrds[iter->first][s*8+p+2];
			unique_vert_ex.insert(recv_back_vrt_ids[iter->first][s*8+p]);
			
			mesh_topo_bl->verts_g2l_ex[recv_back_vrt_ids[iter->first][s*8+p]]=lid_vert;
			mesh_topo_bl->local_ex_verts[lid_vert]=vert;
			lid_vert++;
		}

	    }
        }
    }
}
*/




std::vector<double> Mesh_Topology::ReduceUToVertices(
                        Domain* pDom,
                        std::map<int,double> Uelem)
{
    std::vector<double> Uv;
    
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    
    int im = 0;
    std::map<int,std::vector<int> >::iterator itm;
    for(itm=v2e.begin();itm!=v2e.end();itm++)
    {
        double sum = 0.0;
        double avg = 0.0;
        
        for(int q=0;q<itm->second.size();q++)
        {
            int gEl = itm->second[q];
            sum = sum + Uelem[gEl];
        }
        
        avg = sum/itm->second.size();
        //vval->setVal(im,0,avg);
        Uv.push_back(avg);
        im++;
    }
    
    return Uv;
    
}

std::map<int,double> Mesh_Topology::ReduceFieldToVertices(
                        Domain* pDom,
                        std::map<int,double> Uelem)
{
    std::vector<double> Uv;
    std::map<int,double> Uvm;
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        double sum = 0.0;
        double avg = 0.0;
        int gv = itm->first;
        
        for(int q=0;q<itm->second.size();q++)
        {
            int gEl = itm->second[q];
            sum = sum + Uelem[gEl];
        }
        
        avg = sum/itm->second.size();
        Uv.push_back(avg);
        Uvm[gv]=avg;
        //std::cout <<"hhh "<< gv << " " << gv2  << " " << avg<<std::endl;

        im++;
    }

    return Uvm;

}

std::map<int,Array<double>* > Mesh_Topology::ReduceMetricToVertices(Domain* pDom,
                                                     std::map<int,Array<double>* > Telem)
{
    std::map<int,Array<double>*> Tmet_v;
    
    std::map<int,std::vector<int> > v2e = pDom->vert2elem;
    std::vector<int> glob_part_verts = pDom->glob_part_verts;
    std::map<int,int> lv2gpv = pDom->lv2gpv;
    int im = 0;
    int tel=0;
    std::map<int,std::vector<int> >::iterator itm;
    double sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
    double sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
    double sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
    for(itm=pDom->vert2elem.begin();itm!=pDom->vert2elem.end();itm++)
    {
        sum00 = 0.0,sum01 = 0.0,sum02 = 0.0;
        sum10 = 0.0,sum11 = 0.0,sum12 = 0.0;
        sum20 = 0.0,sum21 = 0.0,sum22 = 0.0;
        double avg = 0.0;
        int gv = itm->first;
        Array<double>* Tv = new Array<double>(3,3);
        for(int q=0;q<itm->second.size();q++)
        {
            int gEl = itm->second[q];
            sum00 = sum00 + Telem[gEl]->getVal(0,0);
            sum01 = sum01 + Telem[gEl]->getVal(0,1);
            sum02 = sum02 + Telem[gEl]->getVal(0,2);
            
            sum10 = sum10 + Telem[gEl]->getVal(1,0);
            sum11 = sum11 + Telem[gEl]->getVal(1,1);
            sum12 = sum12 + Telem[gEl]->getVal(1,2);
            
            sum20 = sum20 + Telem[gEl]->getVal(2,0);
            sum21 = sum21 + Telem[gEl]->getVal(2,1);
            sum22 = sum22 + Telem[gEl]->getVal(2,2);
        }
        
        Tv->setVal(0,0,sum00/itm->second.size());
        Tv->setVal(0,1,sum01/itm->second.size());
        Tv->setVal(0,2,sum02/itm->second.size());
        
        Tv->setVal(1,0,sum10/itm->second.size());
        Tv->setVal(1,1,sum11/itm->second.size());
        Tv->setVal(1,2,sum12/itm->second.size());
        
        Tv->setVal(2,0,sum20/itm->second.size());
        Tv->setVal(2,1,sum21/itm->second.size());
        Tv->setVal(2,2,sum22/itm->second.size());
        
        Tmet_v[gv]=Tv;

        im++;
    }
    
    return Tmet_v;
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
