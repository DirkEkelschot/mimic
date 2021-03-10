#include "adapt_bltopology.h"


BLShellInfo* FindOuterShellBoundaryLayerMesh(int wall_id, int nLayer,
                                             Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g,
                                             ParallelState* xcn_pstate, ParallelState* ien_pstate, std::map<int,std::vector<int> > bnd_face_map, std::map<int,int> vert_ref_map, MPI_Comm comm)
{
    BLShellInfo* BLinfo = new BLShellInfo;
    BLinfo->ShellRef = new Array<int>(xcn_g->getNrow(),1);
    int te1=0;
    int te2=0;
    int te3=0;
    for(int i=0;i<xcn_g->getNrow();i++)
    {
        if(vert_ref_map.find(i)!=vert_ref_map.end())
        {
            BLinfo->ShellRef->setVal(i,0,100+vert_ref_map[i]);
        }
        else
        {
            BLinfo->ShellRef->setVal(i,0,-3);
        }
    }
    
    std::vector<int> outer_shell_faces;
    std::vector<std::vector<int> > outer_shell_elements;
    std::map<int,std::set<int> > outer_shell_Faces2Nodes;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<Vec3D*> dpvec(6);
    //std::cout << "Determining outer shell of BL mesh..." << std::endl;
    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    double* Pijk = new double[8*3];
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    int bvid,opposite_bvid;
    int bvid_b;
    int fv1_b;
    int fv2_b;
    int fv3_b;
    int glob_el_id = 0;
    for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
    {
        std::vector<int> layer;

        int bfaceid = bnd_face_map[wall_id][bf];
        int faceid  = bfaceid;
        int elid0   = ife_g->getVal(faceid,0);
        int elid1   = ife_g->getVal(faceid,1);

        if(elid0<ien_g->getNrow())
        {
            elid_cur = elid0;
        }
        else
        {
            elid_cur = elid1;
        }
        
        layer.push_back(elid_cur);
        std::set<int> local_faces;

        for(int k=0;k<8;k++)
        {
           loc_vid     = ien_g->getVal(elid_cur,k);
           Pijk_id[k]  = loc_vid;
           Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
           Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
           Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
        }
        
        Vert* Vijk = ComputeCenterCoord(Pijk, 8);
        
        Vert* Vface  = new Vert;
        std::vector<Vert*> face;
        std::vector<Vert*> face_turned(4);
        std::vector<Vert*> face_turned2(4);
        std::set<int> conn_bvid;
        for(int r=0;r<4;r++)
        {
            int vid  = ifn_g->getVal(faceid,r);
            if(r==0)
            {
                bvid = vid;
            }
            Vert* V  = new Vert;
            V->x     = xcn_g->getVal(vid,0);
            V->y     = xcn_g->getVal(vid,1);
            V->z     = xcn_g->getVal(vid,2);
            Vface->x = Vface->x+V->x;
            Vface->y = Vface->y+V->y;
            Vface->z = Vface->z+V->z;
            face.push_back(V);
        }
        
        conn_bvid.insert(ifn_g->getVal(faceid,1));
        conn_bvid.insert(ifn_g->getVal(faceid,3));
        bvid_b = bvid;
        fv1_b = ifn_g->getVal(faceid,1);
        fv2_b = ifn_g->getVal(faceid,2);
        fv3_b = ifn_g->getVal(faceid,3);
        
        
        Vface->x = Vface->x/4.0;
        Vface->y = Vface->y/4.0;
        Vface->z = Vface->z/4.0;
                        
        Vec3D* r0 = new Vec3D;
        r0->c0 = (Vface->x-Vijk->x);
        r0->c1 = (Vface->y-Vijk->y);
        r0->c2 = (Vface->z-Vijk->z);
        Vec3D* v0 = new Vec3D;
        v0->c0 = face[1]->x-face[0]->x;
        v0->c1 = face[1]->y-face[0]->y;
        v0->c2 = face[1]->z-face[0]->z;
        Vec3D* v1 = new Vec3D;
        v1->c0 = face[3]->x-face[0]->x;
        v1->c1 = face[3]->y-face[0]->y;
        v1->c2 = face[3]->z-face[0]->z;
        
        Vec3D* nbf     = ComputeSurfaceNormal(v0,v1);
        double orient0 = DotVec3D(r0,nbf);
        
        if(orient0<0.0)
        {
            NegateVec3D(nbf);
        }
        face.clear();
        
        std::vector<std::map<int,std::set<int> > > layer_locN2NEl;
        for(int c=0;c<nLayer;c++)
        {
            for(int k=0;k<8;k++)
            {
               loc_vid     = ien_g->getVal(elid_cur,k);
               Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
               Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
               Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
                
               if(BLinfo->verts_set.find(loc_vid)==BLinfo->verts_set.end())
               {
                   BLinfo->verts_set.insert(loc_vid);
               }
            }
        
            Vert* Vijk = ComputeCenterCoord(Pijk, 8);
            std::vector<std::vector<int> > face_id_stored(6);
            std::vector<std::vector<Vert*> > face_stored(6);
            std::map<int,std::set<int> > local_node2node_element;
            std::vector<std::map<int,std::set<int> > > local_node2node_face(6);
            std::vector<std::map<int,int> > local_node2opponode_face(6);
            for(int k=0;k<6;k++)
            {
                int fid = ief_g->getVal(elid_cur,k);
                Vert* Vface2  = new Vert;
                
                std::vector<int> faceVert_IDs(4);
                std::vector<Vert*> face2;
                for(int r=0;r<4;r++)
                {
                    int vid  = ifn_g->getVal(fid,r);
                    
                    Vert* V  = new Vert;
                    V->x     = xcn_g->getVal(vid,0);
                    V->y     = xcn_g->getVal(vid,1);
                    V->z     = xcn_g->getVal(vid,2);
                    Vface2->x = Vface2->x+V->x;
                    Vface2->y = Vface2->y+V->y;
                    Vface2->z = Vface2->z+V->z;
                    face2.push_back(V);
                    
                    faceVert_IDs[r] = vid;
                }
                
                local_node2node_element[ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,1));
                local_node2node_element[ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,3));
                local_node2node_element[ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,0));
                local_node2node_element[ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,2));
                local_node2node_element[ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,1));
                local_node2node_element[ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,3));
                local_node2node_element[ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,2));
                local_node2node_element[ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,0));
                
                local_node2node_face[k][ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,1));
                local_node2node_face[k][ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,3));
                local_node2node_face[k][ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,0));
                local_node2node_face[k][ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,2));
                local_node2node_face[k][ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,1));
                local_node2node_face[k][ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,3));
                local_node2node_face[k][ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,2));
                local_node2node_face[k][ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,0));
                
                local_node2opponode_face[k][ifn_g->getVal(fid,0)]=ifn_g->getVal(fid,2);
                local_node2opponode_face[k][ifn_g->getVal(fid,1)]=ifn_g->getVal(fid,3);
                local_node2opponode_face[k][ifn_g->getVal(fid,2)]=ifn_g->getVal(fid,0);
                local_node2opponode_face[k][ifn_g->getVal(fid,3)]=ifn_g->getVal(fid,1);
                //layer_locN2NEl.push_back(local_node2node_element);
                Vface2->x = Vface2->x/4.0;
                Vface2->y = Vface2->y/4.0;
                Vface2->z = Vface2->z/4.0;
                                
                Vec3D* r00 = new Vec3D;
                r00->c0 = (Vface2->x-Vijk->x);
                r00->c1 = (Vface2->y-Vijk->y);
                r00->c2 = (Vface2->z-Vijk->z);
                Vec3D* v00 = new Vec3D;
                v00->c0 = face2[1]->x-face2[0]->x;
                v00->c1 = face2[1]->y-face2[0]->y;
                v00->c2 = face2[1]->z-face2[0]->z;
                Vec3D* v11 = new Vec3D;
                v11->c0 = face2[3]->x-face2[0]->x;
                v11->c1 = face2[3]->y-face2[0]->y;
                v11->c2 = face2[3]->z-face2[0]->z;
                
                Vec3D* n00        = ComputeSurfaceNormal(v00,v11);
                double orient00   = DotVec3D(r00,n00);
                
                if(orient00<0.0)
                {
                    NegateVec3D(n00);
                }
                
                dp[k]               =   DotVec3D(nbf,n00);
                dpvec[k]            =   n00;
                face2.clear();
            }
            std::set<int>::iterator its;
            for(its=local_node2node_element[bvid].begin();its!=local_node2node_element[bvid].end();its++)
            {
                if(conn_bvid.find(*its)==conn_bvid.end())
                {
                    opposite_bvid = *its;
                }
            }
        
            int min_index  = std::min_element(dp.begin(),dp.end())-dp.begin();
            //double min_val = *std::min_element(dp.begin(),dp.end());
            int fid_new    = ief_g->getVal(elid_cur,min_index);
            nbf            = dpvec[min_index];
            
            std::map<int,std::set<int> > node2node_face = local_node2node_face[min_index];
            std::set<int>::iterator itu;
            std::vector<int>opposite_tri(3);
            opposite_tri[0] = opposite_bvid;
            int l = 1;
            for(itu=node2node_face[opposite_bvid].begin();itu!=node2node_face[opposite_bvid].end();itu++)
            {
                opposite_tri[l] = *itu;
                l++;
            }

            NegateVec3D(nbf);

            int gEl0=ife_g->getVal(fid_new,0);
            int gEl1=ife_g->getVal(fid_new,1);

            if(gEl0==elid_cur)
            {
                elid_next = gEl1;
            }
            else if(gEl1==elid_cur)
            {
                elid_next = gEl0;
            }
            if(c<nLayer-1)
            {
                layer.push_back(elid_next);
            }
            
            //int changed = ChkHexorient(Pijk,Pijk_id);
            
            if(c==nLayer-1)
            {
                outer_shell_faces.push_back(fid_new);
                outer_shell_Faces2Nodes[fid_new].insert(ifn_g->getVal(fid_new,0));
                outer_shell_Faces2Nodes[fid_new].insert(ifn_g->getVal(fid_new,1));
                outer_shell_Faces2Nodes[fid_new].insert(ifn_g->getVal(fid_new,2));
                outer_shell_Faces2Nodes[fid_new].insert(ifn_g->getVal(fid_new,3));
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,0),0,-1);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,1),0,-1);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,2),0,-1);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,3),0,-1);
                std::set<int> ShellTri0;
                ShellTri0.insert(ifn_g->getVal(fid_new,0));
                ShellTri0.insert(ifn_g->getVal(fid_new,1));
                ShellTri0.insert(ifn_g->getVal(fid_new,3));
                BLinfo->ShellTri2FaceID[ShellTri0] = fid_new;
                std::set<int> ShellTri1;
                ShellTri1.insert(ifn_g->getVal(fid_new,1));
                ShellTri1.insert(ifn_g->getVal(fid_new,2));
                ShellTri1.insert(ifn_g->getVal(fid_new,3));
                BLinfo->ShellTri2FaceID[ShellTri1] = fid_new;
                std::set<int> ShellTri2;
                ShellTri2.insert(ifn_g->getVal(fid_new,0));
                ShellTri2.insert(ifn_g->getVal(fid_new,1));
                ShellTri2.insert(ifn_g->getVal(fid_new,2));
                BLinfo->ShellTri2FaceID[ShellTri2] = fid_new;
                std::set<int> ShellTri3;
                ShellTri3.insert(ifn_g->getVal(fid_new,2));
                ShellTri3.insert(ifn_g->getVal(fid_new,3));
                ShellTri3.insert(ifn_g->getVal(fid_new,0));
                BLinfo->ShellTri2FaceID[ShellTri3] = fid_new;
                
                //            bvid > opposite_bvid
                //            fv1 > opposite_tri[1];
                //            fv3 > opposite_tri[2];
                //            fv2 > local_node2opponode_face[min_index][opposite_bvid]
                std::map<int,int> opposite_verts;
                
//                opposite_verts[bvid_b] = opposite_bvid;
//                opposite_verts[fv1_b]  = opposite_tri[1];
//                opposite_verts[fv3_b]  = opposite_tri[2];
//                opposite_verts[fv2_b]  = local_node2opponode_face[min_index][opposite_bvid];
                
                opposite_verts[opposite_bvid]    = bvid_b;
                opposite_verts[opposite_tri[1]]  = fv1_b;
                opposite_verts[opposite_tri[2]]  = fv3_b;
                opposite_verts[local_node2opponode_face[min_index][opposite_bvid]]  = fv2_b;
                
                BLinfo->ShellFace2ShellVert2OppositeBoundaryVerts[fid_new] = opposite_verts;
                BLinfo->ShellFace2BFace[fid_new] = bfaceid;
                BLinfo->BFace2ShellFace[bfaceid] = fid_new;
                std::vector<int> fe(2);
                fe[0] = elid_cur;
                fe[1] = elid_next;
                outer_shell_elements.push_back(fe);
                opposite_verts.clear();
            }
            
            conn_bvid.clear();
            bvid = opposite_tri[0];
            conn_bvid.insert(opposite_tri[1]);
            conn_bvid.insert(opposite_tri[2]);
            
            elid_cur = elid_next;
            
            BLinfo->BLlayers[bfaceid]=layer;
            //BLinfo->bFace2_locN2NEl[bfaceid]=layer_locN2NEl;
        }
        
        
    }
    delete[] Pijk;
    delete[] Pijk_id;
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing extracting outer shell BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    for(itt=BLinfo->BLlayers.begin();itt!=BLinfo->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            BLinfo->elements_set.insert(itt->second[q]);
        }
    }
//    OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_");
    return BLinfo;
}



Mesh_Topology_BL* ExtractBoundaryLayerMeshFromShell(std::vector<std::vector<int> > u_tris, BLShellInfo* BLshell, int wall_id, int nLayer, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g, ParallelState* xcn_pstate, ParallelState* ien_pstate, std::map<int,std::vector<int> > bnd_face_map, std::map<std::set<int>,int> tria_ref_map, std::map<std::set<int>,int> quad_ref_map,  MPI_Comm comm)
{
    Mesh_Topology_BL* mesh_topology_bl = new Mesh_Topology_BL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<Vec3D*> dpvec(6);

    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    double* Pijk = new double[8*3];
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    //std::vector<Vert*> face_c;
    std::map<int,std::vector<Vert*> > prisms;
    
    std::set<int> tria0;
    std::set<int> tria1;
    std::set<int> quad0;
    std::set<int> quad1;
    std::set<int> quad2;
    
    int or0 = 0;

    int bvid=-1,opposite_bvid=-1;
    std::vector<int> opposite_tri(3);
    std::vector<int> opposite_tri_n(3);
    std::vector<int> opposite_tri1(3);
    std::vector<int> opposite_tri1_n(3);
    std::vector<int> prism0(6);
    std::vector<int> prism1(6);

    std::vector<int> prismStored0(6);
    std::vector<int> prismStored1(6);
    mesh_topology_bl->Nprisms = 0;
    int glob_el_id = 0;
    std::map<int,int> bface2shellface = BLshell->BFace2ShellFace;
    std::map<std::set<int>,int> shelltri2shellface = BLshell->ShellTri2FaceID;
    std::map<int,std::vector<int> > shellfaceID2triID =  BLshell->ShellFaceID2TriID;
    int opposite_p00,opposite_p01,opposite_p02,opposite_p10,opposite_p11,opposite_p12;
    int cnt_turn = 0;
    int fc0 = 0;
    int fc1 = 0;
    int fwrong = 0;
    int fright = 0;
    for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
    {
        
        int bvid=-1,obvid_i=-1,opposite_bvid=-1;
        std::vector<int> layer;
        int bfaceid      = bnd_face_map[wall_id][bf];
        int shell_faceid = BLshell->BFace2ShellFace[bfaceid];
        int triID0 = shellfaceID2triID[shell_faceid][0];
        int triID1 = shellfaceID2triID[shell_faceid][1];
        std::vector<int> tri_shell_0 = u_tris[triID0];
        std::vector<int> tri_shell_1 = u_tris[triID1];
        
//        std::cout << "shell tri0 " << tri_shell_0[0] << " " <<  tri_shell_0[1] << " " <<  tri_shell_0[2] << std::endl;
//        std::cout << "shell tri1 " << tri_shell_1[0] << " " <<  tri_shell_1[1] << " " <<  tri_shell_1[2] << std::endl;
        
        std::vector<int> tri_bound_0(3);
        tri_bound_0[0] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_0[0]];
        tri_bound_0[1] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_0[1]];
        tri_bound_0[2] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_0[2]];
        std::vector<int> tri_bound_1(3);
        tri_bound_1[0] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_1[0]];
        tri_bound_1[1] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_1[1]];
        tri_bound_1[2] = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid][tri_shell_1[2]];
        
//        std::cout << "bound tri0 " << tri_bound_0[0] << " " <<  tri_bound_0[1] << " " <<  tri_bound_0[2] << std::endl;
//        std::cout << "bound tri1 " << tri_bound_1[0] << " " <<  tri_bound_1[1] << " " <<  tri_bound_1[2] << std::endl;
        
        //std::cout << "BFF "<< bf << " " << bfaceid << " " << wall_id << std::endl; 
        
        
        std::set<int> conn_bvid;
        std::vector<int> tri_0n(3);
        std::vector<int> tri_1n(3);
        std::vector<int> tri_0n_tmp(3);
        std::vector<int> tri_1n_tmp(3);
        for(int v=0;v<3;v++)
        {
            int fid0 = tri_bound_0[v];
            if(std::find(tri_bound_1.begin(), tri_bound_1.end(), fid0) != tri_bound_1.end())
            {
                conn_bvid.insert(fid0);
            }
        }
        for(int v=0;v<3;v++)
        {
            int fid0 = tri_bound_0[v];
            int fid1 = tri_bound_1[v];
            if(conn_bvid.find(fid0)==conn_bvid.end())
            {
                bvid = fid0;
            }
            if(conn_bvid.find(fid1)==conn_bvid.end())
            {
                obvid_i = fid1;
            }
            
        }
        
        tri_0n[0] = bvid;
        tri_0n[1] = *next(conn_bvid.begin(),0);
        tri_0n[2] = *next(conn_bvid.begin(),1);
        
        tri_1n[0] = obvid_i;
        tri_1n[1] = *next(conn_bvid.begin(),1);
        tri_1n[2] = *next(conn_bvid.begin(),0);
        
        std::vector<int> facenew(4);
        facenew[0] = bvid;
        facenew[1] = tri_0n[1];
        facenew[2] = obvid_i;
        facenew[3] = tri_0n[2];
        
        int faceid  = bfaceid;
        int elid0   = ife_g->getVal(faceid,0);
        int elid1   = ife_g->getVal(faceid,1);

        if(elid0<ien_g->getNrow())
        {
            elid_cur = elid0;
        }
        else
        {
            elid_cur = elid1;
        }
        layer.push_back(elid_cur);
        
        std::set<int> local_faces;
        
        //std::cout << "Element -> ";
        for(int k=0;k<8;k++)
        {
           loc_vid     = ien_g->getVal(elid_cur,k);
           // std::cout << loc_vid << " ";
           Pijk_id[k]  = loc_vid;
           Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
           Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
           Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
        }
        //std::cout << std::endl;
        //int changed = ChkHexorient(Pijk,Pijk_id);
//
        Vert* Vijk = ComputeCenterCoord(Pijk, 8);
        
        
        Vert* Vface  = new Vert;
        std::vector<Vert*> face;
        std::vector<Vert*> face_turned(4);
        std::vector<Vert*> face_turned2(4);
        for(int r=0;r<4;r++)
        {
            int vid  = ifn_g->getVal(faceid,r);
            
//            if(r==0)
//            {
//                bvid = vid;
//            }
            
            Vert* V  = new Vert;
            V->x     = xcn_g->getVal(vid,0);
            V->y     = xcn_g->getVal(vid,1);
            V->z     = xcn_g->getVal(vid,2);
            Vface->x = Vface->x+V->x;
            Vface->y = Vface->y+V->y;
            Vface->z = Vface->z+V->z;
            
            face.push_back(V);
        }
        std::vector<int> tri0(3);
        std::vector<int> tri1(3);
        
        Vec3D* v_t0 = new Vec3D;
        v_t0->c0 = xcn_g->getVal(tri_0n[1],0)-xcn_g->getVal(tri_0n[0],0);
        v_t0->c1 = xcn_g->getVal(tri_0n[1],1)-xcn_g->getVal(tri_0n[0],1);
        v_t0->c2 = xcn_g->getVal(tri_0n[1],2)-xcn_g->getVal(tri_0n[0],2);
        Vec3D* v_t1 = new Vec3D;
        v_t1->c0 = xcn_g->getVal(tri_0n[2],0)-xcn_g->getVal(tri_0n[0],0);
        v_t1->c1 = xcn_g->getVal(tri_0n[2],1)-xcn_g->getVal(tri_0n[0],1);
        v_t1->c2 = xcn_g->getVal(tri_0n[2],2)-xcn_g->getVal(tri_0n[0],2);
        Vec3D* n_t0        = ComputeSurfaceNormal(v_t0,v_t1);
        
        Vec3D* v_t10 = new Vec3D;
        v_t10->c0 = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10->c1 = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10->c2 = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);
        Vec3D* v_t11 = new Vec3D;
        v_t11->c0 = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11->c1 = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11->c2 = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        Vec3D* n_t10        = ComputeSurfaceNormal(v_t10,v_t11);
        
        tri0[0] = tri_0n[0];
        tri0[1] = tri_0n[1];
        tri0[2] = tri_0n[2];
        mesh_topology_bl->BndFaces.push_back(tri0);

        tri1[0] = tri_1n[0];
        tri1[1] = tri_1n[1];
        tri1[2] = tri_1n[2];
        mesh_topology_bl->BndFaces.push_back(tri1);
        
        Vface->x = Vface->x/4.0;
        Vface->y = Vface->y/4.0;
        Vface->z = Vface->z/4.0;
                        
        Vec3D* r0 = new Vec3D;
        double r0L = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)
                          +(Vface->y-Vijk->y)*(Vface->y-Vijk->y)
                          +(Vface->z-Vijk->z)*(Vface->z-Vijk->z));
        r0->c0 = (Vface->x-Vijk->x)/r0L;
        r0->c1 = (Vface->y-Vijk->y)/r0L;
        r0->c2 = (Vface->z-Vijk->z)/r0L;
        Vec3D* v0 = new Vec3D;
        v0->c0 = face[1]->x-face[0]->x;
        v0->c1 = face[1]->y-face[0]->y;
        v0->c2 = face[1]->z-face[0]->z;
        Vec3D* v1 = new Vec3D;
        v1->c0 = face[3]->x-face[0]->x;
        v1->c1 = face[3]->y-face[0]->y;
        v1->c2 = face[3]->z-face[0]->z;
        
        Vec3D* nbf     = ComputeSurfaceNormal(v0,v1);
        double orient0 = DotVec3D(r0,nbf);
        
        
        double orient_t0 = DotVec3D(r0,n_t0);
        double orient_t1 = DotVec3D(r0,n_t10);
        
        if(orient_t0<0)
        {
            fc0++;
        }
        if(orient_t1<0)
        {
            fc1++;
        }
        //std::cout << "First Check " << orient_t0 << " " << orient_t1 << std::endl;
        
        if(orient_t0<0.0)
        {
            tri_0n_tmp[0] = tri_0n[0];
            tri_0n_tmp[1] = tri_0n[2];
            tri_0n_tmp[2] = tri_0n[1];
            
            tri_0n[0] = tri_0n_tmp[0];
            tri_0n[1] = tri_0n_tmp[1];
            tri_0n[2] = tri_0n_tmp[2];
            
            NegateVec3D(n_t0);
        }
        if(orient_t1<0.0)
        {
            tri_1n_tmp[0] = tri_1n[0];
            tri_1n_tmp[1] = tri_1n[2];
            tri_1n_tmp[2] = tri_1n[1];
            
            tri_1n[0] = tri_1n_tmp[0];
            tri_1n[1] = tri_1n_tmp[1];
            tri_1n[2] = tri_1n_tmp[2];
            
            NegateVec3D(n_t10);
        }
        
        v_t0->c0 = xcn_g->getVal(tri_0n[1],0)-xcn_g->getVal(tri_0n[0],0);
        v_t0->c1 = xcn_g->getVal(tri_0n[1],1)-xcn_g->getVal(tri_0n[0],1);
        v_t0->c2 = xcn_g->getVal(tri_0n[1],2)-xcn_g->getVal(tri_0n[0],2);

        v_t1->c0 = xcn_g->getVal(tri_0n[2],0)-xcn_g->getVal(tri_0n[0],0);
        v_t1->c1 = xcn_g->getVal(tri_0n[2],1)-xcn_g->getVal(tri_0n[0],1);
        v_t1->c2 = xcn_g->getVal(tri_0n[2],2)-xcn_g->getVal(tri_0n[0],2);
        //Vec3D* n_t0_v1 = ComputeSurfaceNormal(v_t0,v_t1);
        //double orient_t0_check = DotVec3D(r0,n_t0_v1);
        
        v_t10->c0 = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10->c1 = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10->c2 = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);

        v_t11->c0 = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11->c1 = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11->c2 = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        //Vec3D* n_t10_v1 = ComputeSurfaceNormal(v_t10,v_t11);
        //double orient_t1_check = DotVec3D(r0,n_t10_v1);
        
        //std::cout << "check = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
        
        prism0[0] = tri_0n[0];
        prism0[1] = tri_0n[1];
        prism0[2] = tri_0n[2];
        
        prism1[0] = tri_1n[0];
        prism1[1] = tri_1n[1];
        prism1[2] = tri_1n[2];
        
        v_t0->c0 = xcn_g->getVal(prism0[1],0)-xcn_g->getVal(prism0[0],0);
        v_t0->c1 = xcn_g->getVal(prism0[1],1)-xcn_g->getVal(prism0[0],1);
        v_t0->c2 = xcn_g->getVal(prism0[1],2)-xcn_g->getVal(prism0[0],2);

        v_t1->c0 = xcn_g->getVal(prism0[2],0)-xcn_g->getVal(prism0[0],0);
        v_t1->c1 = xcn_g->getVal(prism0[2],1)-xcn_g->getVal(prism0[0],1);
        v_t1->c2 = xcn_g->getVal(prism0[2],2)-xcn_g->getVal(prism0[0],2);
        Vec3D* n_t0_v2 = ComputeSurfaceNormal(v_t0,v_t1);
        //orient_t0_check = DotVec3D(r0,n_t0_v2);
//
        //n_t0 = n_t0_v2;
        v_t10->c0 = xcn_g->getVal(prism1[1],0)-xcn_g->getVal(prism1[0],0);
        v_t10->c1 = xcn_g->getVal(prism1[1],1)-xcn_g->getVal(prism1[0],1);
        v_t10->c2 = xcn_g->getVal(prism1[1],2)-xcn_g->getVal(prism1[0],2);

        v_t11->c0 = xcn_g->getVal(prism1[2],0)-xcn_g->getVal(prism1[0],0);
        v_t11->c1 = xcn_g->getVal(prism1[2],1)-xcn_g->getVal(prism1[0],1);
        v_t11->c2 = xcn_g->getVal(prism1[2],2)-xcn_g->getVal(prism1[0],2);
        Vec3D* n_t10_v2 = ComputeSurfaceNormal(v_t10,v_t11);
        //orient_t1_check = DotVec3D(r0,n_t10_v2);
        
        //n_t10 = n_t10_v2;
//
//        if(orient_t0_check > 0 && orient_t1_check > 0)
//        {
//            //std::cout << "Correct = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
//            fwrong++;
//        }
//        if(orient_t0_check < 0 && orient_t1_check < 0)
//        {
//            //std::cout << "Incorrect = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
//            fright++;
//        }
        
        if(orient0<0.0)
        {
            NegateVec3D(nbf);
            face_turned[0] = face[0];
            face_turned[1] = face[3];
            face_turned[2] = face[2];
            face_turned[3] = face[1];
        }
        else
        {
            face_turned[0] = face[0];
            face_turned[1] = face[1];
            face_turned[2] = face[2];
            face_turned[3] = face[3];
        }
        face.clear();
        std::vector<Element*> PElements(nLayer*2);
        std::vector<std::vector<int> > PPrisms(nLayer*2);
        for(int c=0;c<nLayer;c++)
        {
            
            for(int k=0;k<8;k++)
            {
               loc_vid     = ien_g->getVal(elid_cur,k);
               Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
               Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
               Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
            }
            
            //int changed = ChkHexorient(Pijk,Pijk_id);
            
            Vert* Vijk = ComputeCenterCoord(Pijk, 8);
            std::vector<std::vector<int> > face_id_stored(6);
            std::vector<std::vector<Vert*> > face_stored(6);
            map<int,std::set<int> > local_node2node_element;
            std::vector<std::map<int,std::set<int> > > local_node2node_face(6);
            std::vector<std::map<int,int> > local_node2opponode_face(6);
            for(int k=0;k<6;k++)
            {
                int fid = ief_g->getVal(elid_cur,k);
                Vert* Vface2  = new Vert;
                
                std::vector<int> faceVert_IDs(4);
                std::vector<Vert*> face2;
                for(int r=0;r<4;r++)
                {
                    int vid  = ifn_g->getVal(fid,r);
                    
                    Vert* V  = new Vert;
                    V->x     = xcn_g->getVal(vid,0);
                    V->y     = xcn_g->getVal(vid,1);
                    V->z     = xcn_g->getVal(vid,2);
                    Vface2->x = Vface2->x+V->x;
                    Vface2->y = Vface2->y+V->y;
                    Vface2->z = Vface2->z+V->z;
                    face2.push_back(V);
                    
                    faceVert_IDs[r] = vid;
                }
                
                local_node2node_element[ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,1));
                local_node2node_element[ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,3));
                local_node2node_element[ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,0));
                local_node2node_element[ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,2));
                local_node2node_element[ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,1));
                local_node2node_element[ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,3));
                local_node2node_element[ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,2));
                local_node2node_element[ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,0));
                
                local_node2node_face[k][ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,1));
                local_node2node_face[k][ifn_g->getVal(fid,0)].insert(ifn_g->getVal(fid,3));
                local_node2node_face[k][ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,0));
                local_node2node_face[k][ifn_g->getVal(fid,1)].insert(ifn_g->getVal(fid,2));
                local_node2node_face[k][ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,1));
                local_node2node_face[k][ifn_g->getVal(fid,2)].insert(ifn_g->getVal(fid,3));
                local_node2node_face[k][ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,2));
                local_node2node_face[k][ifn_g->getVal(fid,3)].insert(ifn_g->getVal(fid,0));
                
                local_node2opponode_face[k][ifn_g->getVal(fid,0)]=ifn_g->getVal(fid,2);
                local_node2opponode_face[k][ifn_g->getVal(fid,1)]=ifn_g->getVal(fid,3);
                local_node2opponode_face[k][ifn_g->getVal(fid,2)]=ifn_g->getVal(fid,0);
                local_node2opponode_face[k][ifn_g->getVal(fid,3)]=ifn_g->getVal(fid,1);

                Vface2->x = Vface2->x/4.0;
                Vface2->y = Vface2->y/4.0;
                Vface2->z = Vface2->z/4.0;
                                
                Vec3D* r00 = new Vec3D;
                r00->c0 = (Vface2->x-Vijk->x);
                r00->c1 = (Vface2->y-Vijk->y);
                r00->c2 = (Vface2->z-Vijk->z);
                Vec3D* v00 = new Vec3D;
                v00->c0 = face2[1]->x-face2[0]->x;
                v00->c1 = face2[1]->y-face2[0]->y;
                v00->c2 = face2[1]->z-face2[0]->z;
                Vec3D* v11 = new Vec3D;
                v11->c0 = face2[3]->x-face2[0]->x;
                v11->c1 = face2[3]->y-face2[0]->y;
                v11->c2 = face2[3]->z-face2[0]->z;
                
                Vec3D* n00        = ComputeSurfaceNormal(v00,v11);
                double orient00   = DotVec3D(r00,n00);
                
                if(orient00<0.0)
                {
                    NegateVec3D(n00);
                    face_turned2[0] = face2[0];
                    face_turned2[1] = face2[3];
                    face_turned2[2] = face2[2];
                    face_turned2[3] = face2[1];
                }
                else
                {
                    face_turned2[0] = face2[0];
                    face_turned2[1] = face2[1];
                    face_turned2[2] = face2[2];
                    face_turned2[3] = face2[3];
                }
                
                face_stored[k]      =   face_turned2;
                dp[k]               =   DotVec3D(nbf,n00);
                dpvec[k]            =   n00;
                face_id_stored[k]   =   faceVert_IDs;
            }
            
            std::set<int>::iterator itset;
            for(itset=local_node2node_element[bvid].begin();itset!=local_node2node_element[bvid].end();itset++)
            {
                if(conn_bvid.find(*itset)==conn_bvid.end())
                {
                    opposite_bvid = *itset;
                }
            }
        
            int min_index  = std::min_element(dp.begin(),dp.end())-dp.begin();
            //double min_val = *std::min_element(dp.begin(),dp.end());

            int fid_new                  = ief_g->getVal(elid_cur,min_index);
            nbf                          = dpvec[min_index];
            std::vector<int> faceVertIDs = face_id_stored[min_index];
            std::vector<Vert*> faceupdate = face_stored[min_index];
            std::map<int,std::set<int> > node2node_face = local_node2node_face[min_index];
            
            std::set<int>::iterator itu;
            opposite_tri[0] = opposite_bvid;
            int l = 1;
            for(itu=node2node_face[opposite_bvid].begin();itu!=node2node_face[opposite_bvid].end();itu++)
            {
                opposite_tri[l] = *itu;
                l++;
            }
            
            Vec3D* v_toppo0 = new Vec3D;
            v_toppo0->c0 = xcn_g->getVal(opposite_tri[1],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo0->c1 = xcn_g->getVal(opposite_tri[1],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo0->c2 = xcn_g->getVal(opposite_tri[1],2)-xcn_g->getVal(opposite_tri[0],2);
            Vec3D* v_toppo1 = new Vec3D;
            v_toppo1->c0 = xcn_g->getVal(opposite_tri[2],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo1->c1 = xcn_g->getVal(opposite_tri[2],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo1->c2 = xcn_g->getVal(opposite_tri[2],2)-xcn_g->getVal(opposite_tri[0],2);
            
            Vec3D* n_toppo0        = ComputeSurfaceNormal(v_toppo0,v_toppo1);
            double orient0oppo0    = DotVec3D(n_t0_v2 ,n_toppo0 );
            if(orient0oppo0>0)
            {
                opposite_tri_n[0] = opposite_tri[0];
                opposite_tri_n[1] = opposite_tri[2];
                opposite_tri_n[2] = opposite_tri[1];
                
                opposite_tri[0] = opposite_tri_n[0];
                opposite_tri[1] = opposite_tri_n[1];
                opposite_tri[2] = opposite_tri_n[2];
                
                cnt_turn++;
            }

            v_toppo0->c0 = xcn_g->getVal(opposite_tri[1],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo0->c1 = xcn_g->getVal(opposite_tri[1],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo0->c2 = xcn_g->getVal(opposite_tri[1],2)-xcn_g->getVal(opposite_tri[0],2);
            
            v_toppo1->c0 = xcn_g->getVal(opposite_tri[2],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo1->c1 = xcn_g->getVal(opposite_tri[2],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo1->c2 = xcn_g->getVal(opposite_tri[2],2)-xcn_g->getVal(opposite_tri[0],2);
            n_toppo0        = ComputeSurfaceNormal(v_toppo0, v_toppo1);
            
            orient0oppo0    = DotVec3D(n_t0_v2 , n_toppo0 );
            
            
            
            opposite_tri1[0] = local_node2opponode_face[min_index][opposite_bvid];
            opposite_tri1[1] = opposite_tri[2];
            opposite_tri1[2] = opposite_tri[1];
            
            Vec3D* v_toppo10 = new Vec3D;
            v_toppo10->c0 = xcn_g->getVal(opposite_tri1[1],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo10->c1 = xcn_g->getVal(opposite_tri1[1],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo10->c2 = xcn_g->getVal(opposite_tri1[1],2)-xcn_g->getVal(opposite_tri1[0],2);
            Vec3D* v_toppo11 = new Vec3D;
            v_toppo11->c0 = xcn_g->getVal(opposite_tri1[2],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo11->c1 = xcn_g->getVal(opposite_tri1[2],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo11->c2 = xcn_g->getVal(opposite_tri1[2],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            Vec3D* n_toppo10        = ComputeSurfaceNormal(v_toppo10,v_toppo11);
            double orient0oppo10    = DotVec3D(n_t10_v2 , n_toppo10 );

            if(orient0oppo10>0)
            {
                //std::cout << "orientation " << orient0oppo10 << std::endl;
                opposite_tri1_n[0] = opposite_tri1[0];
                opposite_tri1_n[1] = opposite_tri1[2];
                opposite_tri1_n[2] = opposite_tri1[1];
                
                opposite_tri1[0] = opposite_tri1_n[0];
                opposite_tri1[1] = opposite_tri1_n[1];
                opposite_tri1[2] = opposite_tri1_n[2];
                
                cnt_turn++;
            }

            v_toppo10->c0 = xcn_g->getVal(opposite_tri1[1],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo10->c1 = xcn_g->getVal(opposite_tri1[1],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo10->c2 = xcn_g->getVal(opposite_tri1[1],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            v_toppo11->c0 = xcn_g->getVal(opposite_tri1[2],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo11->c1 = xcn_g->getVal(opposite_tri1[2],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo11->c2 = xcn_g->getVal(opposite_tri1[2],2)-xcn_g->getVal(opposite_tri1[0],2);
            n_toppo10        = ComputeSurfaceNormal(v_toppo10, v_toppo11);
            
            orient0oppo10    = DotVec3D(n_t10_v2 , n_toppo10 );
//            if(orient0oppo10>0)
//            {
//                std::cout << " Still not changed "<<std::endl;
//            }
            if(orient0oppo0>0 || orient0oppo10>0)
            {
                or0++;
                //std::cout << "orient0oppoorient0oppo " << orient0oppo0 << " " << orient0oppo10 << std::endl;
            }
            
            
//            for(itset=local_node2node_element[prism0[0]].begin();itset!=local_node2node_element[prism0[0]].end();itset++)
//            {
//
//                if(*itset!=prism0[1] && *itset!=prism0[2])
//                {
//                    opposite_p00 = *itset;
//                }
//            }
//            for(itset=local_node2node_element[prism0[1]].begin();itset!=local_node2node_element[prism0[1]].end();itset++)
//            {
//                if(*itset!=prism0[0] && *itset!=prism0[2])
//                {
//                    opposite_p01 = *itset;
//                }
//            }
//            for(itset=local_node2node_element[prism0[2]].begin();itset!=local_node2node_element[prism0[2]].end();itset++)
//            {
//                if(*itset!=prism0[0] && *itset!=prism0[1])
//                {
//                    opposite_p02 = *itset;
//                }
//            }
//
//
//            for(itset=local_node2node_element[prism1[0]].begin();itset!=local_node2node_element[prism1[0]].end();itset++)
//            {
//                if(*itset!=prism1[1] && *itset!=prism1[2])
//                {
//                    opposite_p10 = *itset;
//                }
//            }
//            for(itset=local_node2node_element[prism1[1]].begin();itset!=local_node2node_element[prism1[1]].end();itset++)
//            {
//                if(*itset!=prism1[0] && *itset!=prism1[2])
//                {
//                    opposite_p11 = *itset;
//                }
//            }
//            for(itset=local_node2node_element[prism1[2]].begin();itset!=local_node2node_element[prism1[2]].end();itset++)
//            {
//                if(*itset!=prism1[0] && *itset!=prism1[1])
//                {
//                    opposite_p12 = *itset;
//                }
//            }
            
            prism0[3] = opposite_tri[0];
            prism0[4] = opposite_tri[1];
            prism0[5] = opposite_tri[2];
//            prism0.push_back(opposite_tri[0]);
//            prism0.push_back(opposite_tri[1]);
//            prism0.push_back(opposite_tri[2]);
            
            prism1[3] = opposite_tri1[0];
            prism1[4] = opposite_tri1[1];
            prism1[5] = opposite_tri1[2];
//            prism1.push_back(opposite_tri1[0]);
//            prism1.push_back(opposite_tri1[1]);
//            prism1.push_back(opposite_tri1[2]);

            NegateVec3D(nbf);

            int gEl0=ife_g->getVal(fid_new,0);
            int gEl1=ife_g->getVal(fid_new,1);

            if(gEl0==elid_cur)
            {
                elid_next = gEl1;
            }
            else if(gEl1==elid_cur)
            {
                elid_next = gEl0;
            }
            
            if(c<nLayer-1)
            {
                layer.push_back(elid_next);
            }
            if(c==nLayer-1)
            {
                mesh_topology_bl->outer_shell_faces.push_back(fid_new);
            }
            
            
            // PRISM 0==================================================================================
            prismStored0[0] = prism0[0];prismStored0[1] = prism0[1];prismStored0[2] = prism0[2];
            prismStored0[3] = prism0[3];prismStored0[4] = prism0[4];prismStored0[5] = prism0[5];
            //std::cout << "prims0 " << glob_el_id << " :: " << prism0[0] << " " << prism0[1] << " " << prism0[2] << " " << prism0[3] << " " << prism0[4] << " " << prism0[5] << std::endl;
            Element* P0     = new Element;
            P0->GlobalNodes = prismStored0;
            P0->globID      = glob_el_id;
            
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[2]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria0.insert(prismStored0[0]);
            tria0.insert(prismStored0[1]);
            tria0.insert(prismStored0[2]);
            if(tria_ref_map.find(tria0)!=tria_ref_map.end())
            {
                int ref0 = tria_ref_map[tria0];
                std::vector<int> bctria(3);
                bctria[0] = prismStored0[0];
                bctria[1] = prismStored0[1];
                bctria[2] = prismStored0[2];
                mesh_topology_bl->bcTria[ref0].push_back(bctria);
            }
            
            
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[5]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria1.insert(prismStored0[3]);
            tria1.insert(prismStored0[4]);
            tria1.insert(prismStored0[5]);
            if(tria_ref_map.find(tria1)!=tria_ref_map.end())
            {
                int ref1 = tria_ref_map[tria1];
                std::vector<int> bctria(3);
                bctria[0] = prismStored0[3];
                bctria[1] = prismStored0[4];
                bctria[2] = prismStored0[5];
                mesh_topology_bl->bcTria[ref1].push_back(bctria);
            }
            
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[3]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad0.insert(prismStored0[0]);
            quad0.insert(prismStored0[2]);
            quad0.insert(prismStored0[4]);
            quad0.insert(prismStored0[3]);
            if(quad_ref_map.find(quad0)!=quad_ref_map.end())
            {
                int ref0 = quad_ref_map[quad0];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[0];
                bcquad[1] = prismStored0[2];
                bcquad[2] = prismStored0[4];
                bcquad[3] = prismStored0[3];
                mesh_topology_bl->bcQuad[ref0].push_back(bcquad);
            }
            
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[2]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad1.insert(prismStored0[1]);
            quad1.insert(prismStored0[5]);
            quad1.insert(prismStored0[4]);
            quad1.insert(prismStored0[2]);
            
            if(quad_ref_map.find(quad1)!=quad_ref_map.end())
            {
                int ref1 = quad_ref_map[quad1];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[1];
                bcquad[1] = prismStored0[5];
                bcquad[2] = prismStored0[4];
                bcquad[3] = prismStored0[2];
                mesh_topology_bl->bcQuad[ref1].push_back(bcquad);

            }
            
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[1]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad2.insert(prismStored0[0]);
            quad2.insert(prismStored0[3]);
            quad2.insert(prismStored0[5]);
            quad2.insert(prismStored0[1]);
            if(quad_ref_map.find(quad2)!=quad_ref_map.end())
            {
                int ref2 = quad_ref_map[quad2];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[0];
                bcquad[1] = prismStored0[3];
                bcquad[2] = prismStored0[5];
                bcquad[3] = prismStored0[1];
                mesh_topology_bl->bcQuad[ref2].push_back(bcquad);
            }
            
            tria0.clear();
            tria1.clear();
            quad0.clear();
            quad1.clear();
            quad2.clear();
            
            glob_el_id = glob_el_id+1;
            
            // PRISM 1================================================================================
            
            prismStored1[0] = prism1[0];prismStored1[1] = prism1[1];prismStored1[2] = prism1[2];
            prismStored1[3] = prism1[3];prismStored1[4] = prism1[4];prismStored1[5] = prism1[5];
            
            Element* P1     = new Element;
            P1->GlobalNodes = prismStored1;
            P1->globID      = glob_el_id;
            
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[2]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria0.insert(prismStored1[0]);
            tria0.insert(prismStored1[1]);
            tria0.insert(prismStored1[2]);
            if(tria_ref_map.find(tria0)!=tria_ref_map.end())
            {
                int ref0 = tria_ref_map[tria0];
                std::vector<int> bctria(3);
                bctria[0] = prismStored1[0];
                bctria[1] = prismStored1[1];
                bctria[2] = prismStored1[2];
                mesh_topology_bl->bcTria[ref0].push_back(bctria);
            }
            
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[5]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria1.insert(prismStored1[3]);
            tria1.insert(prismStored1[4]);
            tria1.insert(prismStored1[5]);
            if(tria_ref_map.find(tria1)!=tria_ref_map.end())
            {
                int ref1 = tria_ref_map[tria1];
                std::vector<int> bctria(3);
                bctria[0] = prismStored1[3];
                bctria[1] = prismStored1[4];
                bctria[2] = prismStored1[5];
                mesh_topology_bl->bcTria[ref1].push_back(bctria);
            }
            
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[2]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[3]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad0.insert(prismStored1[0]);
            quad0.insert(prismStored1[2]);
            quad0.insert(prismStored1[4]);
            quad0.insert(prismStored1[3]);
            if(quad_ref_map.find(quad0)!=quad_ref_map.end())
            {
                int ref0 = quad_ref_map[quad0];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[0];
                bcquad[1] = prismStored1[2];
                bcquad[2] = prismStored1[4];
                bcquad[3] = prismStored1[3];
                mesh_topology_bl->bcQuad[ref0].push_back(bcquad);

            }
            
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[5]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[2]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad1.insert(prismStored1[1]);
            quad1.insert(prismStored1[5]);
            quad1.insert(prismStored1[4]);
            quad1.insert(prismStored1[2]);
            
            if(quad_ref_map.find(quad1)!=quad_ref_map.end())
            {
                int ref1 = quad_ref_map[quad1];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[1];
                bcquad[1] = prismStored1[5];
                bcquad[2] = prismStored1[4];
                bcquad[3] = prismStored1[2];
                mesh_topology_bl->bcQuad[ref1].push_back(bcquad);

            }
            
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[5]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[1]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad2.insert(prismStored1[0]);
            quad2.insert(prismStored1[3]);
            quad2.insert(prismStored1[5]);
            quad2.insert(prismStored1[1]);
            if(quad_ref_map.find(quad2)!=quad_ref_map.end())
            {
                int ref2 = quad_ref_map[quad2];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[0];
                bcquad[1] = prismStored1[3];
                bcquad[2] = prismStored1[5];
                bcquad[3] = prismStored1[1];
                mesh_topology_bl->bcQuad[ref2].push_back(bcquad);

            }
            
            tria0.clear();
            tria1.clear();
            quad0.clear();
            quad1.clear();
            quad2.clear();
            
            
            glob_el_id = glob_el_id+1;
        //    std::cout << prismStored0[0] << " " << prismStored0[1] << " " << prismStored0[2] << " " << prismStored0[3] << " " << prismStored0[4] << " " << prismStored0[5] << std::endl;
         //   std::cout << prismStored1[0] << " " << prismStored1[1] << " " << prismStored1[2] << " " << prismStored1[3] << " " << prismStored1[4] << " " << prismStored1[5] << std::endl;
            PElements[c*2+0]=P0;
            PElements[c*2+1]=P1;
            PPrisms[c*2+0] = prismStored0;
            PPrisms[c*2+1] = prismStored1;

            //prism0.clear();
            //prism1.clear();
                    
            // Store the same initial triangles as the determined opposite triangles.
            // However change orientation of the nodes.
            
            prism0[0] = opposite_tri[0];
            prism0[1] = opposite_tri[2];
            prism0[2] = opposite_tri[1];
            
            prism1[0] = opposite_tri1[0];
            prism1[1] = opposite_tri1[2];
            prism1[2] = opposite_tri1[1];
            
//            prism0.push_back(opposite_tri[0]);
//            prism0.push_back(opposite_tri[1]);
//            prism0.push_back(opposite_tri[2]);
//            prism1.push_back(opposite_tri1[0]);
//            prism1.push_back(opposite_tri1[1]);
//            prism1.push_back(opposite_tri1[2]);
            
            conn_bvid.clear();
            bvid = opposite_tri[0];
            conn_bvid.insert(opposite_tri[1]);
            conn_bvid.insert(opposite_tri[2]);
            
            local_node2node_element.clear();
            local_node2node_face.clear();
            local_node2opponode_face.clear();
            //delete P0;
            //delete P1;
            elid_cur = elid_next;
        }
//        prism0.clear();
//        prism1.clear();
        mesh_topology_bl->BLlayersPrisms[bfaceid]=PPrisms;
        mesh_topology_bl->BLlayersElements[bfaceid]=PElements;
        mesh_topology_bl->BLlayers[bfaceid]=layer;
        mesh_topology_bl->Nprisms = mesh_topology_bl->Nprisms+PPrisms.size();
         
    }
    
    delete[] Pijk_id;
    delete[] Pijk;
    
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing for extracting BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    std::vector<int> elements;
    for(itt=mesh_topology_bl->BLlayers.begin();itt!=mesh_topology_bl->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            elements.push_back(itt->second[q]);
        }
    }
    OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_NEW");
       
    return mesh_topology_bl;
}
