#include "adapt_bltopology.h"


BLShellInfo* GetBoundaryLayerVolume(int wall_id, int nLayer,
                                    Array<double>* xcn_g, Array<int>* ien_g, Array<int>* iee_g,
                                    Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g,
                                    std::map<int,std::vector<int> > bnd_face_map,
                                    std::map<int,int> vert_ref_map, MPI_Comm comm)
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
    
    std::vector<std::vector<int> > outer_shell_elements;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<std::vector<double> > dpvec(6);
    //std::cout << "Determining outer shell of BL mesh..." << std::endl;
    clock_t start;
    start = std::clock();
    //int* Pijk_id = new int[8];
    //double* Pijk = new double[8*3];
    std::vector<double> Pijk(8*3);
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    int bvid,opposite_bvid;
    int bvid_b;
    int fv1_b;
    int fv2_b;
    int fv3_b;
    int glob_el_id = 0;
    std::vector<double> Vface2(3);
//    Vec3D* r00 = new Vec3D;
//    Vec3D* v00 = new Vec3D;
//    Vec3D* v11 = new Vec3D;
    std::vector<double> r00(3);
    std::vector<double> v00(3);
    std::vector<double> v11(3);
  
    std::vector<double> Vface(3);
    std::vector<double> V(3);
  
    int nbface = 0;
    Vec3D* NextFace2OriginalElem = new Vec3D;
    Vec3D* NextElem2OriginalFace = new Vec3D;
   
    std::vector<double> VadjMinAngle(3);
    int elprev = -1;
    int pfid = -1;
    int j;
    int fid;
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
        
        for(int c=0;c<nLayer;c++)
        {
            for(int j=0;j<6;j++)
            {
                int testface_id = ief_g->getVal(elid_cur,j);
                
                if(testface_id == faceid)
                {
                    fid = j;
                }
            }
            
            if(fid==0)
            {
                pfid = 1;
            }
            if(fid==1)
            {
                pfid = 0;
            }
            if(fid==2)
            {
                pfid = 3;
            }
            if(fid==3)
            {
                pfid = 2;
            }
            if(fid==4)
            {
                pfid = 5;
            }
            if(fid==5)
            {
                pfid = 4;
            }
            
            int elnext = iee_g->getVal(elid_cur,pfid);
            
            elid_cur = elnext;
            
            layer.push_back(elid_cur);
        }
        
        
        
        
        
        BLinfo->BLlayers[bfaceid]=layer;
        if(layer.size()!=21)
        {
            std::cout << "layers = " << layer.size() << " :: ";
            for(int t=0;t<layer.size();t++)
            {
                std::cout << layer[t] << " ";
            }
            std::cout << std::endl;
        }
        
        
    }
    
//    delete[] Pijk;
//    delete[] Pijk_id;
//
//    delete r00;
//    delete v00;
//    delete v11;
//    delete Vface2;
//    delete Vface;
//    delete V;
    
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing extracting outer shell BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    
    std::cout << "before " << BLinfo->elements_set.size() << std::endl;

    for(itt=BLinfo->BLlayers.begin();itt!=BLinfo->BLlayers.end();itt++)
    {
        //std::cout << itt->first << " :: " << itt->second.size() <<  " == ";
        
        for(int q=0;q<itt->second.size();q++)
        {
            if(BLinfo->elements_set.find(itt->second[q])==BLinfo->elements_set.end())
            {
                BLinfo->elements_set.insert(itt->second[q]);
                //std::cout << itt->second[q] << ", ";
            }
        }
        //std::cout << endl;
    }

    std::cout << "after " << BLinfo->elements_set.size() << std::endl;
    return BLinfo;
}


void FindOuterShellBoundaryLayerMesh_V2(BLShellInfo* BLinfo, int wall_id, int nLayer,
                            Array<double>* xcn_g, Array<int>* ien_g, Array<int>* iee_g,
                            Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g,
                            std::map<int,std::vector<int> > bnd_face_map,
                            std::map<int,int> vert_ref_map, MPI_Comm comm)
{
//    BLShellInfo* BLinfo = new BLShellInfo;
    BLinfo->ShellRef = new Array<int>(xcn_g->getNrow(),1);
    int te1=0;
    int te2=0;
    int te3=0;
    for(int i=0;i<xcn_g->getNrow();i++)
    {
        
        BLinfo->ShellRef->setVal(i,0,3);

        
//        if(vert_ref_map.find(i)!=vert_ref_map.end())
//        {
//            BLinfo->ShellRef->setVal(i,0,100+vert_ref_map[i]);
//        }
//        else
//        {
//            BLinfo->ShellRef->setVal(i,0,-3);
//        }
    }
    
    std::vector<std::vector<int> > outer_shell_elements;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<std::vector<double> > dpvec(6);
    //std::cout << "Determining outer shell of BL mesh..." << std::endl;
    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    //double* Pijk = new double[8*3];
    std::vector<double> Pijk(8*3);
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    int bvid,opposite_anchor;
    int bvid_b;

    int glob_el_id = 0;
  
    std::vector<double> Vface2(3);
//    Vec3D* r00 = new Vec3D;
//    Vec3D* v00 = new Vec3D;
//    Vec3D* v11 = new Vec3D;
    std::vector<double> r00(3);
    std::vector<double> v00(3);
    std::vector<double> v11(3);
    std::vector<double> Vface(3);
    std::vector<double> V(3);

    int nbface = 0;
    Vec3D* NextFace2OriginalElem = new Vec3D;
    Vec3D* NextElem2OriginalFace = new Vec3D;
    std::vector<double> VadjMinAngle(3);

    int elprev = -1;
    //for(int bf=0;bf<1;bf++)
    int check_f_id;
    int pfid;
    int fv0_b;
    int fv1_b;
    int fv2_b;
    int fv3_b;
    std::vector<int> bface_vertids(4);
    std::cout << "Extracting started  " << bnd_face_map[wall_id].size() << std::endl;
    for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
    {

        int bfaceid = bnd_face_map[wall_id][bf];
        int faceid  = bfaceid;
        int elid0   = ife_g->getVal(faceid,0);
        int elid1   = ife_g->getVal(faceid,1);

        fv1_b = ifn_g->getVal(faceid,1);
        fv2_b = ifn_g->getVal(faceid,2);
        fv3_b = ifn_g->getVal(faceid,3);
        
        if(elid0<ien_g->getNrow())
        {
            elid_cur = elid0;
        }
        else
        {
            elid_cur = elid1;
        }
        
        BLinfo->BLlayers[bfaceid].push_back(elid_cur);
        
       
        std::vector<double> VoriginalElem(3);
        VoriginalElem[0] = 0.0;
        VoriginalElem[1] = 0.0;
        VoriginalElem[2] = 0.0;
        
        for(int k=0;k<8;k++)
        {
            loc_vid     = ien_g->getVal(elid_cur,k);
            VoriginalElem[0] = VoriginalElem[0]+xcn_g->getVal(loc_vid,0);
            VoriginalElem[1] = VoriginalElem[1]+xcn_g->getVal(loc_vid,1);
            VoriginalElem[2] = VoriginalElem[2]+xcn_g->getVal(loc_vid,2);
        }

        VoriginalElem[0] = VoriginalElem[0]/8.0;
        VoriginalElem[1] = VoriginalElem[1]/8.0;
        VoriginalElem[2] = VoriginalElem[2]/8.0;
        
        
        std::vector<double> VoriginalFace(3);
        VoriginalFace[0] = 0.0;
        VoriginalFace[1] = 0.0;
        VoriginalFace[2] = 0.0;
        //std::cout << "face = ";
        std::vector<std::vector<double> > face;
        std::vector<std::vector<double> > face_turned(4);
        for(int k=0;k<4;k++)
        {
            loc_vid = ifn_g->getVal(faceid,k);

            VoriginalFace[0] = VoriginalFace[0]+xcn_g->getVal(loc_vid,0);
            VoriginalFace[1] = VoriginalFace[1]+xcn_g->getVal(loc_vid,1);
            VoriginalFace[2] = VoriginalFace[2]+xcn_g->getVal(loc_vid,2);
            
            std::vector<double> V(3);

            V[0] = xcn_g->getVal(loc_vid,0);
            V[1] = xcn_g->getVal(loc_vid,1);
            V[2] = xcn_g->getVal(loc_vid,2);
            
            face.push_back(V);
        }
        
        VoriginalFace[0] = VoriginalFace[0]/4.0;
        VoriginalFace[1] = VoriginalFace[1]/4.0;
        VoriginalFace[2] = VoriginalFace[2]/4.0;
        int bvid_b       = ifn_g->getVal(faceid,0);
        int anchor_vert  = bvid_b;
        std::set<int> conn_anchor;
        conn_anchor.insert(ifn_g->getVal(faceid,0));
        conn_anchor.insert(ifn_g->getVal(faceid,1));
        conn_anchor.insert(ifn_g->getVal(faceid,3));
        
        bface_vertids[0] = ifn_g->getVal(faceid,0);
        bface_vertids[1] = ifn_g->getVal(faceid,1);
        bface_vertids[2] = ifn_g->getVal(faceid,2);
        bface_vertids[3] = ifn_g->getVal(faceid,3);
        
        std::vector<double> r0(3);
        r0[0] = (VoriginalFace[0]-VoriginalElem[0]);
        r0[1] = (VoriginalFace[1]-VoriginalElem[1]);
        r0[2] = (VoriginalFace[2]-VoriginalElem[2]);
        std::vector<double> v0(3);
        v0[0] = face[1][0]-face[0][0];
        v0[1] = face[1][1]-face[0][1];
        v0[2] = face[1][2]-face[0][2];
        std::vector<double> v1(3);
        v1[0] = face[3][0]-face[0][0];
        v1[1] = face[3][1]-face[0][1];
        v1[2] = face[3][2]-face[0][2];
        std::vector<double> nbf     = ComputeSurfaceNormal(v0,v1);
        double orient0 = DotVec3D(r0,nbf);
        
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
        
        int copy_bface = bfaceid;

        for(int c=0;c<nLayer;c++)
        {
            
            std::map<int,std::set<int> > local_node2node_element;
            std::vector<std::map<int,std::set<int> > > local_node2node_face(6);
            std::vector<std::map<int,int> > local_node2opponode_face(6);
            
            for(int k=0;k<6;k++)
            {
                int el_adj_id  = iee_g->getVal(elid_cur,k);
                
                int face_id = ief_g->getVal(elid_cur,k);
 
                if(face_id == copy_bface)
                {
                    check_f_id = k;
                }
                
                local_node2node_element[ifn_g->getVal(face_id,0)].insert(ifn_g->getVal(face_id,1));
                local_node2node_element[ifn_g->getVal(face_id,0)].insert(ifn_g->getVal(face_id,3));
                local_node2node_element[ifn_g->getVal(face_id,1)].insert(ifn_g->getVal(face_id,0));
                local_node2node_element[ifn_g->getVal(face_id,1)].insert(ifn_g->getVal(face_id,2));
                local_node2node_element[ifn_g->getVal(face_id,2)].insert(ifn_g->getVal(face_id,1));
                local_node2node_element[ifn_g->getVal(face_id,2)].insert(ifn_g->getVal(face_id,3));
                local_node2node_element[ifn_g->getVal(face_id,3)].insert(ifn_g->getVal(face_id,2));
                local_node2node_element[ifn_g->getVal(face_id,3)].insert(ifn_g->getVal(face_id,0));
                
                local_node2node_face[k][ifn_g->getVal(face_id,0)].insert(ifn_g->getVal(face_id,1));
                local_node2node_face[k][ifn_g->getVal(face_id,0)].insert(ifn_g->getVal(face_id,3));
                local_node2node_face[k][ifn_g->getVal(face_id,1)].insert(ifn_g->getVal(face_id,0));
                local_node2node_face[k][ifn_g->getVal(face_id,1)].insert(ifn_g->getVal(face_id,2));
                local_node2node_face[k][ifn_g->getVal(face_id,2)].insert(ifn_g->getVal(face_id,1));
                local_node2node_face[k][ifn_g->getVal(face_id,2)].insert(ifn_g->getVal(face_id,3));
                local_node2node_face[k][ifn_g->getVal(face_id,3)].insert(ifn_g->getVal(face_id,2));
                local_node2node_face[k][ifn_g->getVal(face_id,3)].insert(ifn_g->getVal(face_id,0));
                
                local_node2opponode_face[k][ifn_g->getVal(face_id,0)]=ifn_g->getVal(face_id,2);
                local_node2opponode_face[k][ifn_g->getVal(face_id,1)]=ifn_g->getVal(face_id,3);
                local_node2opponode_face[k][ifn_g->getVal(face_id,2)]=ifn_g->getVal(face_id,0);
                local_node2opponode_face[k][ifn_g->getVal(face_id,3)]=ifn_g->getVal(face_id,1);
            }
            
            std::set<int>::iterator its;
            for(its=local_node2node_element[anchor_vert].begin();
                its!=local_node2node_element[anchor_vert].end();
                its++)
            {
                if(conn_anchor.find(*its)==conn_anchor.end())
                {
                    opposite_anchor = *its;
                }
            }
            
            switch(check_f_id) {
                case 0:
                    pfid = 1;
                    break;
                case 1:
                    pfid = 0;
                    break;
                case 2:
                    pfid = 3;
                    break;
                case 3:
                    pfid = 2;
                    break;
                case 4:
                    pfid = 5;
                    break;
                case 5:
                    pfid = 4;
                    break;
                default :
                    cout << "Error in determining opposite face in FindOuterShellBoundaryLayerMesh_V2()." << endl;
            }
            
            int max_index       = pfid;
            int elpicked        = iee_g->getVal(elid_cur,max_index);
            int fapicked        = ief_g->getVal(elid_cur,max_index);
        
            std::map<int,std::set<int> > node2node_face = local_node2node_face[max_index];
            std::set<int>::iterator itu;
            std::vector<int>opposite_tri(3);
            opposite_tri[0] = opposite_anchor;
            int l = 1;
            int iser = 0;
            if(node2node_face.find(opposite_anchor)!=node2node_face.end())
            {
                iser = 1;
            }
            
            //std::cout << "check c again " << c << " " << node2node_face.size() << " " << iser << " " << max_index << std::endl;

            for(itu=node2node_face[opposite_anchor].begin();
                itu!=node2node_face[opposite_anchor].end();
                itu++)
            {
                opposite_tri[l] = *itu;
                //std::cout << "c - nlayer = opposite_anchor " << opposite_anchor << " " << opposite_tri[l] << " " << l << " c = "<< c << " " << nLayer << " node2node_face " << node2node_face.size() << std::endl;
                
                l++;
            }
            
            if(c<nLayer-1)
            {
                BLinfo->BLlayers[bfaceid].push_back(elpicked);
            }
            
            if(c==nLayer-1)
            {
                for(int u=0;u<4;u++)
                {
                    int vInterFaceID              = ifn_g->getVal(fapicked,u);
                    BLinfo->shellVrts.insert(vInterFaceID);
                }
                
                int fid_new = fapicked;
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,0),0,999);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,1),0,999);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,2),0,999);
                BLinfo->ShellRef->setVal(ifn_g->getVal(fid_new,3),0,999);
                
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
                
                BLinfo->ShellFace2BFace[fid_new] = bfaceid;
                BLinfo->BFace2ShellFace[bfaceid] = fid_new;
                
                std::map<int,int> opposite_verts;
                opposite_verts[opposite_tri[0]]  = bvid_b;
                opposite_verts[opposite_tri[1]]  = fv1_b;
                opposite_verts[opposite_tri[2]]  = fv3_b;
                opposite_verts[local_node2opponode_face[max_index][opposite_anchor]]  = fv2_b;
                
                //std::cout << "opposite_verts (" << opposite_tri[0] << ", " << opposite_verts[opposite_tri[0]] << ") " << " (" << opposite_tri[1] << " " << opposite_verts[opposite_tri[1]] << ") (" << opposite_tri[2] << " " << opposite_verts[opposite_tri[2]] << ") " << std::endl;
                
                BLinfo->ShellFace2ShellVert2OppositeBoundaryVerts[fid_new] = opposite_verts;
                
            }
            
            anchor_vert = opposite_tri[0];
            
            conn_anchor.clear();
            //bvid = opposite_tri[0];
            conn_anchor.insert(opposite_tri[1]);
            conn_anchor.insert(opposite_tri[2]);
            
            local_node2node_element.clear();
            local_node2node_face.clear();
            local_node2opponode_face.clear();
            
            copy_bface = fapicked;
            elid_cur   = elpicked;
        }
    }

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing extracting outer shell BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    
    for(itt=BLinfo->BLlayers.begin();itt!=BLinfo->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            if(BLinfo->elements_set.find(itt->second[q])==BLinfo->elements_set.end())
            {
                BLinfo->elements_set.insert(itt->second[q]);
            }
        }
    }

    //return BLinfo;
}

BLShellInfo* FindOuterShellBoundaryLayerMesh(int wall_id, int nLayer,
                            Array<double>* xcn_g, Array<int>* ien_g, Array<int>* iee_g,
                            Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g,
                            std::map<int,std::vector<int> > bnd_face_map,
                            std::map<int,int> vert_ref_map, MPI_Comm comm)
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
    
    std::vector<std::vector<int> > outer_shell_elements;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<std::vector<double> > dpvec(6);
    //std::cout << "Determining outer shell of BL mesh..." << std::endl;
    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    //double* Pijk = new double[8*3];
    std::vector<double> Pijk(8*3);
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    int bvid,opposite_bvid;
    int bvid_b;
    int fv1_b;
    int fv2_b;
    int fv3_b;
    int glob_el_id = 0;
    std::vector<double> Vface2(3);
//    Vec3D* r00 = new Vec3D;
//    Vec3D* v00 = new Vec3D;
//    Vec3D* v11 = new Vec3D;
    std::vector<double> r00(3);
    std::vector<double> v00(3);
    std::vector<double> v11(3);
    std::vector<double> Vface(3);
    std::vector<double> V(3);

    int nbface = 0;
    Vec3D* NextFace2OriginalElem = new Vec3D;
    Vec3D* NextElem2OriginalFace = new Vec3D;
    
    std::vector<double> VadjMinAngle(3);
    int elprev = -1;
    //for(int bf=0;bf<1;bf++)
    int check_f_id;
    int pfid;
    std::cout << "Extracting started  " << bnd_face_map[wall_id].size() << std::endl;
    for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
    {

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
        
        BLinfo->BLlayers[bfaceid].push_back(elid_cur);
        
        
        std::vector<double> VoriginalElem(3);
        VoriginalElem[0] = 0.0;
        VoriginalElem[1] = 0.0;
        VoriginalElem[2] = 0.0;
        
        for(int k=0;k<8;k++)
        {
            loc_vid     = ien_g->getVal(elid_cur,k);
            VoriginalElem[0] = VoriginalElem[0]+xcn_g->getVal(loc_vid,0);
            VoriginalElem[1] = VoriginalElem[1]+xcn_g->getVal(loc_vid,1);
            VoriginalElem[2] = VoriginalElem[2]+xcn_g->getVal(loc_vid,2);
        }

        VoriginalElem[0] = VoriginalElem[0]/8.0;
        VoriginalElem[1] = VoriginalElem[1]/8.0;
        VoriginalElem[2] = VoriginalElem[2]/8.0;
        
        
        std::vector<double> VoriginalFace(3);
        VoriginalFace[0] = 0.0;
        VoriginalFace[1] = 0.0;
        VoriginalFace[2] = 0.0;
        //std::cout << "face = ";
        std::vector<std::vector<double> > face;
        std::vector<std::vector<double> > face_turned(4);
        for(int k=0;k<4;k++)
        {
            loc_vid = ifn_g->getVal(faceid,k);

            VoriginalFace[0] = VoriginalFace[0]+xcn_g->getVal(loc_vid,0);
            VoriginalFace[1] = VoriginalFace[1]+xcn_g->getVal(loc_vid,1);
            VoriginalFace[2] = VoriginalFace[2]+xcn_g->getVal(loc_vid,2);
            
            
            std::vector<double> V(3);
            V[0] = xcn_g->getVal(loc_vid,0);
            V[1] = xcn_g->getVal(loc_vid,1);
            V[2] = xcn_g->getVal(loc_vid,2);
            
            face.push_back(V);
        }
        

        VoriginalFace[0] = VoriginalFace[0]/4.0;
        VoriginalFace[1] = VoriginalFace[1]/4.0;
        VoriginalFace[2] = VoriginalFace[2]/4.0;
        
        std::vector<double> r0(3);
        r0[0] = (VoriginalFace[0]-VoriginalElem[0]);
        r0[1] = (VoriginalFace[1]-VoriginalElem[1]);
        r0[2] = (VoriginalFace[2]-VoriginalElem[2]);
        std::vector<double> v0(3);
        v0[0] = face[1][0]-face[0][0];
        v0[1] = face[1][1]-face[0][1];
        v0[2] = face[1][2]-face[0][2];
        std::vector<double> v1(3);
        v1[0] = face[3][0]-face[0][0];
        v1[1] = face[3][1]-face[0][1];
        v1[2] = face[3][2]-face[0][2];
        std::vector<double> nbf     = ComputeSurfaceNormal(v0,v1);
        double orient0 = DotVec3D(r0,nbf);
        
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
        
        int copy_bface = bfaceid;

        for(int c=0;c<nLayer;c++)
        {
            std::vector<double> angles;
            std::vector<int> eladjs;
            std::vector<int> faadjs;
            std::vector<int> faLocadjs;

            
            std::vector<std::vector<double> > elVerts;
            std::vector<std::vector<double> > faVerts;
            std::vector<std::vector<double> > normals;
            for(int k=0;k<6;k++)
            {
                int el_adj_id  = iee_g->getVal(elid_cur,k);
                
                int face_id = ief_g->getVal(elid_cur,k);
 
                if(face_id == copy_bface)
                {
                    check_f_id = k;
                }
                
                
                std::vector<double> Vface(3);
                Vface[0] = 0.0;
                Vface[1] = 0.0;
                Vface[2] = 0.0;
                
                std::vector<std::vector<double> > face2;
                
                for(int s=0;s<4;s++)
                {
                    loc_vid  = ifn_g->getVal(face_id,s);
                    Vface[0] = Vface[0] + xcn_g->getVal(loc_vid,0);
                    Vface[1] = Vface[1] + xcn_g->getVal(loc_vid,1);
                    Vface[2] = Vface[2] + xcn_g->getVal(loc_vid,2);
                    
                    
                    std::vector<double> V(3);
                    V[0]     = xcn_g->getVal(loc_vid,0);
                    V[1]     = xcn_g->getVal(loc_vid,1);
                    V[2]     = xcn_g->getVal(loc_vid,2);
                    face2.push_back(V);

                }
                Vface[0] = Vface[0]/4.0;
                Vface[1] = Vface[1]/4.0;
                Vface[2] = Vface[2]/4.0;
                
                std::vector<double> r00(3);
                r00[0] = (Vface[0]-VoriginalElem[0]);
                r00[1] = (Vface[1]-VoriginalElem[1]);
                r00[2] = (Vface[2]-VoriginalElem[2]);
                std::vector<double> v00(3);
                v00[0] = face2[1][0]-face2[0][0];
                v00[1] = face2[1][1]-face2[0][1];
                v00[2] = face2[1][2]-face2[0][2];
                std::vector<double> v11(3);
                v11[0] = face2[3][0]-face2[0][0];
                v11[1] = face2[3][1]-face2[0][1];
                v11[2] = face2[3][2]-face2[0][2];
                
                std::vector<double> n00        = ComputeSurfaceNormal(v00,v11);
                double orient00   = DotVec3D(r00,n00);
                dp[k]             = DotVec3D(nbf,n00);
            
                if(orient00<0.0)
                {
                    NegateVec3D(n00);
                }
            
                angles.push_back(DotVec3D(nbf,n00));
                normals.push_back(n00);
                
                std::vector<double> Velem(3);
                Velem[0] = 0.0;
                Velem[1] = 0.0;
                Velem[2] = 0.0;
                for(int f=0;f<8;f++)
                {
                    loc_vid     = ien_g->getVal(el_adj_id,f);
                    Velem[0]    = Velem[0] + xcn_g->getVal(loc_vid,0);
                    Velem[1]    = Velem[1] + xcn_g->getVal(loc_vid,1);
                    Velem[2]    = Velem[2] + xcn_g->getVal(loc_vid,2);
                }

                Velem[0]        = Velem[0]/8.0;
                Velem[1]        = Velem[1]/8.0;
                Velem[2]        = Velem[2]/8.0;

                eladjs.push_back(el_adj_id);
                faadjs.push_back(face_id);
                faLocadjs.push_back(k);
                elVerts.push_back(Velem);
                faVerts.push_back(Vface);
                face2.clear();
                                
            }
                
            
            if(check_f_id==0)
            {
                pfid = 1;
            }
            if(check_f_id==1)
            {
                pfid = 0;
            }
            if(check_f_id==2)
            {
                pfid = 3;
            }
            if(check_f_id==3)
            {
                pfid = 2;
            }
            if(check_f_id==4)
            {
                pfid = 5;
            }
            if(check_f_id==5)
            {
                pfid = 4;
            }
            
            std::vector<double>::iterator result = std::min_element(angles.begin(),angles.end());
            
            int max_index       = std::distance(angles.begin(), result);
            double min_val      = *std::min_element(angles.begin(),angles.end());
            
            int elpicked        = eladjs[max_index];
            int fapicked        = faadjs[max_index];
            int localFpicked    = faLocadjs[max_index];
            //std::cout << "fapicked " << fapicked << " " << elpicked << " minval " << min_val << std::endl;
//            //std::cout << "fac picked " << fapicked << " --> ";
//            for(int g=0;g<4;g++)
//            {
//                std::cout << ifn_g->getVal(fapicked,g) << " ";
//            }
//            std::cout << std::endl;
//
            VoriginalElem[0]    = elVerts[max_index][0];
            VoriginalElem[1]    = elVerts[max_index][1];
            VoriginalElem[2]    = elVerts[max_index][2];
                
            VoriginalFace[0]    = faVerts[max_index][0];
            VoriginalFace[1]    = faVerts[max_index][1];
            VoriginalFace[2]    = faVerts[max_index][2];
            
            nbf[0] = -normals[max_index][0];
            nbf[1] = -normals[max_index][1];
            nbf[2] = -normals[max_index][2];
            
            if(c<nLayer-1)
            {
                BLinfo->BLlayers[bfaceid].push_back(elpicked);

            }
            std::cout << "local faces " <<  check_f_id << " " << localFpicked << " " << max_index << std::endl;
            
            copy_bface = fapicked;
            
            //elprev   = elid_cur;
            elid_cur = elpicked;
            
            //std::cout << "elid_cur " << elid_cur << std::endl;
            angles.clear();
            eladjs.clear();
            faadjs.clear();
            faLocadjs.clear();
            //elVerts.clear();
            //normals.clear();
            faVerts.clear();
        }
    
        
        
    }
    
//    delete[] Pijk;
//    delete[] Pijk_id;
//
//    delete r00;
//    delete v00;
//    delete v11;
//    delete Vface2;
//    delete Vface;
//    delete V;
    
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing extracting outer shell BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    
    std::cout << "before " << BLinfo->elements_set.size() << std::endl;

    for(itt=BLinfo->BLlayers.begin();itt!=BLinfo->BLlayers.end();itt++)
    {
        //std::cout << itt->first << " :: " << itt->second.size() <<  " == "<<std::endl;
        
        for(int q=0;q<itt->second.size();q++)
        {
            if(BLinfo->elements_set.find(itt->second[q])==BLinfo->elements_set.end())
            {
                BLinfo->elements_set.insert(itt->second[q]);
                //std::cout << itt->second[q] << ", ";
            }
        }
        //std::cout << endl;
    }

    std::cout << "after " << BLinfo->elements_set.size() << std::endl;
    return BLinfo;
}



void ExtractBoundaryLayerMeshFromShell(Mesh_Topology_BL* mesh_topology_bl,std::vector<std::vector<int> > u_tris, BLShellInfo* BLshell, int wall_id, int nLayer, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g, std::map<int,std::vector<int> > bnd_face_map, std::map<std::set<int>,int> tria_ref_map, std::map<std::set<int>,int> quad_ref_map,  MPI_Comm comm)
{
    //Mesh_Topology_BL* mesh_topology_bl = new Mesh_Topology_BL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<double> dp(6);
    std::vector<std::vector<double> > dpvec(6);

    clock_t start;
    start = std::clock();
    //int* Pijk_id = new int[8];
    std::vector<int> Pijk_id(8);
    //double* Pijk = new double[8*3];
    std::vector<double> Pijk(8*3);
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;

    std::map<int,std::vector<std::vector<double> > > prisms;
    
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
    //mesh_topology_bl->Nprisms = 0;
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
    
    
    std::vector<double> Vface2(3);
    
    std::vector<double> v_t0(3);
    std::vector<double> v_t1(3);
    std::vector<double> v_t10(3);
    std::vector<double> v_t11(3);
    

    std::vector<double> Vface(3);
//    Vec3D* r00 = new Vec3D;
//    Vec3D* v00 = new Vec3D;
//    Vec3D* v11 = new Vec3D;
    
    std::vector<double> r00(3);
    std::vector<double> v00(3);
    std::vector<double> v11(3);
    
    int fv1_b,fv2_b,fv3_b;
    
    for(int bf=0;bf<bnd_face_map[wall_id].size();bf++)
    {
        int bvid=-1,obvid_i=-1,opposite_bvid=-1;
        std::vector<int> layer;
        int bfaceid      = bnd_face_map[wall_id][bf];
        int shell_faceid = BLshell->BFace2ShellFace[bfaceid];
        fv1_b = ifn_g->getVal(bfaceid,1);
        fv2_b = ifn_g->getVal(bfaceid,2);
        fv3_b = ifn_g->getVal(bfaceid,3);
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
        
//        std::cout << "tri_shell_0 " << shell_faceid   << " --- " << tri_shell_0[0] << " " << tri_shell_0[0] << " " << tri_shell_0[2] << std::endl;
//        std::cout << "tri_shell_1 " << tri_shell_1[0] << "     " << tri_shell_1[0]  << " " << tri_shell_1[2] << std::endl;
//
//        std::cout << "bound tri0 " << tri_bound_0[0] << " " << tri_bound_0[1] << " " << tri_bound_0[2] << std::endl;
//        std::cout << "bound tri1 " << tri_bound_1[0] << " " << tri_bound_1[1] << " " << tri_bound_1[2] << std::endl;
//
        //std::cout << "BFF "<< bf << " " << bfaceid << " " << wall_id << std::endl; 
        
        
        std::set<int> conn_bvid;
        std::vector<int> tri_0n(3);
        std::vector<int> tri_1n(3);
        std::vector<int> tri_0n_tmp(3);
        std::vector<int> tri_1n_tmp(3);
        for(int v=0;v<3;v++)
        {
            int fid0 = tri_bound_0[v];
            //std::cout << "tri_bound_0 (" << tri_bound_0[v] << " " << tri_shell_0[v] << ") ";
            if(std::find(tri_bound_1.begin(), tri_bound_1.end(), fid0) != tri_bound_1.end())
            {
                conn_bvid.insert(fid0);
            }
            //std::cout << std::endl;
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
        
        std::set<int> local_faces;
        
        //std::cout << "Element -> ";
        for(int k=0;k<8;k++)
        {
           loc_vid     = ien_g->getVal(elid_cur,k);
           Pijk_id[k]  = loc_vid;
           Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
           Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
           Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
        }

        std::vector<double> Vijk = ComputeCenterCoord(Pijk, 8);
        
        Vface[0]=0.0;
        Vface[1]=0.0;
        Vface[2]=0.0;
        std::vector<std::vector<double> > face;
        std::vector<std::vector<double> > face_turned(4);
        std::vector<std::vector<double> > face_turned2(4);
        for(int r=0;r<4;r++)
        {
            int vid  = ifn_g->getVal(faceid,r);
            
            
            std::vector<double> V(3);
            V[0]     = xcn_g->getVal(vid,0);
            V[1]     = xcn_g->getVal(vid,1);
            V[2]     = xcn_g->getVal(vid,2);
            Vface[0] = Vface[0]+V[0];
            Vface[1] = Vface[1]+V[1];
            Vface[2] = Vface[2]+V[2];
            
            face.push_back(V);
        }
        std::vector<int> tri0(3);
        std::vector<int> tri1(3);
        //std::cout << "End loop " << tri_0n[0] << " " << tri_0n[1] << " " << tri_0n[2] << std::endl;
       
        v_t0[0] = xcn_g->getVal(tri_0n[1],0)-xcn_g->getVal(tri_0n[0],0);
        v_t0[1] = xcn_g->getVal(tri_0n[1],1)-xcn_g->getVal(tri_0n[0],1);
        v_t0[2] = xcn_g->getVal(tri_0n[1],2)-xcn_g->getVal(tri_0n[0],2);
        
        v_t1[0] = xcn_g->getVal(tri_0n[2],0)-xcn_g->getVal(tri_0n[0],0);
        v_t1[1] = xcn_g->getVal(tri_0n[2],1)-xcn_g->getVal(tri_0n[0],1);
        v_t1[2] = xcn_g->getVal(tri_0n[2],2)-xcn_g->getVal(tri_0n[0],2);
        std::vector<double> n_t0 = ComputeSurfaceNormal(v_t0,v_t1);
        
        v_t10[0] = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10[1] = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10[2] = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);
        
        v_t11[0] = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11[1] = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11[2] = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        std::vector<double> n_t10 = ComputeSurfaceNormal(v_t10,v_t11);
        
        tri0[0] = tri_0n[0];
        tri0[1] = tri_0n[1];
        tri0[2] = tri_0n[2];
        mesh_topology_bl->BndFaces.push_back(tri0);

        tri1[0] = tri_1n[0];
        tri1[1] = tri_1n[1];
        tri1[2] = tri_1n[2];
        mesh_topology_bl->BndFaces.push_back(tri1);
        
        Vface[0] = Vface[0]/4.0;
        Vface[1] = Vface[1]/4.0;
        Vface[2] = Vface[2]/4.0;
                        
        std::vector<double> r0(3);
        double r0L = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])
                          +(Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])
                          +(Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));
        
        r0[0] = (Vface[0]-Vijk[0])/r0L;
        r0[1] = (Vface[1]-Vijk[1])/r0L;
        r0[2] = (Vface[2]-Vijk[2])/r0L;
        
        std::vector<double> v0(3);
        v0[0] = face[1][0]-face[0][0];
        v0[1] = face[1][1]-face[0][1];
        v0[2] = face[1][2]-face[0][2];
        std::vector<double> v1(3);
        v1[0] = face[3][0]-face[0][0];
        v1[1] = face[3][1]-face[0][1];
        v1[2] = face[3][2]-face[0][2];
        
        std::vector<double> nbf     = ComputeSurfaceNormal(v0,v1);
        
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
        
        v_t0[0] = xcn_g->getVal(tri_0n[1],0)-xcn_g->getVal(tri_0n[0],0);
        v_t0[1] = xcn_g->getVal(tri_0n[1],1)-xcn_g->getVal(tri_0n[0],1);
        v_t0[2] = xcn_g->getVal(tri_0n[1],2)-xcn_g->getVal(tri_0n[0],2);

        v_t1[0] = xcn_g->getVal(tri_0n[2],0)-xcn_g->getVal(tri_0n[0],0);
        v_t1[1] = xcn_g->getVal(tri_0n[2],1)-xcn_g->getVal(tri_0n[0],1);
        v_t1[2] = xcn_g->getVal(tri_0n[2],2)-xcn_g->getVal(tri_0n[0],2);
        
        //Vec3D* n_t0_v1 = ComputeSurfaceNormal(v_t0,v_t1);
        //double orient_t0_check = DotVec3D(r0,n_t0_v1);
        
        v_t10[0] = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10[1] = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10[2] = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);

        v_t11[0] = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11[1] = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11[2] = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        //Vec3D* n_t10_v1 = ComputeSurfaceNormal(v_t10,v_t11);
        //double orient_t1_check = DotVec3D(r0,n_t10_v1);
        
        //std::cout << "check = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
        
        prism0[0] = tri_0n[0];
        prism0[1] = tri_0n[1];
        prism0[2] = tri_0n[2];
        
        //std::cout << "prism0[0] " << prism0[0] << " " << prism0[1] << " " << prism0[2] << std::endl;
        
        prism1[0] = tri_1n[0];
        prism1[1] = tri_1n[1];
        prism1[2] = tri_1n[2];
        
        v_t0[0] = xcn_g->getVal(prism0[1],0)-xcn_g->getVal(prism0[0],0);
        v_t0[1] = xcn_g->getVal(prism0[1],1)-xcn_g->getVal(prism0[0],1);
        v_t0[2] = xcn_g->getVal(prism0[1],2)-xcn_g->getVal(prism0[0],2);

        v_t1[0] = xcn_g->getVal(prism0[2],0)-xcn_g->getVal(prism0[0],0);
        v_t1[1] = xcn_g->getVal(prism0[2],1)-xcn_g->getVal(prism0[0],1);
        v_t1[2] = xcn_g->getVal(prism0[2],2)-xcn_g->getVal(prism0[0],2);
        std::vector<double> n_t0_v2 = ComputeSurfaceNormal(v_t0,v_t1);
        //orient_t0_check = DotVec3D(r0,n_t0_v2);
//
        //n_t0 = n_t0_v2;
        v_t10[0] = xcn_g->getVal(prism1[1],0)-xcn_g->getVal(prism1[0],0);
        v_t10[1] = xcn_g->getVal(prism1[1],1)-xcn_g->getVal(prism1[0],1);
        v_t10[2] = xcn_g->getVal(prism1[1],2)-xcn_g->getVal(prism1[0],2);

        v_t11[0] = xcn_g->getVal(prism1[2],0)-xcn_g->getVal(prism1[0],0);
        v_t11[1] = xcn_g->getVal(prism1[2],1)-xcn_g->getVal(prism1[0],1);
        v_t11[2] = xcn_g->getVal(prism1[2],2)-xcn_g->getVal(prism1[0],2);
        std::vector<double> n_t10_v2 = ComputeSurfaceNormal(v_t10,v_t11);

        
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
        std::vector<std::vector<int> > PPrisms(nLayer*2);
        std::vector<Element*> PElements(nLayer*2);

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
            
            std::vector<double> Vijk = ComputeCenterCoord(Pijk, 8);
            std::vector<std::vector<int> > face_id_stored(6);
            std::vector<std::vector<std::vector<double> > > face_stored(6);
            map<int,std::set<int> > local_node2node_element;
            std::vector<std::map<int,std::set<int> > > local_node2node_face(6);
            std::vector<std::map<int,int> > local_node2opponode_face(6);
            for(int k=0;k<6;k++)
            {
                int fid = ief_g->getVal(elid_cur,k);
                Vface2[0] = 0.0;
                Vface2[1] = 0.0;
                Vface2[2] = 0.0;
                std::vector<int> faceVert_IDs(4);
                std::vector<std::vector<double> > face2;
                for(int r=0;r<4;r++)
                {
                    int vid  = ifn_g->getVal(fid,r);
                    
                    
                    std::vector<double> V(3);
                    V[0]     = xcn_g->getVal(vid,0);
                    V[1]     = xcn_g->getVal(vid,1);
                    V[2]     = xcn_g->getVal(vid,2);
                    Vface2[0] = Vface2[0]+V[0];
                    Vface2[1] = Vface2[1]+V[1];
                    Vface2[2] = Vface2[2]+V[2];
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

                Vface2[0] = Vface2[0]/4.0;
                Vface2[1] = Vface2[1]/4.0;
                Vface2[2] = Vface2[2]/4.0;
                
                r00[0] = (Vface2[0]-Vijk[0]);
                r00[1] = (Vface2[1]-Vijk[1]);
                r00[2] = (Vface2[2]-Vijk[2]);
            
                v00[0] = face2[1][0]-face2[0][0];
                v00[1] = face2[1][1]-face2[0][1];
                v00[2] = face2[1][2]-face2[0][2];
                
                v11[0] = face2[3][0]-face2[0][0];
                v11[1] = face2[3][1]-face2[0][1];
                v11[2] = face2[3][2]-face2[0][2];
                
                std::vector<double> n00 = ComputeSurfaceNormal(v00,v11);
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
            std::vector<std::vector<double> > faceupdate = face_stored[min_index];
            std::map<int,std::set<int> > node2node_face = local_node2node_face[min_index];
            
            std::set<int>::iterator itu;
            opposite_tri[0] = opposite_bvid;
            int l = 1;
            for(itu=node2node_face[opposite_bvid].begin();itu!=node2node_face[opposite_bvid].end();itu++)
            {
                opposite_tri[l] = *itu;
                l++;
            }
            
            std::vector<double> v_toppo0(3);
            v_toppo0[0] = xcn_g->getVal(opposite_tri[1],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo0[1] = xcn_g->getVal(opposite_tri[1],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo0[2] = xcn_g->getVal(opposite_tri[1],2)-xcn_g->getVal(opposite_tri[0],2);
            std::vector<double> v_toppo1(3);
            v_toppo1[0] = xcn_g->getVal(opposite_tri[2],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo1[1] = xcn_g->getVal(opposite_tri[2],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo1[2] = xcn_g->getVal(opposite_tri[2],2)-xcn_g->getVal(opposite_tri[0],2);
            
            std::vector<double> n_toppo0 = ComputeSurfaceNormal(v_toppo0,v_toppo1);
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

           
            v_toppo0[0] = xcn_g->getVal(opposite_tri[1],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo0[1] = xcn_g->getVal(opposite_tri[1],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo0[2] = xcn_g->getVal(opposite_tri[1],2)-xcn_g->getVal(opposite_tri[0],2);
            
            
            v_toppo1[0] = xcn_g->getVal(opposite_tri[2],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo1[1] = xcn_g->getVal(opposite_tri[2],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo1[2] = xcn_g->getVal(opposite_tri[2],2)-xcn_g->getVal(opposite_tri[0],2);
            std::vector<double> n_toppo01        = ComputeSurfaceNormal(v_toppo0, v_toppo1);
            
            orient0oppo0    = DotVec3D(n_t0_v2 , n_toppo01 );
            
            
            
            opposite_tri1[0] = local_node2opponode_face[min_index][opposite_bvid];
            opposite_tri1[1] = opposite_tri[2];
            opposite_tri1[2] = opposite_tri[1];
            
            std::vector<double> v_toppo10(3);
            v_toppo10[0] = xcn_g->getVal(opposite_tri1[1],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo10[1] = xcn_g->getVal(opposite_tri1[1],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo10[2] = xcn_g->getVal(opposite_tri1[1],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            std::vector<double> v_toppo11(3);
            v_toppo11[0] = xcn_g->getVal(opposite_tri1[2],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo11[1] = xcn_g->getVal(opposite_tri1[2],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo11[2] = xcn_g->getVal(opposite_tri1[2],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            std::vector<double> n_toppo10 = ComputeSurfaceNormal(v_toppo10,v_toppo11);
            
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

            v_toppo10[0] = xcn_g->getVal(opposite_tri1[1],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo10[1] = xcn_g->getVal(opposite_tri1[1],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo10[2] = xcn_g->getVal(opposite_tri1[1],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            v_toppo11[0] = xcn_g->getVal(opposite_tri1[2],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo11[1] = xcn_g->getVal(opposite_tri1[2],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo11[2] = xcn_g->getVal(opposite_tri1[2],2)-xcn_g->getVal(opposite_tri1[0],2);
            std::vector<double> n_toppo11 = ComputeSurfaceNormal(v_toppo10, v_toppo11);
            
            orient0oppo10    = DotVec3D(n_t10_v2 , n_toppo11 );

            
            

            
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
            
            
            // PRISM 0==================================================================================
            prismStored0[0] = prism0[0];prismStored0[1] = prism0[1];prismStored0[2] = prism0[2];
            prismStored0[3] = prism0[3];prismStored0[4] = prism0[4];prismStored0[5] = prism0[5];

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
        //mesh_topology_bl->BLlayers[bfaceid]=layer;
        mesh_topology_bl->BLlayersElements[bfaceid]=PElements;
        mesh_topology_bl->Nprisms = mesh_topology_bl->Nprisms+PPrisms.size();
        /* */
        ///delete Vijk;
    }
    
//    delete[] Pijk_id;
//    delete[] Pijk;
//    delete r00;
//    delete v00;
//    delete v11;
//    delete Vface2;
//    delete v_t0;
//    delete v_t1;
//    delete v_t10;
//    delete v_t11;
    
    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Timing for extracting BL mesh = " << duration << std::endl;
//    std::map<int,std::vector<int> >::iterator itt;
//    std::vector<int> elements;
//    for(itt=mesh_topology_bl->BLlayers.begin();itt!=mesh_topology_bl->BLlayers.end();itt++)
//    {
//        for(int q=0;q<itt->second.size();q++)
//        {
//            elements.push_back(itt->second[q]);
//        }
//    }
    //OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_NEW");
       
    //return mesh_topology_bl;
}
