#include "adapt_io.h"
#include "adapt_output.h"


std::vector<double> ReadMetricInputs(const char* fn_metric)
{
    std::ifstream fin;
    fin.open(fn_metric);
    if(!fin.is_open())
    {
        std::cout << "Error:: Make sure there is a metric.inp file in the directory where main.cpp resides. "<<std::endl;
        exit(0);
    }
    
    double v=0.0;
    std::vector<double> metric_inputs;
    int t=0;
    while(fin >> v)
    {
        metric_inputs.push_back(v);
       t++;
    }
    return metric_inputs;
}






void WriteUS3DGridFromMMG_itN(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, US3D* us3d)
{
    std::map<int,std::vector<int> > ref2bface;
    std::map<int,std::vector<int> > ref2bqface;
    std::set<set<int> > bfaces;
    std::set<set<int> > bqfaces;
    std::map<set<int>,int > btfaces_Ref;
    std::map<set<int>,int > bqfaces_Ref;
    std::set<int>face;
    std::set<int> bcrefs;
    int wr = 0;
    
    int nVerts = mmgMesh->np;
    int nTet = mmgMesh->ne;
    int nPrism = mmgMesh->nprism;
    
    for(int i=1;i<=mmgMesh->nt;i++)
    {
        if(mmgMesh->tria[i].ref>0 && mmgMesh->tria[i].ref!=20)// -1 is the tag for internal shell.
        {
            
            ref2bface[mmgMesh->tria[i].ref].push_back(i);
            face.insert(mmgMesh->tria[i].v[0]);
            face.insert(mmgMesh->tria[i].v[1]);
            face.insert(mmgMesh->tria[i].v[2]);
            
            if(btfaces_Ref.find(face)==btfaces_Ref.end())
            {
                btfaces_Ref[face] = mmgMesh->tria[i].ref;
            }
            
            bfaces.insert(face);
            if(bcrefs.find(mmgMesh->tria[i].ref)==bcrefs.end())
            {
                bcrefs.insert(mmgMesh->tria[i].ref);
            }
            face.clear();
        }
    }
    
    for(int i=1;i<=mmgMesh->nquad;i++)
    {
        if(mmgMesh->quadra[i].ref>0 && mmgMesh->quadra[i].ref!=2)// -1 is the tag for internal shell.
        {
            
            ref2bqface[mmgMesh->quadra[i].ref].push_back(i);
            face.insert(mmgMesh->quadra[i].v[0]);
            face.insert(mmgMesh->quadra[i].v[1]);
            face.insert(mmgMesh->quadra[i].v[2]);
            face.insert(mmgMesh->quadra[i].v[3]);
            
            if(bqfaces_Ref.find(face)==bqfaces_Ref.end())
            {
                bqfaces_Ref[face] = mmgMesh->quadra[i].ref;
            }
            
            bqfaces.insert(face);
            if(bcrefs.find(mmgMesh->quadra[i].ref)==bcrefs.end())
            {
                bcrefs.insert(mmgMesh->quadra[i].ref);
            }
            face.clear();
        }
    }
    
    
    //std::map<int,std::vector<int> > bnd_map = bnd_face_map;
    std::map<int,std::vector<int> >::iterator bnd_m;
    std::map<int,int> bnd_Ntri;
    std::map<int,int> bnd_Nquad;
    int i=0;
    std::set<int>::iterator refit;
    
    for(refit=bcrefs.begin();refit!=bcrefs.end();refit++)
    {
        int ref_inq = *refit;
        
        if(ref2bface.find(ref_inq)==ref2bface.end())
        {
            bnd_Ntri[ref_inq]=0;
        }
        else
        {
            bnd_Ntri[ref_inq]=ref2bface[ref_inq].size();
        }
        
        
        if(ref2bqface.find(ref_inq)==ref2bqface.end())
        {
            bnd_Nquad[ref_inq]=0;
        }
        else
        {
            bnd_Nquad[ref_inq]=ref2bqface[ref_inq].size();
        }
    }
    
    std::map<int,std::vector<std::vector<int> > > bctrias;
    std::map<int,std::vector<std::vector<int> > > bcquads;

    Array<double>* xcn_mmg = new Array<double>(mmgMesh->np,3);
    for(int i=0;i<mmgMesh->np;i++)
    {
        xcn_mmg->setVal(i,0,mmgMesh->point[i+1].c[0]);
        xcn_mmg->setVal(i,1,mmgMesh->point[i+1].c[1]);
        xcn_mmg->setVal(i,2,mmgMesh->point[i+1].c[2]);
    }
    
    std::cout<<"-- Writing in HDF5 format..."<<std::endl;
    hid_t ret;
    //Output the new grid.h5 which has the new vertices and ifn map.
    //===================================================================
    //===================================================================
    //===================================================================
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    plist_id               = H5P_DEFAULT;
    //H5Pset_fapl_mpio(plist_id, comm, info);
    hid_t file_id = H5Fcreate("grid_madam.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    hid_t status;
    hid_t att_space;
    hid_t attr_id;
    
    hsize_t dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    hid_t type =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type, 14);
    ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri[] = "US3D Grid File";
    status = H5Awrite(attr_id, type, &stri);
    H5Aclose(attr_id);
    
    hsize_t dimsf_att2 = 1;
    att_space = H5Screate_simple(1, &dimsf_att2, NULL);
    hid_t type2 =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type2, 5);
    ret = H5Tset_strpad(type2,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (file_id, "filevers", type2, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri2[] = "1.1.8";
    status = H5Awrite(attr_id, type2, &stri2);
    H5Aclose(attr_id);
    
    //====================================================================================
    // Add xcn map to the grid.h5 file
    //====================================================================================
    hsize_t     dimsf[2];
    hsize_t    count[2];              // hyperslab selection parameters
    hsize_t    offset[2];
    dimsf[0] = xcn_mmg->getNrow();
    dimsf[1] = xcn_mmg->getNcol();
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);

    hid_t dset_id = H5Dcreate(file_id, "xcn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    hid_t memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_mmg->data);
    delete xcn_mmg;
    //====================================================================================
    
    std::map<std::set<int>, int> qfacemap;
    std::map<std::set<int>, int>  facemap;
    std::set<std::set<int> > faces;
    std::set<std::set<int> > qfaces;
    std::map<int,std::vector<int> > face2node;
    std::set<int> face0;
    std::set<int> face1;
    std::set<int> face2;
    std::set<int> face3;
    std::set<int> qface0;
    std::set<int> qface1;
    std::set<int> qface2;
    
    
    int fid = 0;
    int vid0,vid1,vid2,vid3;
    std::map<int,int> lh;
    std::map<int,int> rh;
    
    std::map<int,int> Nlh;
    std::map<int,int> Nrh;
    
    int of = 0;
    int fset_cnt = 0;
    Array<int>* adapt_iet = new Array<int>(mmgMesh->ne+mmgMesh->nprism,1);
    // local face2vert_map for a tet in mmg  {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
    int bf = 0;
    int bq = 0;
    std::cout << "-- Constructing the new face-2-node and face-2-element map..."<<std::endl;
    int orient0;
    
    for(int i=1;i<=mmgMesh->ne;i++)
    {
        adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
        
        face0.insert(mmgMesh->tetra[i].v[1]);
        face0.insert(mmgMesh->tetra[i].v[2]);
        face0.insert(mmgMesh->tetra[i].v[3]);
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
            facemap[face0]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face0]] = i-1;
            faces.erase(face0);
            facemap.erase(face0);
        }
        
        face1.insert(mmgMesh->tetra[i].v[0]);
        face1.insert(mmgMesh->tetra[i].v[2]);
        face1.insert(mmgMesh->tetra[i].v[3]);
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            facemap[face1]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face1]] = i-1;
            faces.erase(face1);
            facemap.erase(face1);
        }
        
        
        
        
        face2.insert(mmgMesh->tetra[i].v[0]);
        face2.insert(mmgMesh->tetra[i].v[3]);
        face2.insert(mmgMesh->tetra[i].v[1]);
        if( faces.count(face2) != 1)
        {
            faces.insert(face2);
            facemap[face2]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face2]] = i-1;
            faces.erase(face2);
            facemap.erase(face2);
        }

        
        
        face3.insert(mmgMesh->tetra[i].v[0]);
        face3.insert(mmgMesh->tetra[i].v[2]);
        face3.insert(mmgMesh->tetra[i].v[1]);
        if( faces.count(face3) != 1)
        {
            faces.insert(face3);
            facemap[face3]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face3]] = i-1;
            faces.erase(face3);
            facemap.erase(face3);
        }
    
        face0.clear();
        face1.clear();
        face2.clear();
        face3.clear();
        
    }
    

    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
    int qfid  = 0;
    int fnew  = 0;
    int fiold = 0;

    for(int i=1;i<=mmgMesh->nprism;i++)
    {
        adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are dealing with prisms.
  
        face0.insert(mmgMesh->prism[i].v[0]);
        face0.insert(mmgMesh->prism[i].v[2]);
        face0.insert(mmgMesh->prism[i].v[1]);
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
            facemap[face0]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[facemap[face0]] = mmgMesh->ne+i-1;
            faces.erase(face0);
            facemap.erase(face0);
        }
        
        
        face1.insert(mmgMesh->prism[i].v[3]);
        face1.insert(mmgMesh->prism[i].v[4]);
        face1.insert(mmgMesh->prism[i].v[5]);
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            facemap[face1]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[facemap[face1]]  = mmgMesh->ne+i-1;
            faces.erase(face1);
            facemap.erase(face1);
        }
        
        
        
        // Quad faces //
        qface0.insert(mmgMesh->prism[i].v[0]);
        qface0.insert(mmgMesh->prism[i].v[3]);
        qface0.insert(mmgMesh->prism[i].v[4]);
        qface0.insert(mmgMesh->prism[i].v[1]);
        if( qfaces.count(qface0) != 1)
        {
            qfaces.insert(qface0);
            qfacemap[qface0]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[qfacemap[qface0]]  = mmgMesh->ne+i-1;
            qfaces.erase(qface0);
            qfacemap.erase(qface0);
        }

        
        
        
        qface1.insert(mmgMesh->prism[i].v[1]);
        qface1.insert(mmgMesh->prism[i].v[4]);
        qface1.insert(mmgMesh->prism[i].v[5]);
        qface1.insert(mmgMesh->prism[i].v[2]);
        if( qfaces.count(qface1) != 1)
        {
            qfaces.insert(qface1);
            qfacemap[qface1]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[qfacemap[qface1]]  = mmgMesh->ne+i-1;
            qfaces.erase(qface1);
            qfacemap.erase(qface1);
        }
        
        qface2.insert(mmgMesh->prism[i].v[0]);//1
        qface2.insert(mmgMesh->prism[i].v[2]);//2
        qface2.insert(mmgMesh->prism[i].v[5]);//5
        qface2.insert(mmgMesh->prism[i].v[3]);//4
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface2) != 1)
        {
            qfaces.insert(qface2);
            qfacemap[qface2]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[qfacemap[qface2]] = mmgMesh->ne+i-1;
            qfaces.erase(qface2);
            qfacemap.erase(qface2);
        }
        
        face0.clear();
        face1.clear();
        qface0.clear();
        qface1.clear();
        qface2.clear();
         
    }
    
    //====================================================================================
    //====================================================================================
    dimsf[0] = adapt_iet->getNrow();
    dimsf[1] = adapt_iet->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);

    dset_id = H5Dcreate(file_id, "iet", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_iet->data);
    delete adapt_iet;
    //====================================================================================
    //====================================================================================

    std::cout << "AND WE ARE DONE COUNTING " << std::endl;
    
    int ftot = 0;
    std::map<int,int> new2old;
//
    std::map<int,int>::iterator itm;
    int it;
    for(itm=lh.begin();itm!=lh.end();itm++)
    {
        it = itm->first;
        if(rh.find(it)!=rh.end())
        {
            new2old[it] = ftot;
            ftot++;
        }
    }
    
    std::cout << "AND WE ARE DONE DETERMINING THE ORDER" << std::endl;
    std::cout << "Starting creating massive array" << std::endl;

    Array<int>* adapt_ifn = new Array<int>(fid,8);
    std::cout << "Finished creating massive array" << std::endl;

    faces.clear();
    qfaces.clear();
    facemap.clear();
    qfacemap.clear();
    
    fid = 0;
    int idx = 0;
    std::map<std::set<int>,int> bctFace2lh;
    std::map<std::set<int>,int> bcqFace2lh;

    for(int i=1;i<=mmgMesh->ne;i++)
    {
        //adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
        
        face0.insert(mmgMesh->tetra[i].v[1]);
        face0.insert(mmgMesh->tetra[i].v[2]);
        face0.insert(mmgMesh->tetra[i].v[3]);
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
                    
            if(bfaces.find(face0)!=bfaces.end())
            {
                int refe = btfaces_Ref[face0];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->tetra[i].v[1];
                bctria[1] = mmgMesh->tetra[i].v[2];
                bctria[2] = mmgMesh->tetra[i].v[3];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face0]=lh[fid];
                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[1]);
                adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[2]);
                adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[3]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            
            fid++;
        }
        else
        {
            faces.erase(face0);
        }
        
        
        
        
        face1.insert(mmgMesh->tetra[i].v[0]);
        face1.insert(mmgMesh->tetra[i].v[2]);
        face1.insert(mmgMesh->tetra[i].v[3]);
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            if(bfaces.find(face1)!=bfaces.end())
            {
                int refe = btfaces_Ref[face1];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->tetra[i].v[0];
                bctria[1] = mmgMesh->tetra[i].v[3];
                bctria[2] = mmgMesh->tetra[i].v[2];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face1]=lh[fid];
                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[3]);
                adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[2]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            
            fid++;
        }
        else
        {
            faces.erase(face1);
        }
        
        
        
        
        
        face2.insert(mmgMesh->tetra[i].v[0]);
        face2.insert(mmgMesh->tetra[i].v[3]);
        face2.insert(mmgMesh->tetra[i].v[1]);
        if( faces.count(face2) != 1)
        {
            faces.insert(face2);
        
            if(bfaces.find(face2)!=bfaces.end())
            {
                int refe = btfaces_Ref[face2];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->tetra[i].v[0];
                bctria[1] = mmgMesh->tetra[i].v[1];
                bctria[2] = mmgMesh->tetra[i].v[3];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face2]=lh[fid];

                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[1]);
                adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[3]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            
            fid++;
        }
        else
        {
            faces.erase(face2);
        }
        
        
        
        
        

        face3.insert(mmgMesh->tetra[i].v[0]);
        face3.insert(mmgMesh->tetra[i].v[2]);
        face3.insert(mmgMesh->tetra[i].v[1]);
        
        if( faces.count(face3) != 1)
        {
            faces.insert(face3);
            
            if(bfaces.find(face3)!=bfaces.end())
            {
                int refe = btfaces_Ref[face3];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->tetra[i].v[0];
                bctria[1] = mmgMesh->tetra[i].v[2];
                bctria[2] = mmgMesh->tetra[i].v[1];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face3]=lh[fid];

                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[2]);
                adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[1]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            fid++;
        }
        else
        {
            faces.erase(face3);
        }
    
        face0.clear();
        face1.clear();
        face2.clear();
        face3.clear();
        
    }
    
    
    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
    for(int i=1;i<=mmgMesh->nprism;i++)
    {
        //adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are dealing with prisms.
  
        face0.insert(mmgMesh->prism[i].v[0]);
        face0.insert(mmgMesh->prism[i].v[2]);
        face0.insert(mmgMesh->prism[i].v[1]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
            facemap[face0]=fid;
            
            if(bfaces.find(face0)!=bfaces.end())
            {
                int refe = btfaces_Ref[face0];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->prism[i].v[0];
                bctria[1] = mmgMesh->prism[i].v[1];
                bctria[2] = mmgMesh->prism[i].v[2];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face0]=lh[fid];
                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[1]);
                adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[2]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            
            fid++;
            
        }
        else
        {
            faces.erase(face0);
        }
        
        
        
        
        
        face1.insert(mmgMesh->prism[i].v[3]);
        face1.insert(mmgMesh->prism[i].v[4]);
        face1.insert(mmgMesh->prism[i].v[5]);
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            
            facemap[face1]=fid;
            if(bfaces.find(face1)!=bfaces.end())
            {
                int refe = btfaces_Ref[face1];
                std::vector<int> bctria(3);
                bctria[0] = mmgMesh->prism[i].v[3];
                bctria[1] = mmgMesh->prism[i].v[5];
                bctria[2] = mmgMesh->prism[i].v[4];
                bctrias[refe].push_back(bctria);
                bctFace2lh[face1]=lh[fid];
                bf++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,3);
                adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[3]);
                adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[5]);
                adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[4]);
                adapt_ifn->setVal(idx,4,0);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            fid++;
        }
        else
        {
            faces.erase(face1);
        }
        
        
        
        
        
        // Quad faces //
        qface0.insert(mmgMesh->prism[i].v[0]);
        qface0.insert(mmgMesh->prism[i].v[3]);
        qface0.insert(mmgMesh->prism[i].v[4]);
        qface0.insert(mmgMesh->prism[i].v[1]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
        if( qfaces.count(qface0) != 1)
        {
            qfaces.insert(qface0);
            
            if(bqfaces.find(qface0)!=bqfaces.end())
            {
                int refe = bqfaces_Ref[qface0];
                std::vector<int> bcquad(4);
                bcquad[0] = mmgMesh->prism[i].v[0];
                bcquad[1] = mmgMesh->prism[i].v[3];
                bcquad[2] = mmgMesh->prism[i].v[4];
                bcquad[3] = mmgMesh->prism[i].v[1];
                bcquads[refe].push_back(bcquad);
                bcqFace2lh[qface0]=lh[fid];
                bq++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,4);
                adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[3]);
                adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[4]);
                adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[1]);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            
            fid++;
            
        }
        else
        {
            qfaces.erase(qface0);
        }

        
        
        
        qface1.insert(mmgMesh->prism[i].v[1]);
        qface1.insert(mmgMesh->prism[i].v[4]);
        qface1.insert(mmgMesh->prism[i].v[5]);
        qface1.insert(mmgMesh->prism[i].v[2]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface1) != 1)
        {
            qfaces.insert(qface1);
            
            if(bqfaces.find(qface1)!=bqfaces.end())
            {
                int refe = bqfaces_Ref[qface1];
                std::vector<int> bcquad(4);
                bcquad[0] = mmgMesh->prism[i].v[1];
                bcquad[1] = mmgMesh->prism[i].v[4];
                bcquad[2] = mmgMesh->prism[i].v[5];
                bcquad[3] = mmgMesh->prism[i].v[2];
                bcquads[refe].push_back(bcquad);
                bcqFace2lh[qface1]=lh[fid];
                bq++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,4);
                adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[1]);
                adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[4]);
                adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[5]);
                adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[2]);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
           
            fid++;
           
        }
        else
        {
            qfaces.erase(qface1);
        }
        
        qface2.insert(mmgMesh->prism[i].v[0]);//1
        qface2.insert(mmgMesh->prism[i].v[2]);//2
        qface2.insert(mmgMesh->prism[i].v[5]);//5
        qface2.insert(mmgMesh->prism[i].v[3]);//4
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface2) != 1)
        {
            qfaces.insert(qface2);
            qfacemap[qface2]=fid;
            if(bqfaces.find(qface2)!=bqfaces.end())
            {
                int refe = bqfaces_Ref[qface2];
                std::vector<int> bcquad(4);
                bcquad[0] = mmgMesh->prism[i].v[0];
                bcquad[1] = mmgMesh->prism[i].v[2];
                bcquad[2] = mmgMesh->prism[i].v[5];
                bcquad[3] = mmgMesh->prism[i].v[3];
                bcquads[refe].push_back(bcquad);
                bcqFace2lh[qface2]=lh[fid];
                bq++;
            }
            else
            {
                idx = new2old[fid];
                adapt_ifn->setVal(idx,0,4);
                adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[2]);
                adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[5]);
                adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[3]);
                adapt_ifn->setVal(idx,5,rh[fid]+1);
                adapt_ifn->setVal(idx,6,lh[fid]+1);
                adapt_ifn->setVal(idx,7,2);
                //idx++;
            }
            fid++;
        }
        else
        {
            qfaces.erase(qface2);
        }
        
        face0.clear();
        face1.clear();
        qface0.clear();
        qface1.clear();
        qface2.clear();
         
    }
    
    
    
    faces.clear();
    qfaces.clear();
    
    
    
    //====================================================================================
    //====================================================================================
    //MMG3D_Free_allSols(mmgMesh,&mmgSol);
    MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppSols,&mmgSol,
                   MMG5_ARG_end);
    //====================================================================================
    //====================================================================================
    
    
    
    int bc_id = 0;
    int gaa   = 0;
    int typef3 = 0;
    int typef4 = 0;
    
    std::cout << "-- Adding the interior faces to the new ifn array... face2node.size() -> " << face2node.size() << " " << lh.size() << " " << rh.size() <<std::endl;
    int ty=0;
    

    std::map<int,std::vector<int> >::iterator it_bref;
    int faceid;
    std::set<int> iface;
    int nbound = 0;
    int fa=0;
    std::cout << "-- Adding the boundary faces to the new ifn array..."<<std::endl;
    std::map<int,std::vector<std::vector<int> > >::iterator iterbc;
    int t = ftot;
    int lhi;
    for(iterbc=bctrias.begin();iterbc!=bctrias.end();iterbc++)
    {
        int bnd_id     = iterbc->first;
        int Ntris      = iterbc->second.size();
        int Nquads     = bcquads[bnd_id].size();
        std::cout << "bnd ref test = " << bnd_id << " " << Ntris << " " << Nquads << std::endl;
    }
    
    for(iterbc=bctrias.begin();iterbc!=bctrias.end();iterbc++)
    {
        int bnd_id     = iterbc->first;
        int Ntris      = iterbc->second.size();
        int Nquads     = bcquads[bnd_id].size();
                
        for(int q=0;q<Ntris;q++)
        {
            //faceid = ref2bface[bnd_id][q];
            adapt_ifn->setVal(t,0,3);
            adapt_ifn->setVal(t,1,iterbc->second[q][0]);
            adapt_ifn->setVal(t,2,iterbc->second[q][1]);
            adapt_ifn->setVal(t,3,iterbc->second[q][2]);
            adapt_ifn->setVal(t,4,0);
            iface.insert(iterbc->second[q][0]);
            iface.insert(iterbc->second[q][1]);
            iface.insert(iterbc->second[q][2]);
            lhi = bctFace2lh[iface];
            adapt_ifn->setVal(t,5,0);
            adapt_ifn->setVal(t,6,lhi+1);
            adapt_ifn->setVal(t,7,bnd_id);
            iface.clear();
            t++;
        }
        for(int q=0;q<Nquads;q++)
        {
            //faceid = ref2bqface[bnd_id][q];
            adapt_ifn->setVal(t,0,4);
            adapt_ifn->setVal(t,1,bcquads[bnd_id][q][0]);
            adapt_ifn->setVal(t,2,bcquads[bnd_id][q][1]);
            adapt_ifn->setVal(t,3,bcquads[bnd_id][q][2]);
            adapt_ifn->setVal(t,4,bcquads[bnd_id][q][3]);
            iface.insert(bcquads[bnd_id][q][0]);
            iface.insert(bcquads[bnd_id][q][1]);
            iface.insert(bcquads[bnd_id][q][2]);
            iface.insert(bcquads[bnd_id][q][3]);
            lhi = bcqFace2lh[iface];
            adapt_ifn->setVal(t,5,0);
            adapt_ifn->setVal(t,6,lhi+1);
            adapt_ifn->setVal(t,7,bnd_id);
            iface.clear();
            t++;
        }
    }
    
    
    
    int nbo = bcrefs.size();
    std::cout << "-- Constructing the zdefs array..."<<std::endl;
    Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
    std::cout << "faces.size()+qfaces.size()-bfaces.size()-bqfaces.size() " << faces.size() << " " << qfaces.size() << " " << bfaces.size() << " " << bqfaces.size() << std::endl;
    // Collect node data (10) . Starting index-ending index Nodes
    adapt_zdefs->setVal(0,0,10);
    adapt_zdefs->setVal(0,1,-1);
    adapt_zdefs->setVal(0,2,1);
    adapt_zdefs->setVal(0,3,1);
    adapt_zdefs->setVal(0,4,nVerts);
    adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
    adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
    // Collect element data (12) . Starting index-ending index Element
    adapt_zdefs->setVal(1,0,12);
    adapt_zdefs->setVal(1,1,-1);
    adapt_zdefs->setVal(1,2,2);
    adapt_zdefs->setVal(1,3,1);
    adapt_zdefs->setVal(1,4,nTet+nPrism);
    adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
    adapt_zdefs->setVal(1,6,2);
    // Collect internal face data (13) . Starting index-ending index internal face.
    adapt_zdefs->setVal(2,0,13);
    adapt_zdefs->setVal(2,1,-1);
    adapt_zdefs->setVal(2,2, 3);
    adapt_zdefs->setVal(2,3, 1);
    adapt_zdefs->setVal(2,4,lh.size()-bfaces.size()-bqfaces.size());
    adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
    adapt_zdefs->setVal(2,6,2);
    // Collect boundary face data (13) . Starting index-ending index boundary face for each boundary ID.
    int q  = 1;
    int nb = 0;
    int face_start = lh.size()-bfaces.size()-bqfaces.size()+1;
    int face_end;
    std::set<int>::iterator itr;
    for(itr=bcrefs.begin();itr!=bcrefs.end();itr++)
    {
        int bnd_ref = *itr;
        face_end = face_start+bnd_Ntri[bnd_ref]+bnd_Nquad[bnd_ref]-1;
        adapt_zdefs->setVal(3+nb,0,13);
        adapt_zdefs->setVal(3+nb,1,-1);
        adapt_zdefs->setVal(3+nb,2,3+q);
        adapt_zdefs->setVal(3+nb,3,face_start);
        adapt_zdefs->setVal(3+nb,4,face_end);
        adapt_zdefs->setVal(3+nb,5,bnd_ref);
        adapt_zdefs->setVal(3+nb,6,2);
        //std::cout << "us3d->zdefs->getVal(3+nb,5) " << us3d->zdefs->getVal(3+nb,5) << std::endl;
        face_start = face_end+1;
        //std::cout << "nb  = " << nb << " " << ref2bface.size() << " " << ref2bqface.size() << std::endl;
        nb++;
        q++;
    }
    
    //std::cout << "elements = " << " prisms->" << mmgMesh->nprism << " tetrahedra->" << mmgMesh->ne << std::endl;
    std::cout << "lh vs rh = " << " " << lh.size() << " " << rh.size() << std::endl;
    //std::cout << "-- sizingf -> " << lh.size() << " " << rh.size() << " " << t << " " << ty <<std::endl;

    
    
    //====================================================================================
    // Add ifn map to the grid.h5 file
    //====================================================================================
    
    dimsf[0] = adapt_ifn->getNrow();
    dimsf[1] = adapt_ifn->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);
    
    dset_id = H5Dcreate(file_id, "ifn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);             // hyperslab selection parameters
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_ifn->data);
    delete adapt_ifn;

    
    //====================================================================================

    hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hsize_t dimsf_att3 = 1;
    att_space = H5Screate_simple(1, &dimsf_att3, NULL);
    hid_t type3 =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type3, 10);
    ret = H5Tset_strpad(type3,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (group_info_id, "date", type3, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri3[] = "27-05-1987";
    status = H5Awrite(attr_id, type3, &stri3);
    H5Aclose(attr_id);
    
    hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    int value = nTet+nPrism;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = lh.size();
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = bfaces.size()+bqfaces.size();
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = nVerts;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    
    
    
    // Create group;
    //====================================================================================
    hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    // Add attribute to group:
    //====================================================================================
    dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = 3+nbo;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    //====================================================================================
    
    
    // Add dataset to group:
    //====================================================================================
    dimsf[0] = adapt_zdefs->getNrow();
    dimsf[1] = adapt_zdefs->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);
    hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    
    count[0]  = dimsf[0];
    count[1]  = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace  = H5Screate_simple(2, count, NULL);
    filespace = H5Dget_space(dset_zdefs_id);
    
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
    //====================================================================================
    
    dimsf_att = us3d->znames->getNrow();
    filespace = H5Screate_simple(1, &dimsf_att, NULL);
    type =  H5Tcopy (H5T_C_S1);
    ret  = H5Tset_size (type, 20);
    ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
    hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Sclose(filespace);
    
    hsize_t cnt = us3d->znames->getNrow();
    //std::cout << " us3d->znames->getNrow()  " << cnt << std::endl;
    memspace  = H5Screate_simple(1, &cnt, NULL);
    filespace = H5Dget_space(dset_znames_id);
    
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);

    PlotBoundaryData(us3d->znames,adapt_zdefs);
    
//    delete xcn_mmg;
//    delete adapt_zdefs;
//    qfacemap.clear();
//    facemap.clear();
//    faces.clear();
//    qfaces.clear();
//    face2node.clear();
//    lh.clear();
//    rh.clear();
//    Nlh.clear();
//    Nrh.clear();
//    bctrias.clear();
//    bcquads.clear();
}


void WriteUS3DGridFromMMG_it0(MMG5_pMesh mmgMesh,MMG5_pSol mmgSol, US3D* us3d)
{
    std::map<int,std::vector<int> > ref2bface;
    std::map<int,std::vector<int> > ref2bqface;
    std::set<set<int> > bfaces;
    std::set<set<int> > bqfaces;
    std::map<set<int>,int > btfaces_Ref;
    std::map<set<int>,int > bqfaces_Ref;
    std::set<int> face;
    std::set<int> bcrefs;
    int wr = 0;
    
    int nVerts = mmgMesh->np;
    int nTet = mmgMesh->ne;
    int nPrism = mmgMesh->nprism;
    
    for(int i=1;i<=mmgMesh->nt;i++)
    {
        if(mmgMesh->tria[i].ref>0 && mmgMesh->tria[i].ref!=20)// -1 is the tag for internal shell.
        {
            ref2bface[mmgMesh->tria[i].ref].push_back(i);
            face.insert(mmgMesh->tria[i].v[0]);
            face.insert(mmgMesh->tria[i].v[1]);
            face.insert(mmgMesh->tria[i].v[2]);
            
            if(btfaces_Ref.find(face)==btfaces_Ref.end())
            {
                btfaces_Ref[face] = mmgMesh->tria[i].ref;
            }
            
            bfaces.insert(face);
            if(bcrefs.find(mmgMesh->tria[i].ref)==bcrefs.end())
            {
                bcrefs.insert(mmgMesh->tria[i].ref);
            }
            face.clear();
        }
    }
    
    for(int i=1;i<=mmgMesh->nquad;i++)
    {
        if(mmgMesh->quadra[i].ref>0 && mmgMesh->quadra[i].ref!=2)// -1 is the tag for internal shell.
        {
            ref2bqface[mmgMesh->quadra[i].ref].push_back(i);
            
            face.insert(mmgMesh->quadra[i].v[0]);
            face.insert(mmgMesh->quadra[i].v[1]);
            face.insert(mmgMesh->quadra[i].v[2]);
            face.insert(mmgMesh->quadra[i].v[3]);
            
            if(bqfaces_Ref.find(face)==bqfaces_Ref.end())
            {
                bqfaces_Ref[face] = mmgMesh->quadra[i].ref;
            }
            
            
            bqfaces.insert(face);
            if(bcrefs.find(mmgMesh->quadra[i].ref)==bcrefs.end())
            {
                bcrefs.insert(mmgMesh->quadra[i].ref);
            }
            face.clear();
        }
    }
    
    
    std::map<int,int> bnd_Ntri;
    std::map<int,int> bnd_Nquad;
    int i=0;
    std::set<int>::iterator refit;
    
    for(refit=bcrefs.begin();refit!=bcrefs.end();refit++)
    {
        int ref_inq = *refit;
        
        if(ref2bface.find(ref_inq)==ref2bface.end())
        {
            bnd_Ntri[ref_inq]=0;
        }
        else
        {
            bnd_Ntri[ref_inq]=ref2bface[ref_inq].size();
        }
        
        
        if(ref2bqface.find(ref_inq)==ref2bqface.end())
        {
            bnd_Nquad[ref_inq]=0;
        }
        else
        {
            bnd_Nquad[ref_inq]=ref2bqface[ref_inq].size();
        }
    }
    
    std::map<int,std::vector<std::vector<int> > > bctrias;
    std::map<int,std::vector<std::vector<int> > > bcquads;
    
    Array<double>* xcn_mmg = new Array<double>(mmgMesh->np,3);
    for(int i=0;i<mmgMesh->np;i++)
    {
        xcn_mmg->setVal(i,0,mmgMesh->point[i+1].c[0]);
        xcn_mmg->setVal(i,1,mmgMesh->point[i+1].c[1]);
        xcn_mmg->setVal(i,2,mmgMesh->point[i+1].c[2]);
    }
    
    std::cout<<"-- Writing in HDF5 format..."<<std::endl;
    hid_t ret;
    //Output the new grid.h5 which has the new vertices and ifn map.
    //===================================================================
    //===================================================================
    //===================================================================
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    plist_id               = H5P_DEFAULT;
    //H5Pset_fapl_mpio(plist_id, comm, info);
    hid_t file_id = H5Fcreate("grid_madam.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    hid_t status;
    hid_t att_space;
    hid_t attr_id;
    
    hsize_t dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    hid_t type =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type, 14);
    ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri[] = "US3D Grid File";
    status = H5Awrite(attr_id, type, &stri);
    H5Aclose(attr_id);
    
    hsize_t dimsf_att2 = 1;
    att_space = H5Screate_simple(1, &dimsf_att2, NULL);
    hid_t type2 =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type2, 5);
    ret = H5Tset_strpad(type2,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (file_id, "filevers", type2, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri2[] = "1.1.8";
    status = H5Awrite(attr_id, type2, &stri2);
    H5Aclose(attr_id);
    
    //====================================================================================
    // Add xcn map to the grid.h5 file
    //====================================================================================
    hsize_t     dimsf[2];
    hsize_t    count[2];              // hyperslab selection parameters
    hsize_t    offset[2];
    dimsf[0] = xcn_mmg->getNrow();
    dimsf[1] = xcn_mmg->getNcol();
    hid_t filespace = H5Screate_simple(2, dimsf, NULL);

    hid_t dset_id = H5Dcreate(file_id, "xcn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    hid_t memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_mmg->data);
    delete xcn_mmg;
    //====================================================================================
    
    std::map<std::set<int>, int> qfacemap;
    std::map<std::set<int>, int>  facemap;
    std::set<std::set<int> > faces;
    std::set<std::set<int> > qfaces;
    std::set<int> face0;
    std::set<int> face1;
    std::set<int> face2;
    std::set<int> face3;
    std::set<int> qface0;
    std::set<int> qface1;
    std::set<int> qface2;
    int fid = 0;
    int vid0,vid1,vid2,vid3;
    std::map<int,int> lh;
    std::map<int,int> rh;
    
    std::map<int,int> Nlh;
    std::map<int,int> Nrh;
    
    int of = 0;
    int fset_cnt = 0;
    Array<int>* adapt_iet = new Array<int>(nTet+nPrism,1);
    // local face2vert_map for a tet in mmg  {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
    int bf = 0;
    int bq = 0;
    std::cout << "-- Constructing the new face-2-node and face-2-element map..."<<std::endl;

    for(int i=1;i<=mmgMesh->ne;i++)
    {
        adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
        
        face0.insert(mmgMesh->tetra[i].v[1]);
        face0.insert(mmgMesh->tetra[i].v[2]);
        face0.insert(mmgMesh->tetra[i].v[3]);
        
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
            facemap[face0]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face0]] = i-1;
            faces.erase(face0);
            facemap.erase(face0);
        }
        
        face1.insert(mmgMesh->tetra[i].v[0]);
        face1.insert(mmgMesh->tetra[i].v[2]);
        face1.insert(mmgMesh->tetra[i].v[3]);
        
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            facemap[face1]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face1]] = i-1;
            faces.erase(face1);
            facemap.erase(face1);
        }
        
        face2.insert(mmgMesh->tetra[i].v[0]);
        face2.insert(mmgMesh->tetra[i].v[3]);
        face2.insert(mmgMesh->tetra[i].v[1]);
        
        
        
        
        if( faces.count(face2) != 1)
        {
            faces.insert(face2);
            facemap[face2]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face2]] = i-1;
            faces.erase(face2);
            facemap.erase(face2);
        }

        face3.insert(mmgMesh->tetra[i].v[0]);
        face3.insert(mmgMesh->tetra[i].v[2]);
        face3.insert(mmgMesh->tetra[i].v[1]);
        
        
        if( faces.count(face3) != 1)
        {
            faces.insert(face3);
            facemap[face3]=fid;
            lh[fid] = i-1;
            fid++;
        }
        else
        {
            rh[facemap[face3]] = i-1;
            faces.erase(face3);
            facemap.erase(face3);
        }
    
        face0.clear();
        face1.clear();
        face2.clear();
        face3.clear();
        
    }
    
    

    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };


    for(int i=1;i<=mmgMesh->nprism;i++)
    {
        adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are dealing with prisms.
               // std::cout  << "Prism ["<<i<<"]=" << mmgMesh->prism[i].v[0] << " " << mmgMesh->prism[i].v[1] << " " << mmgMesh->prism[i].v[2] << " " << mmgMesh->prism[i].v[3] << " " << mmgMesh->prism[i].v[4] << " " << mmgMesh->prism[i].v[5] << std::endl;
  
        face0.insert(mmgMesh->prism[i].v[0]);
        face0.insert(mmgMesh->prism[i].v[2]);
        face0.insert(mmgMesh->prism[i].v[1]);
        
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
    
        if(faces.count(face0) != 1 )
        {
            faces.insert(face0);
            facemap[face0]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[facemap[face0]] = mmgMesh->ne+i-1;
            faces.erase(face0);
            facemap.erase(face0);
        }
        
        
        face1.insert(mmgMesh->prism[i].v[3]);
        face1.insert(mmgMesh->prism[i].v[4]);
        face1.insert(mmgMesh->prism[i].v[5]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
        
        
        if(faces.count(face1) != 1)
        {
            faces.insert(face1);
            facemap[face1]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[facemap[face1]]  = mmgMesh->ne+i-1;
            faces.erase(face1);
            facemap.erase(face1);
        }
        
        // Quad faces //
        qface0.insert(mmgMesh->prism[i].v[0]);
        qface0.insert(mmgMesh->prism[i].v[2]);
        qface0.insert(mmgMesh->prism[i].v[4]);
        qface0.insert(mmgMesh->prism[i].v[3]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface0) != 1)
        {
            qfaces.insert(qface0);
            qfacemap[qface0]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {
            rh[qfacemap[qface0]]  = mmgMesh->ne+i-1;
            qfaces.erase(qface0);
            qfacemap.erase(qface0);
        }

        qface1.insert(mmgMesh->prism[i].v[1]);
        qface1.insert(mmgMesh->prism[i].v[5]);
        qface1.insert(mmgMesh->prism[i].v[4]);
        qface1.insert(mmgMesh->prism[i].v[2]);
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface1) != 1)
        {
            qfaces.insert(qface1);
            qfacemap[qface1]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {

            rh[qfacemap[qface1]]  = mmgMesh->ne+i-1;
            qfaces.erase(qface1);
            qfacemap.erase(qface1);
        }
        
        qface2.insert(mmgMesh->prism[i].v[0]);//1
        qface2.insert(mmgMesh->prism[i].v[3]);//2
        qface2.insert(mmgMesh->prism[i].v[5]);//5
        qface2.insert(mmgMesh->prism[i].v[1]);//4
        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

        if( qfaces.count(qface2) != 1)
        {
            qfaces.insert(qface2);
            qfacemap[qface2]=fid;
            lh[fid] = mmgMesh->ne+i-1;
            fid++;
        }
        else
        {

            rh[qfacemap[qface2]] = mmgMesh->ne+i-1;
            qfaces.erase(qface2);
            qfacemap.erase(qface2);
        }
        
        face0.clear();
        face1.clear();
        qface0.clear();
        qface1.clear();
        qface2.clear();
         
    }
    
    //====================================================================================
    //====================================================================================
    dimsf[0] = adapt_iet->getNrow();
    dimsf[1] = adapt_iet->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);

    dset_id = H5Dcreate(file_id, "iet", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_iet->data);
    delete adapt_iet;
    //====================================================================================
    //====================================================================================

    std::cout << "AND WE ARE DONE COUNTING " << std::endl;
    
    int ftot = 0;
    std::map<int,int> new2old;
//
    std::map<int,int>::iterator itm;
    int it;
    for(itm=lh.begin();itm!=lh.end();itm++)
    {
        it = itm->first;
        if(rh.find(it)!=rh.end())
        {
            new2old[it] = ftot;
            ftot++;
        }
    }
    
    std::cout << "AND WE ARE DONE DETERMINING THE ORDER" << std::endl;
    std::cout << "Starting creating massive array" << std::endl;

    Array<int>* adapt_ifn = new Array<int>(fid,8);
    std::cout << "Finished creating massive array" << std::endl;

    faces.clear();
    qfaces.clear();
    facemap.clear();
    qfacemap.clear();
    
    fid = 0;
    int idx = 0;
    std::map<std::set<int>,int> bctFace2lh;
    std::map<std::set<int>,int> bcqFace2lh;
     for(int i=1;i<=mmgMesh->ne;i++)
     {
         //adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
         
         face0.insert(mmgMesh->tetra[i].v[1]);
         face0.insert(mmgMesh->tetra[i].v[2]);
         face0.insert(mmgMesh->tetra[i].v[3]);
         
         if(faces.count(face0) != 1 )
         {
             faces.insert(face0);
             //facemap[face0]=fid;

             if(bfaces.find(face0)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face0];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->tetra[i].v[1];
                 bctria[1] = mmgMesh->tetra[i].v[2];
                 bctria[2] = mmgMesh->tetra[i].v[3];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face0]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[1]);
                 adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[2]);
                 adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[3]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             fid++;
         }
         else
         {
             faces.erase(face0);
         }
         
         face1.insert(mmgMesh->tetra[i].v[0]);
         face1.insert(mmgMesh->tetra[i].v[2]);
         face1.insert(mmgMesh->tetra[i].v[3]);
         
         if(faces.count(face1) != 1)
         {
             faces.insert(face1);
             //facemap[face1]=fid;

             if(bfaces.find(face1)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face1];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->tetra[i].v[0];
                 bctria[1] = mmgMesh->tetra[i].v[3];
                 bctria[2] = mmgMesh->tetra[i].v[2];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face1]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[3]);
                 adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[2]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             
             fid++;
         }
         else
         {
             faces.erase(face1);
         }
         
         
         face2.insert(mmgMesh->tetra[i].v[0]);
         face2.insert(mmgMesh->tetra[i].v[3]);
         face2.insert(mmgMesh->tetra[i].v[1]);
         
         
         
         
         if( faces.count(face2) != 1)
         {
             faces.insert(face2);
             //facemap[face2]=fid;

             if(bfaces.find(face2)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face2];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->tetra[i].v[0];
                 bctria[1] = mmgMesh->tetra[i].v[1];
                 bctria[2] = mmgMesh->tetra[i].v[3];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face2]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[1]);
                 adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[3]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
                         
             fid++;
         }
         else
         {
             faces.erase(face2);
         }
         
         face3.insert(mmgMesh->tetra[i].v[0]);
         face3.insert(mmgMesh->tetra[i].v[2]);
         face3.insert(mmgMesh->tetra[i].v[1]);
         
         
         if( faces.count(face3) != 1)
         {
             faces.insert(face3);
             
             //facemap[face3]=fid;
             
             if(bfaces.find(face3)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face3];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->tetra[i].v[0];
                 bctria[1] = mmgMesh->tetra[i].v[2];
                 bctria[2] = mmgMesh->tetra[i].v[1];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face3]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->tetra[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->tetra[i].v[2]);
                 adapt_ifn->setVal(idx,3,mmgMesh->tetra[i].v[1]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             fid++;
         }
         else
         {
             faces.erase(face3);
         }
         
         
         face0.clear();
         face1.clear();
         face2.clear();
         face3.clear();
         
     }
     
     

     for(int i=1;i<=mmgMesh->nprism;i++)
     {
         //adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are dealing with prisms.
                // std::cout  << "Prism ["<<i<<"]=" << mmgMesh->prism[i].v[0] << " " << mmgMesh->prism[i].v[1] << " " << mmgMesh->prism[i].v[2] << " " << mmgMesh->prism[i].v[3] << " " << mmgMesh->prism[i].v[4] << " " << mmgMesh->prism[i].v[5] << std::endl;
   
         face0.insert(mmgMesh->prism[i].v[0]);
         face0.insert(mmgMesh->prism[i].v[2]);
         face0.insert(mmgMesh->prism[i].v[1]);
         
         // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
     
         if(faces.count(face0) != 1 )
         {
             faces.insert(face0);
             //facemap[face0]=fid;
             if(bfaces.find(face0)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face0];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->prism[i].v[0];
                 bctria[1] = mmgMesh->prism[i].v[1];
                 bctria[2] = mmgMesh->prism[i].v[2];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face0]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[1]);
                 adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[2]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             
             fid++;

         }
         else
         {
             faces.erase(face0);
         }
         
         face1.insert(mmgMesh->prism[i].v[3]);
         face1.insert(mmgMesh->prism[i].v[4]);
         face1.insert(mmgMesh->prism[i].v[5]);
         // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
         
         
         if(faces.count(face1) != 1)
         {
             faces.insert(face1);
             //facemap[face1]=fid;
             if(bfaces.find(face1)!=bfaces.end())
             {
                 int refe = btfaces_Ref[face1];
                 std::vector<int> bctria(3);
                 bctria[0] = mmgMesh->prism[i].v[3];
                 bctria[1] = mmgMesh->prism[i].v[5];
                 bctria[2] = mmgMesh->prism[i].v[4];
                 bctrias[refe].push_back(bctria);
                 bctFace2lh[face1]=lh[fid];
                 bf++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,3);
                 adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[3]);
                 adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[4]);
                 adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[5]);
                 adapt_ifn->setVal(idx,4,0);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             
             fid++;
         }
         else
         {
             faces.erase(face1);
         }
         
         // Quad faces //
         qface0.insert(mmgMesh->prism[i].v[0]);
         qface0.insert(mmgMesh->prism[i].v[2]);
         qface0.insert(mmgMesh->prism[i].v[4]);
         qface0.insert(mmgMesh->prism[i].v[3]);
         // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

         if( qfaces.count(qface0) != 1)
         {
             qfaces.insert(qface0);
             //qfacemap[qface0]=fid;
             if(bqfaces.find(qface0)!=bqfaces.end())
             {
                 int refe = bqfaces_Ref[qface0];
                 std::vector<int> bcquad(4);
                 bcquad[0] = mmgMesh->prism[i].v[0];
                 bcquad[1] = mmgMesh->prism[i].v[2];
                 bcquad[2] = mmgMesh->prism[i].v[4];
                 bcquad[3] = mmgMesh->prism[i].v[3];
                 bcquads[refe].push_back(bcquad);
                 bcqFace2lh[qface0]=lh[fid];
                 bq++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,4);
                 adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[2]);
                 adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[4]);
                 adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[3]);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             
             fid++;
         }
         else
         {
             qfaces.erase(qface0);
         }

         qface1.insert(mmgMesh->prism[i].v[1]);
         qface1.insert(mmgMesh->prism[i].v[5]);
         qface1.insert(mmgMesh->prism[i].v[4]);
         qface1.insert(mmgMesh->prism[i].v[2]);
         // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

         if( qfaces.count(qface1) != 1)
         {
             qfaces.insert(qface1);
             //qfacemap[qface1]=fid;
             if(bqfaces.find(qface1)!=bqfaces.end())
             {
                 int refe = bqfaces_Ref[qface1];
                 std::vector<int> bcquad(4);
                 bcquad[0] = mmgMesh->prism[i].v[1];
                 bcquad[1] = mmgMesh->prism[i].v[5];
                 bcquad[2] = mmgMesh->prism[i].v[4];
                 bcquad[3] = mmgMesh->prism[i].v[2];
                 bcquads[refe].push_back(bcquad);
                 bcqFace2lh[qface1]=lh[fid];
                 bq++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,4);
                 adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[1]);
                 adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[5]);
                 adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[4]);
                 adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[2]);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             fid++;

         }
         else
         {
             qfaces.erase(qface1);
         }
         
         
         qface2.insert(mmgMesh->prism[i].v[0]);//1
         qface2.insert(mmgMesh->prism[i].v[3]);//2
         qface2.insert(mmgMesh->prism[i].v[5]);//5
         qface2.insert(mmgMesh->prism[i].v[1]);//4
         // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

         if( qfaces.count(qface2) != 1)
         {
             qfaces.insert(qface2);
             //qfacemap[qface2]=fid;
             if(bqfaces.find(qface2)!=bqfaces.end())
             {
                 int refe = bqfaces_Ref[qface2];
                 std::vector<int> bcquad(4);
                 bcquad[0] = mmgMesh->prism[i].v[0];
                 bcquad[1] = mmgMesh->prism[i].v[3];
                 bcquad[2] = mmgMesh->prism[i].v[5];
                 bcquad[3] = mmgMesh->prism[i].v[1];
                 bcquads[refe].push_back(bcquad);
                 bcqFace2lh[qface2]=lh[fid];
                 bq++;
             }
             else
             {
                 idx = new2old[fid];
                 adapt_ifn->setVal(idx,0,4);
                 adapt_ifn->setVal(idx,1,mmgMesh->prism[i].v[0]);
                 adapt_ifn->setVal(idx,2,mmgMesh->prism[i].v[3]);
                 adapt_ifn->setVal(idx,3,mmgMesh->prism[i].v[5]);
                 adapt_ifn->setVal(idx,4,mmgMesh->prism[i].v[1]);
                 adapt_ifn->setVal(idx,5,rh[fid]+1);
                 adapt_ifn->setVal(idx,6,lh[fid]+1);
                 adapt_ifn->setVal(idx,7,2);
                 //idx++;
             }
             
             fid++;
         }
         else
         {
             qfaces.erase(qface2);
         }
         
         face0.clear();
         face1.clear();
         qface0.clear();
         qface1.clear();
         qface2.clear();
          
     }
     
     
     
     
    faces.clear();
    qfaces.clear();
    
    
    //====================================================================================
    //====================================================================================
    //MMG3D_Free_allSols(mmgMesh,&mmgSol);
    MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppSols,&mmgSol,
                   MMG5_ARG_end);
    //====================================================================================
    //====================================================================================
    Array<int>* zdefs = us3d->zdefs;
    int bc_id = 0;
    int gaa   = 0;
    int typef3 = 0;
    int typef4 = 0;
    
    std::cout << "-- Adding the interior faces to the new ifn array... face2node.size() -> " << " " << lh.size() << " " << rh.size() <<std::endl;
    int ty=0;
    int counte;

    

    std::map<int,std::vector<int> >::iterator it_bref;
    int faceid;
    std::set<int> iface;
    int nbound = 0;
    int fa=0;
    std::cout << "-- Adding the boundary faces to the new ifn array..."<<std::endl;
    std::map<int,std::vector<std::vector<int> > >::iterator iterbc;
    int t = ftot;
    for(iterbc=bctrias.begin();iterbc!=bctrias.end();iterbc++)
    {
        int bnd_id     = iterbc->first;
        int Ntris      = iterbc->second.size();
        int Nquads     = bcquads[bnd_id].size();
        std::cout << "bnd ref test = " << bnd_id << " " << Ntris << " " << Nquads << std::endl;
    }
    int lhi;
    for(iterbc=bctrias.begin();iterbc!=bctrias.end();iterbc++)
    {
        int bnd_id     = iterbc->first;
        int Ntris      = iterbc->second.size();
        int Nquads     = bcquads[bnd_id].size();
        
        //std::cout << "bnd ref = " << bnd_id << " " << Ntris << " " << Nquads << std::endl;
        
        for(int q=0;q<Ntris;q++)
        {
            //faceid = ref2bface[bnd_id][q];
            adapt_ifn->setVal(t,0,3);
            adapt_ifn->setVal(t,1,iterbc->second[q][0]);
            adapt_ifn->setVal(t,2,iterbc->second[q][1]);
            adapt_ifn->setVal(t,3,iterbc->second[q][2]);
            adapt_ifn->setVal(t,4,0);
            iface.insert(iterbc->second[q][0]);
            iface.insert(iterbc->second[q][1]);
            iface.insert(iterbc->second[q][2]);
            
            //std::cout << "3 bc row = " << t << " " << mmgMesh->tria[faceid].v[0] << " " << mmgMesh->tria[faceid].v[1] << " " << mmgMesh->tria[faceid].v[2] << std::endl;

            //fid=facemap[iface];
            lhi = bctFace2lh[iface];
            adapt_ifn->setVal(t,5,0);
            adapt_ifn->setVal(t,6,lhi+1);
            //adapt_ifn->setVal(t,7,us3d->zdefs->getVal(3+bnd_id-1,5));
            adapt_ifn->setVal(t,7,bnd_id);

            iface.clear();
            t++;
        }
        for(int q=0;q<Nquads;q++)
        {
            //faceid = ref2bqface[bnd_id][q];
            adapt_ifn->setVal(t,0,4);
            adapt_ifn->setVal(t,1,bcquads[bnd_id][q][0]);
            adapt_ifn->setVal(t,2,bcquads[bnd_id][q][1]);
            adapt_ifn->setVal(t,3,bcquads[bnd_id][q][2]);
            adapt_ifn->setVal(t,4,bcquads[bnd_id][q][3]);
            iface.insert(bcquads[bnd_id][q][0]);
            iface.insert(bcquads[bnd_id][q][1]);
            iface.insert(bcquads[bnd_id][q][2]);
            iface.insert(bcquads[bnd_id][q][3]);
                        
            //fid=qfacemap[iface];
            lhi = bcqFace2lh[iface];
            adapt_ifn->setVal(t,5,0);
            adapt_ifn->setVal(t,6,lhi+1);
            adapt_ifn->setVal(t,7,bnd_id);
            iface.clear();
            
            t++;
        }
    }
    
    //std::cout << "-- sizing2 -> " << lh.size() << " " << rh.size() << " " << t << " " << ty <<std::endl;
    

    int nbo = bcrefs.size();
    std::cout << "-- Constructing the zdefs array..."<<std::endl;
    Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
    std::cout << "faces.size()+qfaces.size()-bfaces.size()-bqfaces.size() " << faces.size() << " " << qfaces.size() << " " << bfaces.size() << " " << bqfaces.size() << std::endl;
    // Collect node data (10) . Starting index-ending index Nodes
    adapt_zdefs->setVal(0,0,10);
    adapt_zdefs->setVal(0,1,-1);
    adapt_zdefs->setVal(0,2,1);
    adapt_zdefs->setVal(0,3,1);
    adapt_zdefs->setVal(0,4,nVerts);
    adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
    adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
    // Collect element data (12) . Starting index-ending index Element
    adapt_zdefs->setVal(1,0,12);
    adapt_zdefs->setVal(1,1,-1);
    adapt_zdefs->setVal(1,2,2);
    adapt_zdefs->setVal(1,3,1);
    adapt_zdefs->setVal(1,4,nTet+nPrism);
    adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
    adapt_zdefs->setVal(1,6,2);
    // Collect internal face data (13) . Starting index-ending index internal face.
    adapt_zdefs->setVal(2,0,13);
    adapt_zdefs->setVal(2,1,-1);
    adapt_zdefs->setVal(2,2, 3);
    adapt_zdefs->setVal(2,3, 1);
    adapt_zdefs->setVal(2,4,lh.size()-bfaces.size()-bqfaces.size());
    adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
    adapt_zdefs->setVal(2,6,2);
    // Collect boundary face data (13) . Starting index-ending index boundary face for each boundary ID.
    int q  = 1;
    int nb = 0;
    int face_start = lh.size()-bfaces.size()-bqfaces.size()+1;
    int face_end;
    std::set<int>::iterator itr;
    for(itr=bcrefs.begin();itr!=bcrefs.end();itr++)
    {
        int bnd_ref = *itr;
        face_end = face_start+bnd_Ntri[bnd_ref]+bnd_Nquad[bnd_ref]-1;
        adapt_zdefs->setVal(3+nb,0,13);
        adapt_zdefs->setVal(3+nb,1,-1);
        adapt_zdefs->setVal(3+nb,2,3+q);
        adapt_zdefs->setVal(3+nb,3,face_start);
        adapt_zdefs->setVal(3+nb,4,face_end);
        adapt_zdefs->setVal(3+nb,5,bnd_ref);
        adapt_zdefs->setVal(3+nb,6,2);
        //std::cout << "us3d->zdefs->getVal(3+nb,5) " << us3d->zdefs->getVal(3+nb,5) << std::endl;
        face_start = face_end+1;
        //std::cout << "nb  = " << nb << " " << ref2bface.size() << " " << ref2bqface.size() << std::endl;
        nb++;
        q++;
    }
    
    //std::cout << "elements = " << " " << mmgMesh->nprism << " " << mmgMesh->ne << std::endl;
    std::cout << "lh vs rh = " << " " << lh.size() << " " << rh.size() << std::endl;
    //std::cout << "-- sizingf -> " << lh.size() << " " << rh.size() << " " << t << " " << ty <<std::endl;

    
    
    //====================================================================================
    // Add ifn map to the grid.h5 file
    //====================================================================================
    
    dimsf[0] = adapt_ifn->getNrow();
    dimsf[1] = adapt_ifn->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);
    
    dset_id = H5Dcreate(file_id, "ifn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);             // hyperslab selection parameters
    count[0] = dimsf[0];
    count[1] = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_ifn->data);
    
    delete adapt_ifn;
    
    //====================================================================================

    hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hsize_t dimsf_att3 = 1;
    att_space = H5Screate_simple(1, &dimsf_att3, NULL);
    hid_t type3 =  H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type3, 10);
    ret = H5Tset_strpad(type3,H5T_STR_SPACEPAD);
    attr_id   = H5Acreate (group_info_id, "date", type3, att_space, H5P_DEFAULT, H5P_DEFAULT);
    char stri3[] = "27-05-1987";
    status = H5Awrite(attr_id, type3, &stri3);
    H5Aclose(attr_id);
    
    hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    int value = nTet+nPrism;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = lh.size();
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = bfaces.size()+bqfaces.size();
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = nVerts;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    
    
    
    
    // Create group;
    //====================================================================================
    hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    // Add attribute to group:
    //====================================================================================
    dimsf_att = 1;
    att_space = H5Screate_simple(1, &dimsf_att, NULL);
    attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
    value = 3+nbo;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
    H5Aclose(attr_id);
    //====================================================================================
    
    
    // Add dataset to group:
    //====================================================================================
    dimsf[0] = adapt_zdefs->getNrow();
    dimsf[1] = adapt_zdefs->getNcol();
    filespace = H5Screate_simple(2, dimsf, NULL);
    hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    
    count[0]  = dimsf[0];
    count[1]  = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace  = H5Screate_simple(2, count, NULL);
    filespace = H5Dget_space(dset_zdefs_id);
    
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
    //====================================================================================
    
    dimsf_att = us3d->znames->getNrow();
    filespace = H5Screate_simple(1, &dimsf_att, NULL);
    type =  H5Tcopy (H5T_C_S1);
    ret  = H5Tset_size (type, 20);
    ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
    hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Sclose(filespace);
    
    hsize_t cnt = us3d->znames->getNrow();
    memspace  = H5Screate_simple(1, &cnt, NULL);
    filespace = H5Dget_space(dset_znames_id);
    
    
    status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);

    PlotBoundaryData(us3d->znames,adapt_zdefs);
    
    delete adapt_zdefs;
    qfacemap.clear();
    facemap.clear();
    faces.clear();
    qfaces.clear();
    lh.clear();
    rh.clear();
    Nlh.clear();
    Nrh.clear();
    bctrias.clear();
    bcquads.clear();
    
}







//US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, MPI_Comm comm, MPI_Info info)
//{
//    int size;
//    MPI_Comm_size(comm, &size);
//    // Get the rank of the process
//    int rank;
//    MPI_Comm_rank(comm, &rank);
//    US3D* us3d = new US3D;
//    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
//
//    ParArray<int>* ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
//    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
//    ParArray<int>* iee = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);
//
//    Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
//    Array<int>* ife = ReadDataSetFromFile<int>(fn_conn,"ife");
//
//    int Nel = ien->getNglob();
//
//    ParArray<double>* interior  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
//    Array<double>* ghost        = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);
//
//    Array<int>*    zdefs        = ReadDataSetFromGroupFromFile<int>(fn_grid,"zones","zdefs");
//    Array<char>*  znames        = ReadDataSetFromGroupFromFile<char>(fn_grid,"zones","znames");
//    std::map<int,std::vector<int> > bnd_face_map;
//    // Collect boundary data;
//    std::vector<int> bnd_m;
//    int t=0;
//    for(int i=4;i<zdefs->getNrow();i++)
//    {
//        bnd_m.push_back(zdefs->getVal(i,3));
//    }
//    bnd_m.push_back(zdefs->getVal(zdefs->getNrow()-1,4));
//
//    if(rank == 0)
//    {
//       //std::cout << "Rank = " << rank << std::endl;
//       PlotBoundaryData(znames,zdefs);
//    }
//
//    int nBnd = zdefs->getNrow()-3;
//    int* bnd_map = new int[zdefs->getNrow()-3];
//    std::map<int,char*> znames_map;
//    Array<char>* znames_new = new Array<char>(znames->getNrow(),znames->getNcol());
//    for(int i=0;i<zdefs->getNrow();i++)
//    {
//        //bnd_map[i-3] = zdefs->getVal(i,3)-1;
//
//        if(zdefs->getVal(i,5)!=1)
//        {
//            char* name = new char[znames->getNcol()];
//
//            for(int j=0;j<znames->getNcol();j++)
//            {
//               name[j]=znames->getVal(i,j);
//            }
//            znames_map[zdefs->getVal(i,5)] = name;
//        }
//    }
//
//    // number of vertices
//    for(int j=0;j<znames->getNcol();j++)
//    {
//       znames_new->setVal(0,j,znames->getVal(0,j));
//    }
//    std::cout << std::endl;
//    // number of cells
//    for(int j=0;j<znames->getNcol();j++)
//    {
//       znames_new->setVal(1,j,znames->getVal(1,j));
//    }
//
//    std::map<int,char*>::iterator itch;
//    int c=2;
//    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
//    {
//        int bid = itch->first;
//        for(int j=0;j<znames->getNcol();j++)
//        {
//            znames_new->setVal(c,j,znames_map[bid][j]);
//        }
//        c++;
//    }
//
//
//
//    int i,j;
//    int nglob = ien->getNglob();
//    int nrow  = ien->getNrow();
//    int ncol  = ien->getNcol()-1;
//    //
//    ParArray<int>* ien_copy = new ParArray<int>(nglob,ncol,comm);
//    //
//    for(i=0;i<nrow;i++)
//    {
//        for(int j=0;j<ncol;j++)
//        {
//            ien_copy->setVal(i,j,ien->getVal(i,j+1)-1);
//        }
//    }
//    delete ien;
//    //
//    int ncol_ief = ief->getNcol()-1;
//    ParArray<int>* ief_copy = new ParArray<int>(nglob,ncol_ief,comm);
//
//    for(i=0;i<nrow;i++)
//    {
//        for(j=0;j<ncol_ief;j++)
//        {
//            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
//        }
//    }
//    delete ief;
//
//    int ncol_iee = iee->getNcol()-1;
//    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
//
//    for(i=0;i<nrow;i++)
//    {
//        for(j=0;j<ncol_iee;j++)
//        {
//            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
//        }
//    }
//    delete iee;
//
//    int nrow_ifn = ifn->getNrow();
//
//    int ncol_ifn = 4;
//    Array<int>* ifn_copy = new Array<int>(nrow_ifn,ncol_ifn);
//    Array<int>* ifn_ref  = new Array<int>(nrow_ifn,1);
//    int ref;
//    //std::ofstream myfile20;
//    //myfile20.open("ifn_ref.dat");
//    std::map<std::set<int>,int> tria_ref_map;
//    std::map<std::set<int>,int> quad_ref_map;
//    std::map<int,int> vert_ref_map;
//    std::set<int> vert_ref_set;
//
//    std::set<int> tria0;
//    std::set<int> tria00;
//    std::set<int> tria1;
//    std::set<int> tria11;
//
//    std::set<int> quad;
//    int faceid;
//    int nodeid;
//    for(i=0;i<nrow_ifn;i++)
//    {
//        ref = ifn->getVal(i,7);
//        ifn_ref->setVal(i,0,ref);
//        faceid = i;
//        if(ref != 2)
//        {
//            bnd_face_map[ref].push_back(faceid);
//        }
//
//        for(j=0;j<ncol_ifn;j++)
//        {
//            nodeid = ifn->getVal(i,j+1)-1; // This is actually node ID!!!!
//            ifn_copy->setVal(i,j,ifn->getVal(i,j+1)-1);
//
//            if(ref!=2)
//            {
//                if(vert_ref_set.find(nodeid)==vert_ref_set.end())
//                {
//                    vert_ref_set.insert(nodeid);
//                    vert_ref_map[nodeid] = ifn->getVal(i,7);
//                }
//            }
////            else
////            {
////                vert_ref_map[nodeid] = ifn->getVal(i,7);
////            }
//
////            if(vert_ref_set.find(nodeid)==vert_ref_set.end())
////            {
////                vert_ref_set.insert(nodeid);
////                vert_ref_map[nodeid] = ifn_ref->getVal(i,0);
////
//////                if(ifn_ref->getVal(i,0)!=0)
//////                {
//////                    std::cout << "ref not zero " << ifn_ref->getVal(i,0) << std::endl;
//////                }
////            }
//        }
//
//
//        tria0.insert(ifn->getVal(i,0+1)-1);
//        tria0.insert(ifn->getVal(i,1+1)-1);
//        tria0.insert(ifn->getVal(i,2+1)-1);
//        tria00.insert(ifn->getVal(i,0+1)-1);
//        tria00.insert(ifn->getVal(i,2+1)-1);
//        tria00.insert(ifn->getVal(i,3+1)-1);
//
//        tria1.insert(ifn->getVal(i,0+1)-1);
//        tria1.insert(ifn->getVal(i,1+1)-1);
//        tria1.insert(ifn->getVal(i,3+1)-1);
//        tria11.insert(ifn->getVal(i,1+1)-1);
//        tria11.insert(ifn->getVal(i,2+1)-1);
//        tria11.insert(ifn->getVal(i,3+1)-1);
//
//        quad.insert(ifn->getVal(i,0+1)-1);
//        quad.insert(ifn->getVal(i,1+1)-1);
//        quad.insert(ifn->getVal(i,2+1)-1);
//        quad.insert(ifn->getVal(i,3+1)-1);
//
//        if(tria_ref_map.find(tria0)==tria_ref_map.end() && ref!=2)
//        {
//            tria_ref_map[tria0] = ref;
//            tria_ref_map[tria00] = ref;
//        }
//        if(tria_ref_map.find(tria1)==tria_ref_map.end() && ref!=2)
//        {
//            tria_ref_map[tria1] = ref;
//            tria_ref_map[tria11] = ref;
//        }
//
//        if(quad_ref_map.find(quad)==quad_ref_map.end() && ref!=2)
//        {
//            quad_ref_map[quad] = ref;
//        }
//
//        tria0.clear();
//        tria00.clear();
//        tria1.clear();
//        tria11.clear();
//        quad.clear();
//    }
//
////    std::map<int,int>::iterator itje;
////    for(itje=vert_ref_map.begin();itje!=vert_ref_map.end();itje++)
////    {
////        if(itje->second!=0)
////        {
////            std::cout << itje->first << " " << itje->second << std::endl;
////
////        }
////    }
//
//    //std::cout << " vert_ref_map = " << vert_ref_map.size() << std::endl;
//
//    //myfile20.close();
//    if(rank == 0)
//    {
//        std::map<int,std::vector<int> >::iterator its;
//        for(its=bnd_face_map.begin();its!=bnd_face_map.end();its++)
//        {
//            std::cout << "bnd_face_map " << its->first << " " << its->second.size() << std::endl;
//        }
//    }
//
//    delete ifn;
//    int nrow_ife = ife->getNrow();
//    int ncol_ife = 2;
//    Array<int>* ife_copy = new Array<int>(nrow_ife,ncol_ife);
//    for(i=0;i<nrow_ife;i++)
//    {
//        for(j=0;j<ncol_ife;j++)
//        {
//            ife_copy->setVal(i,j,ife->getVal(i,j)-1);
//        }
//    }
//    delete ife;
//
//    us3d->xcn = xcn;
//
//    us3d->ien           = ien_copy;
//    us3d->ief           = ief_copy;
//    us3d->iee           = iee_copy;
//
//    us3d->ifn           = ifn_copy;
//    us3d->ifn_ref       = ifn_copy;
//    us3d->ife           = ife_copy;
//
//    us3d->interior      = interior;
//    us3d->ghost         = ghost;
//
//    us3d->znames        = znames_new;
//    us3d->zdefs         = zdefs;
//    us3d->bnd_m         = bnd_m;
//    us3d->bnd_map       = bnd_map;
//    us3d->bnd_face_map  = bnd_face_map;
//    us3d->nBnd          = nBnd;
//
//    us3d->tria_ref_map  = tria_ref_map;
//    us3d->quad_ref_map  = quad_ref_map;
//    us3d->vert_ref_map  = vert_ref_map;
//
//    return us3d;
//}


US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, int readFromStats, MPI_Comm comm, MPI_Info info)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    US3D* us3d = new US3D;
    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
    //std::cout << "Reading from :: " << fn_conn << std::endl;
    ParArray<int>* ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    ParArray<int>* iee = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);
    ParArray<int>* iet = ReadDataSetFromFileInParallel<int>(fn_grid,"iet",comm,info);
    ParArray<int>* ifn = ReadDataSetFromFileInParallel<int>(fn_grid,"ifn",comm,info);
    ParArray<int>* ife = ReadDataSetFromFileInParallel<int>(fn_conn,"ife",comm,info);

    
    int Nel = ien->getNglob();
    ParArray<double>* interior;
    if(readFromStats==1)
    {
        interior  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","stats-mean",0,Nel,comm,info);
    }
    else{
        interior  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
    }
    
    Array<double>* ghost        = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);

    
    Array<int>*    zdefs        = ReadDataSetFromGroupFromFile<int>(fn_grid,"zones","zdefs");
    Array<char>*  znames        = ReadDataSetFromGroupFromFile<char>(fn_grid,"zones","znames");
    
    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    int t=0;
    for(int i=4;i<zdefs->getNrow();i++)
    {
        bnd_m.push_back(zdefs->getVal(i,3));
    }
    bnd_m.push_back(zdefs->getVal(zdefs->getNrow()-1,4));
    
    if(rank == 0)
    {
       PlotBoundaryData(znames,zdefs);
    }
    
    std::map<int,char*> znames_map;
    Array<char>* znames_new = new Array<char>(znames->getNrow(),znames->getNcol());
    for(int i=0;i<zdefs->getNrow();i++)
    {
        if(zdefs->getVal(i,5)!=1)
        {
            char* name = new char[znames->getNcol()];

            for(int j=0;j<znames->getNcol();j++)
            {
               name[j]=znames->getVal(i,j);
            }
            znames_map[zdefs->getVal(i,5)] = name;
        }
    }
    
    // number of vertices
    for(int j=0;j<znames->getNcol();j++)
    {
       znames_new->setVal(0,j,znames->getVal(0,j));
    }
    // number of cells
    for(int j=0;j<znames->getNcol();j++)
    {
       znames_new->setVal(1,j,znames->getVal(1,j));
    }
    
    std::map<int,char*>::iterator itch;
    int c=2;
    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
    {
        int bid = itch->first;
        for(int j=0;j<znames->getNcol();j++)
        {
            znames_new->setVal(c,j,znames_map[bid][j]);
        }
        c++;
    }
    
    int i,j;
    int nglob = ien->getNglob();
    int nrow  = ien->getNrow();
    int ncol  = ien->getNcol()-1;
    //
    ParArray<int>* ien_copy = new ParArray<int>(nglob,ncol,comm);
    //
    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            ien_copy->setVal(i,j,ien->getVal(i,j+1)-1);
        }
    }
    delete ien;
    //
    int ncol_ief = ief->getNcol()-1;
    ParArray<int>* ief_copy = new ParArray<int>(nglob,ncol_ief,comm);

    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol_ief;j++)
        {
            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
        }
    }
    delete ief;
    
    int ncol_iee = iee->getNcol()-1;
    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
    
    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol_iee;j++)
        {
            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
        }
    }
    
    delete iee;
    
    int nrow_fglob = ifn->getNglob();
    int nrow_floc  = ifn->getNrow();
    int ncol_ifn = 4;
    int ncol_ife = 2;
    ParArray<int>* ifn_copy    = new ParArray<int>(nrow_fglob,ncol_ifn,comm);
    ParArray<int>* ife_copy    = new ParArray<int>(nrow_fglob,ncol_ife,comm);
    ParArray<int>* if_ref_copy = new ParArray<int>(nrow_fglob,1,comm);
    ParArray<int>* if_Nv_copy  = new ParArray<int>(nrow_fglob,1,comm);

    for(i=0;i<nrow_floc;i++)
    {
        if_ref_copy->setVal(i,0,ifn->getVal(i,7));
        if_Nv_copy->setVal(i,0,ifn->getVal(i,0));

        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy->setVal(i,j,ifn->getVal(i,j+1)-1);
        }
        for(j=0;j<ncol_ife;j++)
        {
            ife_copy->setVal(i,j,ife->getVal(i,j)-1);
        }
    }
    
    delete ifn;
    delete ife;
    
    ParArray<int>* ie_Nv    = new ParArray<int>(nglob,1,comm);
    ParArray<int>* ie_Nf    = new ParArray<int>(nglob,1,comm);
    
    int check_hex = 0;
    int check_tet = 0;
    int check_pri = 0;
    
    
    for(int i=0;i<nrow;i++)
    {
        if(iet->getVal(i,0)==2) // Tet
        {
            ie_Nv->setVal(i,0,4);
            ie_Nf->setVal(i,0,4);
            check_tet = 1;
        }
        if(iet->getVal(i,0)==6) // Prism
        {
            ie_Nv->setVal(i,0,6);
            ie_Nf->setVal(i,0,5);
            check_pri = 1;
        }
        if(iet->getVal(i,0)==4) // Hex
        {
            ie_Nv->setVal(i,0,8);
            ie_Nf->setVal(i,0,6);
            check_hex = 1;
        }
    }
    
    Array<int>* elTypes = new Array<int>(3,1);
    elTypes->setVal(0,0,check_tet);
    elTypes->setVal(1,0,check_pri);
    elTypes->setVal(2,0,check_hex);
    
    us3d->xcn           = xcn;
    us3d->elTypes       = elTypes;
    us3d->ien           = ien_copy;
    us3d->ief           = ief_copy;
    us3d->iee           = iee_copy;
    us3d->iet           = iet;
    us3d->ie_Nv         = ie_Nv;
    us3d->ie_Nf         = ie_Nf;
    us3d->if_Nv         = if_Nv_copy;

    us3d->ifn           = ifn_copy;
    us3d->if_ref        = if_ref_copy;
    us3d->ife           = ife_copy;

    us3d->interior      = interior;
    us3d->ghost         = ghost;
    
    us3d->znames        = znames_new;
    us3d->zdefs         = zdefs;

    //delete zdefs;
    delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;
    return us3d;
}

