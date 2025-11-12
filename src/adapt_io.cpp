#include "adapt_io.h"
#include "adapt_output.h"
#include "NekFace.h"

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


double ReadStatisticsTimeFromRunInFileInParallel(const char* file_name, const char* run_name, MPI_Comm comm, MPI_Info info)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    herr_t ret;
    double stime;
    hid_t file_id        = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id       = H5Gopen(file_id,"solution",H5P_DEFAULT);
    hid_t run_id         = H5Gopen(group_id,run_name,H5P_DEFAULT);
    hid_t attr           = H5Aopen(run_id,"stats_time", H5P_DEFAULT);
    ret                  = H5Aread(attr, H5T_NATIVE_DOUBLE, &stime);
        
    return stime;
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
    hid_t file_id = H5Fcreate("grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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
    int fid2 = 0;
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

        face_start = face_end+1;

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
    hid_t file_id = H5Fcreate("grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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
    int fid2 = 0;
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







void WriteUS3DGridFromMMG_it0_NEW(MMG5_pMesh mmgMesh,MMG5_pSol mmgSol, US3D* us3d)
{
    std::map<int,std::vector<int> > ref2bface;
    std::map<int,std::vector<int> > ref2bqface;
    std::set<set<int> > bfaces;
    std::set<set<int> > bqfaces;
    std::map<set<int>,int > btfaces_Ref;
    std::map<set<int>,int > bqfaces_Ref;
    std::set<int> face;
    std::set<int> bcrefs;
    
    FaceSetPointer m_boundTrias_Ref;
    FaceSetPointer m_boundQuads_Ref;
    FaceSetPointer m_boundFaces;

    FaceSetPointer m_Faces;

    
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
            
            std::vector<int> faceVec(3);
            faceVec[0] = mmgMesh->tria[i].v[0];
            faceVec[1] = mmgMesh->tria[i].v[1];
            faceVec[2] = mmgMesh->tria[i].v[2];
            
            FaceSharedPtr BoundFacePointer = std::shared_ptr<NekFace>(new NekFace(faceVec));
            pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_boundTrias_Ref.insert(BoundFacePointer);
            
            if(testInsPointer.second)
            {
                (*testInsPointer.first)->SetFaceRef(mmgMesh->tria[i].ref);
            }
            
            if(btfaces_Ref.find(face)==btfaces_Ref.end())
            {
                btfaces_Ref[face] = mmgMesh->tria[i].ref;
            }
            
            pair<FaceSetPointer::iterator, bool> testPointer;
            testPointer = m_boundFaces.insert(BoundFacePointer);
            
            
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
            
            std::vector<int> faceVec(4);
            faceVec[0] = mmgMesh->quadra[i].v[0];
            faceVec[1] = mmgMesh->quadra[i].v[1];
            faceVec[2] = mmgMesh->quadra[i].v[2];
            faceVec[3] = mmgMesh->quadra[i].v[3];
            
            FaceSharedPtr BoundFacePointer = std::shared_ptr<NekFace>(new NekFace(faceVec));
            
            pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_boundQuads_Ref.insert(BoundFacePointer);
            
            pair<FaceSetPointer::iterator, bool> testPointer;
            testPointer = m_boundFaces.insert(BoundFacePointer);
            
            
            if(testInsPointer.second)
            {
                (*testInsPointer.first)->SetFaceRef(mmgMesh->quadra[i].ref);
            }
            
            
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
    
    std::map<int,std::vector<FaceSharedPtr> > bctrias2;
    std::map<int,std::vector<FaceSharedPtr> > bcquads2;
    
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
    hid_t file_id = H5Fcreate("grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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
    int fid2 = 0;
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

    
    std::vector<std::vector<int> > loc_faces_tet;
    std::vector<int> f0t(3);
    f0t[0] = 1;f0t[1] = 2;f0t[2] = 3;
    loc_faces_tet.push_back(f0t);
    std::vector<int> f1t(3);
    f1t[0] = 0;f1t[1] = 2;f1t[2] = 3;
    loc_faces_tet.push_back(f1t);
    std::vector<int> f2t(3);
    f2t[0] = 0;f2t[1] = 3;f2t[2] = 1;
    loc_faces_tet.push_back(f2t);
    std::vector<int> f3t(3);
    f3t[0] = 0;f3t[1] = 2;f3t[2] = 1;
    loc_faces_tet.push_back(f3t);
    for(int i=1;i<=mmgMesh->ne;i++)
    {
        adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
        
        for(int j=0;j<4;j++)
        {
            std::vector<int> faceVec(3);
            faceVec[0] = mmgMesh->tetra[i].v[loc_faces_tet[j][0]];
            faceVec[1] = mmgMesh->tetra[i].v[loc_faces_tet[j][1]];
            faceVec[2] = mmgMesh->tetra[i].v[loc_faces_tet[j][2]];
            
            pair<FaceSetPointer::iterator, bool> FPointer;
            FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(faceVec));
            
            FPointer = m_Faces.insert(FacePointer);
            
            if(FPointer.second)
            {
                (*FPointer.first)->SetFaceLeftElement(i-1);
                (*FPointer.first)->SetFaceRightElement(0);
                (*FPointer.first)->SetFaceID(fid2);
                fid2++;
            }
            else
            {
                (*FPointer.first)->SetFaceRightElement(i-1);
            }
        }
        
        
        
        //====================================================================
        
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
//            faces.erase(face0);
//            facemap.erase(face0);
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
//            faces.erase(face1);
//            facemap.erase(face1);
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
//            faces.erase(face2);
//            facemap.erase(face2);
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
//            faces.erase(face3);
//            facemap.erase(face3);
        }
    
        face0.clear();
        face1.clear();
        face2.clear();
        face3.clear();
        
    }
    
    
    
    
    
    std::cout << "lhsize - " << lh.size() << " " << m_Faces.size() << std::endl;
    

    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

    
    std::vector<std::vector<int> > loc_faces_prism;
    std::vector<int> f0p(3);
    f0p[0] = 0;
    f0p[1] = 2;
    f0p[2] = 1;
    loc_faces_prism.push_back(f0p);
    std::vector<int> f1p(3);
    f1p[0] = 3;
    f1p[1] = 4;
    f1p[2] = 5;
    loc_faces_prism.push_back(f1p);
    std::vector<int> f2p(4);
    f2p[0] = 0;
    f2p[1] = 2;
    f2p[2] = 4;
    f2p[3] = 3;
    loc_faces_prism.push_back(f2p);
    std::vector<int> f3p(4);
    f3p[0] = 1;
    f3p[1] = 5;
    f3p[2] = 4;
    f3p[3] = 2;
    loc_faces_prism.push_back(f3p);
    std::vector<int> f4p(4);
    f4p[0] = 0;
    f4p[1] = 3;
    f4p[2] = 5;
    f4p[3] = 1;
    loc_faces_prism.push_back(f4p);
    
    
    for(int i=1;i<=mmgMesh->nprism;i++)
    {
        adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are
        
        for(int j=0;j<5;j++)
        {
            int nvrt = loc_faces_prism[j].size();
            
            std::vector<int> faceVec(nvrt);
            
            for(int k=0;k<nvrt;k++)
            {
                faceVec[k] = mmgMesh->prism[i].v[loc_faces_prism[j][k]];
            }
            
            pair<FaceSetPointer::iterator, bool> FPointer;
            FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(faceVec));
            
            FPointer = m_Faces.insert(FacePointer);
            
            if(FPointer.second)
            {
                (*FPointer.first)->SetFaceLeftElement(mmgMesh->ne+i-1);
                (*FPointer.first)->SetFaceRightElement(0);
                (*FPointer.first)->SetFaceID(fid2);
                fid2++;
            }
            else
            {
                (*FPointer.first)->SetFaceRightElement(mmgMesh->ne+i-1);
            }
        }
        
        
        
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
//            faces.erase(face0);
//            facemap.erase(face0);
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
//            faces.erase(face1);
//            facemap.erase(face1);
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
//            qfaces.erase(qface0);
//            qfacemap.erase(qface0);
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
//            qfaces.erase(qface1);
//            qfacemap.erase(qface1);
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
//            qfaces.erase(qface2);
//            qfacemap.erase(qface2);
        }
        
        face0.clear();
        face1.clear();
        qface0.clear();
        qface1.clear();
        qface2.clear();
         
    }
    
    std::cout << "lhsize after - " << lh.size() << " " << m_Faces.size() << std::endl;

    FaceSetPointer::iterator itf;
    int tes = 0;
    for(itf=m_Faces.begin();itf!=m_Faces.end();itf++)
    {
        int fidnew    = (*itf)->GetFaceID();
        int faceLelem = (*itf)->GetFaceLeftElement();
        int faceRelem = (*itf)->GetFaceRightElement();
        
        std::vector<int> fv = (*itf)->GetEdgeIDs();

        std::set<int> fs;
        for(int y=0;y<fv.size();y++)
        {
            fs.insert(fv[y]);
        }
        
        if(fs.size()==3)
        {
            if(facemap.find(fs)!=facemap.end())
            {
                int fsid = facemap[fs];
                if(fsid!=fidnew)
                {
                    std::cout << fsid << " " << fidnew << " HELP " << std::endl;
                    tes++;
                }
            }
        }
        
        if(fs.size()==4)
        {
            if(qfacemap.find(fs)!=qfacemap.end())
            {
                int fsid = qfacemap[fs];
                if(fsid!=fidnew)
                {
                    std::cout << fsid << " " << fidnew << " HELP " << std::endl;
                    tes++;
                }
            }
        }
    }
    
    if(tes>0)
    {
        std::cout << "Test failed" << std::endl;
    }
    else
    {
        std::cout << "Test passed" << std::endl;
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
    
    
    int ftot2 = 0;
    std::map<int,int> new2old2;
    FaceSetPointer::iterator fiter;
    for(fiter=m_Faces.begin();fiter!=m_Faces.end();fiter++)
    {
        int fidnew    = (*fiter)->GetFaceID();
        int faceLelem = (*fiter)->GetFaceLeftElement();
        int faceRelem = (*fiter)->GetFaceRightElement();

        if(faceRelem != 0)
        {
            new2old2[fidnew] = ftot2;
//            (*fiter)->SetFaceID(ftot2);
            ftot2++;
        }
    }
    
    
    
    
    std::cout << "AND WE ARE DONE DETERMINING THE ORDER" << std::endl;
    std::cout << "Starting creating massive array" << std::endl;

    Array<int>* adapt_ifn = new Array<int>(fid,8);
    Array<int>* adapt_ifn2 = new Array<int>(fid2,8);
    std::cout << "Finished creating massive array " << fid << " " << fid2 << " " << new2old.size() << " " << new2old2.size() << std::endl;

    faces.clear();
    qfaces.clear();
    facemap.clear();
    qfacemap.clear();
    
    fid = 0;
    fid2 = 0;
    int idx = 0;
    std::map<std::set<int>,int> bctFace2lh;
    std::map<std::set<int>,int> bcqFace2lh;
    
     for(int i=1;i<=mmgMesh->ne;i++)
     {
         //adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
         
         for(int j=0;j<4;j++)
         {
             std::vector<int> faceVec(3);
             faceVec[0] = mmgMesh->tetra[i].v[loc_faces_tet[j][0]];
             faceVec[1] = mmgMesh->tetra[i].v[loc_faces_tet[j][1]];
             faceVec[2] = mmgMesh->tetra[i].v[loc_faces_tet[j][2]];
             
             FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(faceVec));
             FaceSetPointer::iterator FPointer = m_Faces.find(FacePointer);
             
             int lelem  = (*FPointer)->GetFaceLeftElement();
             int relem  = (*FPointer)->GetFaceRightElement();
             int faceID = (*FPointer)->GetFaceID();
             bool handled = (*FPointer)->GetHandled();

             if(handled==false)
             {
                 if(relem == 0)
                 {
                     FaceSetPointer::iterator BFPointer = m_boundFaces.find(FacePointer);
                     
                     if(BFPointer != m_boundFaces.end())
                     {
                         int bfref = (*BFPointer)->GetFaceRef();
                         bctrias2[bfref].push_back((*BFPointer));
                         (*BFPointer)->SetFaceLeftElement(lelem);
                     }
                     
                     (*FPointer)->SetHandled(true);
                 }
                 else
                 {
                     int idx2 = new2old2[faceID];
                     adapt_ifn2->setVal(idx2,0,3);
                     adapt_ifn2->setVal(idx2,1,faceVec[0]);
                     adapt_ifn2->setVal(idx2,2,faceVec[1]);
                     adapt_ifn2->setVal(idx2,3,faceVec[2]);
                     adapt_ifn2->setVal(idx2,4,0);
                     adapt_ifn2->setVal(idx2,5,relem+1);
                     adapt_ifn2->setVal(idx2,6,lelem+1);
                     adapt_ifn2->setVal(idx2,7,2);
                     std::cout << faceVec[0] << " " << faceVec[1] << " " << faceVec[1] << std::endl;
                     (*FPointer)->SetHandled(true);
                 }
             }
             
         }
         
     }
     
     

     for(int i=1;i<=mmgMesh->nprism;i++)
     {
         //adapt_iet->setVal(mmgMesh->ne+i-1,0,6); // Element type = 6 since we are dealing with prisms.
         
         for(int j=0;j<5;j++)
         {
             int nvrt = loc_faces_prism[j].size();
             
             if(nvrt==3)
             {
                 std::vector<int> faceVec(nvrt);
                 
                 for(int k=0;k<nvrt;k++)
                 {
                     faceVec[k] = mmgMesh->prism[i].v[loc_faces_prism[j][k]];
                 }
                 
                 FaceSharedPtr FacePointer         = std::shared_ptr<NekFace>(new NekFace(faceVec));
                 FaceSetPointer::iterator FPointer = m_Faces.find(FacePointer);
                 
                 
                 int lelem  = (*FPointer)->GetFaceLeftElement();
                 int relem  = (*FPointer)->GetFaceRightElement();
                 int faceID = (*FPointer)->GetFaceID();
                 bool handled = (*FPointer)->GetHandled();
                 
                 if(handled == false)
                 {
                     if(relem == 0)
                     {
                         FaceSetPointer::iterator BFPointer = m_boundFaces.find(FacePointer);
                                      
                         if(BFPointer != m_boundFaces.end())
                         {
                             
                             int bfref = (*BFPointer)->GetFaceRef();
                             bctrias2[bfref].push_back((*BFPointer));
                             (*BFPointer)->SetFaceLeftElement(lelem);
                         }
                         
                         (*FPointer)->SetHandled(true);
                     }
                     else
                     {
                         int idx2 = new2old2[faceID];
                         adapt_ifn2->setVal(idx2,0,3);
                         adapt_ifn2->setVal(idx2,1,faceVec[0]);
                         adapt_ifn2->setVal(idx2,2,faceVec[1]);
                         adapt_ifn2->setVal(idx2,3,faceVec[2]);
                         adapt_ifn2->setVal(idx2,4,0);
                         adapt_ifn2->setVal(idx2,5,relem+1);
                         adapt_ifn2->setVal(idx2,6,lelem+1);
                         adapt_ifn2->setVal(idx2,7,2);
                         
                         (*FPointer)->SetHandled(true);
                     }
                     
                 }
             }
             
             
             if(nvrt==4)
             {
                 std::vector<int> faceVec(nvrt);

                 for(int k=0;k<nvrt;k++)
                 {
                     faceVec[k] = mmgMesh->prism[i].v[loc_faces_prism[j][k]];
                 }
                 
                 FaceSharedPtr FacePointer         = std::shared_ptr<NekFace>(new NekFace(faceVec));
                 FaceSetPointer::iterator FPointer = m_Faces.find(FacePointer);
                 
                 int lelem  = (*FPointer)->GetFaceLeftElement();
                 int relem  = (*FPointer)->GetFaceRightElement();
                 int faceID = (*FPointer)->GetFaceID();
                 bool handled = (*FPointer)->GetHandled();

                 if(handled == false)
                 {
                     if(relem == 0)
                     {
                         FaceSetPointer::iterator BFPointer = m_boundFaces.find(FacePointer);
                         
                         if(BFPointer != m_boundFaces.end())
                         {
                             int bfref = (*BFPointer)->GetFaceRef();
                             bcquads2[bfref].push_back((*BFPointer));
                             (*BFPointer)->SetFaceLeftElement(lelem);
                         }
                         
                         (*FPointer)->SetHandled(true);

                     }
                     else
                     {
                         int idx2 = new2old2[faceID];
                         adapt_ifn2->setVal(idx2,0,4);
                         adapt_ifn2->setVal(idx2,1,faceVec[0]);
                         adapt_ifn2->setVal(idx2,2,faceVec[1]);
                         adapt_ifn2->setVal(idx2,3,faceVec[2]);
                         adapt_ifn2->setVal(idx2,4,faceVec[3]);
                         adapt_ifn2->setVal(idx2,5,relem+1);
                         adapt_ifn2->setVal(idx2,6,lelem+1);
                         adapt_ifn2->setVal(idx2,7,2);
                         
                         (*FPointer)->SetHandled(true);

                     }
                 }
             }
         }
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
    std::map<int,std::vector<FaceSharedPtr> >::iterator iterbc2;
    std::map<int,std::vector<std::vector<int> > >::iterator iterbc;
    
    int t = ftot;
    for(iterbc2=bctrias2.begin();iterbc2!=bctrias2.end();iterbc2++)
    {
        int bnd_id     = iterbc2->first;
        int Ntris      = iterbc2->second.size();
        int Nquads     = bcquads2[bnd_id].size();
        std::cout << "bnd ref test = " << bnd_id << " " << Ntris << " " << Nquads << std::endl;
    }
    int lhi,lhi2;
    
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
    
    t = ftot;
    
    for(iterbc2=bctrias2.begin();iterbc2!=bctrias2.end();iterbc2++)
    {
        int bnd_id     = iterbc2->first;
        int Ntris      = iterbc2->second.size();
        int Nquads     = bcquads2[bnd_id].size();
        
        std::cout << "bnd ref = " << bnd_id << " " << Ntris << " " << Nquads << std::endl;
        
        for(int q=0;q<Ntris;q++)
        {
            //faceid = ref2bface[bnd_id][q];
            adapt_ifn2->setVal(t,0,3);
            adapt_ifn2->setVal(t,1,iterbc2->second[q]->GetEdgeIDs()[0]);
            adapt_ifn2->setVal(t,2,iterbc2->second[q]->GetEdgeIDs()[1]);
            adapt_ifn2->setVal(t,3,iterbc2->second[q]->GetEdgeIDs()[2]);
            adapt_ifn2->setVal(t,4,0);
            
            lhi2 = iterbc2->second[q]->GetFaceLeftElement();
            
            adapt_ifn2->setVal(t,5,0);
            adapt_ifn2->setVal(t,6,lhi2+1);
            adapt_ifn2->setVal(t,7,bnd_id);

            iface.clear();
            t++;
        }
        for(int q=0;q<Nquads;q++)
        {
            //faceid = ref2bqface[bnd_id][q];
            adapt_ifn2->setVal(t,0,4);
            adapt_ifn2->setVal(t,1,bcquads2[bnd_id][q]->GetEdgeIDs()[0]);
            adapt_ifn2->setVal(t,2,bcquads2[bnd_id][q]->GetEdgeIDs()[1]);
            adapt_ifn2->setVal(t,3,bcquads2[bnd_id][q]->GetEdgeIDs()[2]);
            adapt_ifn2->setVal(t,4,bcquads2[bnd_id][q]->GetEdgeIDs()[3]);
            
//          iface.insert(bcquads[bnd_id][q][0]);
//          iface.insert(bcquads[bnd_id][q][1]);
//          iface.insert(bcquads[bnd_id][q][2]);
//          iface.insert(bcquads[bnd_id][q][3]);
                        
            //lhi = bcqFace2lh[iface];
            
            lhi2 = bcquads2[bnd_id][q]->GetFaceLeftElement();
            
            adapt_ifn2->setVal(t,5,0);
            adapt_ifn2->setVal(t,6,lhi2+1);
            adapt_ifn2->setVal(t,7,bnd_id);
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
    count[0]  = dimsf[0];
    count[1]  = dimsf[1];
    offset[0] = 0;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_ifn2->data);
    
    delete adapt_ifn;
    delete adapt_ifn2;
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
//        for(int j=0;j<ncol_ief;j++)
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
//        for(int j=0;j<ncol_iee;j++)
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


int ProvideBoundaryID(int findex, std::map<int,std::vector<int> > ranges)
{
    std::map<int,std::vector<int> >::iterator it;
    int retid = 0;
    for(it=ranges.begin();it!=ranges.end();it++)
    {
        int bndid  = it->first;
        int low    = it->second[0];
        int high   = it->second[1];
        
        if(findex>=low && findex<=high)
        {
            retid = bndid;
            //std::cout << low << " " << findex << " " << high << std::endl;
            break;
        }
    }
    return retid;
}


int ProvideBoundaryRef(int findex, std::map<int,std::vector<int> > ranges)
{
    std::map<int,std::vector<int> >::iterator it;
    int retref;
    for(it=ranges.begin();it!=ranges.end();it++)
    {
        int bndref = it->first;
        int low    = it->second[0];
        int high   = it->second[1];
        
        if(findex>=low && findex<=high)
        {
            retref = bndref;
            break;
        }
    }
    return retref;
}




std::map<int, std::vector<int> > ReadElementsFromPyFRMeshFileInParallel_Lite(const char* file_name, const char* group_name, const char* dataset_name, MPI_Comm comm, MPI_Info info)
{
    std::map<int, std::vector<int> > PAv;
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::cout << rank << " " << size << std::endl;
    
    hid_t acc_tpl1          = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    //herr_t ret            = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    //herr_t ret            = H5Pset_dxpl_mpio(,comm,info);
    acc_tpl1                = H5P_DEFAULT;
    // Open file and data set to get dimensions of array;
    
    hid_t file_id           = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t group_id          = H5Gopen(file_id,group_name,H5P_DEFAULT);
    hid_t dset_id           = H5Dopen(group_id,dataset_name,H5P_DEFAULT);
    hid_t dspace            = H5Dget_space(dset_id);
    // int ndims            = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    size_t N = dims[0]; 
    hsize_t nloc                = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    hsize_t offset              = rank*int(N/size) + MIN(rank, N%size);

    // Detect element order by checking the dataset type
    hid_t dtype = H5Dget_type(dset_id);
    hid_t nodes_member_type = H5Tget_member_type(dtype, 0); // Get the 'nodes' field type
    int nodes_array_rank = H5Tget_array_ndims(nodes_member_type);
    hsize_t nodes_array_dims[1];
    H5Tget_array_dims(nodes_member_type, nodes_array_dims);
    int num_nodes = nodes_array_dims[0];
    
    // Determine order: 4 nodes = first order tet, 6 nodes = first order pri, 10 nodes = second order
    bool is_second_order = (num_nodes == 10);
    int num_corner_nodes = (strcmp(dataset_name, "pri") == 0) ? 6 : 4;
    
    H5Tclose(nodes_member_type);
    H5Tclose(dtype);

    //============================================================================
    if (strcmp(dataset_name, "pri") == 0)
    {
        // Face compound type
        hid_t face_tid = H5Tcreate(H5T_COMPOUND, sizeof(FaceEntry));
        H5Tinsert(face_tid, "cidx", offsetof(FaceEntry, cidx), H5T_NATIVE_INT16);
        H5Tinsert(face_tid, "off", offsetof(FaceEntry, off), H5T_NATIVE_INT64);

        hsize_t offsets[1] = { offset };     // offset for this process
        hsize_t counts[1]  = { nloc };      // number of rows for this process

        H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
        hid_t mspace = H5Screate_simple(1, counts, NULL);

        std::cout << "Reading prism elements, order: " << (is_second_order ? "2" : "1") << ", nloc: " << nloc << std::endl;

        if (is_second_order) {
            // Second order prism: 10 nodes, 4 faces
            hid_t arr_nodes_tid = H5Tarray_create(H5T_NATIVE_INT64, 1, (const hsize_t[]){10});
            hid_t arr_faces_tid = H5Tarray_create(face_tid, 1, (const hsize_t[]){4});
            hid_t prism_tid = H5Tcreate(H5T_COMPOUND, sizeof(PrismEntry_O2));
            H5Tinsert(prism_tid, "nodes", offsetof(PrismEntry_O2, nodes), arr_nodes_tid);
            H5Tinsert(prism_tid, "curved", offsetof(PrismEntry_O2, curved), H5T_NATIVE_INT8);
            H5Tinsert(prism_tid, "faces", offsetof(PrismEntry_O2, faces), arr_faces_tid);

            std::vector<PrismEntry_O2> elements_loc(nloc);
            H5Dread(dset_id, prism_tid, mspace, dspace, H5P_DEFAULT, elements_loc.data());

            int i = 0;
            for (const PrismEntry_O2& element : elements_loc) {
                std::vector<int> row(6, 0);
                int elid = offset + i;
                // Extract only corner nodes (first 6)
                for (int j = 0; j < 6; j++) {
                    row[j] = element.nodes[j];
                }
                PAv[elid] = row;
                i++;
            }

            H5Tclose(arr_faces_tid);
            H5Tclose(arr_nodes_tid);
            H5Tclose(prism_tid);
        } else {
            // First order prism: 6 nodes, 5 faces
            hid_t arr_nodes_tid = H5Tarray_create(H5T_NATIVE_INT64, 1, (const hsize_t[]){6});
            hid_t arr_faces_tid = H5Tarray_create(face_tid, 1, (const hsize_t[]){5});
            hid_t prism_tid = H5Tcreate(H5T_COMPOUND, sizeof(PrismEntry_O1));
            H5Tinsert(prism_tid, "nodes", offsetof(PrismEntry_O1, nodes), arr_nodes_tid);
            H5Tinsert(prism_tid, "curved", offsetof(PrismEntry_O1, curved), H5T_NATIVE_INT8);
            H5Tinsert(prism_tid, "faces", offsetof(PrismEntry_O1, faces), arr_faces_tid);

            std::vector<PrismEntry_O1> elements_loc(nloc);
            H5Dread(dset_id, prism_tid, mspace, dspace, H5P_DEFAULT, elements_loc.data());

            int i = 0;
            for (const PrismEntry_O1& element : elements_loc) {
                std::vector<int> row(6, 0);
                int elid = offset + i;
                for (int j = 0; j < 6; j++) {
                    row[j] = element.nodes[j];
                }
                PAv[elid] = row;
                i++;
            }

            H5Tclose(arr_faces_tid);
            H5Tclose(arr_nodes_tid);
            H5Tclose(prism_tid);
        }

        H5Tclose(face_tid);
        H5Sclose(mspace);
    }


    if (strcmp(dataset_name, "tet") == 0)
    {
        // Face compound type
        hid_t face_tid = H5Tcreate(H5T_COMPOUND, sizeof(FaceEntry));
        H5Tinsert(face_tid, "cidx", offsetof(FaceEntry, cidx), H5T_NATIVE_INT16);
        H5Tinsert(face_tid, "off", offsetof(FaceEntry, off), H5T_NATIVE_INT64);

        hsize_t offsets[1] = { offset };     // offset for this process
        hsize_t counts[1]  = { nloc };      // number of rows for this process

        H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
        hid_t mspace = H5Screate_simple(1, counts, NULL);

        std::cout << "Reading tetrahedron elements, order: " << (is_second_order ? "2" : "1") << ", nloc: " << nloc << std::endl;

        if (is_second_order) {
            // Second order tetrahedron: 10 nodes, 4 faces
            hid_t arr_nodes_tid = H5Tarray_create(H5T_NATIVE_INT64, 1, (const hsize_t[]){10});
            hid_t arr_faces_tid = H5Tarray_create(face_tid, 1, (const hsize_t[]){4});
            hid_t tet_tid = H5Tcreate(H5T_COMPOUND, sizeof(TetEntry_O2));
            H5Tinsert(tet_tid, "nodes", offsetof(TetEntry_O2, nodes), arr_nodes_tid);
            H5Tinsert(tet_tid, "curved", offsetof(TetEntry_O2, curved), H5T_NATIVE_INT8);
            H5Tinsert(tet_tid, "faces", offsetof(TetEntry_O2, faces), arr_faces_tid);

            std::vector<TetEntry_O2> elements_loc(nloc);
            H5Dread(dset_id, tet_tid, mspace, dspace, H5P_DEFAULT, elements_loc.data());

            int i = 0;
            for (const TetEntry_O2& element : elements_loc) {
                std::vector<int> row(4, 0);
                int elid = offset + i;
                // Extract only corner nodes (first 4)
                for (int j = 0; j < 4; j++) {
                    row[j] = element.nodes[j];
                }
                PAv[elid] = row;
                i++;
            }

            H5Tclose(arr_faces_tid);
            H5Tclose(arr_nodes_tid);
            H5Tclose(tet_tid);
        } else {
            // First order tetrahedron: 4 nodes, 4 faces
            hid_t arr_nodes_tid = H5Tarray_create(H5T_NATIVE_INT64, 1, (const hsize_t[]){4});
            hid_t arr_faces_tid = H5Tarray_create(face_tid, 1, (const hsize_t[]){4});
            hid_t tet_tid = H5Tcreate(H5T_COMPOUND, sizeof(TetEntry_O1));
            H5Tinsert(tet_tid, "nodes", offsetof(TetEntry_O1, nodes), arr_nodes_tid);
            H5Tinsert(tet_tid, "curved", offsetof(TetEntry_O1, curved), H5T_NATIVE_INT8);
            H5Tinsert(tet_tid, "faces", offsetof(TetEntry_O1, faces), arr_faces_tid);

            std::vector<TetEntry_O1> elements_loc(nloc);
            H5Dread(dset_id, tet_tid, mspace, dspace, H5P_DEFAULT, elements_loc.data());

            int i = 0;
            for (const TetEntry_O1& element : elements_loc) {
                std::vector<int> row(4, 0);
                int elid = offset + i;
                for (int j = 0; j < 4; j++) {
                    row[j] = element.nodes[j];
                }
                PAv[elid] = row;
                i++;
            }

            H5Tclose(arr_faces_tid);
            H5Tclose(arr_nodes_tid);
            H5Tclose(tet_tid);
        }

        H5Tclose(face_tid);
        H5Sclose(mspace);
    }

    // for (size_t i = 0; i < nloc; ++i)
    // {
    //     std::vector<T> single_prism(elements_loc[i].node_indices,
    //                             elements_loc[i].node_indices + 6);
    //     PAv.push_back(single_prism);
    // }

    // Proper resource management:
    H5Sclose(dspace);
    H5Dclose(dset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
    H5Pclose(acc_tpl1);
    //============================================================================

    
    return PAv;
}


std::vector<std::vector<std::vector<float> > > ReadSolutionFromPyFRFileInParallel_Lite(const char* file_name, const char* element_type, int order, MPI_Comm comm, MPI_Info info)
{
    std::vector<std::vector<std::vector<float> > > solution_data;
    
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Open file and navigate to solution group
    hid_t file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t soln_group_id = H5Gopen(file_id, "soln", H5P_DEFAULT);
    
    // Construct dataset name (e.g., "p1-tet" for first order, "p2-tet" for second order)
    std::string dataset_name;
    if (order == 0) {
        // Auto-detect: try p2 first, then p1
        std::string p2_name = std::string("p2-") + std::string(element_type);
        std::string p1_name = std::string("p1-") + std::string(element_type);
        
        if (H5Lexists(soln_group_id, p2_name.c_str(), H5P_DEFAULT) > 0) {
            dataset_name = p2_name;
            if (rank == 0) std::cout << "Auto-detected second order solution: " << dataset_name << std::endl;
        } else if (H5Lexists(soln_group_id, p1_name.c_str(), H5P_DEFAULT) > 0) {
            dataset_name = p1_name;
            if (rank == 0) std::cout << "Auto-detected first order solution: " << dataset_name << std::endl;
        } else {
            if (rank == 0) {
                std::cerr << "Error: Could not find solution dataset for element type: " << element_type << std::endl;
            }
            H5Gclose(soln_group_id);
            H5Fclose(file_id);
            return solution_data;
        }
    } else {
        dataset_name = std::string("p") + std::to_string(order) + std::string("-") + std::string(element_type);
    }
    
    hid_t dset_id = H5Dopen(soln_group_id, dataset_name.c_str(), H5P_DEFAULT);
    if (dset_id < 0) {
        if (rank == 0) {
            std::cerr << "Error: Could not open dataset: " << dataset_name << std::endl;
        }
        H5Gclose(soln_group_id);
        H5Fclose(file_id);
        return solution_data;
    }
    
    // Get dataset dimensions
    hid_t dspace = H5Dget_space(dset_id);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    
    size_t n_elements = dims[0];  // Number of elements
    size_t n_vars = dims[1];      // Number of variables (typically 5 for Euler: rho, rho*u, rho*v, rho*w, rho*E)
    size_t n_pts = dims[2];       // Number of solution points per element
    
    // Partition elements across processes
    hsize_t nloc = int(n_elements/size) + (rank < n_elements%size);
    hsize_t offset = rank*int(n_elements/size) + MIN(rank, n_elements%size);
    
    std::cout << "Rank " << rank << ": Reading solution data - elements: " << n_elements 
              << ", vars: " << n_vars << ", pts: " << n_pts 
              << ", local: " << nloc << ", offset: " << offset << std::endl;
    
    if (nloc > 0) {
        // Select hyperslab for this process
        hsize_t offsets[3] = {offset, 0, 0};
        hsize_t counts[3] = {nloc, n_vars, n_pts};
        H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
        
        // Create memory space
        hid_t mspace = H5Screate_simple(3, counts, NULL);
        
        // Allocate buffer for local data
        std::vector<float> buffer(nloc * n_vars * n_pts);
        
        // Read data
        H5Dread(dset_id, H5T_NATIVE_FLOAT, mspace, dspace, H5P_DEFAULT, buffer.data());
        
        // Reshape data into nested vector structure: solution[element][variable][solution_point]
        solution_data.resize(nloc);
        for (size_t i = 0; i < nloc; i++) {
            solution_data[i].resize(n_vars);
            for (size_t j = 0; j < n_vars; j++) {
                solution_data[i][j].resize(n_pts);
                for (size_t k = 0; k < n_pts; k++) {
                    solution_data[i][j][k] = buffer[i * n_vars * n_pts + j * n_pts + k];
                }
            }
        }
        
        H5Sclose(mspace);
    }
    
    // Cleanup
    H5Sclose(dspace);
    H5Dclose(dset_id);
    H5Gclose(soln_group_id);
    H5Fclose(file_id);
    
    return solution_data;
}





mesh* ReadUS3DMeshDataWithMetric(const char* fn_conn, const char* fn_grid, const char* fn_data, int readFromStats, int StateVar, int RunNumber, MPI_Comm comm, MPI_Info info)
{
    mesh* mRead = new mesh;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<double> >   xcn = ReadDataSetFromFileInParallel_Lite<double>(fn_grid,"xcn",comm,info);
    std::vector<std::vector<int> >      ien = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ien",comm,info);
    std::vector<std::vector<int> >      ief = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ief",comm,info);
    std::vector<std::vector<int> >      iee = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"iee",comm,info);
    std::vector<std::vector<int> >      iet = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"iet",comm,info);
    std::vector<std::vector<int> >      ifn = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"ifn",comm,info);
    std::vector<std::vector<int> >      ife = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ife",comm,info);
    std::vector<int> vert_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_grid,"xcn",comm,info);
    std::vector<int> elem_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ien",comm,info);
    std::vector<int> face_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ife",comm,info);

    int nElem = elem_rankInfo[0];
    int nFace = face_rankInfo[0];
    int nVert = vert_rankInfo[0];

    Array<char>* zvnames                    = ReadDataSetFromGroupInGroupFromFile<char>(fn_data,"info","solver","svnames");

    std::map<string,int> var_map            = PlotVariableNames(zvnames);
    

    int offset_el = elem_rankInfo[2];

    // if(rank == 0)
    // {
    // std::map<std::string,int>::iterator itv;
    // for(itv=var_map.begin();itv!=var_map.end();itv++)
    // {
    //    std::cout << itv->first << " " << itv->second << std::endl;
    // }

    // }

    int M00_id = var_map["dT2dxx"];
    int M01_id = var_map["dT2dxy"];
    int M02_id = var_map["dT2dxz"];
    int M11_id = var_map["dT2dyy"];
    int M12_id = var_map["dT2dyz"];
    int M22_id = var_map["dT2dzz"];

    int Nel     = elem_rankInfo[0];
    int Nel_loc = elem_rankInfo[1];
    
    std::map<int, std::vector<double> > interior;
    std::map<int, std::vector<double> > ghost_out;

    std::string base_run_name = "run_";

    std::string full_run_name = base_run_name + std::to_string(RunNumber);
    const char* full_run_name_char = full_run_name.c_str();

    if(readFromStats==1)
    {       
        std::cout << "Error: Reading the metric directly when reading statistics is not implemented yet!" << std::endl; 
    }
    
    if(readFromStats==0)
    {
        std::vector<std::vector<double> > readdata  = ReadDataSetFromRunInFileInParallel_Lite<double>(fn_data,full_run_name_char,"interior",0,Nel,comm,info);
        
        double rhoState,uState,vState,wState,TState,VtotState,aState,MState;

        for(int u=0;u<Nel_loc;u++)
        {
            int elid = offset_el + u;
            std::vector<double> row_interior(6,0.0);

            row_interior[0]   = readdata[u][M00_id];
            row_interior[1]   = readdata[u][M01_id];
            row_interior[2]   = readdata[u][M02_id];
            row_interior[3]   = readdata[u][M11_id];
            row_interior[4]   = readdata[u][M12_id];
            row_interior[5]   = readdata[u][M22_id];

            // std::cout << row_interior[0] << " " << row_interior[1] << " " << row_interior[2] << " " << row_interior[3] << " " << row_interior[4] << " " << row_interior[5] << std::endl;

            interior[elid] = row_interior;
        }
    }
    
    std::vector<std::vector<double> > ghost = ReadUS3DGhostCellsFromRun_Lite<double>(fn_data,full_run_name_char,"interior",Nel);
    std::vector<std::vector<int> > zdefs    = ReadDataSetFromGroupFromFile_Lite<int>(fn_grid,"zones","zdefs");
    std::vector<std::vector<char> >  znames = ReadDataSetFromGroupFromFile_Lite<char>(fn_grid,"zones","znames");
    double rhoGhostState,uGhostState,vGhostState,wGhostState,TGhostState,VtotGhostState;
    int nghosts = ghost.size();
    for(int u=0;u<nghosts;u++)
    {
        int elid        = nElem+u;

        std::vector<double> row_ghost(6,0.0);

        row_ghost[0] = ghost[u][M00_id];
        row_ghost[1] = ghost[u][M01_id];
        row_ghost[2] = ghost[u][M02_id];
        row_ghost[3] = ghost[u][M11_id];
        row_ghost[4] = ghost[u][M12_id];
        row_ghost[5] = ghost[u][M22_id];

        ghost_out[elid] = row_ghost;
    }


    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    std::vector<int> low_range;
    std::vector<int> high_range;
    std::vector<int> ref_range;
    
    int t=0;
    int gg=0;
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,std::vector<int> > ranges_ref;
    std::map<int,int> zone2ref;
    
    for(int i=2;i<zdefs.size();i++)
    {
        // std::cout << zdefs[i][2] << " " << zdefs[i][3] << " " << zdefs[i][4] << " " << zdefs[i][5] << std::endl;
        bnd_m.push_back(zdefs[i][5]);
        
        low_range.push_back(zdefs[i][3]);
        high_range.push_back(zdefs[i][4]);
        ref_range.push_back(zdefs[i][5]);
        
        zone2ref[zdefs[i][2]] = zdefs[i][5];

        std::vector<int> ra(2);
        ra[0] = zdefs[i][3]-1;
        ra[1] = zdefs[i][4]-1;
        ranges_id[zdefs[i][2]] = ra;
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][3]-1);
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][4]-1);
        ranges_ref[zdefs[i][5]]=ra;
        
        // gg++;
    }
    
    
    bnd_m.push_back(zdefs[zdefs.size()-1][4]);
   
    if(rank == 0)
    {
       PlotBoundaryData_Lite(znames,zdefs);
    }
    
    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;
    std::map<int,std::vector<int> > bref2zone;
    std::map<int,char*> zone2name;

    std::map<int,int> zone2bcref;
    for(int i=0;i<zdefs.size();i++)
    {
        if(zdefs[i][5]!=1)
        {
            std::string name;
            
            char* namechar = new char[znames[0].size()];
                        
            for(int j=0;j<znames[0].size();j++)
            {
               namechar[j]=znames[i][j];

               char ch = znames[i][j];
               char chref = ' ';
               if(ch!=chref)
               {
                   name.push_back(ch);
               }
            }
            if(i>2)
            {
                bref2zone[zdefs[i][5]].push_back(i);
                zone2bcref[zdefs[i][2]] = zdefs[i][5];
                znames_map[zdefs[i][2]] = name;
                znames_map_inv[name]    = zdefs[i][2];
                zone2name[zdefs[i][2]]  = namechar;
            }
            
        }
    }
    
    // number of vertices

    std::vector<std::vector<char> > znames_new(znames.size());
    
    std::vector<char> znames_new_row0(znames[0].size());

    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row0[j] = znames[0][j];
    }
    znames_new[0] = znames_new_row0;

    std::vector<char> znames_new_row1(znames[0].size());
    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row1[j] = znames[1][j];
    }
    znames_new[1] = znames_new_row1;

    
    std::map<int,std::string>::iterator itch;
    int c=2;
    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
    {
        std::vector<char> znames_new_row(znames[0].size());
        int bid = itch->first;
        for(int j=0;j<itch->second.size();j++)
        {
            znames_new_row[j] = itch->second[j];
        }
        znames_new[c] = znames_new_row;
        c++;
    }




    //==========================================================================
    //==========================================================================
    //==========================================================================    
    int i,j;

    int nglob     = elem_rankInfo[0];
    int nrow      = elem_rankInfo[1];
    //

    std::map<int,std::vector<int> > ien_copy;
    std::map<int,std::vector<int> > ief_copy;
    std::map<int,std::vector<int> > iee_copy;
    std::map<int,int> iet_copy;
    std::map<int,int> ie_Nv_copy;
    std::map<int,int> ie_Nf_copy;


    int ncol_ien = ien[0].size()-1;
    int ncol_ief = ief[0].size()-1;
    int ncol_iee = iee[0].size()-1;
    //
    std::vector<int> element2rank_loc(nglob,0);
    //std::vector<int> iet_loc(nglob,0);
    int check_hex = 0;
    int check_tet = 0;
    int check_pyr = 0;
    int check_pri = 0;

    int tetCount  = 0;
    int hexCount  = 0;
    int pyrCount  = 0;
    int priCount  = 0;
        
    for(i=0;i<nrow;i++)
    {
        int elid = offset_el+i;

        element2rank_loc[elid] = rank;
        iet_copy[elid]         = iet[i][0];

        if(iet[i][0]==2) // Tet
        {
            ie_Nv_copy[elid] = 4;
            ie_Nf_copy[elid] = 4;
            check_tet        = check_tet+1;
            ncol_ien         = 4;
            ncol_iee         = 4;
            ncol_ief         = 4;
            tetCount++;
        }
        if(iet[i][0]==4) // Hex
        {
            ie_Nv_copy[elid] = 8;
            ie_Nf_copy[elid] = 6;
            ncol_ien         = 8;
            ncol_ief         = 6;
            ncol_iee         = 6;
            check_hex        = check_hex+1;
            hexCount++;
        }
        if(iet[i][0]==5) // Pyramid
        {
           ie_Nv_copy[elid]  = 5;
           ie_Nf_copy[elid]  = 5;
           ncol_ien          = 5;
           ncol_ief          = 5;
           ncol_iee          = 5;
           check_pyr         = check_pyr+1;
           pyrCount++;
        }
        if(iet[i][0]==6) // Prism
        {
            ie_Nv_copy[elid] = 6;
            ie_Nf_copy[elid] = 5;
            ncol_ien         = 6;
            ncol_ief         = 5;
            ncol_iee         = 5;
            check_pri        = check_pri+1;
            priCount++;
        }

        std::vector<int> ien_copy_row(ncol_ien);
        for(int j=0;j<ncol_ien;j++)
        {
            ien_copy_row[j] = ien[i][j+1]-1;

            if(ien[i][j+1]-1 > vert_rankInfo[0])
            {
                std::cout << "Raar while reading" << std::endl;
            }
        }
        
        std::vector<int> ief_copy_row(ncol_ief);
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy_row[j] = fabs(ief[i][j+1])-1;
        }

        std::vector<int> iee_copy_row(ncol_iee);
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy_row[j] = iee[i][j+1]-1;
        }

        ien_copy[elid] = ien_copy_row;
        ief_copy[elid] = ief_copy_row;
        iee_copy[elid] = iee_copy_row;

        // if(iet[i][0]!=2 && iet[i][0]!=4 && iet[i][0]!=6)
        // {
        //     std::cout << "Warning: this mesh has pyramids! " << iet[i][0] << std::endl;
        // }

    }

    ien.clear();
    iee.clear();
    iet.clear();
    ief.clear();

    std::vector<int> element2rank_glob(nglob,0);
    MPI_Allreduce(element2rank_loc.data(), element2rank_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    element2rank_loc.clear();

     int* colTetCount    = new int[size];
    int* RedcolTetCount = new int[size];
    int* OffcolTetCount = new int[size];

    for(int i=0;i<size;i++)
    {
        colTetCount[i]    = 0;
        RedcolTetCount[i] = 0;
        if(i==rank)
        {
            colTetCount[i] = tetCount;
        }
    }
    
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    
    int offset_tetC = 0;
    for(int i=0;i<size;i++)
    {
        OffcolTetCount[i] = offset_tetC;
        offset_tetC = offset_tetC+RedcolTetCount[i];
    }
    
    std::vector<int> ie_tetCnt(nrow,0);

    int tett=0;
    int pris=0;
    for(int i=0;i<nrow;i++)
    {
        ie_tetCnt[i] = -1;
        
        if(iet[i][0]==2) // Tet
        {
            ie_tetCnt[i] = OffcolTetCount[rank]+tett;
            //ie_tetCnt.setVal(i,0,OffcolTetCount[rank]+tett);
            tett++;
        }
        else
        {
            pris++;
        }
    }
    //std::cout << "before partitioning rank = " << rank << " #tets = " << tett << " #prisms " << pris << std::endl;
    std::vector<int> elTypes(4);
    int check_tet_glob, check_pri_glob, check_pyr_glob, check_hex_glob;
    MPI_Allreduce(&check_tet, &check_tet_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pri, &check_pri_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pyr, &check_pyr_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_hex, &check_hex_glob, 1, MPI_INT, MPI_SUM, comm);

    elTypes[0] = check_tet_glob;
    elTypes[1] = check_pri_glob;
    elTypes[2] = check_pyr_glob;
    elTypes[3] = check_hex_glob;
    
    delete[] OffcolTetCount;

    // std::vector<int> iet_glob(nglob,0);
    // MPI_Allreduce(iet_loc.data(), iet_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    // iet_loc.clear();

    //==========================================================================
    //==========================================================================
    //==========================================================================

    int nrow_fglob  = face_rankInfo[0];
    int nrow_floc   = face_rankInfo[1];
    int offset_floc = face_rankInfo[2];
    // int ncol_ifn    = 4;
    int ncol_ife    = 2;
    
    // ParMatrix<int> ifn_copy        = new ParMatrix<int>(nrow_fglob,ncol_ifn,comm);
    // ParMatrix<int> ife_copy        = new ParMatrix<int>(nrow_fglob,ncol_ife,comm);
    // ParMatrix<int> if_ref_copy     = new ParMatrix<int>(nrow_fglob,1,comm);
    // ParMatrix<int> if_Nv_copy      = new ParMatrix<int>(nrow_fglob,1,comm);
    
    std::map<int,std::vector<int> > ifn_copy;
    std::map<int,std::vector<int> > ife_copy;
    std::map<int,std::vector<int> > if_ref_copy;
    std::map<int,std::vector<int> > if_Nv_copy;    
    int fref2;
    int index_range = 0;
    
    for(i=0;i<nrow_floc;i++)
    {
        int fid        = offset_floc+i; 
        int fref       = ifn[i][7];
        int index      = i + offset_floc;
        int fzone      = ProvideBoundaryID(index,ranges_id);
        int fref2      = zone2ref[fzone];
        int ncol_ifn   = ifn[i][0];
        if((fref2!=fref))
        {
            std::cout << "mapping ranges are wrong"  << std::endl;
        }

        std::vector<int> ifn_copy_row(ncol_ifn);
        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy_row[j] = ifn[i][j+1]-1;
        }

        std::vector<int> ife_copy_row(ncol_ife);
        for(j=0;j<ncol_ife;j++)
        {
            ife_copy_row[j] = ife[i][j]-1;
        }

        if_ref_copy[fid].push_back(ifn[i][7]);
        if_Nv_copy[fid].push_back(ifn[i][0]);

        ifn_copy[fid] = ifn_copy_row;
        ife_copy[fid] = ife_copy_row;
    }
    ife.clear();
    ifn.clear();

    
    std::map<int,std::vector<double> > xcn_copy;
    int nloc_verts   = vert_rankInfo[1];
    int offset_verts = vert_rankInfo[2];
    for(int i=0;i<nloc_verts;i++)
    {
        int vid = offset_verts+i;
        //std::cout << "vud " << vid << " " << rank << std::endl;
        
        std::vector<double> coords(3);
        for(int j=0;j<3;j++)
        {
            coords[j] = xcn[i][j];
        }
        xcn_copy[vid] = coords;
    }
    xcn.clear();
    

    
    mRead->nElem          = nElem;
    mRead->nFace          = nFace;
    mRead->nVert          = nVert;
    mRead->xcn            = xcn_copy;
    mRead->ien            = ien_copy;
    mRead->ief            = ief_copy;
    mRead->iee            = iee_copy;
    mRead->iet            = iet_copy;
    mRead->elTypes        = elTypes;
    mRead->ie_Nv          = ie_Nv_copy;
    mRead->ie_Nf          = ie_Nf_copy;
    mRead->if_Nv          = if_Nv_copy;
    mRead->ifn            = ifn_copy;
    mRead->if_ref         = if_ref_copy;
    mRead->ife            = ife_copy;
    mRead->ie_tetCnt      = ie_tetCnt;
    mRead->interior       = interior;
    mRead->ghost          = ghost_out;
    mRead->zone2bcref     = zone2bcref;
    mRead->element2rank   = element2rank_glob;
    mRead->bref2zone      = bref2zone;
    mRead->znames_map     = znames_map;
    mRead->znames_map_inv = znames_map_inv;  
    mRead->zone2name      = zone2name;
    mRead->znames         = znames;
    mRead->zdefs          = zdefs;
    mRead->ranges_id      = ranges_id;
    mRead->ranges_ref     = ranges_ref;
    mRead->ntetra         = tetCount;
    mRead->nprism         = priCount;
    mRead->nhexahedral    = hexCount;
    mRead->npyramid       = pyrCount;

    /**/
    //delete zdefs;
    //delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;

    return mRead;

    
}








mesh* ReadUS3DMeshData(const char* fn_conn, const char* fn_grid, const char* fn_data, int readFromStats, int StateVar, int RunNumber, MPI_Comm comm, MPI_Info info)
{
    mesh* mRead = new mesh;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<double> >   xcn = ReadDataSetFromFileInParallel_Lite<double>(fn_grid,"xcn",comm,info);
    std::vector<std::vector<int> >      ien = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ien",comm,info);
    std::vector<std::vector<int> >      ief = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ief",comm,info);
    std::vector<std::vector<int> >      iee = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"iee",comm,info);
    std::vector<std::vector<int> >      iet = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"iet",comm,info);
    std::vector<std::vector<int> >      ifn = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"ifn",comm,info);
    std::vector<std::vector<int> >      ife = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ife",comm,info);
    std::vector<int> vert_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_grid,"xcn",comm,info);
    std::vector<int> elem_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ien",comm,info);
    std::vector<int> face_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ife",comm,info);

    int nElem = elem_rankInfo[0];
    int nFace = face_rankInfo[0];
    int nVert = vert_rankInfo[0];

    Array<char>* zvnames                    = ReadDataSetFromGroupInGroupFromFile<char>(fn_data,"info","solver","svnames");

    std::map<string,int> var_map            = PlotVariableNames(zvnames);
    

    int offset_el = elem_rankInfo[2];

    // if(rank == 0)
    // {
    // std::map<std::string,int>::iterator itv;
    // for(itv=var_map.begin();itv!=var_map.end();itv++)
    // {
    //    std::cout << itv->first << " " << itv->second << std::endl;
    // }

    // }

    int uid = var_map["u"];
    int vid = var_map["v"];
    int wid = var_map["w"];
    int Tid = var_map["T"];
  

    int Nel     = elem_rankInfo[0];
    int Nel_loc = elem_rankInfo[1];
    
    std::map<int, std::vector<double> > interior;
    std::map<int, std::vector<double> > ghost_out;

    std::string base_run_name = "run_";

    std::string full_run_name = base_run_name + std::to_string(RunNumber);
    const char* full_run_name_char = full_run_name.c_str();

    if(readFromStats==1)
    {       
        double time_stats = ReadStatisticsTimeFromRunInFileInParallel(fn_data,full_run_name_char,comm,info);

        std::vector<std::vector<double> > mean  = ReadDataSetFromRunInFileInParallel_Lite<double>(fn_data,full_run_name_char,"stats-mean",0,Nel,comm,info);

        double rhoState,uState,vState,wState,TState,VtotState,aState,MState;

        if(rank == 0)
        {
            std::cout << "Statistics time = " << time_stats << std::endl;
        }
        std::vector<std::vector<double> > vel_mean(Nel_loc);

        for(int u=0;u<Nel_loc;u++)
        {
            int elid = offset_el + u;
            rhoState = mean[u][0]/time_stats;
            uState   = mean[u][uid]/time_stats;
            vState   = mean[u][vid]/time_stats;
            wState   = mean[u][wid]/time_stats;
            TState   = mean[u][Tid]/time_stats;

            std::vector<double> vel_mean_row(3);
            vel_mean_row[0] = uState;
            vel_mean_row[1] = vState;
            vel_mean_row[2] = wState;

            vel_mean[u] = vel_mean_row;

            VtotState   = sqrt(uState*uState+vState*vState+wState*wState);

            MState      = VtotState;
            
            std::vector<double> row_interior(2,0.0);

            row_interior[0] = 0.0;
            row_interior[1] = 0.0;

            if(StateVar==0)
            {
                row_interior[1] = MState;
            }
            if(StateVar==1)
            {
                row_interior[1] = TState;
            }

            interior[elid] = row_interior;
        }
        
        std::vector<std::vector<double> > stats  = ReadDataSetFromRunInFileInParallel_Lite<double>(fn_data,full_run_name_char,"stats-stat",0,Nel,comm,info);

        double upup,vpvp,wpwp,tke;
        std::vector<double> tkevec(Nel_loc);
        for(int u=0;u<Nel_loc;u++)
        {
            upup = stats[u][0]/time_stats-vel_mean[u][0]*vel_mean[u][0];
            vpvp = stats[u][1]/time_stats-vel_mean[u][1]*vel_mean[u][1];
            wpwp = stats[u][2]/time_stats-vel_mean[u][2]*vel_mean[u][2];
            tke = 0.5*(upup+vpvp+wpwp);
            tkevec[u] = tke;
        }

        
        double tkeMax = *std::max_element(tkevec.begin(),tkevec.end());
        double tkeMax_glob = 0.0;
        MPI_Allreduce(&tkeMax, &tkeMax_glob, 1, MPI_DOUBLE, MPI_MAX, comm);
        //std::cout << "tkeMax_glob " << tkeMax_glob << std::endl;

        for(int u=0;u<Nel_loc;u++)
        {
            int elid = offset_el + u;
            interior[elid][0] = tkevec[u]/tkeMax_glob;
        }
    }
    
    if(readFromStats==0)
    {
        std::vector<std::vector<double> > readdata  = ReadDataSetFromRunInFileInParallel_Lite<double>(fn_data,full_run_name_char,"interior",0,Nel,comm,info);
        
        double rhoState,uState,vState,wState,TState,VtotState,aState,MState;

        for(int u=0;u<Nel_loc;u++)
        {
            int elid = offset_el + u;
            rhoState  = readdata[u][0];
            uState    = readdata[u][uid];
            vState    = readdata[u][vid];
            wState    = readdata[u][wid];
            TState    = readdata[u][Tid];
            VtotState = sqrt(uState*uState+vState*vState+wState*wState);

            std::vector<double> row_interior(2,0.0);

            row_interior[0] = VtotState;
            row_interior[1] = TState;

            interior[elid]  = row_interior;
        }
    }
    
    std::vector<std::vector<double> > ghost = ReadUS3DGhostCellsFromRun_Lite<double>(fn_data,full_run_name_char,"interior",Nel);
    std::vector<std::vector<int> > zdefs    = ReadDataSetFromGroupFromFile_Lite<int>(fn_grid,"zones","zdefs");
    std::vector<std::vector<char> >  znames = ReadDataSetFromGroupFromFile_Lite<char>(fn_grid,"zones","znames");
    double rhoGhostState,uGhostState,vGhostState,wGhostState,TGhostState,VtotGhostState;
    int nghosts = ghost.size();
    for(int u=0;u<nghosts;u++)
    {
        int elid        = nElem+u;
        rhoGhostState       = ghost[u][0];
        uGhostState          = ghost[u][uid];
        vGhostState          = ghost[u][vid];
        wGhostState          = ghost[u][wid];
        TGhostState          = ghost[u][Tid];
        VtotGhostState       = sqrt(uGhostState*uGhostState
                                    +vGhostState*vGhostState
                                    +wGhostState*wGhostState);

        std::vector<double> row_ghost(2,0.0);

        row_ghost[0] = VtotGhostState;
        row_ghost[1] = TGhostState;

        ghost_out[elid] = row_ghost;
    }


    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    std::vector<int> low_range;
    std::vector<int> high_range;
    std::vector<int> ref_range;
    
    int t=0;
    int gg=0;
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,std::vector<int> > ranges_ref;
    std::map<int,int> zone2ref;
    
    for(int i=2;i<zdefs.size();i++)
    {
        // std::cout << zdefs[i][2] << " " << zdefs[i][3] << " " << zdefs[i][4] << " " << zdefs[i][5] << std::endl;
        bnd_m.push_back(zdefs[i][5]);
        
        low_range.push_back(zdefs[i][3]);
        high_range.push_back(zdefs[i][4]);
        ref_range.push_back(zdefs[i][5]);
        
        zone2ref[zdefs[i][2]] = zdefs[i][5];

        std::vector<int> ra(2);
        ra[0] = zdefs[i][3]-1;
        ra[1] = zdefs[i][4]-1;
        ranges_id[zdefs[i][2]] = ra;
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][3]-1);
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][4]-1);
        ranges_ref[zdefs[i][5]]=ra;
        
        // gg++;
    }
    
    
    bnd_m.push_back(zdefs[zdefs.size()-1][4]);
   
    if(rank == 0)
    {
       PlotBoundaryData_Lite(znames,zdefs);
    }
    
    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;
    std::map<int,std::vector<int> > bref2zone;
    std::map<int,char*> zone2name;

    std::map<int,int> zone2bcref;
    for(int i=0;i<zdefs.size();i++)
    {
        if(zdefs[i][5]!=1)
        {
            std::string name;
            
            char* namechar = new char[znames[0].size()];
                        
            for(int j=0;j<znames[0].size();j++)
            {
               namechar[j]=znames[i][j];

               char ch = znames[i][j];
               char chref = ' ';
               if(ch!=chref)
               {
                   name.push_back(ch);
               }
            }
            if(i>2)
            {
                bref2zone[zdefs[i][5]].push_back(i);
                zone2bcref[zdefs[i][2]] = zdefs[i][5];
                znames_map[zdefs[i][2]] = name;
                znames_map_inv[name]    = zdefs[i][2];
                zone2name[zdefs[i][2]]  = namechar;
            }
            
        }
    }
    
    // number of vertices

    std::vector<std::vector<char> > znames_new(znames.size());
    
    std::vector<char> znames_new_row0(znames[0].size());

    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row0[j] = znames[0][j];
    }
    znames_new[0] = znames_new_row0;

    std::vector<char> znames_new_row1(znames[0].size());
    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row1[j] = znames[1][j];
    }
    znames_new[1] = znames_new_row1;

    
    std::map<int,std::string>::iterator itch;
    int c=2;
    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
    {
        std::vector<char> znames_new_row(znames[0].size());
        int bid = itch->first;
        for(int j=0;j<itch->second.size();j++)
        {
            znames_new_row[j] = itch->second[j];
        }
        znames_new[c] = znames_new_row;
        c++;
    }




    //==========================================================================
    //==========================================================================
    //==========================================================================    
    int i,j;

    int nglob     = elem_rankInfo[0];
    int nrow      = elem_rankInfo[1];
    //

    std::map<int,std::vector<int> > ien_copy;
    std::map<int,std::vector<int> > ief_copy;
    std::map<int,std::vector<int> > iee_copy;
    std::map<int,int> iet_copy;
    std::map<int,int> ie_Nv_copy;
    std::map<int,int> ie_Nf_copy;


    int ncol_ien = ien[0].size()-1;
    int ncol_ief = ief[0].size()-1;
    int ncol_iee = iee[0].size()-1;
    //
    std::vector<int> element2rank_loc(nglob,0);
    //std::vector<int> iet_loc(nglob,0);
    int check_hex = 0;
    int check_tet = 0;
    int check_pyr = 0;
    int check_pri = 0;

    int tetCount  = 0;
    int hexCount  = 0;
    int pyrCount  = 0;
    int priCount  = 0;
        
    for(i=0;i<nrow;i++)
    {
        int elid = offset_el+i;

        element2rank_loc[elid] = rank;
        iet_copy[elid]         = iet[i][0];

        if(iet[i][0]==2) // Tet
        {
            ie_Nv_copy[elid] = 4;
            ie_Nf_copy[elid] = 4;
            check_tet        = check_tet+1;
            ncol_ien         = 4;
            ncol_iee         = 4;
            ncol_ief         = 4;
            tetCount++;
        }
        if(iet[i][0]==4) // Hex
        {
            ie_Nv_copy[elid] = 8;
            ie_Nf_copy[elid] = 6;
            ncol_ien         = 8;
            ncol_ief         = 6;
            ncol_iee         = 6;
            check_hex        = check_hex+1;
            hexCount++;
        }
        if(iet[i][0]==5) // Pyramid
        {
           ie_Nv_copy[elid]  = 5;
           ie_Nf_copy[elid]  = 5;
           ncol_ien          = 5;
           ncol_ief          = 5;
           ncol_iee          = 5;
           check_pyr         = check_pyr+1;
           pyrCount++;
        }
        if(iet[i][0]==6) // Prism
        {
            ie_Nv_copy[elid] = 6;
            ie_Nf_copy[elid] = 5;
            ncol_ien         = 6;
            ncol_ief         = 5;
            ncol_iee         = 5;
            check_pri        = check_pri+1;
            priCount++;
        }

        std::vector<int> ien_copy_row(ncol_ien);
        for(int j=0;j<ncol_ien;j++)
        {
            ien_copy_row[j] = ien[i][j+1]-1;

            if(ien[i][j+1]-1 > vert_rankInfo[0])
            {
                std::cout << "Raar while reading" << std::endl;
            }
        }
        
        std::vector<int> ief_copy_row(ncol_ief);
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy_row[j] = fabs(ief[i][j+1])-1;
        }

        std::vector<int> iee_copy_row(ncol_iee);
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy_row[j] = iee[i][j+1]-1;
        }

        ien_copy[elid] = ien_copy_row;
        ief_copy[elid] = ief_copy_row;
        iee_copy[elid] = iee_copy_row;

        // if(iet[i][0]!=2 && iet[i][0]!=4 && iet[i][0]!=6)
        // {
        //     std::cout << "Warning: this mesh has pyramids! " << iet[i][0] << std::endl;
        // }

    }

    ien.clear();
    iee.clear();
    iet.clear();
    ief.clear();

    std::vector<int> element2rank_glob(nglob,0);
    MPI_Allreduce(element2rank_loc.data(), element2rank_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    element2rank_loc.clear();

     int* colTetCount    = new int[size];
    int* RedcolTetCount = new int[size];
    int* OffcolTetCount = new int[size];

    for(int i=0;i<size;i++)
    {
        colTetCount[i]    = 0;
        RedcolTetCount[i] = 0;
        if(i==rank)
        {
            colTetCount[i] = tetCount;
        }
    }
    
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    
    int offset_tetC = 0;
    for(int i=0;i<size;i++)
    {
        OffcolTetCount[i] = offset_tetC;
        offset_tetC = offset_tetC+RedcolTetCount[i];
    }
    
    std::vector<int> ie_tetCnt(nrow,0);

    int tett=0;
    int pris=0;
    for(int i=0;i<nrow;i++)
    {
        ie_tetCnt[i] = -1;
        
        if(iet[i][0]==2) // Tet
        {
            ie_tetCnt[i] = OffcolTetCount[rank]+tett;
            //ie_tetCnt.setVal(i,0,OffcolTetCount[rank]+tett);
            tett++;
        }
        else
        {
            pris++;
        }
    }
    //std::cout << "before partitioning rank = " << rank << " #tets = " << tett << " #prisms " << pris << std::endl;
    std::vector<int> elTypes(4);
    int check_tet_glob, check_pri_glob, check_pyr_glob, check_hex_glob;
    MPI_Allreduce(&check_tet, &check_tet_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pri, &check_pri_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pyr, &check_pyr_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_hex, &check_hex_glob, 1, MPI_INT, MPI_SUM, comm);

    elTypes[0] = check_tet_glob;
    elTypes[1] = check_pri_glob;
    elTypes[2] = check_pyr_glob;
    elTypes[3] = check_hex_glob;
    
    delete[] OffcolTetCount;

    // std::vector<int> iet_glob(nglob,0);
    // MPI_Allreduce(iet_loc.data(), iet_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    // iet_loc.clear();

    //==========================================================================
    //==========================================================================
    //==========================================================================

    int nrow_fglob  = face_rankInfo[0];
    int nrow_floc   = face_rankInfo[1];
    int offset_floc = face_rankInfo[2];
    // int ncol_ifn    = 4;
    int ncol_ife    = 2;
    
    // ParMatrix<int> ifn_copy        = new ParMatrix<int>(nrow_fglob,ncol_ifn,comm);
    // ParMatrix<int> ife_copy        = new ParMatrix<int>(nrow_fglob,ncol_ife,comm);
    // ParMatrix<int> if_ref_copy     = new ParMatrix<int>(nrow_fglob,1,comm);
    // ParMatrix<int> if_Nv_copy      = new ParMatrix<int>(nrow_fglob,1,comm);
    
    std::map<int,std::vector<int> > ifn_copy;
    std::map<int,std::vector<int> > ife_copy;
    std::map<int,std::vector<int> > if_ref_copy;
    std::map<int,std::vector<int> > if_Nv_copy;    
    int fref2;
    int index_range = 0;
    
    for(i=0;i<nrow_floc;i++)
    {
        int fid        = offset_floc+i; 
        int fref       = ifn[i][7];
        int index      = i + offset_floc;
        int fzone      = ProvideBoundaryID(index,ranges_id);
        int fref2      = zone2ref[fzone];
        int ncol_ifn   = ifn[i][0];
        if((fref2!=fref))
        {
            std::cout << "mapping ranges are wrong"  << std::endl;
        }

        std::vector<int> ifn_copy_row(ncol_ifn);
        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy_row[j] = ifn[i][j+1]-1;
        }

        std::vector<int> ife_copy_row(ncol_ife);
        for(j=0;j<ncol_ife;j++)
        {
            ife_copy_row[j] = ife[i][j]-1;
        }

        if_ref_copy[fid].push_back(ifn[i][7]);
        if_Nv_copy[fid].push_back(ifn[i][0]);

        ifn_copy[fid] = ifn_copy_row;
        ife_copy[fid] = ife_copy_row;
    }
    ife.clear();
    ifn.clear();

    
    std::map<int,std::vector<double> > xcn_copy;
    int nloc_verts   = vert_rankInfo[1];
    int offset_verts = vert_rankInfo[2];
    for(int i=0;i<nloc_verts;i++)
    {
        int vid = offset_verts+i;
        //std::cout << "vud " << vid << " " << rank << std::endl;
        
        std::vector<double> coords(3);
        for(int j=0;j<3;j++)
        {
            coords[j] = xcn[i][j];
        }
        xcn_copy[vid] = coords;
    }
    xcn.clear();
    

    
    mRead->nElem          = nElem;
    mRead->nFace          = nFace;
    mRead->nVert          = nVert;
    mRead->xcn            = xcn_copy;
    mRead->ien            = ien_copy;
    mRead->ief            = ief_copy;
    mRead->iee            = iee_copy;
    mRead->iet            = iet_copy;
    mRead->elTypes        = elTypes;
    mRead->ie_Nv          = ie_Nv_copy;
    mRead->ie_Nf          = ie_Nf_copy;
    mRead->if_Nv          = if_Nv_copy;
    mRead->ifn            = ifn_copy;
    mRead->if_ref         = if_ref_copy;
    mRead->ife            = ife_copy;
    mRead->ie_tetCnt      = ie_tetCnt;
    mRead->interior       = interior;
    mRead->ghost          = ghost_out;
    mRead->zone2bcref     = zone2bcref;
    mRead->element2rank   = element2rank_glob;
    mRead->bref2zone      = bref2zone;
    mRead->znames_map     = znames_map;
    mRead->znames_map_inv = znames_map_inv;  
    mRead->zone2name      = zone2name;
    mRead->znames         = znames;
    mRead->zdefs          = zdefs;
    mRead->ranges_id      = ranges_id;
    mRead->ranges_ref     = ranges_ref;
    mRead->ntetra         = tetCount;
    mRead->nprism         = priCount;
    mRead->nhexahedral    = hexCount;
    mRead->npyramid       = pyrCount;

    //std::cout << "ReadUS3DMeshData " << tetCount << std::endl;
    /**/
    //delete zdefs;
    //delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;

    return mRead;

    
}






mesh* ReadUS3DMesh(const char* fn_conn, const char* fn_grid, int readFromStats, int StateVar, MPI_Comm comm, MPI_Info info)
{
    mesh* mRead = new mesh;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<double> >   xcn = ReadDataSetFromFileInParallel_Lite<double>(fn_grid,"xcn",comm,info);
    std::vector<std::vector<int> >      ien = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ien",comm,info);
    std::vector<std::vector<int> >      ief = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ief",comm,info);
    std::vector<std::vector<int> >      iee = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"iee",comm,info);
    std::vector<std::vector<int> >      iet = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"iet",comm,info);
    std::vector<std::vector<int> >      ifn = ReadDataSetFromFileInParallel_Lite<int>(fn_grid,"ifn",comm,info);
    std::vector<std::vector<int> >      ife = ReadDataSetFromFileInParallel_Lite<int>(fn_conn,"ife",comm,info);
    std::vector<int> vert_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_grid,"xcn",comm,info);
    std::vector<int> elem_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ien",comm,info);
    std::vector<int> face_rankInfo          = ReadDataSetSizeFromFileInParallel(fn_conn,"ife",comm,info);

    int nElem = elem_rankInfo[0];
    int nFace = face_rankInfo[0];
    int nVert = vert_rankInfo[0];
    

    int offset_el = elem_rankInfo[2];

    // if(rank == 0)
    // {
    // std::map<std::string,int>::iterator itv;
    // for(itv=var_map.begin();itv!=var_map.end();itv++)
    // {
    //    std::cout << itv->first << " " << itv->second << std::endl;
    // }

    // }

  

    int Nel     = elem_rankInfo[0];
    int Nel_loc = elem_rankInfo[1];
    
    std::map<int, std::vector<double> > interior;
    std::map<int, std::vector<double> > ghost_out;

    std::vector<std::vector<int> > zdefs    = ReadDataSetFromGroupFromFile_Lite<int>(fn_grid,"zones","zdefs");
    std::vector<std::vector<char> >  znames = ReadDataSetFromGroupFromFile_Lite<char>(fn_grid,"zones","znames");
    double rhoGhostState,uGhostState,vGhostState,wGhostState,TGhostState,VtotGhostState;
    


    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    std::vector<int> low_range;
    std::vector<int> high_range;
    std::vector<int> ref_range;
    
    int t=0;
    int gg=0;
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,std::vector<int> > ranges_ref;
    std::map<int,int> zone2ref;
    
    for(int i=2;i<zdefs.size();i++)
    {
        // std::cout << zdefs[i][2] << " " << zdefs[i][3] << " " << zdefs[i][4] << " " << zdefs[i][5] << std::endl;
        bnd_m.push_back(zdefs[i][5]);
        
        low_range.push_back(zdefs[i][3]);
        high_range.push_back(zdefs[i][4]);
        ref_range.push_back(zdefs[i][5]);
        
        zone2ref[zdefs[i][2]] = zdefs[i][5];

        std::vector<int> ra(2);
        ra[0] = zdefs[i][3]-1;
        ra[1] = zdefs[i][4]-1;
        ranges_id[zdefs[i][2]] = ra;
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][3]-1);
        // ranges_id[zdefs[i][2]].push_back(zdefs[i][4]-1);
        ranges_ref[zdefs[i][5]]=ra;
        
        // gg++;
    }
    
    
    bnd_m.push_back(zdefs[zdefs.size()-1][4]);
   
    if(rank == 0)
    {
       PlotBoundaryData_Lite(znames,zdefs);
    }
    
    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;
    std::map<int,std::vector<int> > bref2zone;
    std::map<int,char*> zone2name;

    std::map<int,int> zone2bcref;
    for(int i=0;i<zdefs.size();i++)
    {
        if(zdefs[i][5]!=1)
        {
            std::string name;
            
            char* namechar = new char[znames[0].size()];
                        
            for(int j=0;j<znames[0].size();j++)
            {
               namechar[j]=znames[i][j];

               char ch = znames[i][j];
               char chref = ' ';
               if(ch!=chref)
               {
                   name.push_back(ch);
               }
            }
            if(i>2)
            {
                bref2zone[zdefs[i][5]].push_back(i);
                zone2bcref[zdefs[i][2]] = zdefs[i][5];
                znames_map[zdefs[i][2]] = name;
                znames_map_inv[name]    = zdefs[i][2];
                zone2name[zdefs[i][2]]  = namechar;
            }
            
        }
    }
    
    // number of vertices

    std::vector<std::vector<char> > znames_new(znames.size());
    
    std::vector<char> znames_new_row0(znames[0].size());

    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row0[j] = znames[0][j];
    }
    znames_new[0] = znames_new_row0;

    std::vector<char> znames_new_row1(znames[0].size());
    for(int j=0;j<znames[0].size();j++)
    {
        znames_new_row1[j] = znames[1][j];
    }
    znames_new[1] = znames_new_row1;

    
    std::map<int,std::string>::iterator itch;
    int c=2;
    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
    {
        std::vector<char> znames_new_row(znames[0].size());
        int bid = itch->first;
        for(int j=0;j<itch->second.size();j++)
        {
            znames_new_row[j] = itch->second[j];
        }
        znames_new[c] = znames_new_row;
        c++;
    }




    //==========================================================================
    //==========================================================================
    //==========================================================================    
    int i,j;

    int nglob     = elem_rankInfo[0];
    int nrow      = elem_rankInfo[1];
    //

    std::map<int,std::vector<int> > ien_copy;
    std::map<int,std::vector<int> > ief_copy;
    std::map<int,std::vector<int> > iee_copy;
    std::map<int,int> iet_copy;
    std::map<int,int> ie_Nv_copy;
    std::map<int,int> ie_Nf_copy;


    int ncol_ien = ien[0].size()-1;
    int ncol_ief = ief[0].size()-1;
    int ncol_iee = iee[0].size()-1;
    //
    std::vector<int> element2rank_loc(nglob,0);
    //std::vector<int> iet_loc(nglob,0);
    int check_hex = 0;
    int check_tet = 0;
    int check_pyr = 0;
    int check_pri = 0;

    int tetCount  = 0;
    int hexCount  = 0;
    int pyrCount  = 0;
    int priCount  = 0;
        
    for(i=0;i<nrow;i++)
    {
        int elid = offset_el+i;

        element2rank_loc[elid] = rank;
        iet_copy[elid]         = iet[i][0];

        if(iet[i][0]==2) // Tet
        {
            ie_Nv_copy[elid] = 4;
            ie_Nf_copy[elid] = 4;
            check_tet        = check_tet+1;
            ncol_ien         = 4;
            ncol_iee         = 4;
            ncol_ief         = 4;
            tetCount++;
        }
        if(iet[i][0]==4) // Hex
        {
            ie_Nv_copy[elid] = 8;
            ie_Nf_copy[elid] = 6;
            ncol_ien         = 8;
            ncol_ief         = 6;
            ncol_iee         = 6;
            check_hex        = check_hex+1;
            hexCount++;
        }
        if(iet[i][0]==5) // Pyramid
        {
           ie_Nv_copy[elid]  = 5;
           ie_Nf_copy[elid]  = 5;
           ncol_ien          = 5;
           ncol_ief          = 5;
           ncol_iee          = 5;
           check_pyr         = check_pyr+1;
           pyrCount++;
        }
        if(iet[i][0]==6) // Prism
        {
            ie_Nv_copy[elid] = 6;
            ie_Nf_copy[elid] = 5;
            ncol_ien         = 6;
            ncol_ief         = 5;
            ncol_iee         = 5;
            check_pri        = check_pri+1;
            priCount++;
        }

        std::vector<int> ien_copy_row(ncol_ien);
        for(int j=0;j<ncol_ien;j++)
        {
            ien_copy_row[j] = ien[i][j+1]-1;

            if(ien[i][j+1]-1 > vert_rankInfo[0])
            {
                std::cout << "Raar while reading" << std::endl;
            }
        }
        
        std::vector<int> ief_copy_row(ncol_ief);
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy_row[j] = fabs(ief[i][j+1])-1;
        }

        std::vector<int> iee_copy_row(ncol_iee);
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy_row[j] = iee[i][j+1]-1;
        }

        ien_copy[elid] = ien_copy_row;
        ief_copy[elid] = ief_copy_row;
        iee_copy[elid] = iee_copy_row;

        // if(iet[i][0]!=2 && iet[i][0]!=4 && iet[i][0]!=6)
        // {
        //     std::cout << "Warning: this mesh has pyramids! " << iet[i][0] << std::endl;
        // }

    }

    ien.clear();
    iee.clear();
    iet.clear();
    ief.clear();

    std::vector<int> element2rank_glob(nglob,0);
    MPI_Allreduce(element2rank_loc.data(), element2rank_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    element2rank_loc.clear();

     int* colTetCount    = new int[size];
    int* RedcolTetCount = new int[size];
    int* OffcolTetCount = new int[size];

    for(int i=0;i<size;i++)
    {
        colTetCount[i]    = 0;
        RedcolTetCount[i] = 0;
        if(i==rank)
        {
            colTetCount[i] = tetCount;
        }
    }
    
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    
    int offset_tetC = 0;
    for(int i=0;i<size;i++)
    {
        OffcolTetCount[i] = offset_tetC;
        offset_tetC = offset_tetC+RedcolTetCount[i];
    }
    
    std::vector<int> ie_tetCnt(nrow,0);

    int tett=0;
    int pris=0;
    for(int i=0;i<nrow;i++)
    {
        ie_tetCnt[i] = -1;
        
        if(iet[i][0]==2) // Tet
        {
            ie_tetCnt[i] = OffcolTetCount[rank]+tett;
            //ie_tetCnt.setVal(i,0,OffcolTetCount[rank]+tett);
            tett++;
        }
        else
        {
            pris++;
        }
    }
    //std::cout << "before partitioning rank = " << rank << " #tets = " << tett << " #prisms " << pris << std::endl;
    std::vector<int> elTypes(4);
    int check_tet_glob, check_pri_glob, check_pyr_glob, check_hex_glob;
    MPI_Allreduce(&check_tet, &check_tet_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pri, &check_pri_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_pyr, &check_pyr_glob, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&check_hex, &check_hex_glob, 1, MPI_INT, MPI_SUM, comm);

    elTypes[0] = check_tet_glob;
    elTypes[1] = check_pri_glob;
    elTypes[2] = check_pyr_glob;
    elTypes[3] = check_hex_glob;
    
    delete[] OffcolTetCount;

    // std::vector<int> iet_glob(nglob,0);
    // MPI_Allreduce(iet_loc.data(), iet_glob.data(), nglob, MPI_INT, MPI_SUM, comm);
    // iet_loc.clear();

    //==========================================================================
    //==========================================================================
    //==========================================================================

    int nrow_fglob  = face_rankInfo[0];
    int nrow_floc   = face_rankInfo[1];
    int offset_floc = face_rankInfo[2];
    // int ncol_ifn    = 4;
    int ncol_ife    = 2;
    
    // ParMatrix<int> ifn_copy        = new ParMatrix<int>(nrow_fglob,ncol_ifn,comm);
    // ParMatrix<int> ife_copy        = new ParMatrix<int>(nrow_fglob,ncol_ife,comm);
    // ParMatrix<int> if_ref_copy     = new ParMatrix<int>(nrow_fglob,1,comm);
    // ParMatrix<int> if_Nv_copy      = new ParMatrix<int>(nrow_fglob,1,comm);
    
    std::map<int,std::vector<int> > ifn_copy;
    std::map<int,std::vector<int> > ife_copy;
    std::map<int,std::vector<int> > if_ref_copy;
    std::map<int,std::vector<int> > if_Nv_copy;    
    int fref2;
    int index_range = 0;
    
    for(i=0;i<nrow_floc;i++)
    {
        int fid        = offset_floc+i; 
        int fref       = ifn[i][7];
        int index      = i + offset_floc;
        int fzone      = ProvideBoundaryID(index,ranges_id);
        int fref2      = zone2ref[fzone];
        int ncol_ifn   = ifn[i][0];
        if((fref2!=fref))
        {
            std::cout << "mapping ranges are wrong"  << std::endl;
        }

        std::vector<int> ifn_copy_row(ncol_ifn);
        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy_row[j] = ifn[i][j+1]-1;
        }

        std::vector<int> ife_copy_row(ncol_ife);
        for(j=0;j<ncol_ife;j++)
        {
            ife_copy_row[j] = ife[i][j]-1;
        }

        if_ref_copy[fid].push_back(ifn[i][7]);
        if_Nv_copy[fid].push_back(ifn[i][0]);

        ifn_copy[fid] = ifn_copy_row;
        ife_copy[fid] = ife_copy_row;
    }
    ife.clear();
    ifn.clear();

    
    std::map<int,std::vector<double> > xcn_copy;
    int nloc_verts   = vert_rankInfo[1];
    int offset_verts = vert_rankInfo[2];
    for(int i=0;i<nloc_verts;i++)
    {
        int vid = offset_verts+i;
        //std::cout << "vud " << vid << " " << rank << std::endl;
        
        std::vector<double> coords(3);
        for(int j=0;j<3;j++)
        {
            coords[j] = xcn[i][j];
        }
        xcn_copy[vid] = coords;
    }
    xcn.clear();
    

    
    mRead->nElem          = nElem;
    mRead->nFace          = nFace;
    mRead->nVert          = nVert;
    mRead->xcn            = xcn_copy;
    mRead->ien            = ien_copy;
    mRead->ief            = ief_copy;
    mRead->iee            = iee_copy;
    mRead->iet            = iet_copy;
    mRead->elTypes        = elTypes;
    mRead->ie_Nv          = ie_Nv_copy;
    mRead->ie_Nf          = ie_Nf_copy;
    mRead->if_Nv          = if_Nv_copy;
    mRead->ifn            = ifn_copy;
    mRead->if_ref         = if_ref_copy;
    mRead->ife            = ife_copy;
    mRead->ie_tetCnt      = ie_tetCnt;
    mRead->interior       = interior;
    mRead->ghost          = ghost_out;
    mRead->zone2bcref     = zone2bcref;
    mRead->element2rank   = element2rank_glob;
    mRead->bref2zone      = bref2zone;
    mRead->znames_map     = znames_map;
    mRead->znames_map_inv = znames_map_inv;  
    mRead->zone2name      = zone2name;
    mRead->znames         = znames;
    mRead->zdefs          = zdefs;
    mRead->ranges_id      = ranges_id;
    mRead->ranges_ref     = ranges_ref;
    mRead->ntetra         = tetCount;
    mRead->nprism         = priCount;
    mRead->nhexahedral    = hexCount;
    mRead->npyramid       = pyrCount;

    /**/
    //delete zdefs;
    //delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;

    return mRead;

    
}





std::vector<int> ReadDataSetSizeFromFileInParallel(const char* file_name, const char* dataset_name, MPI_Comm comm, MPI_Info info)
{
    
    // Get the size of the process;
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //std::cout << rank << " " << size << std::endl;
    
    hid_t acc_tpl1          = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret;
    //herr_t ret            = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    //herr_t ret            = H5Pset_dxpl_mpio(,comm,info);
    acc_tpl1                = H5P_DEFAULT;
    // Open file and data set to get dimensions of array;
    
    hid_t file_id           = H5Fopen(file_name, H5F_ACC_RDONLY,H5P_DEFAULT);
    hid_t dset_id           = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
    hid_t dspace            = H5Dget_space(dset_id);
    int ndims               = H5Sget_simple_extent_ndims(dspace);
    
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int nrow                = dims[0];
    int ncol                = dims[1];
    int N                   = nrow;
    int nloc                = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset              = rank*int(N/size) + MIN(rank, N%size);



    H5Sclose(dspace);

    H5Dclose(dset_id);
    H5Fclose(file_id);

    std::vector<int> comm_info_entity(5);
    comm_info_entity[0] = N;
    comm_info_entity[1] = nloc;
    comm_info_entity[2] = offset;
    comm_info_entity[3] = size;
    comm_info_entity[4] = rank;

    return comm_info_entity;
}








US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, int readFromStats, int StateVar, MPI_Comm comm, MPI_Info info)
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

    Array<char>* zvnames = ReadDataSetFromGroupInGroupFromFile<char>(fn_data,"info","solver","svnames");

    
    std::map<string,int> var_map = PlotVariableNames(zvnames);
    
    if(rank == 0)
    {
  	std::map<std::string,int>::iterator itv;
	for(itv=var_map.begin();itv!=var_map.end();itv++)
	{
	   std::cout << itv->first << " " << itv->second << std::endl;
	}

    }

    int uid = var_map["u"];
    int vid = var_map["v"];
    int wid = var_map["w"];
    int Tid = var_map["T"];
  
    if(rank == 0)
    { 
    
    std::cout <<"variable IDs " << var_map.size() << " " << uid << " " << vid << " " << wid << " " << Tid << std::endl;
    }
    int Nel_loc = ien->getNrow();

    int Nel = ien->getNglob();
    ParArray<double>* interior;

    if(readFromStats==1)
    {
        ParArray<double>* mean;
        ParArray<double>* stats;
        
        double time_stats = ReadStatisticsTimeFromRunInFileInParallel(fn_data,"run_1",comm,info);
        
        interior = new ParArray<double>(Nel,2,comm);
        mean  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","stats-mean",0,Nel,comm,info);

        double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
        
        if(rank == 0)
        {
            std::cout << "Statistics time = " << time_stats << std::endl;
        }
        Array<double>* vel_mean = new Array<double>(Nel_loc,3);
        for(int u=0;u<Nel_loc;u++)
        {
            
            rhoState = mean->getVal(u,0)/time_stats;
            uState   = mean->getVal(u,uid)/time_stats;
            vState   = mean->getVal(u,vid)/time_stats;
            wState   = mean->getVal(u,wid)/time_stats;
            TState   = mean->getVal(u,Tid)/time_stats;
            vel_mean->setVal(u,0,uState);
            vel_mean->setVal(u,1,vState);
            vel_mean->setVal(u,2,wState);
            VtotState = sqrt(uState*uState+vState*vState+wState*wState);
            //aState   = sqrt(1.4*287.05*TState);
            //aState   = sqrt(1.29*188.92*TState);
            //MState = VtotState/aState;
            MState = TState;
            
            
 	    if(StateVar==0)
	    {
            	interior->setVal(u,1,MState);
	    }
	    if(StateVar==1)
	    {
		interior->setVal(u,1,TState);
	    }
	    //std::cout << "rhoState" << rhoState << " uState " << uState << " vState " << vState << " wState " << wState << " TState " << TState << " MState " << MState << std::endl;   
        }
        
        stats  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","stats-stat",0,Nel,comm,info);
        //std::cout << "stats size " << stats->getNrow() << " " << stats->getNcol() << std::endl;
        double upup,vpvp,wpwp,tke;
        std::vector<double> tkevec(Nel_loc);
        for(int u=0;u<Nel_loc;u++)
        {
            upup = stats->getVal(u,0)/time_stats-vel_mean->getVal(u,0)*vel_mean->getVal(u,0);
            vpvp = stats->getVal(u,1)/time_stats-vel_mean->getVal(u,1)*vel_mean->getVal(u,1);
            wpwp = stats->getVal(u,2)/time_stats-vel_mean->getVal(u,2)*vel_mean->getVal(u,2);
            tke = 0.5*(upup+vpvp+wpwp);
            tkevec[u] = tke;
        }
        delete vel_mean;
        
        double tkeMax = *std::max_element(tkevec.begin(),tkevec.end());
        double tkeMax_glob = 0.0;
        MPI_Allreduce(&tkeMax, &tkeMax_glob, 1, MPI_DOUBLE, MPI_MAX, comm);
        //std::cout << "tkeMax_glob " << tkeMax_glob << std::endl;
        for(int u=0;u<Nel_loc;u++)
        {
            interior->setVal(u,0,tkevec[u]/tkeMax_glob);
        }
        
        
        delete mean;
        delete stats;
        
    }
    if(readFromStats==0)
    {
        Array<double>* readdata  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
        
        interior = new ParArray<double>(Nel,1,comm);
        double rhoState,uState,vState,wState,TState,VtotState,aState,MState;

        for(int u=0;u<Nel_loc;u++)
        {
            rhoState = readdata->getVal(u,0);
            uState   = readdata->getVal(u,uid);
            vState   = readdata->getVal(u,vid);
            wState   = readdata->getVal(u,wid);
            TState   = readdata->getVal(u,Tid);
            VtotState = sqrt(uState*uState+vState*vState+wState*wState);
//            aState   = sqrt(1.4*287.05*TState);
//            MState  = VtotState/aState;
            interior->setVal(u,0,TState);
        }
        
        delete readdata;
    }
    
    

    
    
    Array<double>* ghost        = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);

    
    Array<int>*    zdefs        = ReadDataSetFromGroupFromFile<int>(fn_grid,"zones","zdefs");
    Array<char>*  znames        = ReadDataSetFromGroupFromFile<char>(fn_grid,"zones","znames");
    
    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    std::vector<int> low_range;
    std::vector<int> high_range;
    std::vector<int> ref_range;
    
    int t=0;
    int gg=0;
    std::map<int,std::vector<int> > ranges_id;
    std::map<int,int> zone2ref;
    for(int i=2;i<zdefs->getNrow();i++)
    {
        bnd_m.push_back(zdefs->getVal(i,5));
        
        low_range.push_back(zdefs->getVal(i,3));
        high_range.push_back(zdefs->getVal(i,4));
        ref_range.push_back(zdefs->getVal(i,5));
        
        zone2ref[zdefs->getVal(i,2)] = zdefs->getVal(i,5);
        ranges_id[zdefs->getVal(i,2)].push_back(zdefs->getVal(i,3)-1);
        ranges_id[zdefs->getVal(i,2)].push_back(zdefs->getVal(i,4)-1);
        
        if(rank == 0)
        {
        	std::cout << zdefs->getVal(i,2) << " range " << zdefs->getVal(i,3)-1 << ", " << zdefs->getVal(i,4)-1 << std::endl;
        }
        
        gg++;
    }
    bnd_m.push_back(zdefs->getVal(zdefs->getNrow()-1,4));
    
    if(rank == 0)
    {
       PlotBoundaryData(znames,zdefs);
    }
    
    std::map<int,std::string> znames_map;
    std::map<std::string,int> znames_map_inv;
    std::map<int,std::vector<int> > bref2zone;
    std::map<int,char*> zone2name;
    Array<char>* znames_new = new Array<char>(znames->getNrow(),znames->getNcol());
    std::map<int,int> zone2bcref;
    for(int i=0;i<zdefs->getNrow();i++)
    {
        if(zdefs->getVal(i,5)!=1)
        {
        	std::string name;
        	
        	char* namechar = new char[znames->getNcol()];
			            
            for(int j=0;j<znames->getNcol();j++)
            {
 			   namechar[j]=znames->getVal(i,j);

               char ch = znames->getVal(i,j);
               char chref = ' ';
               if(ch!=chref)
               {
            	   name.push_back(ch);
               }
            }
            if(i>2)
            {
            	bref2zone[zdefs->getVal(i,5)].push_back(i);
            	zone2bcref[zdefs->getVal(i,2)] = zdefs->getVal(i,5);
            	znames_map[zdefs->getVal(i,2)] = name;
            	znames_map_inv[name] = zdefs->getVal(i,2);
            	zone2name[zdefs->getVal(i,2)]= namechar;
            }
            
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
    
    std::map<int,std::string>::iterator itch;
    int c=2;
    for(itch=znames_map.begin();itch!=znames_map.end();itch++)
    {
        int bid = itch->first;
        for(int j=0;j<itch->second.size();j++)
        {
            znames_new->setVal(c,j,itch->second[j]);
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
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
        }
    }
    delete ief;
    
    int ncol_iee = iee->getNcol()-1;
    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
    
    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
        }
    }
    
    delete iee;
    
    int nrow_fglob  = ifn->getNglob();
    int nrow_floc   = ifn->getNrow();
    int ncol_ifn    = 4;
    int ncol_ife    = 2;
    
    ParArray<int>* ifn_copy        = new ParArray<int>(nrow_fglob,ncol_ifn,comm);
    ParArray<int>* ife_copy        = new ParArray<int>(nrow_fglob,ncol_ife,comm);
    ParArray<int>* if_ref_copy     = new ParArray<int>(nrow_fglob,1,comm);
    ParArray<int>* if_Nv_copy      = new ParArray<int>(nrow_fglob,1,comm);
    
    int foffset                    = if_ref_copy->getOffset(rank);
    int fref2;
    
    int index_range = 0;
    
    for(i=0;i<nrow_floc;i++)
    {
        int fref       = ifn->getVal(i,7);
        int index      = i + foffset;
        int fzone      = ProvideBoundaryID(index,ranges_id);
//        if(fzone == 3)
//        {
//        	//std::cout << "Its actually " << fzone << std::endl;
//        	
//        	fzone = 2;
//        }
        int fref2	   = zone2ref[fzone];
        
        if((fref2!=fref))
        {
        	std::cout << "mapping ranges are wrong"  << std::endl;
        }
        
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
        
        //if(fref!=3 && fref!=13 && fref!=36 && fref!=7&&fref!=2)
        //{
        //    std::cout << "its wrong here already " << fref << std::endl;
        //}
    }
    
    
//    for(int q=0;q<if_ref_copy->getNrow();q++)
//    {
//        if(if_ref_copy->getVal(q,0)!=3 && if_ref_copy->getVal(q,0)!=7 && if_ref_copy->getVal(q,0)!=10 && if_ref_copy->getVal(q,0)!=36 && if_ref_copy->getVal(q,0)!=2)
//        {
//            std::cout<<"while reading TEFFIE " << if_ref_copy->getVal(q,0) << std::endl;
//        }
//    }
    
    
    delete ifn;
    delete ife;
    
    ParArray<int>* ie_Nv    = new ParArray<int>(nglob,1,comm);
    ParArray<int>* ie_Nf    = new ParArray<int>(nglob,1,comm);
    
    int check_hex = 0;
    int check_tet = 0;
    int check_pyr = 0;
    int check_pri = 0;
    
    int tetCount=0;
    for(int i=0;i<nrow;i++)
    {
        if(iet->getVal(i,0)==2) // Tet
        {
            ie_Nv->setVal(i,0,4);
            ie_Nf->setVal(i,0,4);
            check_tet = 1;
            tetCount++;
        }
        if(iet->getVal(i,0)==4) // Hex
		{
			ie_Nv->setVal(i,0,8);
			ie_Nf->setVal(i,0,6);
			check_hex = 1;
		}
        if(iet->getVal(i,0)==5) // Pyramid
	    {
		   ie_Nv->setVal(i,0,5);
		   ie_Nf->setVal(i,0,5);
		   check_pyr = 1;
	    }
        if(iet->getVal(i,0)==6) // Prism
        {
            ie_Nv->setVal(i,0,6);
            ie_Nf->setVal(i,0,5);
            check_pri = 1;
        }
        
        // if(iet->getVal(i,0)!=2 && iet->getVal(i,0)!=4 && iet->getVal(i,0)!=6)
        // {
        //     std::cout << "Warning: this mesh has pyramids! " << iet->getVal(i,0) << std::endl;
        // }
    }
    
    int* colTetCount    = new int[size];
    int* RedcolTetCount = new int[size];
    int* OffcolTetCount = new int[size];

    for(int i=0;i<size;i++)
    {
        colTetCount[i]    = 0;
        RedcolTetCount[i] = 0;
        if(i==rank)
        {
            colTetCount[i] = tetCount;
        }
    }
    
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    
    int offset_tetC = 0;
    for(int i=0;i<size;i++)
    {
        OffcolTetCount[i] = offset_tetC;
        offset_tetC = offset_tetC+RedcolTetCount[i];
    }
    
    Array<int>* ie_tetCnt    = new Array<int>(nrow,1);

    int tett=0;
    int pris=0;
    for(int i=0;i<nrow;i++)
    {
        ie_tetCnt->setVal(i,0,-1);
        
        if(iet->getVal(i,0)==2) // Tet
        {
            ie_tetCnt->setVal(i,0,OffcolTetCount[rank]+tett);
            tett++;
        }
        else
        {
            pris++;
        }
    }
    //std::cout << "before partitioning rank = " << rank << " #tets = " << tett << " #prisms " << pris << std::endl;
    Array<int>* elTypes = new Array<int>(4,1);
    elTypes->setVal(0,0,check_tet);
    elTypes->setVal(1,0,check_pri);
    elTypes->setVal(2,0,check_pyr);
    elTypes->setVal(3,0,check_hex);
    
    delete[] OffcolTetCount;
    
    us3d->xcn            = xcn;
    us3d->elTypes        = elTypes;
    us3d->ien            = ien_copy;
    us3d->ief            = ief_copy;
    us3d->iee            = iee_copy;
    us3d->iet            = iet;
    us3d->ie_Nv          = ie_Nv;
    us3d->ie_Nf          = ie_Nf;
    us3d->if_Nv          = if_Nv_copy;

    us3d->ifn            = ifn_copy;
    us3d->if_ref         = if_ref_copy;
    us3d->ife            = ife_copy;
    us3d->ie_tetCnt      = ie_tetCnt;
    us3d->interior       = interior;
    us3d->ghost          = ghost;
    us3d->zone2bcref	 = zone2bcref;
    
    us3d->bref2zone		 = bref2zone;
    us3d->znames_map	 = znames_map;
    us3d->znames_map_inv = znames_map_inv;
                	
    us3d->zone2name		 = zone2name;
    us3d->znames         = znames;
    us3d->zdefs          = zdefs;
    us3d->ranges_id		 = ranges_id;
    //delete zdefs;
    //delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;
    return us3d;
}







US3D* ReadUS3DGrid(const char* fn_conn, const char* fn_grid, int readFromStats, MPI_Comm comm, MPI_Info info)
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

    int Nel_loc = ien->getNrow();

    int Nel = ien->getNglob();
    ParArray<double>* interior = new ParArray<double>(1,1,comm);
    Array<double>* ghost = new Array<double>(1,1);

    
    Array<int>*    zdefs        = ReadDataSetFromGroupFromFile<int>(fn_grid,"zones","zdefs");
    Array<char>*  znames        = ReadDataSetFromGroupFromFile<char>(fn_grid,"zones","znames");
    
    
    
    std::map<int,std::vector<int> > bnd_face_map;
    // Collect boundary data;
    std::vector<int> bnd_m;
    std::vector<int> low_range;
    std::vector<int> high_range;
    std::vector<int> ref_range;
    
    int t=0;
    int gg=0;
    std::map<int,std::vector<int> > ranges;
    for(int i=2;i<zdefs->getNrow();i++)
    {
        bnd_m.push_back(zdefs->getVal(i,5));
        
        low_range.push_back(zdefs->getVal(i,3)-1);
        high_range.push_back(zdefs->getVal(i,4)-1);
        ref_range.push_back(zdefs->getVal(i,5));
        ranges[zdefs->getVal(i,5)].push_back(zdefs->getVal(i,3));
        ranges[zdefs->getVal(i,5)].push_back(zdefs->getVal(i,4));
//        if(rank == 0)
//        {
//            std::cout << zdefs->getVal(i,5) << " " << zdefs->getVal(i,3)-1 << " " << zdefs->getVal(i,4)-1 << std::endl;
//        }
        gg++;
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
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
        }
    }
    delete ief;
    
    int ncol_iee = iee->getNcol()-1;
    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
    
    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
        }
    }
    
    delete iee;
    
    int nrow_fglob  = ifn->getNglob();
    int nrow_floc   = ifn->getNrow();
    int ncol_ifn    = 4;
    int ncol_ife    = 2;
    ParArray<int>* ifn_copy    = new ParArray<int>(nrow_fglob,ncol_ifn,comm);
    ParArray<int>* ife_copy    = new ParArray<int>(nrow_fglob,ncol_ife,comm);
    ParArray<int>* if_ref_copy = new ParArray<int>(nrow_fglob,1,comm);
    ParArray<int>* if_Nv_copy  = new ParArray<int>(nrow_fglob,1,comm);
    int foffset                = if_ref_copy->getOffset(rank);
    int fref2;
    
    int index_range = 0;
    
    for(i=0;i<nrow_floc;i++)
    {
        int fref   = ifn->getVal(i,7);
        int index  = i + foffset;
        //int fref2  = ProvideBoundaryRef(index,ranges,fref,rank);
        
//        if(fref!=fref2)
//        {
//            fref=fref2;
//        }
        
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
        //if(fref!=3 && fref!=13 && fref!=36 && fref!=7&&fref!=2)
        //{
        //    std::cout << "its wrong here already " << fref << std::endl;
        //}
    }
    
    
//    for(int q=0;q<if_ref_copy->getNrow();q++)
//    {
//        if(if_ref_copy->getVal(q,0)!=3 && if_ref_copy->getVal(q,0)!=7 && if_ref_copy->getVal(q,0)!=10 && if_ref_copy->getVal(q,0)!=36 && if_ref_copy->getVal(q,0)!=2)
//        {
//            std::cout<<"while reading TEFFIE " << if_ref_copy->getVal(q,0) << std::endl;
//        }
//    }
    
    
    delete ifn;
    delete ife;
    
    ParArray<int>* ie_Nv    = new ParArray<int>(nglob,1,comm);
    ParArray<int>* ie_Nf    = new ParArray<int>(nglob,1,comm);
    
    int check_hex = 0;
    int check_tet = 0;
    int check_pri = 0;
    
    int tetCount=0;
    for(int i=0;i<nrow;i++)
    {
        if(iet->getVal(i,0)==2) // Tet
        {
            ie_Nv->setVal(i,0,4);
            ie_Nf->setVal(i,0,4);
            check_tet = 1;
            tetCount++;
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
        if(iet->getVal(i,0)!=2 && iet->getVal(i,0)!=4 && iet->getVal(i,0)!=6)
        {
            std::cout << "What is this type " << iet->getVal(i,0) << std::endl;
        }
    }
    
    int* colTetCount = new int[size];
    int* RedcolTetCount = new int[size];
    int* OffcolTetCount = new int[size];

    for(int i=0;i<size;i++)
    {
        colTetCount[i]    = 0;
        RedcolTetCount[i] = 0;
        if(i==rank)
        {
            colTetCount[i] = tetCount;
        }
    }
    
    
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(colTetCount,  RedcolTetCount,  size, MPI_INT, MPI_SUM, comm);
    
    int offset_tetC = 0;
    for(int i=0;i<size;i++)
    {
        OffcolTetCount[i] = offset_tetC;
        offset_tetC = offset_tetC+RedcolTetCount[i];
    }
    
    Array<int>* ie_tetCnt    = new Array<int>(nrow,1);

    int tett=0;
    int pris=0;
    for(int i=0;i<nrow;i++)
    {
        ie_tetCnt->setVal(i,0,-1);
        
        if(iet->getVal(i,0)==2) // Tet
        {
            ie_tetCnt->setVal(i,0,OffcolTetCount[rank]+tett);
            tett++;
        }
        else
        {
            pris++;
        }
    }
    //std::cout << "before partitioning rank = " << rank << " #tets = " << tett << " #prisms " << pris << std::endl;
    Array<int>* elTypes = new Array<int>(3,1);
    elTypes->setVal(0,0,check_tet);
    elTypes->setVal(1,0,check_pri);
    elTypes->setVal(2,0,check_hex);
    
    delete[] OffcolTetCount;
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
    us3d->ie_tetCnt     = ie_tetCnt;
    us3d->interior      = interior;
    us3d->ghost         = ghost;
    
    us3d->znames        = znames_new;
    us3d->zdefs         = zdefs;

    //delete zdefs;
    delete znames;
    
    //std::cout << interior->getNrow() << " " << interior->getNcol() << std::endl;
    return us3d;
}





std::map<int, std::vector<double> > ReadHyperSolveHessianData(const char* hessian_data, int Nel, MPI_Comm comm, MPI_Info info)
{
    std::map<int, std::vector<double> > T_hessian;
     int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    T_hessian = ReadDataSetFromHyperSolveFileInParallel_Lite<double>(hessian_data,
                                                                    "temperature", 
                                                                    Nel, comm, info);

    return T_hessian;
                                                                                                      
}

