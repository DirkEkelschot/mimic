#include "adapt_io.h"
#include "adapt_recongrad.h"
#include "adapt_recongrad2.h"
#include "adapt_output.h"
#include "adapt_geometry.h"

#include "hex2tet.h"
int mpi_size, mpi_rank;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void OutputAdapt_Grid()
{
    
    std::ifstream fin;
    fin.open("xcn_mmg.dat");

    // Read the file row by row
    std::vector<double> row(3);
    std::vector<std::vector<double> > xcn_vec;
    int t=0;
    while(fin >> row[0] >> row[1] >> row[2])
    {
       xcn_vec.push_back(row);
       t++;
    }
    int L = xcn_vec.size();
    std::ifstream finhex;
    finhex.open("ien_mmg.dat");

    // Read the file row by row
    std::vector<int> rowTet(4);
    std::vector<std::vector<int> > arrTet;

    while(finhex >> rowTet[0] >> rowTet[1] >> rowTet[2] >> rowTet[3] )
    {
       arrTet.push_back(rowTet);
    }
    
    std::map<std::set<int>, int> face2id;
    std::set<std::set<int> > faces;
    std::map<int,std::vector<int> > face2node;
    std::map<int,std::vector<int> > element2face;
    std::map<int,std::vector<int> > face2element;
    int fid = 0;
    std::set<int> face0;
    std::set<int> face1;
    std::set<int> face2;
    std::set<int> face3;
    for(int i=1;i<=arrTet.size();i++)
    { face0.insert(arrTet[i-1][0]);face0.insert(arrTet[i-1][1]);face0.insert(arrTet[i-1][2]);
        if( faces.count(face0) != 1 )
        {
            faces.insert(face0);
            face2id[face0]=fid;
            face2node[fid].push_back(arrTet[i-1][0]);
            face2node[fid].push_back(arrTet[i-1][1]);
            face2node[fid].push_back(arrTet[i-1][2]);
            
            element2face[i-1].push_back(fid);
            face2element[fid].push_back(i-1);
            
            fid++;
        }
        else{
            element2face[i-1].push_back(face2id[face0]);
            face2element[face2id[face0]].push_back(i-1);
        }
    face1.insert(arrTet[i-1][1]);face1.insert(arrTet[i-1][2]);face1.insert(arrTet[i-1][3]);
        if( faces.count(face1) != 1 )
        {
            faces.insert(face1);
            face2id[face1]=fid;
            
            face2node[fid].push_back(arrTet[i-1][1]);
            face2node[fid].push_back(arrTet[i-1][2]);
            face2node[fid].push_back(arrTet[i-1][3]);
            
            element2face[i-1].push_back(fid);
            face2element[fid].push_back(i-1);
            
            fid++;
        }
        else{
            element2face[i-1].push_back(face2id[face1]);
            face2element[face2id[face1]].push_back(i-1);
        }
        
    face2.insert(arrTet[i-1][2]);face2.insert(arrTet[i-1][3]);face2.insert(arrTet[i-1][0]);
        if( faces.count(face2) != 1 )
        {
            faces.insert(face2);
            face2id[face2]=fid;
            face2node[fid].push_back(arrTet[i-1][2]);
            face2node[fid].push_back(arrTet[i-1][3]);
            face2node[fid].push_back(arrTet[i-1][0]);
            
            element2face[i-1].push_back(fid);
            face2element[fid].push_back(i-1);
            
            fid++;
        }
        else{
            element2face[i-1].push_back(face2id[face2]);
            face2element[face2id[face2]].push_back(i-1);
        }
    face3.insert(arrTet[i-1][3]);face3.insert(arrTet[i-1][0]);face3.insert(arrTet[i-1][1]);
        if( faces.count(face3) != 1 )
        {
            faces.insert(face3);
            face2id[face3]=fid;
            face2node[fid].push_back(arrTet[i-1][3]);
            face2node[fid].push_back(arrTet[i-1][0]);
            face2node[fid].push_back(arrTet[i-1][1]);
            
            element2face[i-1].push_back(fid);
            face2element[fid].push_back(i-1);
            
            fid++;
        }
        else{
            element2face[i-1].push_back(face2id[face3]);
            face2element[face2id[face3]].push_back(i-1);
        }
        
        face0.clear();
        face1.clear();
        face2.clear();
        face3.clear();
    }
    

    std::cout << "number of unique faces is: " << faces.size() << " " << face2node.size() << std::endl;
    
    std::map<int,std::vector<int> >::iterator it;
    
    for(it=face2element.begin();it!=face2element.end();it++)
    {
        std::cout << it->first << " -> ";
        std::vector<int>::iterator it2;
        for(it2=it->second.begin();it2!=it->second.end();it2++)
        {
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
    }
    
//    int ufaces = faces.size();
//    Array<int>* ifn_mmg = new Array<int>(faces.size(),3);
//    for(int i=0;i<ufaces;i++)
//    {
//        ifn_mmg->setVal(i,0,face2node[i][0]);
//        ifn_mmg->setVal(i,1,face2node[i][1]);
//        ifn_mmg->setVal(i,2,face2node[i][2]);
//    }
    
    //Output the new grid.h5 which has the new vertices and ifn map.
//    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
//    plist_id               = H5P_DEFAULT;
//    //H5Pset_fapl_mpio(plist_id, comm, info);
//    hid_t file_id = H5Fcreate("adapt_grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//    H5Pclose(plist_id);
//    hsize_t     dimsf[2];
//    dimsf[0] = arrTet;
//    dimsf[1] = arrTet;
//    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
//
//    hid_t dset_id = H5Dcreate(file_id, "xcn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(filespace);
//
//    hsize_t    count[2];              /* hyperslab selection parameters */
//    hsize_t    offset[2];
//    count[0] = dimsf[0];
//    count[1] = dimsf[1];
//    offset[0] = 0;
//    offset[1] = 0;
//    hid_t memspace = H5Screate_simple(2, count, NULL);
//
//    filespace = H5Dget_space(dset_id);
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    hid_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_mmg->data);
}


MMG5_pMesh ReadMMG_pMesh(US3D* us3d, MPI_Comm comm, MPI_Info info)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol mmgSol   = NULL;
    
    MMG3D_Init_mesh(MMG5_ARG_start,
    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
    MMG5_ARG_end);
    
    std::ifstream fin;
    fin.open("metric_restart.dat");

    // Read the file row by row
    std::vector<double> row(9);
    std::vector<std::vector<double> > Vmetric;
    int t=0;
    while(fin >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5] >> row[6] >> row[7] >> row[8])
    {
       Vmetric.push_back(row);
       t++;
    }
    int L = Vmetric.size();
    std::ifstream finhex;
    finhex.open("elements_restart.dat");

    // Read the file row by row
    std::vector<int> rowHex(8);
    std::vector<std::vector<int> > arrHex;

    while(finhex >> rowHex[0] >> rowHex[1] >> rowHex[2] >> rowHex[3] >> rowHex[4] >> rowHex[5] >> rowHex[6] >> rowHex[7])
    {
       arrHex.push_back(rowHex);
    }
    
    int nbHex      = arrHex.size();
    int nbVertices = Vmetric.size();
    int  nbTriangles = us3d->tria_ref_map.size()/2;
    if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,nbTriangles,0,0) != 1 )  exit(EXIT_FAILURE);
    
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
    
    for(int i=0;i<nbVertices;i++)
    {
        mmgMesh->point[i+1].c[0] = Vmetric[i][0];
        mmgMesh->point[i+1].c[1] = Vmetric[i][1];
        mmgMesh->point[i+1].c[2] = Vmetric[i][2];
        mmgMesh->point[i+1].ref  = 1;
        
        double m11 = Vmetric[i][3];
        double m12 = Vmetric[i][4];
        double m13 = Vmetric[i][5];
        double m22 = Vmetric[i][6];
        double m23 = Vmetric[i][7];
        double m33 = Vmetric[i][8];
        
        if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
    }
    
    int* hexTab = new int[9*(nbHex+1)];
    int ref = 0;
    for(int i=0;i<nbHex;i++)
    {
        int hexTabPosition = 9*(i+1);
        for(int j=0;j<8;j++)
        {
            //int val = ien->getVal(i,j+1);
            int val = arrHex[i][j];
            hexTab[hexTabPosition+j] = val;
        }
        hexTab[hexTabPosition+8] = ref;
    }
    
    int num = H2T_chkorient(mmgMesh,hexTab,nbHex);
    
    int* adjahex = NULL;
    adjahex = (int*)calloc(6*nbHex+7,sizeof(int));
    assert(adjahex);
    //
    if(!H2T_hashHexa(hexTab,adjahex,nbHex))
    {
        std::cout << "Error :: setting up the new adjacency for the hexes after reorientation." << std::endl;
    }

    Hedge        hed2;
    hed2.size  = 6*nbHex;
    hed2.hnxt  = 6*nbHex;
    hed2.nhmax = (int)(16*6*nbHex);
    hed2.item  = NULL;
    hed2.item  = (hedge*)calloc(hed2.nhmax+1,sizeof(hedge));

    for (int k=6*nbHex; k<hed2.nhmax; k++)
    {
        hed2.item[k].nxt = k+1;
    }
    int ret = H2T_cuthex(mmgMesh, &hed2, hexTab, adjahex, nbHex);
    // allocate boundary ids:
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Begin of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::set<int> tria0;
    std::set<int> tria1;
    std::set<int> tria2;
    std::set<int> tria3;
    int ref0,ref1,ref2,ref3;
    int tel = 0;
    std::set<std::set<int> > tria_unique;
    int offset_NE = (int)mmgMesh->ne/2;
    t = 1;
    for(int i=1;i<=offset_NE;i++)
    {
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
        if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
        {
            ref0 = us3d->tria_ref_map[tria0];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[1];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[2];
            mmgMesh->tria[t].ref  = ref0;
            t++;
        }
        
        tria1.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
        tria1.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
        tria1.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
        if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
        {
            ref1 = us3d->tria_ref_map[tria1];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[1];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
            mmgMesh->tria[t].ref  = ref1;
            t++;
        }
        
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        if(us3d->tria_ref_map.find(tria2)!=us3d->tria_ref_map.end())
        {
            ref2 = us3d->tria_ref_map[tria2];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[2];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].ref  = ref2;
            t++;
        }
        
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
        if(us3d->tria_ref_map.find(tria3)!=us3d->tria_ref_map.end())
        {
            ref3 = us3d->tria_ref_map[tria3];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[3];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[1];
            mmgMesh->tria[t].ref  = ref3;
            t++;
        }
        
        tria0.clear();
        tria1.clear();
        tria2.clear();
        tria3.clear();
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //End of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    MMG3D_Set_handGivenMesh(mmgMesh);
    
//    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 1.5) != 1 )
//    exit(EXIT_FAILURE);

    int ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);
    
    OutputMesh_MMG(mmgMesh);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,1);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,2);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,3);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,4);
    WriteUS3DGridFromMMG(mmgMesh, us3d);
//    std::map<int,std::vector<int> > ref2bface;
//    std::set<set<int> > bfaces;
//    std::set<int>face;
//    for(int i=1;i<=mmgMesh->nt;i++)
//    {
//        ref2bface[mmgMesh->tria[i].ref].push_back(i);
//        face.insert(mmgMesh->tria[i].v[0]);
//        face.insert(mmgMesh->tria[i].v[1]);
//        face.insert(mmgMesh->tria[i].v[2]);
//
//        bfaces.insert(face);
//
//        face.clear();
//    }
//
//
//    Array<double>* xcn_mmg = new Array<double>(mmgMesh->np,3);
//    for(int i=0;i<mmgMesh->np;i++)
//    {
//        xcn_mmg->setVal(i,0,mmgMesh->point[i+1].c[0]);
//        xcn_mmg->setVal(i,1,mmgMesh->point[i+1].c[1]);
//        xcn_mmg->setVal(i,2,mmgMesh->point[i+1].c[2]);
//    }
//
//    std::map<std::set<int>, int> facemap;
//    std::set<std::set<int> > faces;
//    std::map<int,std::vector<int> > element2face;
//    std::map<int,std::vector<int> > face2element;
//    std::map<int,std::vector<int> > face2node;
//    std::set<int> face0;
//    std::set<int> face1;
//    std::set<int> face2;
//    std::set<int> face3;
//    int fid = 0;
//    int vid0,vid1,vid2,vid3;
//    std::map<int,int> lh;
//    std::map<int,int> rh;
//    int of = 0;
//    Array<int>* adapt_iet = new Array<int>(mmgMesh->ne,1);
//    for(int i=1;i<=mmgMesh->ne;i++)
//    {
//        adapt_iet->setVal(i-1,0,2); // Element type = 2 since we are dealing with tetrahedra.
//
//        face0.insert(mmgMesh->tetra[i].v[1]);
//        face0.insert(mmgMesh->tetra[i].v[2]);
//        face0.insert(mmgMesh->tetra[i].v[3]);
//
//        if(faces.count(face0) != 1 )
//        {
//            faces.insert(face0);
//            facemap[face0]=fid;
//            face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
//            element2face[i-1].push_back(fid);
//            face2element[fid].push_back(i-1);
//            lh[fid] = i-1;
//            fid++;
//        }
//        else
//        {
//            rh[facemap[face0]] = i-1;
//        }
//
//        face1.insert(mmgMesh->tetra[i].v[0]);
//        face1.insert(mmgMesh->tetra[i].v[3]);
//        face1.insert(mmgMesh->tetra[i].v[2]);
//
//        if(faces.count(face1) != 1)
//        {
//            faces.insert(face1);
//            facemap[face1]=fid;
//            face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
//            element2face[i-1].push_back(fid);
//            face2element[fid].push_back(i-1);
//            lh[fid] = i-1;
//            fid++;
//        }
//        else
//        {
//            rh[facemap[face1]] = i-1;
//        }
//
//        face2.insert(mmgMesh->tetra[i].v[0]);
//        face2.insert(mmgMesh->tetra[i].v[1]);
//        face2.insert(mmgMesh->tetra[i].v[3]);
//
//        if( faces.count(face2) != 1)
//        {
//            faces.insert(face2);
//            facemap[face2]=fid;
//            face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
//            element2face[i-1].push_back(fid);
//            face2element[fid].push_back(i-1);
//            lh[fid] = i-1;
//            fid++;
//        }
//        else
//        {
//            rh[facemap[face2]] = i-1;
//        }
//
//        face3.insert(mmgMesh->tetra[i].v[0]);
//        face3.insert(mmgMesh->tetra[i].v[2]);
//        face3.insert(mmgMesh->tetra[i].v[1]);
//
//        if( faces.count(face3) != 1)
//        {
//            faces.insert(face3);
//            facemap[face3]=fid;
//            face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
//            face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
//            element2face[i-1].push_back(fid);
//            face2element[fid].push_back(i-1);
//            lh[fid] = i-1;
//            fid++;
//        }
//        else
//        {
//            rh[facemap[face3]] = i-1;
//        }
//
//        face0.clear();
//        face1.clear();
//        face2.clear();
//        face3.clear();
//    }
//
//
//    std::map<int,int>::iterator itm;
//    int it;
//    Array<int>* adapt_ifn = new Array<int>(face2node.size(),8);
//    t = 0;
//    Array<int>* zdefs = us3d->zdefs;
//    int bc_id = 0;
//    for(itm=lh.begin();itm!=lh.end();itm++)
//    {
//        it = itm->first;
//
//        if(rh.find(it)==rh.end())
//        {
//            rh[it] = 0;
//        }
//        else
//        {
//            adapt_ifn->setVal(t,0,3);
//            adapt_ifn->setVal(t,1,face2node[it][0]);
//            adapt_ifn->setVal(t,2,face2node[it][1]);
//            adapt_ifn->setVal(t,3,face2node[it][2]);
//            adapt_ifn->setVal(t,4,0);
//            adapt_ifn->setVal(t,5,rh[it]+1);
//            adapt_ifn->setVal(t,6,itm->second+1);
//            adapt_ifn->setVal(t,7,us3d->zdefs->getVal(3+bc_id-1,5));
//            t++;
//
//        }
//    }
//    std::map<int,std::vector<int> >::iterator it_bref;
//    int faceid;
//    std::set<int> iface;
//    int nbound = 0;
//    for(it_bref=ref2bface.begin();it_bref!=ref2bface.end();it_bref++)
//    {
//        bc_id = it_bref->first;
//        for(int q=0;q<it_bref->second.size();q++)
//        {
//            faceid = it_bref->second[q];
//            adapt_ifn->setVal(t,0,3);
//            adapt_ifn->setVal(t,1,mmgMesh->tria[faceid].v[0]);
//            adapt_ifn->setVal(t,2,mmgMesh->tria[faceid].v[1]);
//            adapt_ifn->setVal(t,3,mmgMesh->tria[faceid].v[2]);
//            adapt_ifn->setVal(t,4,0);
//            iface.insert(mmgMesh->tria[faceid].v[0]);
//            iface.insert(mmgMesh->tria[faceid].v[1]);
//            iface.insert(mmgMesh->tria[faceid].v[2]);
//
//            fid=facemap[iface];
//
//            adapt_ifn->setVal(t,5,rh[fid]);
//            adapt_ifn->setVal(t,6,lh[fid]+1);
//            adapt_ifn->setVal(t,7,us3d->zdefs->getVal(3+bc_id-1,5));
//            iface.clear();
//            t++;
//        }
//        nbound = nbound + it_bref->second.size();
//    }
//
//    int nbo = ref2bface.size();
//    Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
//
//    // Collect node data (10) . Starting index-ending index Nodes
//    adapt_zdefs->setVal(0,0,10);
//    adapt_zdefs->setVal(0,1,-1);
//    adapt_zdefs->setVal(0,2,1);
//    adapt_zdefs->setVal(0,3,1);
//    adapt_zdefs->setVal(0,4,mmgMesh->np);
//    adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
//    adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
//    // Collect element data (12) . Starting index-ending index Element
//    adapt_zdefs->setVal(1,0,12);
//    adapt_zdefs->setVal(1,1,-1);
//    adapt_zdefs->setVal(1,2,2);
//    adapt_zdefs->setVal(1,3,1);
//    adapt_zdefs->setVal(1,4,mmgMesh->ne);
//    adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
//    adapt_zdefs->setVal(1,6,2);
//
//    // Collect internal face data (13) . Starting index-ending index internal face.
//    adapt_zdefs->setVal(2,0,13);
//    adapt_zdefs->setVal(2,1,-1);
//    adapt_zdefs->setVal(2,2, 3);
//    adapt_zdefs->setVal(2,3, 1);
//    adapt_zdefs->setVal(2,4,faces.size()-bfaces.size());
//    adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
//    adapt_zdefs->setVal(2,6,2);
//    // Collect boundary face data (13) . Starting index-ending index boundary face for each boundary ID.
//    int q  = 1;
//    int nb = 0;
//    int face_start = faces.size()-bfaces.size()+1;
//    int face_end;
//    std::map<int,std::vector<int> >::iterator itr;
//    for(itr=ref2bface.begin();itr!=ref2bface.end();itr++)
//    {
//        face_end = face_start+itr->second.size()-1;
//        adapt_zdefs->setVal(3+nb,0,13);
//        adapt_zdefs->setVal(3+nb,1,-1);
//        adapt_zdefs->setVal(3+nb,2,3+q);
//        adapt_zdefs->setVal(3+nb,3,face_start);
//        adapt_zdefs->setVal(3+nb,4,face_end);
//        adapt_zdefs->setVal(3+nb,5,us3d->zdefs->getVal(3+nb,5));
//        adapt_zdefs->setVal(3+nb,6,2);
//        face_start = face_end+1;
//        nb++;
//        q++;
//    }
//
//    //Output the new grid.h5 which has the new vertices and ifn map.
//    //===================================================================
//    //===================================================================
//    //===================================================================
//    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
//    plist_id               = H5P_DEFAULT;
//    //H5Pset_fapl_mpio(plist_id, comm, info);
//    hid_t file_id = H5Fcreate("adapt_grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//    H5Pclose(plist_id);
//    hid_t status;
//    hid_t att_space;
//    hid_t attr_id;
//
//    hsize_t dimsf_att = 1;
//    att_space = H5Screate_simple(1, &dimsf_att, NULL);
//    hid_t type =  H5Tcopy (H5T_C_S1);
//    ret = H5Tset_size (type, 14);
//    ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
//    attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    char stri[] = "US3D Grid File";
//    status = H5Awrite(attr_id, type, &stri);
//    H5Aclose(attr_id);
//
//    //====================================================================================
//    hsize_t     dimsf[2];
//    dimsf[0] = adapt_iet->getNrow();
//    dimsf[1] = adapt_iet->getNcol();
//    hid_t filespace = H5Screate_simple(2, dimsf, NULL);
//
//    hid_t dset_id = H5Dcreate(file_id, "iet", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(filespace);
//
//    hsize_t    count[2];              // hyperslab selection parameters
//    hsize_t    offset[2];
//    count[0] = dimsf[0];
//    count[1] = dimsf[1];
//    offset[0] = 0;
//    offset[1] = 0;
//    hid_t memspace = H5Screate_simple(2, count, NULL);
//
//    filespace = H5Dget_space(dset_id);
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_iet->data);
//
//    //====================================================================================
//    // Add ifn map to the grid.h5 file
//    //====================================================================================
//
//    dimsf[0] = adapt_ifn->getNrow();
//    dimsf[1] = adapt_ifn->getNcol();
//    filespace = H5Screate_simple(2, dimsf, NULL);
//
//    dset_id = H5Dcreate(file_id, "ifn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(filespace);             // hyperslab selection parameters
//    count[0] = dimsf[0];
//    count[1] = dimsf[1];
//    offset[0] = 0;
//    offset[1] = 0;
//    memspace = H5Screate_simple(2, count, NULL);
//
//    filespace = H5Dget_space(dset_id);
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_ifn->data);
//    //====================================================================================
//
//    hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//
//
//    hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//    dimsf_att = 1;
//    att_space = H5Screate_simple(1, &dimsf_att, NULL);
//    attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    int value = mmgMesh->ne;
//    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
//    H5Aclose(attr_id);
//    attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    value = faces.size();
//    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
//    H5Aclose(attr_id);
//    attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    value = mmgMesh->nt;
//    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
//    H5Aclose(attr_id);
//    attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    value = mmgMesh->np;
//    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
//    H5Aclose(attr_id);
//
//
//
//    //====================================================================================
//    // Add iet map to the grid.h5 file
//    //====================================================================================
//
//    dimsf[0] = xcn_mmg->getNrow();
//    dimsf[1] = xcn_mmg->getNcol();
//    filespace = H5Screate_simple(2, dimsf, NULL);
//
//    dset_id = H5Dcreate(file_id, "xcn", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(filespace);
//    count[0] = dimsf[0];
//    count[1] = dimsf[1];
//    offset[0] = 0;
//    offset[1] = 0;
//    memspace = H5Screate_simple(2, count, NULL);
//
//    filespace = H5Dget_space(dset_id);
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_mmg->data);
//
//
//    // Create group;
//    //====================================================================================
//    hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//    // Add attribute to group:
//    //====================================================================================
//    dimsf_att = 1;
//    att_space = H5Screate_simple(1, &dimsf_att, NULL);
//    attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
//    value = 3+nbo;
//    status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
//    H5Aclose(attr_id);
//    //====================================================================================
//
//
//    // Add dataset to group:
//    //====================================================================================
//    dimsf[0] = adapt_zdefs->getNrow();
//    dimsf[1] = adapt_zdefs->getNcol();
//    filespace = H5Screate_simple(2, dimsf, NULL);
//    hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//    H5Sclose(filespace);
//
//    count[0]  = dimsf[0];
//    count[1]  = dimsf[1];
//    offset[0] = 0;
//    offset[1] = 0;
//    memspace  = H5Screate_simple(2, count, NULL);
//    filespace = H5Dget_space(dset_zdefs_id);
//
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
//    //====================================================================================
//
//    dimsf_att = us3d->znames->getNrow();
//    filespace = H5Screate_simple(1, &dimsf_att, NULL);
//    type =  H5Tcopy (H5T_C_S1);
//    ret  = H5Tset_size (type, 20);
//    ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
//    hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//
//    H5Sclose(filespace);
//
//    hsize_t cnt = us3d->znames->getNrow();
//
//    memspace  = H5Screate_simple(1, &cnt, NULL);
//    filespace = H5Dget_space(dset_znames_id);
//
//    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
//
//    status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);

    return mmgMesh;
}



MMG5_pMesh GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, Array<double>* Hv, Array<double>* Mv, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol mmgSol   = NULL;
    
    int nlElem = us3d->ien->getNrow();
    int nElem = us3d->ien->getNglob();
    int nvg = us3d->xcn->getNglob();
    
    
    ParallelState* xcn_pstate = P->getXcnParallelState();
    ParallelState* ien_pstate = P->getIenParallelState();
    
    std::map<int,int> gV2lV = P->getGlobalVert2LocalVert();

    std::vector<double> Uvloc;
    std::vector<int> vids;
    std::vector<double> Uids;
    std::vector<double> Hids;
    int gid,lid;
    int nval = 6;
    std::set<int> vdone;
    
    for(i=0;i<nlElem;i++)
    {
       for(j=0;j<8;j++)
       {
           gid = us3d->ien->getVal(i,j);
           lid = gV2lV[gid];
           int lrange = xcn_pstate->getOffsets()[world_rank];
           int hrange = xcn_pstate->getOffsets()[world_rank]+xcn_pstate->getNlocs()[world_rank];
           if(vdone.find(lid)==vdone.end() && (lrange<gid<hrange))
           {
               vdone.insert(lid);
               vids.push_back(gid);

               Uids.push_back(Mv->getVal(lid,0));
               Uids.push_back(Mv->getVal(lid,1));
               Uids.push_back(Mv->getVal(lid,2));
               Uids.push_back(Mv->getVal(lid,4));
               Uids.push_back(Mv->getVal(lid,5));
               Uids.push_back(Mv->getVal(lid,8));
               
               Hids.push_back(Hv->getVal(lid,0));
               Hids.push_back(Hv->getVal(lid,1));
               Hids.push_back(Hv->getVal(lid,2));
               
               Hids.push_back(Hv->getVal(lid,3));
               Hids.push_back(Hv->getVal(lid,4));
               Hids.push_back(Hv->getVal(lid,5));
               
               Hids.push_back(Hv->getVal(lid,6));
               Hids.push_back(Hv->getVal(lid,7));
               Hids.push_back(Hv->getVal(lid,8));
           }
       }
    }
   
    int* vglob_nloc     = new int[world_size];
    int* vglob_nlocs    = new int[world_size];
    int* vglob_offsets  = new int[world_size];
   
    int* Uvglob_nloc    = new int[world_size];
    int* Uvglob_nlocs   = new int[world_size];
    int* Uvglob_offsets = new int[world_size];
    
    int* Hvglob_nloc    = new int[world_size];
    int* Hvglob_nlocs   = new int[world_size];
    int* Hvglob_offsets = new int[world_size];
   

    for(i=0;i<world_size;i++)
    {
       if(i == world_rank)
       {
          vglob_nloc[i] = vids.size();
       }
       else
       {
           vglob_nloc[i] = 0;
       }
    }
   
    MPI_Allreduce(vglob_nloc,  vglob_nlocs, world_size, MPI_INT, MPI_SUM, comm);
   
    int offset  = 0;
    int vglobt  = 0;
    int Uvglobt = 0;
    int Hvglobt = 0;
    int Uoffset = 0;
    int Hoffset = 0;
    for(i=0;i<world_size;i++)
    {
        vglob_offsets[i]=offset;
        offset=offset+vglob_nlocs[i];
        vglobt = vglobt+vglob_nlocs[i];
       
        Uvglob_nlocs[i] = vglob_nlocs[i]*nval;
       
        Uvglob_offsets[i]=Uoffset;
        Uoffset=Uoffset+Uvglob_nlocs[i];
        Uvglobt = Uvglobt+Uvglob_nlocs[i];
        
        Hvglob_nlocs[i] = vglob_nlocs[i]*9;
        
        Hvglob_offsets[i]=Hoffset;
        Hoffset=Hoffset+Hvglob_nlocs[i];
        Hvglobt = Hvglobt+Hvglob_nlocs[i];
    }
   
    int* vids_t = 0;
    double* Uvids_t = 0;
    double* Hvids_t = 0;
    if(world_rank == 0)
    {
       vids_t  = new int[vglobt];
       Uvids_t = new double[Uvglobt];
       Hvids_t = new double[Hvglobt];
    }

    MPI_Gatherv(&vids[0],
                vids.size(),
                MPI_INT,
                &vids_t[0],
                vglob_nlocs,
                vglob_offsets,
                MPI_INT, 0, comm);
   
    MPI_Gatherv(&Uids[0],
                Uids.size(),
                MPI_DOUBLE,
                &Uvids_t[0],
                Uvglob_nlocs,
                Uvglob_offsets,
                MPI_DOUBLE, 0, comm);
    
    MPI_Gatherv(&Hids[0],
                Hids.size(),
                MPI_DOUBLE,
                &Hvids_t[0],
                Hvglob_nlocs,
                Hvglob_offsets,
                MPI_DOUBLE, 0, comm);
   
    Array<double>* Ug;
    Array<double>* Hg;
    if(world_rank == 0)
    {
        Ug = new Array<double>(us3d->xcn->getNglob(),nval);
        Hg = new Array<double>(us3d->xcn->getNglob(),9);
        int cid=0;
        for(i=0;i<world_size;i++)
        {
            for(j=0;j<vglob_nlocs[i];j++)
            {
                cid = vids_t[vglob_offsets[i]+j];
            
                for(int k=0;k<nval;k++)
                {
                    Ug->setVal(cid,k,Uvids_t[(vglob_offsets[i]+j)*nval+k]);
                }
                for(int k=0;k<nval;k++)
                {
                    Hg->setVal(cid,k,Hvids_t[(vglob_offsets[i]+j)*9+k]);
                }
            }
        }
    }
    else
    {
        Ug = new Array<double>(1,1);
        Hg = new Array<double>(1,1);
    }
   
    Array<double>* xcn_g;
    Array<int>* ien_g;
   
    if(world_rank == 0)
    {
        xcn_g = new Array<double>(nvg,3);
        ien_g = new Array<int>(nElem,8);
    }
    else
    {
        xcn_g = new Array<double>(1,1);
        ien_g = new Array<int>(1,1);
    }
   
    int* ien_nlocs      = new int[world_size];
    int* ien_offsets    = new int[world_size];
    int* xcn_nlocs      = new int[world_size];
    int* xcn_offsets    = new int[world_size];
   
    for(int i=0;i<world_size;i++)
    {
        
        xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]*3;
        xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;
       
        ien_nlocs[i]   = ien_pstate->getNlocs()[i]*8;
        ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
    }
   
    MPI_Gatherv(&us3d->xcn->data[0],
                us3d->xcn->getNrow()*3,
                MPI_DOUBLE,
                &xcn_g->data[0],
                xcn_nlocs,
                xcn_offsets,
                MPI_DOUBLE, 0, comm);
   

    MPI_Gatherv(&us3d->ien->data[0],
                us3d->ien->getNrow()*8,
                MPI_INT,
                &ien_g->data[0],
                ien_nlocs,
                ien_offsets,
                MPI_INT, 0, comm);

    if(world_rank == 0)
    {
        int nbHex = nElem;
        int nbVertices = nvg;
        int nbTriangles = us3d->tria_ref_map.size();
        
        std::cout << " number of tets = " << nbHex*6 << std::endl;
        
        std::cout << " number of triangle = " << nbTriangles<<  " " << us3d->tria_ref.size() << std::endl;
        
        MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
        MMG5_ARG_end);
        std::ofstream myfile4;
        myfile4.open("hexElements.dat");
        for(int i=0;i<nbHex;i++)
        {
            for(int j=0;j<8;j++)
            {
                myfile4 << ien_g->getVal(i,j) << " ";
            }
            myfile4 << std::endl;
        }
        myfile4.close();
        std::ofstream myfile5;
        myfile5.open("hexVertices.dat");
        for(int i=0;i<nbVertices;i++)
        {
            for(int j=0;j<3;j++)
            {
                myfile5 << xcn_g->getVal(i,j) << " ";
            }
            myfile5 << std::endl;
        }
        myfile5.close();
        
        if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
        
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
        
        std::ofstream myfile10;
        myfile10.open("metric.dat");
        std::ofstream myfile11;
        myfile11.open("hessian.dat");
        for(int i=0;i<xcn_g->getNrow();i++)
        {
            mmgMesh->point[i+1].c[0] = xcn_g->getVal(i,0);
            mmgMesh->point[i+1].c[1] = xcn_g->getVal(i,1);
            mmgMesh->point[i+1].c[2] = xcn_g->getVal(i,2);
            mmgMesh->point[i+1].ref = 1;
            
            double m11 = Ug->getVal(i,0);
            double m12 = Ug->getVal(i,1);
            double m13 = Ug->getVal(i,2);
            double m22 = Ug->getVal(i,3);
            double m23 = Ug->getVal(i,4);
            double m33 = Ug->getVal(i,5);
            
            myfile10 << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2) << Ug->getVal(i,0) << " " << Ug->getVal(i,1) << " " << Ug->getVal(i,2) << " " << Ug->getVal(i,3) << " " << Ug->getVal(i,4) << " " << Ug->getVal(i,5) << std::endl;
            
            myfile11 << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2) << Hg->getVal(i,0) << " " << Hg->getVal(i,1) << " " << Hg->getVal(i,2) << " " << Hg->getVal(i,3) << " " << Hg->getVal(i,4) << " " << Hg->getVal(i,5) << " " << Hg->getVal(i,6) << " " << Hg->getVal(i,7) << " " << Hg->getVal(i,8) << std::endl;
            
            if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
        }
        myfile10.close();
        int ref = 0;
        int* hexTab = new int[9*(nbHex+1)];
        for(int i=0;i<nbHex;i++)
        {
            int hexTabPosition = 9*(i+1);
            for(int j=0;j<8;j++)
            {
                //int val = ien->getVal(i,j+1);
                int val = ien_g->getVal(i,j)+1;
                hexTab[hexTabPosition+j] = val;
            }
            hexTab[hexTabPosition+8] = ref;
        }
        
        int num = H2T_chkorient(mmgMesh,hexTab,nbHex);
        
        int* adjahex = NULL;
        adjahex = (int*)calloc(6*nbHex+7,sizeof(int));
        assert(adjahex);
        //
        if(!H2T_hashHexa(hexTab,adjahex,nbHex))
        {
            std::cout << "Error :: setting up the new adjacency for the hexes after reorientation." << std::endl;
        }
        //
        //    //for(int i=0;i<6*nbHex+7;i++)
        //
        Hedge        hed2;
        hed2.size  = 6*nbHex;
        hed2.hnxt  = 6*nbHex;
        hed2.nhmax = (int)(16*6*nbHex);
        hed2.item  = NULL;
        hed2.item  = (hedge*)calloc(hed2.nhmax+1,sizeof(hedge));
        ////
        for (int k=6*nbHex; k<hed2.nhmax; k++)
        {
            hed2.item[k].nxt = k+1;
        }
        
        int ret = H2T_cuthex(mmgMesh, &hed2, hexTab, adjahex, nbHex);

        std::cout << "check faces triangles boundary "  << us3d->tria_ref.size() << std::endl;
        
        
        
        std::set<int> tria0;
        std::set<int> tria1;
        std::set<int> tria2;
        std::set<int> tria3;
        int ref0,ref1,ref2,ref3;
        int tel = 0;
        std::set<std::set<int> > tria_unique;
        int offset_NE = (int)mmgMesh->ne/2;
        int t = 1;
        for(int i=1;i<=mmgMesh->ne;i++)
        {
            tria0.insert(mmgMesh->tetra[offset_NE+i].v[0]);
            tria0.insert(mmgMesh->tetra[offset_NE+i].v[1]);
            tria0.insert(mmgMesh->tetra[offset_NE+i].v[2]);
            std::cout << std::endl;
            std::cout << mmgMesh->tetra[offset_NE+i].v[0] <<
            " " << mmgMesh->tetra[offset_NE+i].v[0] <<
            " " << mmgMesh->tetra[offset_NE+i].v[0] << std::endl;
            
            if(tria_unique.find(tria0)==tria_unique.end())
            {
                tria_unique.insert(tria0);
                ref0 = us3d->tria_ref_map[tria0];
//                mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
//                mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[1];
//                mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[2];
//                mmgMesh->tria[t].ref  = ref0;
                t++;
            }
            tria1.insert(mmgMesh->tetra[offset_NE+i].v[1]);
            tria1.insert(mmgMesh->tetra[offset_NE+i].v[2]);
            tria1.insert(mmgMesh->tetra[offset_NE+i].v[3]);
            
            std::cout << mmgMesh->tetra[offset_NE+i].v[1] <<
            " " << mmgMesh->tetra[offset_NE+i].v[2] <<
            " " << mmgMesh->tetra[offset_NE+i].v[3] << std::endl;
            
            if(tria_unique.find(tria1)==tria_unique.end())
            {
                tria_unique.insert(tria1);
                ref1 = us3d->tria_ref_map[tria1];
//                mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[1];
//                mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
//                mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
//                mmgMesh->tria[t].ref  = ref1;
                t++;
            }
            tria2.insert(mmgMesh->tetra[offset_NE+i].v[2]);
            tria2.insert(mmgMesh->tetra[offset_NE+i].v[3]);
            tria2.insert(mmgMesh->tetra[offset_NE+i].v[0]);

            std::cout << mmgMesh->tetra[offset_NE+i].v[2] <<
            " " << mmgMesh->tetra[offset_NE+i].v[3] <<
            " " << mmgMesh->tetra[offset_NE+i].v[0] << std::endl;
            

            
            if(tria_unique.find(tria2)==tria_unique.end())
            {
                tria_unique.insert(tria2);
                ref2 = us3d->tria_ref_map[tria2];
//                mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[2];
//                mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
//                mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[0];
//                mmgMesh->tria[t].ref  = ref2;
                t++;
            }
            tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]);
            tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]);
            tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]);
            
            std::cout << mmgMesh->tetra[offset_NE+i].v[3] <<
            " " << mmgMesh->tetra[offset_NE+i].v[0] <<
            " " << mmgMesh->tetra[offset_NE+i].v[1] << std::endl;
            
            if(tria_unique.find(tria3)==tria_unique.end())
            {
                tria_unique.insert(tria3);
                ref3 = us3d->tria_ref_map[tria3];
//                mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[3];
//                mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[0];
//                mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[1];
//                mmgMesh->tria[t].ref  = ref3;
                t++;
            }
//            if(ref0!=0 || ref1!=0 || ref2!=0 || ref3!=0)
//            {
//                std::cout << "tel = " << tel << " :: " << ref0 << " " << ref1 << " " << ref2 << " " << ref3 << std::endl;
//                tel++;
//            }
            tria0.clear();
            tria1.clear();
            tria2.clear();
            tria3.clear();
        }
         
//            std::cout << mmgMesh->tetra[offset_NE+i].v[0] << " " << mmgMesh->tetra[offset_NE+i].v[1] << " " << mmgMesh->tetra[offset_NE+i].v[2] << " " << mmgMesh->tetra[offset_NE+i].v[3] << " " << mmgMesh->tetra[offset_NE+i].ref <<  std::endl;
//            std::cout << mmgMesh->tetra[offset_NE+i+1].v[0] << " " << mmgMesh->tetra[offset_NE+i+1].v[1] << " " << mmgMesh->tetra[offset_NE+i+1].v[2] << " " << mmgMesh->tetra[offset_NE+i+1].v[3] << " " << mmgMesh->tetra[offset_NE+i+1].ref <<  std::endl;
//            std::cout << mmgMesh->tetra[offset_NE+i+2].v[0] << " " << mmgMesh->tetra[offset_NE+i+2].v[1] << " " << mmgMesh->tetra[offset_NE+i+2].v[2] << " " << mmgMesh->tetra[offset_NE+i+2].v[3] << " " << mmgMesh->tetra[offset_NE+i+2].ref <<  std::endl;
//            std::cout << mmgMesh->tetra[offset_NE+i+3].v[0] << " " << mmgMesh->tetra[offset_NE+i+3].v[1] << " " << mmgMesh->tetra[offset_NE+i+3].v[2] << " " << mmgMesh->tetra[offset_NE+i+3].v[3] << " " << mmgMesh->tetra[offset_NE+i+3].ref <<  std::endl;
//            std::cout << mmgMesh->tetra[offset_NE+i+4].v[0] << " " << mmgMesh->tetra[offset_NE+i+4].v[1] << " " << mmgMesh->tetra[offset_NE+i+4].v[2] << " " << mmgMesh->tetra[offset_NE+i+4].v[3] << " " << mmgMesh->tetra[offset_NE+i+4].ref << std::endl;
//            std::cout << mmgMesh->tetra[offset_NE+i+5].v[0] << " " << mmgMesh->tetra[offset_NE+i+5].v[1] << " " << mmgMesh->tetra[offset_NE+i+5].v[2] << " " << mmgMesh->tetra[offset_NE+i+5].v[3] << " " << mmgMesh->tetra[offset_NE+i+5].ref << std::endl;
//
        
        
        
        
        
        MMG3D_Set_handGivenMesh(mmgMesh);
            
        if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 3.5) != 1 )
        exit(EXIT_FAILURE);
        
        
        int ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);
        std::ofstream myfile2;
        myfile2.open("mmgMesh_v2.dat");
        myfile2 << "TITLE=\"new_volume.tec\"" << std::endl;
        myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        //std::cout << " verts check " << LVerts.size() << " " << hx.size() << std::endl;
        myfile2 <<"ZONE N = " << mmgMesh->np << ", E = " << mmgMesh->ne << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
        
        Array<double>* xcn_mmg = new Array<double>(mmgMesh->np,3);
        for(int i=0;i<mmgMesh->np;i++)
        {
            myfile2 << mmgMesh->point[i+1].c[0] << " " <<mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] <<  std::endl;
            
            xcn_mmg->setVal(i,0,mmgMesh->point[i+1].c[0]);
            xcn_mmg->setVal(i,1,mmgMesh->point[i+1].c[1]);
            xcn_mmg->setVal(i,2,mmgMesh->point[i+1].c[2]);
        }
        std::set<std::set<int> > faces;
        std::map<int,std::vector<int> > face2node;
        int fid = 0;
        std::set<int> face0;
        std::set<int> face1;
        std::set<int> face2;
        std::set<int> face3;
        for(int i=1;i<=mmgMesh->ne;i++)
        {
            myfile2 << mmgMesh->tetra[i].v[0] << " " << mmgMesh->tetra[i].v[1] << " " << mmgMesh->tetra[i].v[2] << " " << mmgMesh->tetra[i].v[3] << std::endl;
        face0.insert(mmgMesh->tetra[i].v[0]);face0.insert(mmgMesh->tetra[i].v[1]);face0.insert(mmgMesh->tetra[i].v[2]);
            if( faces.count(face0) != 1 )
            {
                faces.insert(face0);
                face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
                fid++;
            }
        face1.insert(mmgMesh->tetra[i].v[1]);face1.insert(mmgMesh->tetra[i].v[2]);face1.insert(mmgMesh->tetra[i].v[3]);
            if( faces.count(face1) != 1 )
            {
                faces.insert(face1);
                face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
                fid++;
            }
        face2.insert(mmgMesh->tetra[i].v[2]);face2.insert(mmgMesh->tetra[i].v[3]);face2.insert(mmgMesh->tetra[i].v[0]);
            if( faces.count(face2) != 1 )
            {
                faces.insert(face2);
                face2node[fid].push_back(mmgMesh->tetra[i].v[2]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
                fid++;
            }
        face3.insert(mmgMesh->tetra[i].v[3]);face3.insert(mmgMesh->tetra[i].v[0]);face3.insert(mmgMesh->tetra[i].v[1]);
            if( faces.count(face3) != 1 )
            {
                faces.insert(face3);
                face2node[fid].push_back(mmgMesh->tetra[i].v[3]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[0]);
                face2node[fid].push_back(mmgMesh->tetra[i].v[1]);
                fid++;
            }
            face0.clear();
            face1.clear();
            face2.clear();
            face3.clear();
        }
        myfile2.close();
        
        std::map<int, std::vector<int> >ref2face;
        for(int i=0;i<mmgMesh->nt;i++)
        {
            if(mmgMesh->tria[i].ref==2)
            {
              std::cout << "outcome2 = " << mmgMesh->tria[i].v[0] << " " << mmgMesh->tria[i].v[1] << " " << mmgMesh->tria[i].v[2] << " " << mmgMesh->tria[i].ref << std::endl;
            }
        }
        /*
        std::cout << "number of unique faces is: " << faces.size() << " " << face2node.size() << std::endl;
        int ufaces = faces.size();
        Array<int>* ifn_mmg = new Array<int>(faces.size(),3);
        for(int i=0;i<ufaces;i++)
        {
            ifn_mmg->setVal(i,0,face2node[i][0]);
            ifn_mmg->setVal(i,1,face2node[i][1]);
            ifn_mmg->setVal(i,2,face2node[i][2]);
        }
         */
    }
    //==================OUTPUT ORIGINAL MESH=======================
    //==================OUTPUT ORIGINAL MESH=======================
    /*
    if(world_rank == 0)
    {
        string filename = "d2UdX2i_global.dat";
        ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"gxx\", \"gxy\", \"gxz\", \"gyy\", \"gyz\", \"gzz\"" << std::endl;
        myfile <<"ZONE N = " << nvg << ", E = " << nElem << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
       
        string filename2 = "metric_restart.dat";
        ofstream myfile2;
        myfile2.open(filename2);
        
        string filename3 = "elements_restart.dat";
        ofstream myfile3;
        myfile3.open(filename3);
        
        for(int i=0;i<nvg;i++)
        {
            myfile << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
            << " " << Ug->getVal(i,0) << " " << Ug->getVal(i,1) << " " << Ug->getVal(i,2)
            << " " << Ug->getVal(i,3) << " " << Ug->getVal(i,4) << " " << Ug->getVal(i,5) << std::endl;
            
            myfile2 << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
            << " " << Ug->getVal(i,0) << " " << Ug->getVal(i,1) << " " << Ug->getVal(i,2)
            << " " << Ug->getVal(i,3) << " " << Ug->getVal(i,4) << " " << Ug->getVal(i,5) << std::endl;
        }
        for(int i=0;i<nElem;i++)
        {
            myfile << ien_g->getVal(i,0)+1 << " " <<
                      ien_g->getVal(i,1)+1 << " " <<
                      ien_g->getVal(i,2)+1 << " " <<
                      ien_g->getVal(i,3)+1 << " " <<
                      ien_g->getVal(i,4)+1 << " " <<
                      ien_g->getVal(i,5)+1 << " " <<
                      ien_g->getVal(i,6)+1 << " " <<
                      ien_g->getVal(i,7)+1 << std::endl;
            
            myfile3 << ien_g->getVal(i,0)+1 << " " <<
                       ien_g->getVal(i,1)+1 << " " <<
                       ien_g->getVal(i,2)+1 << " " <<
                       ien_g->getVal(i,3)+1 << " " <<
                       ien_g->getVal(i,4)+1 << " " <<
                       ien_g->getVal(i,5)+1 << " " <<
                       ien_g->getVal(i,6)+1 << " " <<
                       ien_g->getVal(i,7)+1 << std::endl;
        }
        myfile.close();
        myfile2.close();
        myfile3.close();
    }
    //==================OUTPUT ORIGINAL MESH=======================
    //==================OUTPUT ORIGINAL MESH=======================
    
    delete[] vglob_nloc;
    delete[] vglob_nlocs;
    delete[] vglob_offsets;
    
    delete[] Uvglob_nloc;
    delete[] Uvglob_nlocs;
    delete[] Uvglob_offsets;
    delete xcn_g;
    delete ien_g;
    */
    return mmgMesh;
}

int main(int argc, char** argv) {
    
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    //======================================================================================
    //====================Parsing the inputs from the command line==========================
    //======================================================================================
    std::string fn_grid_t;
    std::string fn_conn_t;
    std::string fn_data_t;
    std::string destination;
    
    char buffer[20];
    
    const char* fn_grid;
    const char* fn_conn;
    const char* fn_data;
    std::map<int,const char*> fnames;
    for(int i = 1;i<argc;i++)
    {
        std::string str = argv[i];
        std::size_t l = str.copy(buffer,6,0);
        buffer[l]='\0';
        
        if (str.compare(0,6,"--grid") == 0)
        {
            char fn_grid_n[20];
            std::size_t ln = str.copy(fn_grid_n,str.size(),7);
            fn_grid_n[ln]='\0';
            fn_grid = fn_grid_n;
            fnames[i]=fn_grid;
        }
        else if (str.compare(0,6,"--conn") == 0)
        {
            char fn_conn_n[20];
            std::size_t ln = str.copy(fn_conn_n,str.size(),7);
            fn_conn_n[ln]='\0';
            fn_conn = fn_conn_n;
            fnames[i]=fn_conn;
        }
        else if (str.compare(0,6,"--data") == 0)
        {
            char fn_data_n[20];
            std::size_t ln = str.copy(fn_data_n,str.size(),7);
            fn_data_n[ln]='\0';
            fn_data = fn_data_n;
            fnames[i]=fn_data;
        }
        
        
    }
    if(fnames.size()!=3)
    {
        if(world_rank == 0)
        {
            std::cout << "For metric determination, the required command line inputs are : --grid=\"grid_name\".h5 --conn=\"conn_name\".h5 --data=\"data_name\".h5." << std::endl;
            std::cout << "For now, reading existing metric data."  << std::endl;
            
            //OutputAdapt_Grid();
            
            
        }
        MPI_Finalize();
    }
    else
    {
        //========================================================================
        //========================================================================
        //========================================================================
        
        US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);
//
        MMG5_pMesh mmgMesh = ReadMMG_pMesh(us3d,comm,info);
//
        /*
        int Nel_part = us3d->ien->getNrow();
        
        ParallelState* ien_pstate = new ParallelState(us3d->ien->getNglob(),comm);

        ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);

        ParallelState* xcn_pstate = new ParallelState(us3d->xcn->getNglob(),comm);
        Array<double>* Uivar = new Array<double>(Nel_part,1);
        for(int i=0;i<Nel_part;i++)
        {
            Uivar->setVal(i,0,us3d->interior->getVal(i,0));
        }
        
        delete us3d->interior;
        
        Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, parmetis_pstate,ien_pstate, us3d->xcn,xcn_pstate,Uivar,comm);
        
        //Array<double>* dUdXi_v2 = ComputedUdx(P, pstate, iee_copy, iee_loc, ief_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm, ife_copy);
        
        std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uivar,comm);
        
        Mesh_Topology* meshTopo = new Mesh_Topology(P,us3d->ifn,us3d->ghost,UauxNew,us3d->bnd_map,us3d->nBnd,comm);
        //OutputBoundaryID(P, meshTopo, us3d, 0);
//        OutputBoundaryID(P, meshTopo, us3d, 1);
//        OutputBoundaryID(P, meshTopo, us3d, 2);
//        OutputBoundaryID(P, meshTopo, us3d, 3);
//        OutputBoundaryID(P, meshTopo, us3d, 4);
        
        //std::map<int,int> faceref = meshTopo->getFaceRef();

        //Gradients* dudxObj = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","MMG",comm);
        
        //Array<double>* dUdXi = dudxObj->getdUdXi();
        //clock_t t;
        double tmax = 0.0;
        double tn = 0.0;
        //t = clock();
        Array<double>* dUdXi   = ComputedUdx_MGG(P,UauxNew,meshTopo,us3d->ghost,comm);
//        Gradients* dudxObj   = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","LSQ",comm);
//        Array<double>* dUdXi = dudxObj->getdUdXi();
        
        //Array<double>* dUdXi = dudxObj->getdUdXi();

        Array<double>* dUidxi = new Array<double>(dUdXi->getNrow(),1);
        Array<double>* dUidyi = new Array<double>(dUdXi->getNrow(),1);
        Array<double>* dUidzi = new Array<double>(dUdXi->getNrow(),1);

        for(int i=0;i<dUdXi->getNrow();i++)
        {
            dUidxi->setVal(i,0,dUdXi->getVal(i,0));
            dUidyi->setVal(i,0,dUdXi->getVal(i,1));
            dUidzi->setVal(i,0,dUdXi->getVal(i,2));
        }
        
        //std::map<int,std::vector<double> > dUdxiauxNew = P->CommunicateAdjacentDataUS3DNew(dUdXi,comm);
    //    std::map<int,std::vector<double> > dUdxauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidxi,comm);
    //    std::map<int,std::vector<double> > dUdyauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidyi,comm);
    //    std::map<int,std::vector<double> > dUdzauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidzi,comm);
        
        std::map<int,double > dUdxauxNew  = P->CommunicateAdjacentDataUS3D(dUidxi,comm);
        std::map<int,double > dUdyauxNew  = P->CommunicateAdjacentDataUS3D(dUidyi,comm);
        std::map<int,double > dUdzauxNew  = P->CommunicateAdjacentDataUS3D(dUidzi,comm);
        
    //    std::map<int,std::vector<double> >::iterator itmv;
    //    int t=0;
    //    for(itmv=dUdxiauxNew.begin();itmv!=dUdxiauxNew.end();itmv++)
    //    {
    //        std::cout << itmv->second[0]-dUdxauxNew[itmv->first][0] << " " <<                  itmv->second[1]-dUdyauxNew[itmv->first][0]  << " " << itmv->second[2]-dUdzauxNew[itmv->first][0]  << std::endl;
    //        t++;
    //    }
        delete dUdXi;
        
        Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,us3d->ghost,comm);
        Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,us3d->ghost,comm);
        Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,us3d->ghost,comm);
        
        Array<double>* d2udx2 = new Array<double>(dU2dXi2->getNrow(),1);
        Array<double>* d2udxy = new Array<double>(dU2dXi2->getNrow(),1);
        Array<double>* d2udxz = new Array<double>(dU2dXi2->getNrow(),1);
        
        Array<double>* d2udyx = new Array<double>(dU2dYi2->getNrow(),1);
        Array<double>* d2udy2 = new Array<double>(dU2dYi2->getNrow(),1);
        Array<double>* d2udyz = new Array<double>(dU2dYi2->getNrow(),1);

        Array<double>* d2udzx = new Array<double>(dU2dZi2->getNrow(),1);
        Array<double>* d2udzy = new Array<double>(dU2dZi2->getNrow(),1);
        Array<double>* d2udz2 = new Array<double>(dU2dZi2->getNrow(),1);

        for(int i=0;i<dU2dXi2->getNrow();i++)
        {
            d2udx2->setVal(i,0,dU2dXi2->getVal(i,0));
            d2udxy->setVal(i,0,dU2dXi2->getVal(i,1));
            d2udxz->setVal(i,0,dU2dXi2->getVal(i,2));
            
            d2udyx->setVal(i,0,dU2dYi2->getVal(i,0));
            d2udy2->setVal(i,0,dU2dYi2->getVal(i,1));
            d2udyz->setVal(i,0,dU2dYi2->getVal(i,2));
            
            d2udzx->setVal(i,0,dU2dZi2->getVal(i,0));
            d2udzy->setVal(i,0,dU2dZi2->getVal(i,1));
            d2udz2->setVal(i,0,dU2dZi2->getVal(i,2));
        }
        
        std::vector<double> u_v    = meshTopo->ReduceUToVertices(Uivar);
        
        std::vector<double> dudx_v = meshTopo->ReduceUToVertices(dUidxi);
        std::vector<double> dudy_v = meshTopo->ReduceUToVertices(dUidyi);
        std::vector<double> dudz_v = meshTopo->ReduceUToVertices(dUidzi);
        
        std::vector<double> d2udx2_v = meshTopo->ReduceUToVertices(d2udx2);
        std::vector<double> d2udxy_v = meshTopo->ReduceUToVertices(d2udxy);
        std::vector<double> d2udxz_v = meshTopo->ReduceUToVertices(d2udxz);
        
        std::vector<double> d2udyx_v = meshTopo->ReduceUToVertices(d2udyx);
        std::vector<double> d2udy2_v = meshTopo->ReduceUToVertices(d2udy2);
        std::vector<double> d2udyz_v = meshTopo->ReduceUToVertices(d2udyz);
        
        std::vector<double> d2udzx_v = meshTopo->ReduceUToVertices(d2udzx);
        std::vector<double> d2udzy_v = meshTopo->ReduceUToVertices(d2udzy);
        std::vector<double> d2udz2_v = meshTopo->ReduceUToVertices(d2udz2);

        std::vector<Vert> Verts =  P->getLocalVerts();
        int nVerts = Verts.size();
        Array<double>* hessian = new Array<double>(nVerts,9);
        Array<double>* grad    = new Array<double>(nVerts,3);
        for(int i=0;i<nVerts;i++)
        {
            grad->setVal(i,0,dudx_v[i]);grad->setVal(i,1,dudy_v[i]);grad->setVal(i,2,dudz_v[i]);
            hessian->setVal(i,0,d2udx2_v[i]); hessian->setVal(i,1,d2udxy_v[i]); hessian->setVal(i,2,d2udxz_v[i]);
            hessian->setVal(i,3,d2udyx_v[i]); hessian->setVal(i,4,d2udy2_v[i]); hessian->setVal(i,5,d2udyz_v[i]);
            hessian->setVal(i,6,d2udzx_v[i]); hessian->setVal(i,7,d2udzy_v[i]); hessian->setVal(i,8,d2udz2_v[i]);
        }
        
        //Array<double>* UgRoot = GetSolutionOnRoot(P,us3d,hessian,comm);
        
        std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
        double max_v = *std::max_element(d2udx2_v.begin(), d2udx2_v.end());
        Array<double>* metric = ComputeMetric(Verts,grad,hessian,max_v);
        
//=================================================================
        //==================Output the data in Tecplot format==============
        //=================================================================
        
        string filename11 = "metric_rank_" + std::to_string(world_rank) + ".dat";
        ofstream myfile11;
        myfile11.open(filename11);
        string filename12 = "elements_rank_" + std::to_string(world_rank) + ".dat";
        ofstream myfile12;
        myfile11.open(filename12);
        string filename = "d2UdX2i_" + std::to_string(world_rank) + ".dat";
        ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
        
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"dx\", \"dy\", \"dz\", \"d2udx2\", \"d2udxy\", \"d2udxz\", \"d2udyx\", \"d2udy2\", \"d2udyz\", \"d2udzx\", \"d2udzy\", \"d2udz2\", \"M00\", \"M01\", \"M02\", \"M10\", \"M11\", \"M12\", \"M20\", \"M21\", \"M22\"" << std::endl;
        int nvert = Verts.size();
        myfile <<"ZONE N = " << nvert << ", E = " << us3d->ien->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
        
        std::map<int,int> gV2lV = P->getGlobalVert2LocalVert();
        
        for(int i=0;i<nvert;i++)
        {
            myfile<< Verts[i].x << " " << Verts[i].y << " " << Verts[i].z <<
            " " << u_v[i] << " "<< dudx_v[i] << " " << dudy_v[i] << " " << dudz_v[i] <<
            " " << hessian->getVal(i,0) << " " << hessian->getVal(i,1) << " " << hessian->getVal(i,2) <<
            " " << hessian->getVal(i,3) << " " << hessian->getVal(i,4) << " " << hessian->getVal(i,5) <<
            " " << hessian->getVal(i,6) << " " << hessian->getVal(i,7) << " " << hessian->getVal(i,8) <<
            " " << metric->getVal(i,0) << " " << metric->getVal(i,1) << " " << metric->getVal(i,2) <<
            " " << metric->getVal(i,3) << " " << metric->getVal(i,4) << " " << metric->getVal(i,5) <<
            " " << metric->getVal(i,6) << " " << metric->getVal(i,7) << " " << metric->getVal(i,8) << std::endl;
            
            myfile11 << Verts[i].x << " " << Verts[i].y << " " << Verts[i].z << " " << metric->getVal(i,0) << " " << metric->getVal(i,1) << " " << metric->getVal(i,2) << " " << metric->getVal(i,3) << " " << metric->getVal(i,4) << " " << metric->getVal(i,5) << " " << metric->getVal(i,6) << " " << metric->getVal(i,7) << " " << metric->getVal(i,8) << std::endl;
        }

        for(int i=0;i<us3d->ien->getNrow();i++)
        {
           myfile << loc_elem2verts_loc[i][0]+1 << " " <<
                     loc_elem2verts_loc[i][1]+1 << " " <<
                     loc_elem2verts_loc[i][2]+1 << " " <<
                     loc_elem2verts_loc[i][3]+1 << " " <<
                     loc_elem2verts_loc[i][4]+1 << " " <<
                     loc_elem2verts_loc[i][5]+1 << " " <<
                     loc_elem2verts_loc[i][6]+1 << " " <<
                     loc_elem2verts_loc[i][7]+1 << std::endl;
            
           myfile12 << loc_elem2verts_loc[i][0]+1 << " " <<
                       loc_elem2verts_loc[i][1]+1 << " " <<
                       loc_elem2verts_loc[i][2]+1 << " " <<
                       loc_elem2verts_loc[i][3]+1 << " " <<
                       loc_elem2verts_loc[i][4]+1 << " " <<
                       loc_elem2verts_loc[i][5]+1 << " " <<
                       loc_elem2verts_loc[i][6]+1 << " " <<
                       loc_elem2verts_loc[i][7]+1 << std::endl;
        }
        
        myfile.close();
        myfile11.close();
        myfile12.close();
        
        MMG5_pMesh mmgMesh = GetOptimizedMMG3DMeshOnRoot(P, us3d, hessian, metric, comm);
    
        delete d2udx2;
        delete d2udxy;
        delete d2udxz;

        delete d2udyx;
        delete d2udy2;
        delete d2udyz;

        delete d2udzx;
        delete d2udzy;
        delete d2udz2;


        d2udx2_v.erase(d2udx2_v.begin(),d2udx2_v.end());
        d2udxy_v.clear();
        d2udxz_v.clear();

        d2udyx_v.clear();
        d2udy2_v.clear();
        d2udyz_v.clear();

        d2udzx_v.clear();
        d2udzy_v.clear();
        d2udz2_v.clear();
        delete us3d->ifn;
        delete us3d->ien;
        delete us3d->iee;
        
        delete hessian;
        delete grad;
        delete ien_pstate;
        delete parmetis_pstate;
        //delete UgRoot;
        */
        MPI_Finalize();
        
    }
     
    return 0;
}
