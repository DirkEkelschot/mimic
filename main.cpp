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
    fin.open("restart/metric_restart.dat");

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
    finhex.open("restart/elements_restart.dat");

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


    return mmgMesh;
}



MMG_Mesh* GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, Array<double>* Hv, Array<double>* Mv, MPI_Comm comm)
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
        
        Hvglob_nlocs[i] = vglob_nlocs[i]*6;
        
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
        Hg = new Array<double>(us3d->xcn->getNglob(),6);
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
                    Hg->setVal(cid,k,Hvids_t[(vglob_offsets[i]+j)*6+k]);
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
        
        MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
        MMG5_ARG_end);
        if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,nbTriangles,0,0) != 1 )  exit(EXIT_FAILURE);
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
        string filename = "metric_restart.dat";
        ofstream myfile;
        myfile.open(filename);
        std::cout <<  xcn_g->getNrow() << std::endl; 
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
            
            myfile << mmgMesh->point[i+1].c[0] << " " << mmgMesh->point[i+1].c[0] << " " << mmgMesh->point[i+1].c[0] << " " << m11 << " " << m12 << " " << m13 << " " << m22 << " " << m23 << " " << m33 << std::endl;
            
            if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
        
        }
        myfile.close();
        
        int ref = 0;
        int* hexTab = new int[9*(nbHex+1)];
        for(int i=0;i<nbHex;i++)
        {
            int hexTabPosition = 9*(i+1);
            for(int j=0;j<8;j++)
            {
                int val = ien_g->getVal(i,j)+1;
                hexTab[hexTabPosition+j] = val;
            }
            hexTab[hexTabPosition+8] = ref;
        }
            
        int num = H2T_chkorient(mmgMesh,hexTab,nbHex);
        
        int* adjahex = NULL;
        adjahex = (int*)calloc(6*nbHex+7,sizeof(int));
        assert(adjahex);
        
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
        int t = 1;
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
        
//        Array<double>* xcn_mmg = new Array<double>(mmgMesh->np,3);
//        for(int i=0;i<mmgMesh->np;i++)
//        {
//            xcn_mmg->setVal(i,0,mmgMesh->point[i+1].c[0]);
//            xcn_mmg->setVal(i,1,mmgMesh->point[i+1].c[1]);
//            xcn_mmg->setVal(i,2,mmgMesh->point[i+1].c[2]);
//        }
    }
    
    //==================OUTPUT ORIGINAL MESH=======================
    //==================OUTPUT ORIGINAL MESH=======================
  
    
    if(world_rank == 0)
    {
        string filename2 = "metric_restart.dat";
        ofstream myfile2;
        myfile2.open(filename2);
        
        string filename3 = "elements_restart.dat";
        ofstream myfile3;
        myfile3.open(filename3);
        
        for(int i=0;i<nvg;i++)
        {
            myfile2 << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
            << " " << Ug->getVal(i,0) << " " << Ug->getVal(i,1) << " " << Ug->getVal(i,2)
            << " " << Ug->getVal(i,3) << " " << Ug->getVal(i,4) << " " << Ug->getVal(i,5) << std::endl;
        }
        for(int i=0;i<nElem;i++)
        {
            myfile3 << ien_g->getVal(i,0)+1 << " " <<
                       ien_g->getVal(i,1)+1 << " " <<
                       ien_g->getVal(i,2)+1 << " " <<
                       ien_g->getVal(i,3)+1 << " " <<
                       ien_g->getVal(i,4)+1 << " " <<
                       ien_g->getVal(i,5)+1 << " " <<
                       ien_g->getVal(i,6)+1 << " " <<
                       ien_g->getVal(i,7)+1 << std::endl;
        }
        myfile2.close();
        myfile3.close();
    }
    //==================OUTPUT ORIGINAL MESH=======================
    //==================OUTPUT ORIGINAL MESH=======================
    
    
    MMG_Mesh* mmg = new MMG_Mesh;
    
    mmg->mmgMesh = mmgMesh;
    mmg->mmgSol  = mmgSol;
    
    return mmg;
}

//ien = 6;
//ief = 6;
//ifn = 6*4;

Vec3D* ComputeOutwardNormal(Vert* Vijk, Vert* Vface, std::vector<Vert*> face)
{
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
    
    Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
    double orient0   = DotVec3D(r0,n0);
    
    if(orient0<0.0)
    {
        NegateVec3D(n0);
    }
    return n0;
}


double ComputeQuickVol(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];

  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];

  vol =   (c2[0]-c1[0]) * (ay*bz - az*by) \
        + (c2[1]-c1[1]) * (az*bx - ax*bz) \
        + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return vol;
}


int ChkHexorient(double* P, int* Pid) {
  int     changed,k,i;
  double  volref,volhex;
  int* Pnew = new int[8];
  volref = 1;
  double* c1 = new double[3];
  double* c2 = new double[3];
  double* c3 = new double[3];
  double* c4 = new double[3];
  c1[0] = P[0*8+0];c1[1] = P[0*8+1];c1[2] = P[0*8+2];
  c2[0] = P[1*8+0];c1[1] = P[1*8+1];c1[2] = P[1*8+2];
  c3[0] = P[3*8+0];c1[1] = P[3*8+1];c1[2] = P[3*8+2];
  c4[0] = P[4*8+0];c1[1] = P[4*8+1];c1[2] = P[4*8+2];
  /** check the orientability of the hexahedra : vol of tet p0 p1 p3 p4 */
  volhex = ComputeQuickVol(c1,c2,c3,c4);
  changed = 0;
  if ( volref*volhex < 0 )
  {
      changed = 1;
      Pnew[0] = Pid[0];
      Pnew[1] = Pid[3];
      Pnew[2] = Pid[2];
      Pnew[3] = Pid[1];
      Pnew[4] = Pid[4];
      Pnew[5] = Pid[7];
      Pnew[6] = Pid[6];
      Pnew[7] = Pid[5];
      
      Pid[0] = Pnew[0];
      Pid[1] = Pnew[1];
      Pid[2] = Pnew[2];
      Pid[3] = Pnew[3];
      Pid[4] = Pnew[4];
      Pid[5] = Pnew[5];
      Pid[6] = Pnew[6];
      Pid[7] = Pnew[7];
  }

  return changed;
}

std::map<int,std::set<int> > FindOuterShellBoundaryLayerMesh(int wall_id, int nLayer, US3D* us3d, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, ParallelState* xcn_pstate, ParallelState* ien_pstate, MPI_Comm comm)
{
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
    std::cout << "Determining outer shell of BL mesh..." << std::endl;
    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    double* Pijk = new double[8*3];
    int nb = us3d->bnd_face_map[wall_id].size();
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    //std::vector<Vert*> face_c;
    int local_face_id;
    
    int bvid,opposite_bvid;
    
    int glob_el_id = 0;
    
    for(int bf=0;bf<us3d->bnd_face_map[wall_id].size();bf++)
    {
        
        std::vector<int> layer;
        int bfaceid = us3d->bnd_face_map[wall_id][bf];
        int faceid  = bfaceid;
        int elid0   = us3d->ife->getVal(faceid,0);
        int elid1   = us3d->ife->getVal(faceid,1);

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
        for(int k=0;k<6;k++)
        {
            if(ief_g->getVal(elid_cur,k)==faceid)
            {
                local_face_id = k;
            }
        }

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
            int vid  = us3d->ifn->getVal(faceid,r);
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
        
        for(int c=0;c<nLayer;c++)
        {
            for(int k=0;k<8;k++)
            {
               loc_vid     = ien_g->getVal(elid_cur,k);
               Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
               Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
               Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
            }
        
            Vert* Vijk = ComputeCenterCoord(Pijk, 8);
            
            for(int k=0;k<6;k++)
            {
                int fid = ief_g->getVal(elid_cur,k);
                Vert* Vface2  = new Vert;
                
                std::vector<int> faceVert_IDs(4);
                std::vector<Vert*> face2;
                for(int r=0;r<4;r++)
                {
                    int vid  = us3d->ifn->getVal(fid,r);
                    
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
        
            int min_index  = std::min_element(dp.begin(),dp.end())-dp.begin();
            double min_val = *std::min_element(dp.begin(),dp.end());

            int fid_new    = ief_g->getVal(elid_cur,min_index);
            nbf            = dpvec[min_index];

            NegateVec3D(nbf);

            int gEl0=us3d->ife->getVal(fid_new,0);
            int gEl1=us3d->ife->getVal(fid_new,1);

            if(gEl0==elid_cur)
            {
                elid_next = gEl1;
            }
            else if(gEl1==elid_cur)
            {
                elid_next = gEl0;
            }
            layer.push_back(elid_next);
            
            int changed = ChkHexorient(Pijk,Pijk_id);
            
            if(c==nLayer-1)
            {
                outer_shell_faces.push_back(fid_new);
                
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,0));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,1));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,2));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,3));
                
                std::vector<int> fe(2);
                fe[0] = elid_cur;
                fe[1] = elid_next;
                
                outer_shell_elements.push_back(fe);
                if(changed!=0)
                {
                    std::cout << "elid_cur " << elid_cur << std::endl;
                }
            }
            
            elid_cur = elid_next;
        }
    }

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << " extracting outer shell BL mesh = " << duration << std::endl;
//    std::map<int,std::vector<int> >::iterator itt;
//    std::vector<int> elements;
//    for(itt=mesh_topology_bl->BLlayers.begin();itt!=mesh_topology_bl->BLlayers.end();itt++)
//    {
//        for(int q=0;q<itt->second.size();q++)
//        {
//            elements.push_back(itt->second[q]);
//        }
//    }
//    OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_");
    return outer_shell_Faces2Nodes;
}

Mesh_Topology_BL* ExtractBoundaryLayerMesh(int wall_id, int nLayer, US3D* us3d, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, ParallelState* xcn_pstate, ParallelState* ien_pstate, MPI_Comm comm)
{
    Mesh_Topology_BL* mesh_topology_bl = new Mesh_Topology_BL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    std::vector<double> dp(6);
    std::vector<Vec3D*> dpvec(6);

    std::cout << "Extracting BL mesh" << std::endl;
    clock_t start;
    start = std::clock();
    int* Pijk_id = new int[8];
    double* Pijk = new double[8*3];
    int nb = us3d->bnd_face_map[wall_id].size();
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    //std::vector<Vert*> face_c;
    int local_face_id;
    std::map<int,std::vector<Vert*> > prisms;
    
    Vec3D* cut_dir_face0 = new Vec3D;
    Vec3D* cut_dir_facet0 = new Vec3D;
    Vec3D* cut_dir_facet1 = new Vec3D;
    int bvid,opposite_bvid;
    std::vector<int> opposite_tri(3);
    
    std::vector<int> prism0;
    std::vector<int> prism1;

    std::vector<int> prismStored0(6);
    std::vector<int> prismStored1(6);
    mesh_topology_bl->Nprisms = 0;
    int glob_el_id = 0;
    
    for(int bf=0;bf<us3d->bnd_face_map[wall_id].size();bf++)
    {
        std::vector<int> layer;
        int bfaceid = us3d->bnd_face_map[wall_id][bf];
        int faceid  = bfaceid;
        int elid0   = us3d->ife->getVal(faceid,0);
        int elid1   = us3d->ife->getVal(faceid,1);

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
        for(int k=0;k<6;k++)
        {
            if(ief_g->getVal(elid_cur,k)==faceid)
            {
                local_face_id = k;
            }
        }

        for(int k=0;k<8;k++)
        {
           loc_vid     = ien_g->getVal(elid_cur,k);
           Pijk_id[k]  = loc_vid;
           Pijk[k*3+0] = xcn_g->getVal(loc_vid,0);
           Pijk[k*3+1] = xcn_g->getVal(loc_vid,1);
           Pijk[k*3+2] = xcn_g->getVal(loc_vid,2);
        }
        
        int changed = ChkHexorient(Pijk,Pijk_id);
//
        Vert* Vijk = ComputeCenterCoord(Pijk, 8);
        
        Vert* Vface  = new Vert;
        std::vector<Vert*> face;
        std::vector<Vert*> face_turned(4);
        std::vector<Vert*> face_turned2(4);
        std::set<int> conn_bvid;
        for(int r=0;r<4;r++)
        {
            int vid  = us3d->ifn->getVal(faceid,r);
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
        std::vector<int> tri0(3);
        std::vector<int> tri1(3);
        prism0.push_back(bvid);
        prism0.push_back(us3d->ifn->getVal(faceid,1));
        prism0.push_back(us3d->ifn->getVal(faceid,3));
        tri0[0] = bvid;
        tri0[1] = us3d->ifn->getVal(faceid,1);
        tri0[2] = us3d->ifn->getVal(faceid,3);
        mesh_topology_bl->BndFaces.push_back(tri0);
        prism1.push_back(us3d->ifn->getVal(faceid,2));
        prism1.push_back(us3d->ifn->getVal(faceid,3));
        prism1.push_back(us3d->ifn->getVal(faceid,1));
        tri1[0] = us3d->ifn->getVal(faceid,2);
        tri1[1] = us3d->ifn->getVal(faceid,3);
        tri1[2] = us3d->ifn->getVal(faceid,1);
        mesh_topology_bl->BndFaces.push_back(tri1);
        
        conn_bvid.insert(us3d->ifn->getVal(faceid,1));
        conn_bvid.insert(us3d->ifn->getVal(faceid,3));
        
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
            
            int changed = ChkHexorient(Pijk,Pijk_id);
            
            Vert* Vijk = ComputeCenterCoord(Pijk, 8);
            std::vector<std::vector<int> > face_id_stored(6);
            std::vector<std::vector<Vert*> > face_stored(6);
            map<int,std::set<int> > local_node2node_element;
            std::vector<map<int,std::set<int> > > local_node2node_face(6);
            std::vector<map<int,int> > local_node2opponode_face(6);
            for(int k=0;k<6;k++)
            {
                int fid = ief_g->getVal(elid_cur,k);
                Vert* Vface2  = new Vert;
                
                std::vector<int> faceVert_IDs(4);
                std::vector<Vert*> face2;
                for(int r=0;r<4;r++)
                {
                    int vid  = us3d->ifn->getVal(fid,r);
                    
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
                
                local_node2node_element[us3d->ifn->getVal(fid,0)].insert(us3d->ifn->getVal(fid,1));
                local_node2node_element[us3d->ifn->getVal(fid,0)].insert(us3d->ifn->getVal(fid,3));
                local_node2node_element[us3d->ifn->getVal(fid,1)].insert(us3d->ifn->getVal(fid,0));
                local_node2node_element[us3d->ifn->getVal(fid,1)].insert(us3d->ifn->getVal(fid,2));
                local_node2node_element[us3d->ifn->getVal(fid,2)].insert(us3d->ifn->getVal(fid,1));
                local_node2node_element[us3d->ifn->getVal(fid,2)].insert(us3d->ifn->getVal(fid,3));
                local_node2node_element[us3d->ifn->getVal(fid,3)].insert(us3d->ifn->getVal(fid,2));
                local_node2node_element[us3d->ifn->getVal(fid,3)].insert(us3d->ifn->getVal(fid,0));
                
                local_node2node_face[k][us3d->ifn->getVal(fid,0)].insert(us3d->ifn->getVal(fid,1));
                local_node2node_face[k][us3d->ifn->getVal(fid,0)].insert(us3d->ifn->getVal(fid,3));
                local_node2node_face[k][us3d->ifn->getVal(fid,1)].insert(us3d->ifn->getVal(fid,0));
                local_node2node_face[k][us3d->ifn->getVal(fid,1)].insert(us3d->ifn->getVal(fid,2));
                local_node2node_face[k][us3d->ifn->getVal(fid,2)].insert(us3d->ifn->getVal(fid,1));
                local_node2node_face[k][us3d->ifn->getVal(fid,2)].insert(us3d->ifn->getVal(fid,3));
                local_node2node_face[k][us3d->ifn->getVal(fid,3)].insert(us3d->ifn->getVal(fid,2));
                local_node2node_face[k][us3d->ifn->getVal(fid,3)].insert(us3d->ifn->getVal(fid,0));
                
                local_node2opponode_face[k][us3d->ifn->getVal(fid,0)]=us3d->ifn->getVal(fid,2);
                local_node2opponode_face[k][us3d->ifn->getVal(fid,1)]=us3d->ifn->getVal(fid,3);
                local_node2opponode_face[k][us3d->ifn->getVal(fid,2)]=us3d->ifn->getVal(fid,0);
                local_node2opponode_face[k][us3d->ifn->getVal(fid,3)]=us3d->ifn->getVal(fid,1);

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
            
            std::set<int>::iterator its;
            for(its=local_node2node_element[bvid].begin();its!=local_node2node_element[bvid].end();its++)
            {
                if(conn_bvid.find(*its)==conn_bvid.end())
                {
                    opposite_bvid = *its;
                }
            }
        
            int min_index  = std::min_element(dp.begin(),dp.end())-dp.begin();
            double min_val = *std::min_element(dp.begin(),dp.end());

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

            prism0.push_back(opposite_tri[0]);
            prism0.push_back(opposite_tri[1]);
            prism0.push_back(opposite_tri[2]);
            
            prism1.push_back(local_node2opponode_face[min_index][opposite_bvid]);
            prism1.push_back(opposite_tri[2]);
            prism1.push_back(opposite_tri[1]);
            
            NegateVec3D(nbf);

            int gEl0=us3d->ife->getVal(fid_new,0);
            int gEl1=us3d->ife->getVal(fid_new,1);

            if(gEl0==elid_cur)
            {
                elid_next = gEl1;
            }
            else if(gEl1==elid_cur)
            {
                elid_next = gEl0;
            }
            layer.push_back(elid_next);
            
            if(c==nLayer-1)
            {
                mesh_topology_bl->outer_shell_faces.push_back(fid_new);
            }
            
            prismStored0[0] = prism0[0];prismStored0[1] = prism0[1];prismStored0[2] = prism0[2];
            prismStored0[3] = prism0[3];prismStored0[4] = prism0[4];prismStored0[5] = prism0[5];
            //std::cout << "prims0 " << glob_el_id << " :: " << prism0[0] << " " << prism0[1] << " " << prism0[2] << " " << prism0[3] << " " << prism0[4] << " " << prism0[5] << std::endl;
            Element* P0     = new Element;
            P0->GlobalNodes = prismStored0;
            P0->globID      = glob_el_id;
            
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[2]);
            
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[5]);
            
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[3]);
            
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[4]);
            
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[5]);
            glob_el_id = glob_el_id+1;
            
            
            prismStored1[0] = prism1[0];prismStored1[1] = prism1[1];prismStored1[2] = prism1[2];
            prismStored1[3] = prism1[3];prismStored1[4] = prism1[4];prismStored1[5] = prism1[5];
            //std::cout << "prims1 " << glob_el_id << " :: " <<prism1[0] << " " << prism1[1] << " " << prism1[2] << " " << prism1[3] << " " << prism1[4] << " " << prism1[5] << std::endl;
            
            Element* P1     = new Element;
            P1->GlobalNodes = prismStored1;
            P1->globID      = glob_el_id;
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[0].push_back(prismStored1[2]);
            
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[5]);
            
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[2]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[3]);
            
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[2]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[5]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[4]);
            
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[5]);
            glob_el_id = glob_el_id+1;
        //    std::cout << prismStored0[0] << " " << prismStored0[1] << " " << prismStored0[2] << " " << prismStored0[3] << " " << prismStored0[4] << " " << prismStored0[5] << std::endl;
         //   std::cout << prismStored1[0] << " " << prismStored1[1] << " " << prismStored1[2] << " " << prismStored1[3] << " " << prismStored1[4] << " " << prismStored1[5] << std::endl;
            PElements[c*2+0]=P0;
            PElements[c*2+1]=P1;
            PPrisms[c*2+0] = prismStored0;
            PPrisms[c*2+1] = prismStored1;

            prism0.clear();
            prism1.clear();
                    
            // Store the same initial triangles as the determined opposite triangles.
            // However change orientation of the nodes.
            prism0.push_back(opposite_tri[0]);
            prism0.push_back(opposite_tri[2]);
            prism0.push_back(opposite_tri[1]);
            prism1.push_back(local_node2opponode_face[min_index][opposite_bvid]);
            prism1.push_back(opposite_tri[1]);
            prism1.push_back(opposite_tri[2]);
            
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
        prism0.clear();
        prism1.clear();
        mesh_topology_bl->BLlayersPrisms[bfaceid]=PPrisms;
        mesh_topology_bl->BLlayersElements[bfaceid]=PElements;
        mesh_topology_bl->BLlayers[bfaceid]=layer;
        mesh_topology_bl->Nprisms = mesh_topology_bl->Nprisms+PPrisms.size();
    }

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << " extracting BL mesh = " << duration << std::endl;
    std::map<int,std::vector<int> >::iterator itt;
    std::vector<int> elements;
    for(itt=mesh_topology_bl->BLlayers.begin();itt!=mesh_topology_bl->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            elements.push_back(itt->second[q]);
        }
    }
    OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_");
       
    return mesh_topology_bl;
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
    
    if(world_rank == 0)
    {
        UnitTestEigenDecomp();
    }
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
        int varia = 0;
        US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);

        //MMG5_pMesh mmgMesh = ReadMMG_pMesh(us3d,comm,info);
        
        int Nel_part = us3d->ien->getNrow();
        
        ParallelState* ien_pstate = new ParallelState(us3d->ien->getNglob(),comm);

        ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);

        ParallelState* xcn_pstate = new ParallelState(us3d->xcn->getNglob(),comm);
        Array<double>* Uivar = new Array<double>(Nel_part,1);
        for(int i=0;i<Nel_part;i++)
        {
            Uivar->setVal(i,0,us3d->interior->getVal(i,varia));
            std::cout << us3d->interior->getVal(i,varia) << std::endl;
        }
        
        delete us3d->interior;
        
//        Array<double>* xcn_g;
//        Array<int>* ief_g;
//        Array<int>* ien_g;
//        if(world_rank == 0)
//        {
//            xcn_g = new Array<double>(us3d->xcn->getNglob(),3);
//            ief_g = new Array<int>(us3d->ief->getNglob(),6);
//            ien_g = new Array<int>(us3d->ien->getNglob(),8);
//        }
//        else
//        {
//            xcn_g = new Array<double>(1,1);
//            ief_g = new Array<int>(1,1);
//            ien_g = new Array<int>(1,1);
//        }
//
//        int* ien_nlocs      = new int[world_size];
//        int* ien_offsets    = new int[world_size];
//        int* ief_nlocs      = new int[world_size];
//        int* ief_offsets    = new int[world_size];
//        int* xcn_nlocs      = new int[world_size];
//        int* xcn_offsets    = new int[world_size];
//        for(int i=0;i<world_size;i++)
//        {
//            xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]*3;
//            xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;
//
//            ief_nlocs[i]   = ien_pstate->getNlocs()[i]*6;
//            ief_offsets[i] = ien_pstate->getOffsets()[i]*6;
//
//            ien_nlocs[i]   = ien_pstate->getNlocs()[i]*8;
//            ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
//        }
//
//        MPI_Gatherv(&us3d->xcn->data[0],
//                    us3d->xcn->getNrow()*3,
//                    MPI_DOUBLE,
//                    &xcn_g->data[0],
//                    xcn_nlocs,
//                    xcn_offsets,
//                    MPI_DOUBLE, 0, comm);
//
//        MPI_Gatherv(&us3d->ief->data[0],
//                    us3d->ief->getNrow()*6,
//                    MPI_INT,
//                    &ief_g->data[0],
//                    ief_nlocs,
//                    ief_offsets,
//                    MPI_INT, 0, comm);
//
//        MPI_Gatherv(&us3d->ien->data[0],
//                    us3d->ien->getNrow()*8,
//                    MPI_INT,
//                    &ien_g->data[0],
//                    ien_nlocs,
//                    ien_offsets,
//                    MPI_INT, 0, comm);
//
//        if(world_rank == 0)
//        {
//            int wall_id = 4;
//            int nLayer = 10;
//            std::map<int,std::set<int> > shell = FindOuterShellBoundaryLayerMesh(wall_id, nLayer, us3d,xcn_g,ien_g,ief_g,xcn_pstate,ien_pstate,comm);
//            std::map<int,std::set<int> >::iterator itt;
//
//            for(itt=shell.begin();itt!=shell.end();itt++)
//            {
//                std::cout << itt->first << " -> ";
//                std::set<int>::iterator itt2;
//                for(itt2=itt->second.begin();itt2!=itt->second.end();itt2++)
//                {
//                    std::cout << *itt2 << " ";
//                }
//                std::cout << std::endl;
//            }
//
////            for(int i=0;i<shell.size();i++)
////            {
////                std::cout << shell[i] << std::endl;
////            }
////
////            Mesh_Topology_BL* BLmesh = ExtractBoundaryLayerMesh(wall_id, nLayer, us3d,xcn_g,ien_g,ief_g,xcn_pstate,ien_pstate,comm);
////
////            OutputBoundaryLayerPrisms(xcn_g, BLmesh, comm);
//        }
        
//
//        (/
        
        Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, parmetis_pstate,ien_pstate, us3d->xcn,xcn_pstate,Uivar,comm);

        //Array<double>* dUdXi_v2 = ComputedUdx(P, pstate, iee_copy, iee_loc, ief_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm, ife_copy);


        std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uivar,comm);

        Mesh_Topology* meshTopo = new Mesh_Topology(P, us3d->ifn, us3d->ife,UauxNew,us3d->bnd_map,us3d->bnd_face_map,us3d->nBnd,comm);
        
        
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
        Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
        for(int i=0;i<us3d->ghost->getNrow();i++)
        {
            gB->setVal(i,0,us3d->ghost->getVal(i,varia));
        }
//        Array<double>* dUdXi = ComputedUdx_LSQ_US3D_v3(P,UauxNew,meshTopo,gB,comm);
        Array<double>* dUdXi   = ComputedUdx_MGG(P,UauxNew,meshTopo,gB,comm);
        //Gradients* dudxObj   = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","LSQ",comm);
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
        
//        Array<double>* dU2dXi2 = ComputedUdx_LSQ_US3D_v3(P,dUdxauxNew,meshTopo,gB,comm);
//        Array<double>* dU2dYi2 = ComputedUdx_LSQ_US3D_v3(P,dUdyauxNew,meshTopo,gB,comm);
//        Array<double>* dU2dZi2 = ComputedUdx_LSQ_US3D_v3(P,dUdzauxNew,meshTopo,gB,comm);

        Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,gB,comm);
        Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,gB,comm);
        Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,gB,comm);
        
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
        Array<double>* hessian = new Array<double>(nVerts,6);
        Array<double>* grad    = new Array<double>(nVerts,3);
        for(int i=0;i<nVerts;i++)
        {
            grad->setVal(i,0,dudx_v[i]);grad->setVal(i,1,dudy_v[i]);grad->setVal(i,2,dudz_v[i]);
            hessian->setVal(i,0,d2udx2_v[i]); hessian->setVal(i,1,d2udxy_v[i]); hessian->setVal(i,2,d2udxz_v[i]);
            hessian->setVal(i,3,d2udy2_v[i]); hessian->setVal(i,4,d2udyz_v[i]);
            hessian->setVal(i,5,d2udz2_v[i]);
        }
        
        //Array<double>* UgRoot = GetSolutionOnRoot(P,us3d,hessian,comm);
        
        std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
        double max_v = *std::max_element(d2udx2_v.begin(), d2udx2_v.end());
        int dim = 3;
        
        Array<double>* metric = ComputeMetric(Verts,grad,hessian,max_v,loc_elem2verts_loc,us3d->ien->getNrow(),comm,dim);
        
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
            myfile << Verts[i].x << " " << Verts[i].y << " " << Verts[i].z <<
            " " << u_v[i] << " "<< dudx_v[i] << " " << dudy_v[i] << " " << dudz_v[i] <<
            " " << hessian->getVal(i,0) << " " << hessian->getVal(i,1) << " " << hessian->getVal(i,2) <<
            " " << hessian->getVal(i,1) << " " << hessian->getVal(i,3) << " " << hessian->getVal(i,4) <<
            " " << hessian->getVal(i,2) << " " << hessian->getVal(i,4) << " " << hessian->getVal(i,5) <<
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
        
        MMG_Mesh* mmg = GetOptimizedMMG3DMeshOnRoot(P, us3d, hessian, metric, comm);
        
        if (world_rank == 0)
        {
            MMG3D_Set_handGivenMesh(mmg->mmgMesh);
            if ( MMG3D_Set_dparameter(mmg->mmgMesh,mmg->mmgSol,MMG3D_DPARAM_hgrad, 2.0) != 1 )    exit(EXIT_FAILURE);
            
            int ier = MMG3D_mmg3dlib(mmg->mmgMesh,mmg->mmgSol);
            
            OutputMesh_MMG(mmg->mmgMesh);
//            OutputBoundaryID_MMG(mmgMesh,ref2bface,1);
//            OutputBoundaryID_MMG(mmgMesh,ref2bface,2);
//            OutputBoundaryID_MMG(mmgMesh,ref2bface,3);
//            OutputBoundaryID_MMG(mmgMesh,ref2bface,4);
            WriteUS3DGridFromMMG(mmg->mmgMesh, us3d);
        }
        
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
            
        MPI_Finalize();
        
    }
     
    return 0;
}
