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

struct Mdata{
    std::vector<std::vector<double> > Vmetric;
    std::vector<std::vector<int> > arrHex;
};

Mdata* ReadMetricData()
{
    std::ifstream fin;
    fin.open("metric_restart_latest_msl.dat");

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
    finhex.open("elements_restart_latest_msl.dat");

    // Read the file row by row
    std::vector<int> rowHex(8);
    std::vector<std::vector<int> > arrHex;

    while(finhex >> rowHex[0] >> rowHex[1] >> rowHex[2] >> rowHex[3] >> rowHex[4] >> rowHex[5] >> rowHex[6] >> rowHex[7])
    {
       arrHex.push_back(rowHex);
    }
    Mdata* md = new Mdata;
    md->arrHex = arrHex;
    md->Vmetric = Vmetric;
    return md;
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
    
    OutputMesh_MMG(mmgMesh,0,mmgMesh->ne,"MMgOutput.dat");
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,1);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,2);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,3);
    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,4);
    //    WriteUS3DGridFromMMG(mmgMesh, us3d);


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


void MatchBoundaryTags(US3D* us3d, MMG5_pMesh mmgMesh,int offset_NE, int Nel)
{
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
    int t = 1;
    // local face2vert_map for a tet in mmg  {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}

    for(int i=1;i<=offset_NE;i++)
    {
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
        tria0.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
        if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
        {
            ref0 = us3d->tria_ref_map[tria0];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[1];
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
        
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
        tria2.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
        if(us3d->tria_ref_map.find(tria2)!=us3d->tria_ref_map.end())
        {
            ref2 = us3d->tria_ref_map[tria2];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[2];
            mmgMesh->tria[t].ref  = ref2;
            t++;
        }
        
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
        tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
        if(us3d->tria_ref_map.find(tria3)!=us3d->tria_ref_map.end())
        {
            ref3 = us3d->tria_ref_map[tria3];
            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[1];
            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
            mmgMesh->tria[t].ref  = ref3;
            t++;
        }
        tria0.clear();
        tria1.clear();
        tria2.clear();
        tria3.clear();
    }
    
    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
    int nPrism = mmgMesh->nprism;
    
    if (nPrism != 0)
    {
        std::set<int> quad0;
        std::set<int> quad1;
        std::set<int> quad2;
        int tq = 1;
        for(int i=0;i<nPrism;i++)
        {
            int v0 = mmgMesh->prism[i+1].v[0];
            int v1 = mmgMesh->prism[i+1].v[1];
            int v2 = mmgMesh->prism[i+1].v[2];
            int v3 = mmgMesh->prism[i+1].v[3];
            int v4 = mmgMesh->prism[i+1].v[4];
            int v5 = mmgMesh->prism[i+1].v[5];
            
            tria0.insert(v0);
            tria0.insert(v1);
            tria0.insert(v2);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

            if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
            {
                ref0 = us3d->tria_ref_map[tria0];
                mmgMesh->tria[t].v[0] = v0;
                mmgMesh->tria[t].v[1] = v1;
                mmgMesh->tria[t].v[2] = v2;
                mmgMesh->tria[t].ref  = ref0;
                t++;
            }
            
            tria1.insert(v3);
            tria1.insert(v4);
            tria1.insert(v5);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

            if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
            {
                ref1 = us3d->tria_ref_map[tria1];
                mmgMesh->tria[t].v[0] = v3;
                mmgMesh->tria[t].v[1] = v5;
                mmgMesh->tria[t].v[2] = v4;
                mmgMesh->tria[t].ref  = ref1;
                t++;
            }

            quad0.insert(v0);
            quad0.insert(v3);
            quad0.insert(v4);
            quad0.insert(v1);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

            if(us3d->quad_ref_map.find(quad0)!=us3d->quad_ref_map.end())
            {
                ref0 = us3d->quad_ref_map[quad0];
                mmgMesh->quadra[tq].v[0] = v0;
                mmgMesh->quadra[tq].v[1] = v3;
                mmgMesh->quadra[tq].v[2] = v4;
                mmgMesh->quadra[tq].v[3] = v1;
                mmgMesh->quadra[tq].ref  = ref0;
                std::cout << "ref0 " << ref1 << std::endl;

                tq++;
            }

            quad1.insert(v1);
            quad1.insert(v4);
            quad1.insert(v5);
            quad1.insert(v2);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

            if(us3d->quad_ref_map.find(quad1)!=us3d->quad_ref_map.end())
            {
                ref1 = us3d->quad_ref_map[quad1];
                mmgMesh->quadra[tq].v[0] = v1;
                mmgMesh->quadra[tq].v[1] = v4;
                mmgMesh->quadra[tq].v[2] = v5;
                mmgMesh->quadra[tq].v[3] = v2;
                mmgMesh->quadra[tq].ref  = ref1;
                std::cout << "ref1 " << ref1 << std::endl;
                tq++;
            }

            quad2.insert(v2);
            quad2.insert(v0);
            quad2.insert(v3);
            quad2.insert(v5);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };

            if(us3d->quad_ref_map.find(quad2)!=us3d->quad_ref_map.end())
            {
                ref2 = us3d->quad_ref_map[quad2];
                mmgMesh->quadra[tq].v[0] = v0;
                mmgMesh->quadra[tq].v[1] = v2;
                mmgMesh->quadra[tq].v[2] = v5;
                mmgMesh->quadra[tq].v[3] = v3;
                mmgMesh->quadra[tq].ref  = ref2;
                std::cout << "ref2 " << ref2 << std::endl;

                tq++;
            }
            tria0.clear();
            tria1.clear();
            quad0.clear();
            quad1.clear();
            quad2.clear();
        
        }
    }
    
    
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //End of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}



struct BLShellInfo{
  
    std::map<int,int> ShellFace2BFace;
    std::map<int,int> BFace2ShellFace;
    std::map<std::set<int>,int> ShellTri2FaceID;
    std::map<int,std::vector<int> > ShellFaceID2TriID;
    std::map<int,int> FaceID2TopoType;
    std::map<int,std::map<int,int> > ShellFace2ShellVert2OppositeBoundaryVerts;
    Array<int>* ShellRef;
    std::map<int,std::vector<int> > BLlayers;
    std::map<int, std::vector<map<int,std::set<int> > > > bFace2_locN2NEl;
    std::set<int> elements_set;
    std::set<int> verts_set;
};


BLShellInfo* FindOuterShellBoundaryLayerMesh(int wall_id, int nLayer, US3D* us3d,
                                             Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g,
                                             ParallelState* xcn_pstate, ParallelState* ien_pstate, MPI_Comm comm)
{
    BLShellInfo* BLinfo = new BLShellInfo;
    BLinfo->ShellRef = new Array<int>(xcn_g->getNrow(),1);
    int te1=0;
    int te2=0;
    int te3=0;
    for(int i=0;i<xcn_g->getNrow();i++)
    {
        if(us3d->vert_ref_map[i]!=0)
        {
            BLinfo->ShellRef->setVal(i,0,100+us3d->vert_ref_map[i]);
            //std::cout << "fqwk " << 100+us3d->vert_ref_map[i] << std::endl;
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
    int bvid_b;
    int fv1_b;
    int fv2_b;
    int fv3_b;
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
        
        conn_bvid.insert(us3d->ifn->getVal(faceid,1));
        conn_bvid.insert(us3d->ifn->getVal(faceid,3));
        bvid_b = bvid;
        fv1_b = us3d->ifn->getVal(faceid,1);
        fv2_b = us3d->ifn->getVal(faceid,2);
        fv3_b = us3d->ifn->getVal(faceid,3);
        
        
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
            double min_val = *std::min_element(dp.begin(),dp.end());
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
            if(c<nLayer-1)
            {
                layer.push_back(elid_next);
            }
            
            int changed = ChkHexorient(Pijk,Pijk_id);
            
            if(c==nLayer-1)
            {
                outer_shell_faces.push_back(fid_new);
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,0));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,1));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,2));
                outer_shell_Faces2Nodes[fid_new].insert(us3d->ifn->getVal(fid_new,3));
                BLinfo->ShellRef->setVal(us3d->ifn->getVal(fid_new,0),0,-1);
                BLinfo->ShellRef->setVal(us3d->ifn->getVal(fid_new,1),0,-1);
                BLinfo->ShellRef->setVal(us3d->ifn->getVal(fid_new,2),0,-1);
                BLinfo->ShellRef->setVal(us3d->ifn->getVal(fid_new,3),0,-1);
                std::set<int> ShellTri0;
                ShellTri0.insert(us3d->ifn->getVal(fid_new,0));
                ShellTri0.insert(us3d->ifn->getVal(fid_new,1));
                ShellTri0.insert(us3d->ifn->getVal(fid_new,3));
                BLinfo->ShellTri2FaceID[ShellTri0] = fid_new;
                std::set<int> ShellTri1;
                ShellTri1.insert(us3d->ifn->getVal(fid_new,1));
                ShellTri1.insert(us3d->ifn->getVal(fid_new,2));
                ShellTri1.insert(us3d->ifn->getVal(fid_new,3));
                BLinfo->ShellTri2FaceID[ShellTri1] = fid_new;
                std::set<int> ShellTri2;
                ShellTri2.insert(us3d->ifn->getVal(fid_new,0));
                ShellTri2.insert(us3d->ifn->getVal(fid_new,1));
                ShellTri2.insert(us3d->ifn->getVal(fid_new,2));
                BLinfo->ShellTri2FaceID[ShellTri2] = fid_new;
                std::set<int> ShellTri3;
                ShellTri3.insert(us3d->ifn->getVal(fid_new,2));
                ShellTri3.insert(us3d->ifn->getVal(fid_new,3));
                ShellTri3.insert(us3d->ifn->getVal(fid_new,0));
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
    std::cout << "Shell quad faces = " << outer_shell_faces.size() << std::endl;

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << " extracting outer shell BL mesh = " << duration << std::endl;
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
    std::vector<int> opposite_tri1(3);
    
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
            
            std::set<int>::iterator itset;
            for(itset=local_node2node_element[bvid].begin();itset!=local_node2node_element[bvid].end();itset++)
            {
                if(conn_bvid.find(*itset)==conn_bvid.end())
                {
                    opposite_bvid = *itset;
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
            
            if(c<nLayer-1)
            {
                layer.push_back(elid_next);
            }
            
            if(c==nLayer-1)
            {
                mesh_topology_bl->outer_shell_faces.push_back(fid_new);
            }
            
            prismStored0[0] = prism0[0];prismStored0[1] = prism0[1];prismStored0[2] = prism0[2];
            prismStored0[3] = prism0[3];prismStored0[4] = prism0[4];prismStored0[5] = prism0[5];

            Element* P0     = new Element;
            P0->GlobalNodes = prismStored0;
            P0->globID      = glob_el_id;
            
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[0].push_back(prismStored0[2]);
            
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[4]);
            
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[3]);
            
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[5]);
            
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[5]);
            glob_el_id = glob_el_id+1;
            
            prismStored1[0] = prism1[0];prismStored1[1] = prism1[1];prismStored1[2] = prism1[2];
            prismStored1[3] = prism1[3];prismStored1[4] = prism1[4];prismStored1[5] = prism1[5];
            
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
     for(itt=mesh_topology_bl->BLlayers.begin();itt!=mesh_topology_bl->BLlayers.end();itt++)
    {
        for(int q=0;q<itt->second.size();q++)
        {
            mesh_topology_bl->elements.push_back(itt->second[q]);
        }
    }
    OutputBLElementsOnRoot(xcn_g,ien_g,mesh_topology_bl->elements,comm,"BL_Root_");
       
    return mesh_topology_bl;
}



//            struct BLShellInfo
//            {
//                std::map<int,int> ShellFace2BFace;
//                std::map<std::set<int>,int> ShellTri2FaceID;
//                std::map<std::set<int>,std::vector<int> > ShellFaceID2TriID
//                std::map<int,int> FaceID2TopoType;
//                std::map<int,std::map<int,int> > ShellFace2ShellVert2OppositeBoundaryVerts;
//                Array<int>* ShellRef;
//            };


Mesh_Topology_BL* ExtractBoundaryLayerMeshFromShell(std::vector<std::vector<int> > u_tris, BLShellInfo* BLshell, int wall_id, int nLayer, US3D* us3d, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, ParallelState* xcn_pstate, ParallelState* ien_pstate, MPI_Comm comm)
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
    int nb = us3d->bnd_face_map[wall_id].size();
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    //std::vector<Vert*> face_c;
    int local_face_id;
    std::map<int,std::vector<Vert*> > prisms;
    
    std::set<int> tria0;
    std::set<int> tria1;
    std::set<int> quad0;
    std::set<int> quad1;
    std::set<int> quad2;
    
    
    Vec3D* cut_dir_face0 = new Vec3D;
    Vec3D* cut_dir_facet0 = new Vec3D;
    Vec3D* cut_dir_facet1 = new Vec3D;
    int bvid,opposite_bvid;
    std::vector<int> opposite_tri(3);
    std::vector<int> opposite_tri_n(3);
    std::vector<int> opposite_tri1(3);
    std::vector<int> opposite_tri1_n(3);
    std::vector<int> prism0;
    std::vector<int> prism1;

    std::vector<int> prismStored0(6);
    std::vector<int> prismStored1(6);
    mesh_topology_bl->Nprisms = 0;
    int glob_el_id = 0;
    std::map<int,int> bface2shellface = BLshell->BFace2ShellFace;
    std::map<std::set<int>,int> shelltri2shellface = BLshell->ShellTri2FaceID;
    std::map<int,std::vector<int> > shellfaceID2triID =  BLshell->ShellFaceID2TriID;
    int opposite_p00,opposite_p01,opposite_p02,opposite_p10,opposite_p11,opposite_p12;
    int cnt_turn = 0;
    for(int bf=0;bf<us3d->bnd_face_map[wall_id].size();bf++)
    {
        
        int bvid,obvid_i,opposite_bvid;
        std::vector<int> layer;
        int bfaceid      = us3d->bnd_face_map[wall_id][bf];
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
            int vid  = us3d->ifn->getVal(faceid,r);
            
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
        
        
        double orient_t0 = DotVec3D(r0,n_t0);
        double orient_t1 = DotVec3D(r0,n_t10);
        
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
        n_t0 = ComputeSurfaceNormal(v_t0,v_t1);
        double orient_t0_check = DotVec3D(r0,n_t0);
        
        v_t10->c0 = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10->c1 = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10->c2 = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);

        v_t11->c0 = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11->c1 = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11->c2 = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        n_t10 = ComputeSurfaceNormal(v_t10,v_t11);
        double orient_t1_check = DotVec3D(r0,n_t10);
        //std::cout << "check = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
        
        prism0.push_back(tri_0n[0]);
        prism0.push_back(tri_0n[2]);
        prism0.push_back(tri_0n[1]);
        
        prism1.push_back(tri_1n[0]);
        prism1.push_back(tri_1n[2]);
        prism1.push_back(tri_1n[1]);
        
        
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
            
            std::set<int>::iterator itset;
            for(itset=local_node2node_element[bvid].begin();itset!=local_node2node_element[bvid].end();itset++)
            {
                if(conn_bvid.find(*itset)==conn_bvid.end())
                {
                    opposite_bvid = *itset;
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
            
            opposite_tri1[0] = local_node2opponode_face[min_index][opposite_bvid];
            opposite_tri1[1] = opposite_tri[2];
            opposite_tri1[2] = opposite_tri[1];
            
            Vec3D* v_toppo0 = new Vec3D;
            v_toppo0->c0 = xcn_g->getVal(opposite_tri[1],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo0->c1 = xcn_g->getVal(opposite_tri[1],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo0->c2 = xcn_g->getVal(opposite_tri[1],2)-xcn_g->getVal(opposite_tri[0],2);
            Vec3D* v_toppo1 = new Vec3D;
            v_toppo1->c0 = xcn_g->getVal(opposite_tri[2],0)-xcn_g->getVal(opposite_tri[0],0);
            v_toppo1->c1 = xcn_g->getVal(opposite_tri[2],1)-xcn_g->getVal(opposite_tri[0],1);
            v_toppo1->c2 = xcn_g->getVal(opposite_tri[2],2)-xcn_g->getVal(opposite_tri[0],2);
            
            Vec3D* n_toppo0        = ComputeSurfaceNormal(v_toppo0,v_toppo1);
            double orient0oppo0   = DotVec3D(n_t0 ,n_toppo0 );
            
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

            Vec3D* v_toppo10 = new Vec3D;
            v_toppo10->c0 = xcn_g->getVal(opposite_tri1[1],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo10->c1 = xcn_g->getVal(opposite_tri1[1],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo10->c2 = xcn_g->getVal(opposite_tri1[1],2)-xcn_g->getVal(opposite_tri1[0],2);
            Vec3D* v_toppo11 = new Vec3D;
            v_toppo11->c0 = xcn_g->getVal(opposite_tri1[2],0)-xcn_g->getVal(opposite_tri1[0],0);
            v_toppo11->c1 = xcn_g->getVal(opposite_tri1[2],1)-xcn_g->getVal(opposite_tri1[0],1);
            v_toppo11->c2 = xcn_g->getVal(opposite_tri1[2],2)-xcn_g->getVal(opposite_tri1[0],2);
            
            Vec3D* n_toppo10        = ComputeSurfaceNormal(v_toppo10,v_toppo11);
            double orient0oppo10    = DotVec3D(n_t10 , n_toppo10 );
            
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
            n_toppo10        = ComputeSurfaceNormal(v_toppo10,v_toppo11);
            orient0oppo10    = DotVec3D(n_t10 , n_toppo10 );
            if(orient0oppo10>0)
            {
                std::cout << " Still not changed "<<std::endl;
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
            

            prism0.push_back(opposite_tri[0]);
            prism0.push_back(opposite_tri[1]);
            prism0.push_back(opposite_tri[2]);
            
            prism1.push_back(opposite_tri1[0]);
            prism1.push_back(opposite_tri1[1]);
            prism1.push_back(opposite_tri1[2]);

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
            if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
            {
                int ref0 = us3d->tria_ref_map[tria0];
                std::vector<int> bctria(3);
                bctria[0] = prismStored0[0];
                bctria[1] = prismStored0[1];
                bctria[2] = prismStored0[2];
                mesh_topology_bl->bcTria[ref0].push_back(bctria);
            }
            
            
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[1].push_back(prismStored0[4]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria1.insert(prismStored0[3]);
            tria1.insert(prismStored0[4]);
            tria1.insert(prismStored0[5]);
            if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
            {
                int ref1 = us3d->tria_ref_map[tria1];
                std::vector<int> bctria(3);
                bctria[0] = prismStored0[3];
                bctria[1] = prismStored0[4];
                bctria[2] = prismStored0[5];
                mesh_topology_bl->bcTria[ref1].push_back(bctria);
            }
            
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[4]);
            P0->LocalFace2GlobalNode[2].push_back(prismStored0[3]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad0.insert(prismStored0[0]);
            quad0.insert(prismStored0[1]);
            quad0.insert(prismStored0[4]);
            quad0.insert(prismStored0[3]);
            if(us3d->quad_ref_map.find(quad0)!=us3d->quad_ref_map.end())
            {
                int ref0 = us3d->quad_ref_map[quad0];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[0];
                bcquad[1] = prismStored0[1];
                bcquad[2] = prismStored0[4];
                bcquad[3] = prismStored0[3];
                mesh_topology_bl->bcQuad[ref0].push_back(bcquad);
            }
            
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[1]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[5]);
            P0->LocalFace2GlobalNode[3].push_back(prismStored0[4]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad1.insert(prismStored0[1]);
            quad1.insert(prismStored0[2]);
            quad1.insert(prismStored0[5]);
            quad1.insert(prismStored0[4]);
            
            if(us3d->quad_ref_map.find(quad1)!=us3d->quad_ref_map.end())
            {
                int ref1 = us3d->quad_ref_map[quad1];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[1];
                bcquad[1] = prismStored0[2];
                bcquad[2] = prismStored0[4];
                bcquad[3] = prismStored0[5];
                mesh_topology_bl->bcQuad[ref1].push_back(bcquad);

            }
            
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[2]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[0]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[3]);
            P0->LocalFace2GlobalNode[4].push_back(prismStored0[5]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad2.insert(prismStored0[2]);
            quad2.insert(prismStored0[0]);
            quad2.insert(prismStored0[3]);
            quad2.insert(prismStored0[5]);
            if(us3d->quad_ref_map.find(quad2)!=us3d->quad_ref_map.end())
            {
                int ref2 = us3d->quad_ref_map[quad2];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored0[2];
                bcquad[1] = prismStored0[0];
                bcquad[2] = prismStored0[3];
                bcquad[3] = prismStored0[5];
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
            if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
            {
                int ref0 = us3d->tria_ref_map[tria0];
                std::vector<int> bctria(3);
                bctria[0] = prismStored1[0];
                bctria[1] = prismStored1[1];
                bctria[2] = prismStored1[2];
                mesh_topology_bl->bcTria[ref0].push_back(bctria);
            }
            
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[5]);
            P1->LocalFace2GlobalNode[1].push_back(prismStored1[4]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            tria1.insert(prismStored1[3]);
            tria1.insert(prismStored1[5]);
            tria1.insert(prismStored1[4]);
            if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
            {
                int ref1 = us3d->tria_ref_map[tria1];
                std::vector<int> bctria(3);
                bctria[0] = prismStored1[3];
                bctria[1] = prismStored1[5];
                bctria[2] = prismStored1[4];
                mesh_topology_bl->bcTria[ref1].push_back(bctria);
            }
            
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[2].push_back(prismStored1[3]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad0.insert(prismStored1[0]);
            quad0.insert(prismStored1[1]);
            quad0.insert(prismStored1[4]);
            quad0.insert(prismStored1[3]);
            if(us3d->quad_ref_map.find(quad0)!=us3d->quad_ref_map.end())
            {
                int ref0 = us3d->quad_ref_map[quad0];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[0];
                bcquad[1] = prismStored1[1];
                bcquad[2] = prismStored1[4];
                bcquad[3] = prismStored1[3];
                mesh_topology_bl->bcQuad[ref0].push_back(bcquad);

            }
            
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[1]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[2]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[4]);
            P1->LocalFace2GlobalNode[3].push_back(prismStored1[5]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad1.insert(prismStored1[1]);
            quad1.insert(prismStored1[2]);
            quad1.insert(prismStored1[4]);
            quad1.insert(prismStored1[5]);
            
            if(us3d->quad_ref_map.find(quad1)!=us3d->quad_ref_map.end())
            {
                int ref1 = us3d->quad_ref_map[quad1];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[1];
                bcquad[1] = prismStored1[2];
                bcquad[2] = prismStored1[4];
                bcquad[3] = prismStored1[5];
                mesh_topology_bl->bcQuad[ref1].push_back(bcquad);

            }
            
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[2]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[0]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[3]);
            P1->LocalFace2GlobalNode[4].push_back(prismStored1[5]);
            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
            quad2.insert(prismStored1[2]);
            quad2.insert(prismStored1[0]);
            quad2.insert(prismStored1[3]);
            quad2.insert(prismStored1[5]);
            if(us3d->quad_ref_map.find(quad2)!=us3d->quad_ref_map.end())
            {
                int ref2 = us3d->quad_ref_map[quad2];
                std::vector<int> bcquad(4);
                bcquad[0] = prismStored1[2];
                bcquad[1] = prismStored1[0];
                bcquad[2] = prismStored1[3];
                bcquad[3] = prismStored1[5];
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

            prism0.clear();
            prism1.clear();
                    
            // Store the same initial triangles as the determined opposite triangles.
            // However change orientation of the nodes.
            prism0.push_back(opposite_tri[0]);
            prism0.push_back(opposite_tri[1]);
            prism0.push_back(opposite_tri[2]);
            prism1.push_back(opposite_tri1[0]);
            prism1.push_back(opposite_tri1[1]);
            prism1.push_back(opposite_tri1[2]);
            
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
    std::cout << "cnt_turn " << cnt_turn << std::endl;
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
    OutputBLElementsOnRoot(xcn_g,ien_g,elements,comm,"BL_Root_NEW");
       
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
        }
        
        delete us3d->interior;
        
        Array<double>* xcn_g;
        Array<int>* ief_g;
        Array<int>* ien_g;
        if(world_rank == 0)
        {
            xcn_g = new Array<double>(us3d->xcn->getNglob(),3);
            ief_g = new Array<int>(us3d->ief->getNglob(),6);
            ien_g = new Array<int>(us3d->ien->getNglob(),8);
        }
        else
        {
            xcn_g = new Array<double>(1,1);
            ief_g = new Array<int>(1,1);
            ien_g = new Array<int>(1,1);
        }

        int* ien_nlocs      = new int[world_size];
        int* ien_offsets    = new int[world_size];
        int* ief_nlocs      = new int[world_size];
        int* ief_offsets    = new int[world_size];
        int* xcn_nlocs      = new int[world_size];
        int* xcn_offsets    = new int[world_size];
        
        for(int i=0;i<world_size;i++)
        {
            xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]*3;
            xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;

            ief_nlocs[i]   = ien_pstate->getNlocs()[i]*6;
            ief_offsets[i] = ien_pstate->getOffsets()[i]*6;

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

        MPI_Gatherv(&us3d->ief->data[0],
                    us3d->ief->getNrow()*6,
                    MPI_INT,
                    &ief_g->data[0],
                    ief_nlocs,
                    ief_offsets,
                    MPI_INT, 0, comm);

        MPI_Gatherv(&us3d->ien->data[0],
                    us3d->ien->getNrow()*8,
                    MPI_INT,
                    &ien_g->data[0],
                    ien_nlocs,
                    ien_offsets,
                    MPI_INT, 0, comm);

        if(world_rank == 0)
        {
            int wall_id = 3;
            int nLayer  = 10;
            
            if(nLayer>0)
            {
                int counter = 0;
                Mdata* Md = ReadMetricData();
                BLShellInfo* BLshell = FindOuterShellBoundaryLayerMesh(wall_id, nLayer, us3d,xcn_g,ien_g,ief_g,xcn_pstate,ien_pstate,comm);
                            
                
                //==============================================================================
                //==============================================================================
                //            struct BLShellInfo
                //            {
                //                std::map<int,int> ShellFace2BFace;
                //                std::map<std::set<int>,int> ShellTri2FaceID;
                //                std::map<std::set<int>,std::vector<int> > ShellFaceID2TriID
                //                std::map<int,int> FaceID2TopoType;
                //                std::map<int,std::map<int,int> > ShellFace2ShellVert2OppositeBoundaryVerts;
                //                Array<int>* ShellRef;
                //            };
                //==============================================================================
                //==============================================================================
                int nbHex = Md->arrHex.size();
                int nbPrisms    =  us3d->bnd_face_map[wall_id].size()*(nLayer)*2;
                int nbTets      = (nbHex-us3d->bnd_face_map[wall_id].size()*(nLayer))*6;
                int nbHexsNew   = (nbHex-us3d->bnd_face_map[wall_id].size()*(nLayer));
                
                std::cout << " Initial number of prims " << nbPrisms << std::endl;
                std::cout << " Initial number of tets "  << nbTets << std::endl;
                std::cout << " Initial number of verts " << xcn_g->getNrow() << std::endl;
                
                int ith = 0;
                std::set<int> u_tet_vert;
                std::map<int,int> lv2gv_tet_mesh;
                std::map<int,int> gv2lv_tet_mesh;
                std::map<int,double*> metric_hex2tet;
                Array<int>* ien_hex2tet = new Array<int>(nbHexsNew,8);
                int sv = 0;
                std::vector<int> locTet_verts;
                for(int i=0;i<ien_g->getNrow();i++)
                {
                    if(BLshell->elements_set.find(i)==BLshell->elements_set.end())
                    {
                        for(int j=0;j<8;j++)
                        {
                            int val = ien_g->getVal(i,j);
                            
                            if(u_tet_vert.find(val)==u_tet_vert.end())
                            {
                                u_tet_vert.insert(val);
                                locTet_verts.push_back(val);
                                gv2lv_tet_mesh[val]=sv;
                                lv2gv_tet_mesh[sv]=val;
                                ien_hex2tet->setVal(ith,j,sv);
                                double* met = new double[6];
                                met[0] = Md->Vmetric[val][3];
                                met[1] = Md->Vmetric[val][4];
                                met[2] = Md->Vmetric[val][5];
                                met[3] = Md->Vmetric[val][6];
                                met[4] = Md->Vmetric[val][7];
                                met[5] = Md->Vmetric[val][8];
                                metric_hex2tet[val]=met;
                                sv++;
                            }
                            else
                            {
                                int lv = gv2lv_tet_mesh[val];
                                ien_hex2tet->setVal(ith,j,lv);
                            }
                        }
                        ith++;
                    }
                }
                std::cout << "lv2gv_tet_mesh.SIIIIIZE " << lv2gv_tet_mesh.size() << " " << sv  << std::endl;
                MMG5_pMesh mmgMesh_TET = NULL;
                MMG5_pSol mmgSol_TET   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppMet,&mmgSol_TET,
                MMG5_ARG_end);
                
                int nbVerts_TET=locTet_verts.size();
                std::cout << "nbVerts_TET " << nbVerts_TET << std::endl;
                
                if ( MMG3D_Set_meshSize(mmgMesh_TET,nbVerts_TET,nbHexsNew*6,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh_TET,mmgSol_TET,MMG5_Vertex,mmgMesh_TET->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                
                int cshell = 0;
                int nshell = 0;
                int reffie0 = 0;
                int reffie1 = 0;
                int reffie2 = 0;
                int reffie3 = 0;
                int reffiemindrie = 0;
                int reffiemineen  = 0;
                int reffiezeros   = 0;
                for(int i=0;i<nbVerts_TET;i++)
                {
                    int gv=lv2gv_tet_mesh[i];
                    
                    mmgMesh_TET->point[i+1].c[0] = xcn_g->getVal(gv,0);
                    mmgMesh_TET->point[i+1].c[1] = xcn_g->getVal(gv,1);
                    mmgMesh_TET->point[i+1].c[2] = xcn_g->getVal(gv,2);
                    mmgMesh_TET->point[i+1].ref  = BLshell->ShellRef->getVal(gv,0);
                    if(BLshell->ShellRef->getVal(gv,0)==0)
                    {
                        std::cout << i << " " << gv << " " << xcn_g->getNrow() << " zero here already" <<std::endl;
                    }
                    if(BLshell->ShellRef->getVal(gv,0)==-1)
                    {
                        //std::cout << i << " " << gv << " " << xcn_g->getNrow() << " -1 here already" <<std::endl;

                        cshell++;
                    }
                    if(BLshell->ShellRef->getVal(gv,0)==-3)
                    {
                        //std::cout << i << " " << gv << " " << xcn_g->getNrow() << " -3 here already" <<std::endl;

                        nshell++;
                    }
                    
                    double m11 = Md->Vmetric[gv][3];
                    double m12 = Md->Vmetric[gv][4];
                    double m13 = Md->Vmetric[gv][5];
                    double m22 = Md->Vmetric[gv][6];
                    double m23 = Md->Vmetric[gv][7];
                    double m33 = Md->Vmetric[gv][8];
                    
                    if ( MMG3D_Set_tensorSol(mmgSol_TET, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
                }
                std::cout <<  "CSHELL = " << cshell << std::endl;
                int ref = 20;
                
                int* hexTabNew = new int[9*(nbHexsNew+1)];
                for(int i=0;i<nbHexsNew;i++)
                {
                    int hexTabPosition = 9*(i+1);
                    for(int j=0;j<8;j++)
                    {
                        int val = ien_hex2tet->getVal(i,j)+1;
                        hexTabNew[hexTabPosition+j] = val;
                    }
                    hexTabNew[hexTabPosition+8] = ref;
                }
                     
                int num = H2T_chkorient(mmgMesh_TET,hexTabNew,nbHexsNew);
                std::cout << "ORIENTATION " << num << std::endl;
                int* adjahexNew = NULL;
                adjahexNew = (int*)calloc(6*nbHexsNew+7,sizeof(int));
                assert(adjahexNew);
                 
                if(!H2T_hashHexa(hexTabNew,adjahexNew,nbHexsNew))
                {
                    std::cout << "Error :: setting up the new adjacency for the hexes after reorientation." << std::endl;
                }
                
                Hedge       hed22;
                hed22.size  = 6*nbHexsNew;
                hed22.hnxt  = 6*nbHexsNew;
                hed22.nhmax = (int)(16*6*nbHexsNew);
                hed22.item  = NULL;
                hed22.item  = (hedge*)calloc(hed22.nhmax+1,sizeof(hedge));

                for (int k=6*nbHexsNew; k<hed22.nhmax; k++)
                {
                    hed22.item[k].nxt = k+1;
                }
                
                int nPoints_before_split  = mmgMesh_TET->np;

                std::cout << "Initial number of vertices outside BL mesh: " << nPoints_before_split << std::endl;
                
                int nTets_before_split  = mmgMesh_TET->ne;
                std::cout << "Estimated number of tetrahedra outside BL mesh based on the initial number of hexes: " << nTets_before_split << " = " << nbHexsNew << " x 6"  << std::endl;
                
                int ret = H2T_cuthex(mmgMesh_TET, &hed22, hexTabNew, adjahexNew, nbHexsNew);
                int nel_tets = mmgMesh_TET->ne;
                std::map<int,std::vector<int> > unique_shell_tri_map;

                std::set<std::set<int> > unique_shell_tris;
                int shell_T_id=0;
                int shell_T_id2 = 0;
                int tellertOr = 0;
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref == 20)
                    {
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                            
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        if( mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1 )
                        {
                            std::set<int> shell_tri;
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            shell_tri.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            std::vector<int> tri(3);
                            tri[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                            tri[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
                            tri[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
                            if(unique_shell_tris.find(shell_tri)==unique_shell_tris.end())
                            {
                                unique_shell_tri_map[shell_T_id]=tri;
                                unique_shell_tris.insert(shell_tri);
                                shell_T_id++;
                            }
                            shell_T_id2++;
                            shell_tri.clear();
                        }
                        
                        tellertOr++;
                    }
                }
                
                // {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
                std::cout << "compare sizes " << unique_shell_tris.size() << " " << us3d->bnd_face_map[wall_id].size()*2  << " " << counter << " " << BLshell->ShellTri2FaceID.size()  << " " << shell_T_id << " " << shell_T_id2 << " " << shell_T_id << std::endl;
                
                std::map<std::set<int>,int > shelltri2fid=BLshell->ShellTri2FaceID;
                std::map<int,int> shellFace2bFace=BLshell->ShellFace2BFace;
                std::map<int,int> bFace2shellFace=BLshell->BFace2ShellFace;
                std::set<std::set<int> >::iterator itset;
                std::map<int,std::vector<int> > TriID2ShellFaceID;
                std::vector<std::vector<int> > u_tris(unique_shell_tris.size());
                int teller=0;
                std::set<int> u_shell_verts;
                for(itset=unique_shell_tris.begin();itset!=unique_shell_tris.end();itset++)
                {
                    std::set<int> shell_tri2 = *itset;
                    std::set<int>::iterator itsh;
                    std::vector<int> u_tri_vec(3);
                    int bb = 0;
                    for(itsh=shell_tri2.begin();itsh!=shell_tri2.end();itsh++)
                    {
                        u_tri_vec[bb]=*itsh;
                        //std::cout << *itsh << " ";
                        if(u_shell_verts.find(*itsh)==u_shell_verts.end())
                        {
                            u_shell_verts.insert(*itsh);
                        }
                        bb++;
                    }
                    //std::cout << std::endl;
                    
                    u_tris[teller] = u_tri_vec;
                    int shell_faceid        = shelltri2fid[shell_tri2];
                    int bfaceID             = shellFace2bFace[shell_faceid];
                    int shell_faceid2       = bFace2shellFace[bfaceID];
                    std::map<int,int> v2v   = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid];
                    //=================================================================
    //                std::cout << " shell face " << shell_faceid << " -> " << std::endl;
    //                std::map<int,int>::iterator itm;
    //                for(itm=v2v.begin();itm!=v2v.end();itm++)
    //                {
    //                    std::cout << itm->first << " " << itm->second << std::endl;
    //                }
    //                std::cout << std::endl;
                    //=================================================================
                    TriID2ShellFaceID[teller].push_back(shell_faceid);
                    BLshell->ShellFaceID2TriID[shell_faceid].push_back(teller);
                    teller++;
                }
                
                std::cout << "teller !!!" << teller << " " << BLshell->ShellFaceID2TriID.size() << " p-p "<< BLshell->ShellTri2FaceID.size() << std::endl;
                Mesh_Topology_BL* mesh_topo_bl2 =  ExtractBoundaryLayerMeshFromShell(u_tris, BLshell, wall_id, nLayer, us3d, xcn_g, ien_g, ief_g, xcn_pstate, ien_pstate, comm);
                
                int nTriangles_BL  = 0;
                int nQuads_BL      = 0;
                
                std::map<int,std::vector<std::vector<int> > >::iterator itrr;
                for(itrr=mesh_topo_bl2->bcTria.begin();itrr!=mesh_topo_bl2->bcTria.end();itrr++)
                {
                    nTriangles_BL  = nTriangles_BL+itrr->second.size();
                    //std::cout << " itrr->second.size()  tri " << itrr->second.size() << std::endl;

                }
                for(itrr=mesh_topo_bl2->bcQuad.begin();itrr!=mesh_topo_bl2->bcQuad.end();itrr++)
                {
                    nQuads_BL  = nQuads_BL+itrr->second.size();
                    std::cout << " itrr->second.size()  quad " << itrr->second.size() << std::endl;
                }
                
                std::map<int,int> loc2glob_final_verts;
                std::map<int,int> glob2loc_final_verts;
                std::map<int,std::vector<Element* > >::iterator iter;
                std::set<int> unique_prism_verts;
                int locp = 0;
                for(iter=mesh_topo_bl2->BLlayersElements.begin();
                   iter!=mesh_topo_bl2->BLlayersElements.end();iter++)
                {
                    int numit=iter->second.size();
                    
                    for(int p=0;p<numit;p++)
                    {
                        std::vector<int> prism = iter->second[p]->GlobalNodes;
                        for(int q=0;q<prism.size();q++)
                        {
                            int pp = prism[q];
                            if(unique_prism_verts.find(pp)==unique_prism_verts.end())
                            {
                                unique_prism_verts.insert(pp);
                                loc2glob_final_verts[locp]=pp;
                                glob2loc_final_verts[pp]=locp;
                                locp++;
                            }
                        }
                        
                    }
                }
                
                //OutputMesh_MMG(mmgMesh_TET,offset_el,offset_el,"OuterVolume_TET.dat");
                
                int refer;
                int tellert = 0;
                int tellert2 = 0;
                std::map<int,std::vector< int* > > bound_tet;
                double m11,m12,m13,m22,m23,m33;
                std::set<int> unique_new_verts;
                std::map<int,int> newVID2locID;
                std::map<int,int> locID2newVID;
                int ut = 0;
                int wtel = 0;
                std::map<int,int> loc2glob_hyb;
                std::map<int,int> glob2loc_hyb;
                
                int bndv = 0;
                int tellie = 0;
                int weight = 0;
                std::set<int> uv_or;
                int ytel = 0;
                std::map<int,std::set<int> > newvert2elem;
                std::map<int,std::set<int> > newvert2vert;
                std::map<int,int*> bndtrisVol;
                std::map<int,int> bndtrisVolRef;
                int tra = 0;
                int st = 0;
                int stt = 0;
                int gnt = 0;
                std::cout << "lv2gv_tet_mesh before " << lv2gv_tet_mesh.size() << " " << mmgMesh_TET->np << std::endl;
                int yte = 0;
                int missing = 0;
                std::set<int> setones;
                
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref != 0)
                    {
                        // Determine the vertices that are newly introduced by the tesselation of the hexes into tets.
                        int weight = 0;
                        std::vector<int> which;
                        std::vector<int> indexes;
                        for(int s=0;s<4;s++)
                        {
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == 0)
                            {
                                if (unique_new_verts.find(mmgMesh_TET->tetra[i].v[s])==unique_new_verts.end())
                                {
                                    unique_new_verts.insert(mmgMesh_TET->tetra[i].v[s]);
                                    locID2newVID[mmgMesh_TET->tetra[i].v[s]] = ut;
                                    newVID2locID[ut] = mmgMesh_TET->tetra[i].v[s];
                                    
                                    ut++;
                                    weight++;
                                    indexes.push_back(s);
                                }
                                
                                newvert2elem[mmgMesh_TET->tetra[i].v[s]].insert(i);
                                
                                for(int t=0;t<4;t++)
                                {
                                    if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[t]].ref != 0
                                    && mmgMesh_TET->tetra[i].v[t]!=mmgMesh_TET->tetra[i].v[s])
                                    {
                                        newvert2vert[mmgMesh_TET->tetra[i].v[s]].insert(mmgMesh_TET->tetra[i].v[t]);
                                    }
                                }
                            }
                            else
                            {
                                which.push_back(s);
                            }
                            
                            
                            
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == -1)
                            {
                                if(setones.find(mmgMesh_TET->tetra[i].v[s])==setones.end())
                                {
                                    setones.insert(mmgMesh_TET->tetra[i].v[s]);
                                }
                            }
                        }
                        if(weight != 0)
                        {
                            wtel++;
                        }
                    
                        //===========================================================================
//                        int v0 = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
//                        int v1 = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
//                        int v2 = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
//                        int v3 = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                        //===========================================================================
//                        int* tetras_or= new int[4];
//                        tetras_or[0] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1];
//                        tetras_or[1] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1];
//                        tetras_or[2] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1];
//                        tetras_or[3] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1];
                        int* tetras= new int[4];
                        for(int s=0;s<4;s++)
                        {
                            if(lv2gv_tet_mesh.find(mmgMesh_TET->tetra[i].v[s]-1)!=lv2gv_tet_mesh.end())
                            {
                                tetras[s] = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[s]-1];
                            }
                            else
                            {
                                tetras[s] = -(s+1);
                                std::cout << "mmgMesh_TET->tetra[i].v[s]-1 " << mmgMesh_TET->tetra[i].v[s]-1 << " "  << std::endl;
                            }
                        }

                        int v0=tetras[0];
                        int v1=tetras[1];
                        int v2=tetras[2];
                        int v3=tetras[3];
                        
                        bndv = 0;
                        std::set<int> bface;
                        
                        std::set<int> tria0;
                        
                        tria0.insert(v0);
                        tria0.insert(v2);
                        tria0.insert(v1);
                        //
                        
                        if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
                        {
                            refer = us3d->tria_ref_map[tria0];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v2+1;
                            tria[2] = v1+1;
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        std::set<int> tria1;
                            
                        tria1.insert(v1);
                        tria1.insert(v2);
                        tria1.insert(v3);
                            
                        if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
                        {
                            refer = us3d->tria_ref_map[tria1];
                            int* tria = new int[3];
                            tria[0] = v1+1;
                            tria[1] = v2+1;
                            tria[2] = v3+1;
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        
                        std::set<int> tria2;
                        
                        tria2.insert(v0);
                        tria2.insert(v3);
                        tria2.insert(v2);
                        
                        if(us3d->tria_ref_map.find(tria2)!=us3d->tria_ref_map.end())
                        {
                            refer = us3d->tria_ref_map[tria2];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v3+1;
                            tria[2] = v2+1;
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra] = tria;
                            bndtrisVolRef[tra] = refer;
                            tra++;
                        }
                        
                        std::set<int> tria3;
                                                   
                        tria3.insert(v0);
                        tria3.insert(v1);
                        tria3.insert(v3);
                       
                        if(us3d->tria_ref_map.find(tria3)!=us3d->tria_ref_map.end())
                        {
                           refer = us3d->tria_ref_map[tria3];
                           int* tria = new int[3];
                           tria[0] = v0+1;
                           tria[1] = v1+1;
                           tria[2] = v3+1;
                           bound_tet[refer].push_back(tria);
                           bndtrisVol[tra] = tria;
                           bndtrisVolRef[tra] = refer;
                           tra++;
                        }
                        
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1 )
                        {
                            std::set<int> tria00;
                            tria00.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            tria00.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            tria00.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            if(unique_shell_tris.find(tria00)!=unique_shell_tris.end())
                            {
                                ytel++;
                            }
                            tria00.clear();
                            gnt++;
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1)
                        {
                            std::set<int> tria11;
                            tria11.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            tria11.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            tria11.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            
                            if(unique_shell_tris.find(tria11)!=unique_shell_tris.end())
                            {
                                ytel++;
                            }
                            gnt++;
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==-1)
                        {
                            std::set<int> tria22;
                            tria22.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            tria22.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            tria22.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[2]-1]);
                            
                            if(unique_shell_tris.find(tria22)!=unique_shell_tris.end())
                            {
                                ytel++;
                            }
                            gnt++;
                        }
                        if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==-1
                           && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==-1)
                        {
                            std::set<int> tria33;
                            tria33.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[0]-1]);
                            tria33.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[1]-1]);
                            tria33.insert(lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[3]-1]);
                            
                            if(unique_shell_tris.find(tria33)!=unique_shell_tris.end())
                            {
                                ytel++;
                            }
                            tria33.clear();
                            gnt++;
                        }
                        
                        tria0.clear();
                        tria1.clear();
                        tria2.clear();
                        tria3.clear();
                                            
                        indexes.clear();
                        tellert++;
                    }
                }
                std::cout << "AFTER@!@@@ " << lv2gv_tet_mesh.size() << std::endl;
                int nTriangles_Vol = 0;
                std::map<int,std::vector< int* > >::iterator itrrb;
                for(itrrb=bound_tet.begin();itrrb!=bound_tet.end();itrrb++)
                {
                    nTriangles_Vol  = nTriangles_Vol+itrrb->second.size();
                    std::cout << itrrb->first << " itrrb->second.size() " << itrrb->second.size() << std::endl;
                }
                
                std::cout << "unique_new_verts " << unique_new_verts.size() << std::endl;
                std::cout << "tellert = " << tellert << " nTets_before_split = " << nTets_before_split << std::endl;
                std::cout << "newly introduced vertices " << unique_new_verts.size() << std::endl;
                std::cout << "initial number of vertices outside of BL mesh " << u_tet_vert.size() << std::endl;
                std::cout << "After number of vertices outside BL mesh " << mmgMesh_TET->np << std::endl;
                std::cout << "initial number of estimated tetrahedra outside of BL mesh " << nbHexsNew*6 << std::endl;
                std::cout << "After number of tetrahedra outside BL mesh " << tellert2 << std::endl;
                std::cout << "Check:" << std::endl;
                std::cout << "=======================================" << std::endl;
                std::cout << "Number of newly introduced vertices are accounted for when --> " << (mmgMesh_TET->np-nPoints_before_split) << " == " << unique_new_verts.size() << std::endl;
                std::cout << "ytel = " << ytel << " " << wtel << std::endl;
                std::cout << nTriangles_BL<<" "<<nTriangles_Vol<<std::endl;
                
//
//                if((double)(tellert-nTets_before_split)/(mmgMesh_TET->np-nPoints_before_split)!=6)
//                {
//                    std::cout << "ERROR :: Number of newly introduced tetrahedra are not accounted since  -->  (double)(tellert-nTets_before_split)/(mmgMesh_TET->np-nPoints_before_split) " << "!=" << 6 << std::endl;
//                    exit(-1);
//                }
//                else
//                {
//                    std::cout << "SUCCES :: Number of newly introduced tetrahedra are accounted for! -> Number of new tetrahedra Nnt = " << (double)(tellert-nTets_before_split) << " and number of new vertices Nnv = " << (mmgMesh_TET->np-nPoints_before_split) << " which results in Nnt/Nnv " << (double)(tellert-nTets_before_split)/(mmgMesh_TET->np-nPoints_before_split) << std::endl;
//                }
//
//                if((nTriangles_BL+nTriangles_Vol)/2+nQuads_BL!=Initial_nBndFaces)
//                {
//                    std::cout << "ERROR :: Boundary triangles are not acounted for since -->  (nTriangles_BL+nTriangles_Vol)/2+nQuads_BL != << Initial_nBndFaces " << std::endl;
//                    exit(-1);
//                }
//                else
//                {
//                    std::cout << "SUCCES :: Boundary triangles are acounted for since -->  (nTriangles_BL+nTriangles_Vol)/2+nQuads_BL = "<<(nTriangles_BL+nTriangles_Vol)/2+nQuads_BL<< " and Initial_nBndFaces = " << Initial_nBndFaces << " and they are are matching!" << std::endl;
//                }
//
                
                std::map<int,std::set<int> >::iterator nve;
                double m11n=0.0;double m12n=0.0;double m13n=0.0;
                double m22n=0.0;double m23n=0.0;
                double m33n=0.0;
                
                std::map<int,double*> newvert2metric;
                for(nve=newvert2vert.begin();nve!=newvert2vert.end();nve++)
                {
                    m11n=0.0;m12n=0.0;m13n=0.0;
                    m22n=0.0;m23n=0.0;
                    m33n=0.0;
                    
                    std::set<int>::iterator nves;
                    for(nves=nve->second.begin();nves!=nve->second.end();nves++)
                    {
                        m11n = m11n+Md->Vmetric[*nves][3];
                        m12n = m12n+Md->Vmetric[*nves][4];
                        m13n = m13n+Md->Vmetric[*nves][5];
                        m22n = m22n+Md->Vmetric[*nves][6];
                        m23n = m23n+Md->Vmetric[*nves][7];
                        m33n = m33n+Md->Vmetric[*nves][8];
                    }
                    
                    double* newM = new double[6];
                    newM[0] = m11n;newM[1] = m12n;newM[2] = m13n;
                    newM[3] = m22n;newM[4] = m23n;newM[5] = m33n;
                    newvert2metric[nve->first]=newM;
                }
                
                //====================================================================
                //====================================================================
                //====================================================================
                //====================================================================
                //              Begin Create hybrid mesh
                //====================================================================
                //====================================================================
                //====================================================================
                //====================================================================
                
                MMG5_pMesh mmgMesh_hyb = NULL;
                MMG5_pSol mmgSol_hyb   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
                MMG5_ARG_end);
                
                int nVertices_New  = mmgMesh_TET->np+unique_prism_verts.size()-cshell;
                int nbTets_New     = tellert;
                
                std::cout << "vertices hyb = " << mmgMesh_TET->np << " " << unique_prism_verts.size() << " " << cshell << std::endl;
                if ( MMG3D_Set_meshSize(mmgMesh_hyb,nVertices_New,nbTets_New,nbPrisms,nTriangles_BL+nTriangles_Vol,nQuads_BL,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh_hyb,mmgSol_hyb,MMG5_Vertex,mmgMesh_hyb->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                std::cout << "nTriangles_BL+nTriangles_Vol " << nTriangles_BL+nTriangles_Vol << std::endl;
                std::cout << "nVert = " << nVertices_New << " nTet = " << nbTets_New << " " << " " << nbPrisms << std::endl;
                std::cout << "vertices hyb compare = " << mmgMesh_TET->np+unique_prism_verts.size()-cshell << " == " << xcn_g->getNrow()+(mmgMesh_TET->np-nPoints_before_split) << " " << bndtrisVol.size() << " " << nTriangles_Vol << " " << mmgMesh_hyb->nt << std::endl;
                
//                for(int i=0;i<nbVertices;i++)
//                {
//                    mmgMesh_hyb->point[i+1].c[0] = Md->Vmetric[i][0];
//                    mmgMesh_hyb->point[i+1].c[1] = Md->Vmetric[i][1];
//                    mmgMesh_hyb->point[i+1].c[2] = Md->Vmetric[i][2];
//
//                    mmgMesh_hyb->point[i+1].ref = BLshell->ShellRef->getVal(i,0);
//
//                    double m11 = Md->Vmetric[i][3];
//                    double m12 = Md->Vmetric[i][4];
//                    double m13 = Md->Vmetric[i][5];
//                    double m22 = Md->Vmetric[i][6];
//                    double m23 = Md->Vmetric[i][7];
//                    double m33 = Md->Vmetric[i][8];
//
//                    if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
//                }
                
                
                int i   = 0;
                int tt  = 1;
                int qt  = 1;
                int qt0 = 0;
                int fnt = 0;
                int qfid = 0;
                int qfidL = 0;
                int qfidR = 0;
                int tfid = 0;
                
                std::map<std::set<int>, int> tfacesmap;
                std::set<std::set<int> >tfaces;
                std::map<int,int> tlh;
                std::map<int,int> trh;
                
                std::map<std::set<int>, int> qfacesmap;
                std::set<std::set<int> >qfaces;
                std::map<int,int> qlh;
                std::map<int,int> qrh;
                for(iter=mesh_topo_bl2->BLlayersElements.begin();
                   iter!=mesh_topo_bl2->BLlayersElements.end();iter++)
                {
                    int numit=iter->second.size();
                    
                    for(int p=0;p<numit;p++)
                    {
                        std::vector<int> prism = iter->second[p]->GlobalNodes;

                        mmgMesh_hyb->prism[i+1].v[0] = prism[0]+1;
                        m11 = Md->Vmetric[prism[0]][3];
                        m12 = Md->Vmetric[prism[0]][4];
                        m13 = Md->Vmetric[prism[0]][5];
                        m22 = Md->Vmetric[prism[0]][6];
                        m23 = Md->Vmetric[prism[0]][7];
                        m33 = Md->Vmetric[prism[0]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[0]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[1] = prism[1]+1;
                        m11 = Md->Vmetric[prism[1]][3];
                        m12 = Md->Vmetric[prism[1]][4];
                        m13 = Md->Vmetric[prism[1]][5];
                        m22 = Md->Vmetric[prism[1]][6];
                        m23 = Md->Vmetric[prism[1]][7];
                        m33 = Md->Vmetric[prism[1]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[1]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[2] = prism[2]+1;
                        m11 = Md->Vmetric[prism[2]][3];
                        m12 = Md->Vmetric[prism[2]][4];
                        m13 = Md->Vmetric[prism[2]][5];
                        m22 = Md->Vmetric[prism[2]][6];
                        m23 = Md->Vmetric[prism[2]][7];
                        m33 = Md->Vmetric[prism[2]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[2]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[3] = prism[3]+1;
                        m11 = Md->Vmetric[prism[3]][3];
                        m12 = Md->Vmetric[prism[3]][4];
                        m13 = Md->Vmetric[prism[3]][5];
                        m22 = Md->Vmetric[prism[3]][6];
                        m23 = Md->Vmetric[prism[3]][7];
                        m33 = Md->Vmetric[prism[3]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[3]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[4] = prism[4]+1;
                        m11 = Md->Vmetric[prism[4]][3];
                        m12 = Md->Vmetric[prism[4]][4];
                        m13 = Md->Vmetric[prism[4]][5];
                        m22 = Md->Vmetric[prism[4]][6];
                        m23 = Md->Vmetric[prism[4]][7];
                        m33 = Md->Vmetric[prism[4]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[4]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].v[5] = prism[5]+1;
                        m11 = Md->Vmetric[prism[5]][3];
                        m12 = Md->Vmetric[prism[5]][4];
                        m13 = Md->Vmetric[prism[5]][5];
                        m22 = Md->Vmetric[prism[5]][6];
                        m23 = Md->Vmetric[prism[5]][7];
                        m33 = Md->Vmetric[prism[5]][8];
                        if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,prism[5]+1) != 1 ) exit(EXIT_FAILURE);
                        
                        mmgMesh_hyb->prism[i+1].ref  = 0;
                        // std::cout << prism[0] << " " << prism[1] << " " << prism[2] << " " << prism[3] << " " << prism[4] << " " << prism[5] << std::endl;
                        std::set<int> tria0;
                        std::set<int> tria1;
                        std::set<int> quad0;
                        std::set<int> quad1;
                        std::set<int> quad2;
                        
                        tria0.insert(prism[0]);
                        tria0.insert(prism[1]);
                        tria0.insert(prism[2]);
                        
                        if(tfaces.find(tria0)==tfaces.end())
                        {
                            tfaces.insert(tria0);
                            tfacesmap[tria0] = tfid;
                            tlh[tfid]=i;
                            tfid++;
                        }
                        else
                        {
                            trh[tfacesmap[tria0]]=i;
                        }
                        
                        if(unique_shell_tris.find(tria0)!=unique_shell_tris.end())
                        {
                            fnt++;
                        }
                        
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(us3d->tria_ref_map.find(tria0)!=us3d->tria_ref_map.end())
                        {
                            refer = us3d->tria_ref_map[tria0];
                            mmgMesh_hyb->tria[tt].v[0] = prism[0]+1;
                            mmgMesh_hyb->tria[tt].v[1] = prism[1]+1;
                            mmgMesh_hyb->tria[tt].v[2] = prism[2]+1;
                            mmgMesh_hyb->tria[tt].ref  = refer;
                            tt++;
                        }
                        tria1.insert(prism[3]);
                        tria1.insert(prism[4]);
                        tria1.insert(prism[5]);
                        if(tfaces.find(tria1)==tfaces.end())
                        {
                            tfaces.insert(tria1);
                            tfacesmap[tria1] = tfid;
                            tlh[tfid]=i;
                            tfid++;
                        }
                        else
                        {
                            trh[tfacesmap[tria0]]=i;
                        }
                        if(unique_shell_tris.find(tria1)!=unique_shell_tris.end())
                        {
                            fnt++;
                        }
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(us3d->tria_ref_map.find(tria1)!=us3d->tria_ref_map.end())
                        {
                            refer = us3d->tria_ref_map[tria1];
                            mmgMesh_hyb->tria[tt].v[0] = prism[3]+1;
                            mmgMesh_hyb->tria[tt].v[1] = prism[4]+1;
                            mmgMesh_hyb->tria[tt].v[2] = prism[5]+1;
                            mmgMesh_hyb->tria[tt].ref  = refer;
                            tt++;
                        }
                        
                        
                        
                        quad0.insert(prism[0]);
                        quad0.insert(prism[1]);
                        quad0.insert(prism[4]);
                        quad0.insert(prism[3]);
                        
                        if(qfaces.find(quad0)==qfaces.end())
                        {
                            qfaces.insert(quad0);
//                            qfacesmap[quad0] = qfid;
//                            qlh[qfid]=i;
//                            qfid++;
                            qfidL++;
                        }
                        else
                        {
                            //qrh[qfacesmap[quad0]]=i;
                            qfidR++;
                        }
                        
                        
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(us3d->quad_ref_map.find(quad0)!=us3d->quad_ref_map.end())
                        {
                            refer = us3d->quad_ref_map[quad0];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[1]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                            qt0++;
                        }
                        quad1.insert(prism[1]);
                        quad1.insert(prism[2]);
                        quad1.insert(prism[5]);
                        quad1.insert(prism[4]);
                        
                        if(qfaces.find(quad1)==qfaces.end())
                        {
                            qfaces.insert(quad1);
//                            qfacesmap[quad1] = qfid;
//                            qlh[qfid]=i;
//                            qfid++;
                            qfidL++;
                        }
                        else
                        {
                            //qrh[qfacesmap[quad1]]=i;
                            qfidR++;
                        }
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(us3d->quad_ref_map.find(quad1)!=us3d->quad_ref_map.end())
                        {
                            refer = us3d->quad_ref_map[quad1];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[1]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                            qt0++;
                        }
                        quad2.insert(prism[0]);
                        quad2.insert(prism[3]);
                        quad2.insert(prism[5]);
                        quad2.insert(prism[2]);
                        if(qfaces.find(quad2)==qfaces.end())
                        {
                            qfaces.insert(quad2);
                            //qfacesmap[quad2] = qfid;
                            //qlh[qfid]=i;
                            qfidL++;
                        }
                        else
                        {
                            //qrh[qfacesmap[quad2]]=i;
                            qfidR++;
                        }
                        // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
                        if(us3d->quad_ref_map.find(quad2)!=us3d->quad_ref_map.end())
                        {
                            refer = us3d->quad_ref_map[quad2];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                            qt0++;
                        }
                        
                        tria0.clear();
                        tria1.clear();
                        quad0.clear();
                        quad1.clear();
                        quad2.clear();
                        i++;
                    }
                }
                
                std::cout << "FNT = " << fnt << " " << qfidL << " " << qfidR << " " << tlh.size() << " " << trh.size() << " qt " << qt << " " << qt0 << " " << qfid << " " << i << std::endl;
                
                for(int i=0;i<bndtrisVol.size();i++)
                {
                    mmgMesh_hyb->tria[tt].v[0] = bndtrisVol[i][0];
                    mmgMesh_hyb->tria[tt].v[1] = bndtrisVol[i][1];
                    mmgMesh_hyb->tria[tt].v[2] = bndtrisVol[i][2];
                    
                    mmgMesh_hyb->tria[tt].ref = bndtrisVolRef[i];
                    tt++;
                }
                
                
                
                
                // tets first then prisms.
                // verts building up the tets first then the verts building up the prisms.
             
                 
                int tet = 0;
                int off_v=0;
                int cntr0 = 0;
                int cntr1 = 0;
                std::set<int> cntr00set;
                std::set<int> cntr0set;
                std::set<int> cntr1set;
                int www = 0;
                int nbVertices = xcn_g->getNrow();
                for(int i=0;i<nbVertices;i++)
                {
                    mmgMesh_hyb->point[i+1].c[0] = xcn_g->getVal(i,0);
                    mmgMesh_hyb->point[i+1].c[1] = xcn_g->getVal(i,1);
                    mmgMesh_hyb->point[i+1].c[2] = xcn_g->getVal(i,2);
                    
                    m11 = Md->Vmetric[i][3];
                    m12 = Md->Vmetric[i][4];
                    m13 = Md->Vmetric[i][5];
                    m22 = Md->Vmetric[i][6];
                    m23 = Md->Vmetric[i][7];
                    m33 = Md->Vmetric[i][8];
                    
                    if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
                }
                
                std::map<int,double*>::iterator v2m;
                std::map<int,int> locNew2globNew;
                std::map<int,int> globNew2locNew;
                int p = 0;
                for(v2m=newvert2metric.begin();v2m!=newvert2metric.end();v2m++)
                {
                    mmgMesh_hyb->point[nbVertices+p+1].c[0] = mmgMesh_TET->point[v2m->first].c[0];
                    mmgMesh_hyb->point[nbVertices+p+1].c[1] = mmgMesh_TET->point[v2m->first].c[1];
                    mmgMesh_hyb->point[nbVertices+p+1].c[2] = mmgMesh_TET->point[v2m->first].c[2];
                    
                    m11 = v2m->second[0];
                    m12 = v2m->second[1];
                    m13 = v2m->second[2];
                    m22 = v2m->second[3];
                    m23 = v2m->second[4];
                    m33 = v2m->second[5];
                    
                    locNew2globNew[nbVertices+p+1] = v2m->first;
                    globNew2locNew[v2m->first]     = nbVertices+p+1;
                    if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,nbVertices+p+1) != 1 ) exit(EXIT_FAILURE);

                    p++;
                }
                
                int here=0;
                int nel_tets4real = 0;
                int fset_cnt = 0;
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref == 20)
                    {
//                        int v0 = mmgMesh_TET->tetra[i].v[0];
//                        int v1 = mmgMesh_TET->tetra[i].v[1];
//                        int v2 = mmgMesh_TET->tetra[i].v[2];
//                        int v3 = mmgMesh_TET->tetra[i].v[3];
//
//                        mmgMesh_hyb->tetra[tet].v[0] = lv2gv_tet_mesh[v0-1];
//                        mmgMesh_hyb->tetra[tet].v[1] = lv2gv_tet_mesh[v1-1];
//                        mmgMesh_hyb->tetra[tet].v[2] = lv2gv_tet_mesh[v2-1];
//                        mmgMesh_hyb->tetra[tet].v[3] = lv2gv_tet_mesh[v3-1];
                        
                        tet++;
                        int* tetra = new int[4];
                        
                        std::set<int> face00;
                        std::set<int> face11;
                        std::set<int> face22;
                        std::set<int> face33;
                        
                        for(int s=0;s<4;s++)
                        {
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == 0)
                            {
                                m11 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][0];
                                m12 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][1];
                                m13 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][2];
                                m22 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][3];
                                m23 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][4];
                                m33 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][5];
                                
                                int locNew = globNew2locNew[mmgMesh_TET->tetra[i].v[s]];
                                tetra[s]   = locNew;
                                mmgMesh_hyb->tetra[tet].v[s] = locNew;

                                if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,locNew) != 1 ) exit(EXIT_FAILURE);
                                here++;
                            }
                            else
                            {
                                int vg = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[s]-1];
                                mmgMesh_hyb->tetra[tet].v[s] = vg+1;
                                tetra[s]   = vg;
                                cntr1++;
                                cntr1set.insert(mmgMesh_TET->tetra[i].v[s]);
                                m11 = Md->Vmetric[vg][3];
                                m12 = Md->Vmetric[vg][4];
                                m13 = Md->Vmetric[vg][5];
                                m22 = Md->Vmetric[vg][6];
                                m23 = Md->Vmetric[vg][7];
                                m33 = Md->Vmetric[vg][8];
                                if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,vg+1) != 1 ) exit(EXIT_FAILURE);
                            }
                            
                        }
                        face00.insert(tetra[1]);
                        face00.insert(tetra[2]);
                        face00.insert(tetra[3]);
                        if(unique_shell_tris.find(face00)!=unique_shell_tris.end())
                        {
                            fset_cnt++;
                        }
                        face11.insert(tetra[0]);
                        face11.insert(tetra[2]);
                        face11.insert(tetra[3]);
                        if(unique_shell_tris.find(face11)!=unique_shell_tris.end())
                        {
                            fset_cnt++;
                        }
                        face22.insert(tetra[0]);
                        face22.insert(tetra[1]);
                        face22.insert(tetra[3]);
                        if(unique_shell_tris.find(face22)!=unique_shell_tris.end())
                        {
                            fset_cnt++;
                        }
                        face33.insert(tetra[0]);
                        face33.insert(tetra[1]);
                        face33.insert(tetra[2]);
                        if(unique_shell_tris.find(face33)!=unique_shell_tris.end())
                        {
                            fset_cnt++;
                        }
                        
                        nel_tets4real++;
                    }
                }
                std::cout << "here = " << here << " " << nel_tets << " " << nel_tets4real << " " << fset_cnt << std::endl;
                
//                int rrr = 0;
//                int uuu = 0;
//                double r11,r12,r13,r22,r23,r33;
//                for (i=0; i<mmgMesh_hyb->np; ++i )
//                {
//                    if ( MMG3D_Get_tensorSol(mmgSol_hyb,&r11,&r12,&r13,&r22,&r23,&r33) != 1) exit(EXIT_FAILURE);
//                }
                
//                int iet = MMG3D_Get_tensorSol(mmgSol_hyb,s11,s12,s13,s22,s23,s33);
//                for(int i=0;i<mmgMesh_hyb->np;i++)
//                {
//                    std::cout << s11[i] << " " << s12[i]  << " " << s13[i]  << " " << s22[i]  << " " << s23[i] << " " << s33[i] << std::endl;
//                }
                //======================================================================================================
                //======================================================================================================
                //======================================================================================================
                //======================================================================================================
                //             End Create hybrid mesh
                //======================================================================================================
                //======================================================================================================
                //======================================================================================================
                //======================================================================================================
                
                //MatchBoundaryTags(us3d,mmgMesh_hyb,0,nel_tets);
                
                
                
                
                std::cout << "print the triangles " << " " << mmgMesh_hyb->nt << " " << nTriangles_BL+nTriangles_Vol << std::endl;
    //            for(int i=0;i<mmgMesh_hyb->nt;i++)
    //            {
    //                if(mmgMesh_hyb->tria[i].ref==4)
    //                {
    //                    std::cout   << i << " :: " << mmgMesh_hyb->tria[i+1].v[0] << " "
    //                    << mmgMesh_hyb->tria[i+1].v[1] << " "
    //                    << mmgMesh_hyb->tria[i+1].v[2]  << " "
    //                    << mmgMesh_hyb->tria[i].ref << std::endl;
    //                }
    //
    //            }
                
//                std::ofstream myfile;
//                myfile.open("OuterVolume_TET_Metric.dat");
//                myfile << "TITLE=\"new_volume.tec\"" << std::endl;
//                myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"m00\", \"m01\", \"m02\", \"m11\", \"m12\", \"m22\"" << std::endl;
//                myfile <<"ZONE N = " << mmgMesh_TET->np << ", E = " << nel_tets << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
//
//                for(int i=0;i<mmgMesh_hyb->np;i++)
//                {
//                    myfile << mmgMesh_hyb->point[i+1].c[0] << " " <<mmgMesh_hyb->point[i+1].c[1] << " " << mmgMesh_hyb->point[i+1].c[2] << " " << Md->Vmetric[i][3] << " " << Md->Vmetric[i][4] << " " << Md->Vmetric[i][5] << " " << Md->Vmetric[i][6]<< " " << Md->Vmetric[i][7] << " " << Md->Vmetric[i][8] <<  std::endl;
//                }
//                for(int i=1;i<=nel_tets;i++)
//                {
//                    myfile << mmgMesh_hyb->tetra[i].v[0]
//                    << " " << mmgMesh_hyb->tetra[i].v[1]
//                    << " " << mmgMesh_hyb->tetra[i].v[2]
//                    << " " << mmgMesh_hyb->tetra[i].v[3] << std::endl;
//                }
//                myfile.close();
                
                //=========================================================================================
                //=========================================================================================
                //=========================================================================================
                //=========================================================================================
                //=========================================================================================
                //=========================================================================================
                
                //MMG3D_Set_handGivenMesh(mmgMesh_hyb);
                if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, 4.0) != 1 )    exit(EXIT_FAILURE);

                //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
                MMG3D_Set_dparameter( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_DPARAM_hgradreq , -1 );
                
                //int ier = MMG3D_mmg3dlib(mmgMesh_hyb,mmgSol_hyb);

                std::cout << " Final number of prims " << mmgMesh_hyb->nprism << std::endl;
                std::cout << " Final number of tets "  << mmgMesh_hyb->ne << std::endl;
                std::cout << " Final number of verts " << mmgMesh_hyb->np << std::endl;
                
                MMG5_pMesh mmgMesh_TETCOPY = NULL;
                MMG5_pSol mmgSol_TETCOPY   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppMet,&mmgSol_TETCOPY,
                MMG5_ARG_end);
      
                if ( MMG3D_Set_meshSize(mmgMesh_TETCOPY,mmgMesh_hyb->np,mmgMesh_hyb->ne,0,0,0,0) != 1 )  exit(EXIT_FAILURE);
                
                //if ( MMG3D_Set_solSize(mmgMesh_TETCOPY,mmgSol_TETCOPY,MMG5_Vertex,mmgMesh_TETCOPY->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                for(int i=0;i<mmgMesh_hyb->np;i++)
                {
                    mmgMesh_TETCOPY->point[i+1].c[0] = mmgMesh_hyb->point[i+1].c[0];
                    mmgMesh_TETCOPY->point[i+1].c[1] = mmgMesh_hyb->point[i+1].c[1];
                    mmgMesh_TETCOPY->point[i+1].c[2] = mmgMesh_hyb->point[i+1].c[2];
                    
                    mmgMesh_TETCOPY->point[i+1].ref  = 0;
                }
                
                for(int i=0;i<mmgMesh_hyb->ne;i++)
                {
                    mmgMesh_TETCOPY->tetra[i+1].v[0] = mmgMesh_hyb->tetra[i+1].v[0];
                    mmgMesh_TETCOPY->tetra[i+1].v[1] = mmgMesh_hyb->tetra[i+1].v[1];
                    mmgMesh_TETCOPY->tetra[i+1].v[2] = mmgMesh_hyb->tetra[i+1].v[2];
                    mmgMesh_TETCOPY->tetra[i+1].v[3] = mmgMesh_hyb->tetra[i+1].v[3];
                    mmgMesh_TETCOPY->tetra[i+1].ref  = 0;
                }
                
                //OutputMesh_MMG(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
                
                WriteUS3DGridFromMMG(mmgMesh_hyb, us3d, unique_shell_tris);
                
                //
            }
            else
            {
                Mdata* Md = ReadMetricData();
                
                MMG5_pMesh mmgMesh = NULL;
                MMG5_pSol mmgSol   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                MMG5_ARG_end);
                
                int nbHex      = ien_g->getNrow();
                int nbVertices = xcn_g->getNrow();
                int nbTriangles = us3d->tria_ref_map.size();
                
                
                if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,nbTriangles,0,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
                
                for(int i=0;i<nbVertices;i++)
                {
                    mmgMesh->point[i+1].c[0] = xcn_g->getVal(i,0);
                    mmgMesh->point[i+1].c[1] = xcn_g->getVal(i,1);
                    mmgMesh->point[i+1].c[2] = xcn_g->getVal(i,2);
                    
                    mmgMesh->point[i+1].ref  = 0;
                    
                    double m11 = Md->Vmetric[i][3];
                    double m12 = Md->Vmetric[i][4];
                    double m13 = Md->Vmetric[i][5];
                    double m22 = Md->Vmetric[i][6];
                    double m23 = Md->Vmetric[i][7];
                    double m33 = Md->Vmetric[i][8];
                    
                    if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
                }
                
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
                std::cout << "mmgMesh->ne " << mmgMesh->nt << std::endl;

                int ref0,ref1,ref2,ref3;
                int tel = 0;
                std::set<std::set<int> > tria_unique;
                int offset_NE = (int)mmgMesh->ne/2;
                int t = 1;
                for(int i=1;i<=offset_NE;i++)
                {
                    std::set<int> tria0;
                    std::set<int> tria1;
                    std::set<int> tria2;
                    std::set<int> tria3;
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
//
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
//
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
//
                    tria0.clear();
                    tria1.clear();
                    tria2.clear();
                    tria3.clear();
                }
                
                MMG3D_Set_handGivenMesh(mmgMesh);
                
                if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 4.0) != 1 )    exit(EXIT_FAILURE);

                //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
                MMG3D_Set_dparameter( mmgMesh,  mmgSol,  MMG3D_DPARAM_hgradreq , -1 );
                
                int ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);
                
                OutputMesh_MMG(mmgMesh,0,mmgMesh->ne,"OuterVolumeFull.dat");
                
                //WriteUS3DGridFromMMG(mmgMesh, us3d);
            }
            

        //=================================================================================================
        //=================================================================================================
        //=================================================================================================
        //=================================================================================================
        //=================================================================================================
        //=================================================================================================
            
            
        }
        
        

        
        
        
        
        
        
        
        
        
        
        
        /*
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
        Array<double>* dUdXi = ComputedUdx_LSQ_US3D_v3(P,UauxNew,meshTopo,gB,comm);
//        Array<double>* dUdXi   = ComputedUdx_MGG(P,UauxNew,meshTopo,gB,comm);
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
        
        Array<double>* dU2dXi2 = ComputedUdx_LSQ_US3D_v3(P,dUdxauxNew,meshTopo,gB,comm);
        Array<double>* dU2dYi2 = ComputedUdx_LSQ_US3D_v3(P,dUdyauxNew,meshTopo,gB,comm);
        Array<double>* dU2dZi2 = ComputedUdx_LSQ_US3D_v3(P,dUdzauxNew,meshTopo,gB,comm);

//        Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,gB,comm);
//        Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,gB,comm);
//        Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,gB,comm);
        
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
        */
        MPI_Finalize();
        
    }
     
    return 0;
}
