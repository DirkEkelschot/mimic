#include "src/adapt_io.h"
#include "src/adapt_recongrad.h"
#include "src/adapt_output.h"
#include "src/adapt_geometry.h"
#include "src/hex2tet.h"
#include <iomanip>


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
        else
        {
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
    fin.open("metric.dat");

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
    finhex.open("elements.dat");

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

//MMG5_pMesh ReadMMG_pMesh(US3D* us3d, MPI_Comm comm, MPI_Info info)
//{
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//    int i,j;
//
//    MMG5_pMesh mmgMesh = NULL;
//    MMG5_pSol mmgSol   = NULL;
//
//    MMG3D_Init_mesh(MMG5_ARG_start,
//    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
//    MMG5_ARG_end);
//
//    std::ifstream fin;
//    fin.open("restart/metric_restart.dat");
//
//    // Read the file row by row
//    std::vector<double> row(9);
//    std::vector<std::vector<double> > Vmetric;
//    int t=0;
//    while(fin >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5] >> row[6] >> row[7] >> row[8])
//    {
//       Vmetric.push_back(row);
//       t++;
//    }
//    int L = Vmetric.size();
//    std::ifstream finhex;
//    finhex.open("restart/elements_restart.dat");
//
//    // Read the file row by row
//    std::vector<int> rowHex(8);
//    std::vector<std::vector<int> > arrHex;
//
//    while(finhex >> rowHex[0] >> rowHex[1] >> rowHex[2] >> rowHex[3] >> rowHex[4] >> rowHex[5] >> rowHex[6] >> rowHex[7])
//    {
//       arrHex.push_back(rowHex);
//    }
//
//    int nbHex      = arrHex.size();
//    int nbVertices = Vmetric.size();
//    int  nbTriangles = tria_ref_map.size()/2;
//    if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,nbTriangles,0,0) != 1 )  exit(EXIT_FAILURE);
//
//    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
//
//    for(int i=0;i<nbVertices;i++)
//    {
//        mmgMesh->point[i+1].c[0] = Vmetric[i][0];
//        mmgMesh->point[i+1].c[1] = Vmetric[i][1];
//        mmgMesh->point[i+1].c[2] = Vmetric[i][2];
//        mmgMesh->point[i+1].ref  = 1;
//
//        double m11 = Vmetric[i][3];
//        double m12 = Vmetric[i][4];
//        double m13 = Vmetric[i][5];
//        double m22 = Vmetric[i][6];
//        double m23 = Vmetric[i][7];
//        double m33 = Vmetric[i][8];
//        if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
//    }
//
//    int* hexTab = new int[9*(nbHex+1)];
//    int ref = 0;
//    for(int i=0;i<nbHex;i++)
//    {
//        int hexTabPosition = 9*(i+1);
//        for(int j=0;j<8;j++)
//        {
//            //int val = ien->getVal(i,j+1);
//            int val = arrHex[i][j];
//            hexTab[hexTabPosition+j] = val;
//        }
//        hexTab[hexTabPosition+8] = ref;
//    }
//
//    int num = H2T_chkorient(mmgMesh,hexTab,nbHex);
//
//    int* adjahex = NULL;
//    adjahex = (int*)calloc(6*nbHex+7,sizeof(int));
//    assert(adjahex);
//    //
//    if(!H2T_hashHexa(hexTab,adjahex,nbHex))
//    {
//        std::cout << "Error :: setting up the new adjacency for the hexes after reorientation." << std::endl;
//    }
//
//    Hedge        hed2;
//    hed2.size  = 6*nbHex;
//    hed2.hnxt  = 6*nbHex;
//    hed2.nhmax = (int)(16*6*nbHex);
//    hed2.item  = NULL;
//    hed2.item  = (hedge*)calloc(hed2.nhmax+1,sizeof(hedge));
//
//    for (int k=6*nbHex; k<hed2.nhmax; k++)
//    {
//        hed2.item[k].nxt = k+1;
//    }
//    int ret = H2T_cuthex(mmgMesh, &hed2, hexTab, adjahex, nbHex);
//    // allocate boundary ids:
//
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //Begin of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    std::set<int> tria0;
//    std::set<int> tria1;
//    std::set<int> tria2;
//    std::set<int> tria3;
//    int ref0,ref1,ref2,ref3;
//    int tel = 0;
//    std::set<std::set<int> > tria_unique;
//    int offset_NE = (int)mmgMesh->ne/2;
//    t = 1;
//    for(int i=1;i<=offset_NE;i++)
//    {
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
//        {
//            ref0 = tria_ref_map[tria0];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].ref  = ref0;
//            t++;
//        }
//
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
//        {
//            ref1 = tria_ref_map[tria1];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].ref  = ref1;
//            t++;
//        }
//
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        if(tria_ref_map.find(tria2)!=tria_ref_map.end())
//        {
//            ref2 = tria_ref_map[tria2];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].ref  = ref2;
//            t++;
//        }
//
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        if(tria_ref_map.find(tria3)!=tria_ref_map.end())
//        {
//            ref3 = tria_ref_map[tria3];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].ref  = ref3;
//            t++;
//        }
//
//        tria0.clear();
//        tria1.clear();
//        tria2.clear();
//        tria3.clear();
//    }
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //End of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//
//    MMG3D_Set_handGivenMesh(mmgMesh);
//
////    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 1.5) != 1 )
////    exit(EXIT_FAILURE);
//
//    int ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);
//
//    OutputMesh_MMG(mmgMesh,0,mmgMesh->ne,"MMgOutput.dat");
//    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,1);
//    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,2);
//    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,3);
//    //    OutputBoundaryID_MMG(mmgMesh,ref2bface,4);
//    //    WriteUS3DGridFromMMG(mmgMesh, us3d);
//
//
//    return mmgMesh;
//}



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

    
    
    
    //==================OUTPUT ORIGINAL MESH=======================
    //==================OUTPUT ORIGINAL MESH=======================
  
    
    if(world_rank == 0)
    {
        string filename2 = "metric.dat";
        ofstream myfile2;
        myfile2.open(filename2);
        
        string filename3 = "elements.dat";
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


//void MatchBoundaryTags(US3D* us3d, MMG5_pMesh mmgMesh,int offset_NE, int Nel)
//{
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //Begin of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//    std::set<int> tria0;
//    std::set<int> tria1;
//    std::set<int> tria2;
//    std::set<int> tria3;
//    int ref0,ref1,ref2,ref3;
//    int tel = 0;
//    std::set<std::set<int> > tria_unique;
//    int t = 1;
//    // local face2vert_map for a tet in mmg  {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
//
//    for(int i=1;i<=offset_NE;i++)
//    {
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        tria0.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
//        {
//            ref0 = tria_ref_map[tria0];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].ref  = ref0;
//            t++;
//        }
//
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        tria1.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
//        {
//            ref1 = tria_ref_map[tria1];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].ref  = ref1;
//            t++;
//        }
//
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        tria2.insert(mmgMesh->tetra[offset_NE+i].v[2]-1);
//        if(tria_ref_map.find(tria2)!=tria_ref_map.end())
//        {
//            ref2 = tria_ref_map[tria2];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[2];
//            mmgMesh->tria[t].ref  = ref2;
//            t++;
//        }
//
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
//        tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
//        if(tria_ref_map.find(tria3)!=tria_ref_map.end())
//        {
//            ref3 = tria_ref_map[tria3];
//            mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[0];
//            mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[1];
//            mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[3];
//            mmgMesh->tria[t].ref  = ref3;
//            t++;
//        }
//        tria0.clear();
//        tria1.clear();
//        tria2.clear();
//        tria3.clear();
//    }
//
//    // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//    int nPrism = mmgMesh->nprism;
//
//    if (nPrism != 0)
//    {
//        std::set<int> quad0;
//        std::set<int> quad1;
//        std::set<int> quad2;
//        int tq = 1;
//        for(int i=0;i<nPrism;i++)
//        {
//            int v0 = mmgMesh->prism[i+1].v[0];
//            int v1 = mmgMesh->prism[i+1].v[1];
//            int v2 = mmgMesh->prism[i+1].v[2];
//            int v3 = mmgMesh->prism[i+1].v[3];
//            int v4 = mmgMesh->prism[i+1].v[4];
//            int v5 = mmgMesh->prism[i+1].v[5];
//
//            tria0.insert(v0);
//            tria0.insert(v1);
//            tria0.insert(v2);
//            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//
//            if(tria_ref_map.find(tria0)!=tria_ref_map.end())
//            {
//                ref0 = tria_ref_map[tria0];
//                mmgMesh->tria[t].v[0] = v0;
//                mmgMesh->tria[t].v[1] = v1;
//                mmgMesh->tria[t].v[2] = v2;
//                mmgMesh->tria[t].ref  = ref0;
//                t++;
//            }
//
//            tria1.insert(v3);
//            tria1.insert(v4);
//            tria1.insert(v5);
//            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//
//            if(tria_ref_map.find(tria1)!=tria_ref_map.end())
//            {
//                ref1 = tria_ref_map[tria1];
//                mmgMesh->tria[t].v[0] = v3;
//                mmgMesh->tria[t].v[1] = v5;
//                mmgMesh->tria[t].v[2] = v4;
//                mmgMesh->tria[t].ref  = ref1;
//                t++;
//            }
//
//            quad0.insert(v0);
//            quad0.insert(v3);
//            quad0.insert(v4);
//            quad0.insert(v1);
//            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//
//            if(quad_ref_map.find(quad0)!=quad_ref_map.end())
//            {
//                ref0 = quad_ref_map[quad0];
//                mmgMesh->quadra[tq].v[0] = v0;
//                mmgMesh->quadra[tq].v[1] = v3;
//                mmgMesh->quadra[tq].v[2] = v4;
//                mmgMesh->quadra[tq].v[3] = v1;
//                mmgMesh->quadra[tq].ref  = ref0;
//                std::cout << "ref0 " << ref1 << std::endl;
//
//                tq++;
//            }
//
//            quad1.insert(v1);
//            quad1.insert(v4);
//            quad1.insert(v5);
//            quad1.insert(v2);
//            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//
//            if(quad_ref_map.find(quad1)!=quad_ref_map.end())
//            {
//                ref1 = quad_ref_map[quad1];
//                mmgMesh->quadra[tq].v[0] = v1;
//                mmgMesh->quadra[tq].v[1] = v4;
//                mmgMesh->quadra[tq].v[2] = v5;
//                mmgMesh->quadra[tq].v[3] = v2;
//                mmgMesh->quadra[tq].ref  = ref1;
//                std::cout << "ref1 " << ref1 << std::endl;
//                tq++;
//            }
//
//            quad2.insert(v2);
//            quad2.insert(v0);
//            quad2.insert(v3);
//            quad2.insert(v5);
//            // local face2vert_map for a prism in mmg {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
//
//            if(quad_ref_map.find(quad2)!=quad_ref_map.end())
//            {
//                ref2 = quad_ref_map[quad2];
//                mmgMesh->quadra[tq].v[0] = v0;
//                mmgMesh->quadra[tq].v[1] = v2;
//                mmgMesh->quadra[tq].v[2] = v5;
//                mmgMesh->quadra[tq].v[3] = v3;
//                mmgMesh->quadra[tq].ref  = ref2;
//                std::cout << "ref2 " << ref2 << std::endl;
//
//                tq++;
//            }
//            tria0.clear();
//            tria1.clear();
//            quad0.clear();
//            quad1.clear();
//            quad2.clear();
//
//        }
//    }
//
//
//
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //End of a hack to match the boundary condition tags of the hexahedral mesh onto the tetrahedral mesh;
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//}



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
    int nb = bnd_face_map[wall_id].size();
    int elid_cur,elid_next;
    int t=0;
    int loc_vid;
    int local_face_id;
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
            
            int changed = ChkHexorient(Pijk,Pijk_id);
            
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



Mesh_Topology_BL* ExtractBoundaryLayerMeshFromShell(std::vector<std::vector<int> > u_tris, BLShellInfo* BLshell, int wall_id, int nLayer, US3D* us3d, Array<double>* xcn_g, Array<int>* ien_g, Array<int>* ief_g, Array<int>* ife_g, Array<int>* ifn_g, ParallelState* xcn_pstate, ParallelState* ien_pstate, std::map<int,std::vector<int> > bnd_face_map, std::map<std::set<int>,int> tria_ref_map, std::map<std::set<int>,int> quad_ref_map,  MPI_Comm comm)
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
    int nb = bnd_face_map[wall_id].size();
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
    
    int or0 = 0;
    Vec3D* cut_dir_face0 = new Vec3D;
    Vec3D* cut_dir_facet0 = new Vec3D;
    Vec3D* cut_dir_facet1 = new Vec3D;
    int bvid,opposite_bvid;
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
        
        int bvid,obvid_i,opposite_bvid;
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
        Vec3D* n_t0_v1 = ComputeSurfaceNormal(v_t0,v_t1);
        double orient_t0_check = DotVec3D(r0,n_t0_v1);
        
        v_t10->c0 = xcn_g->getVal(tri_1n[1],0)-xcn_g->getVal(tri_1n[0],0);
        v_t10->c1 = xcn_g->getVal(tri_1n[1],1)-xcn_g->getVal(tri_1n[0],1);
        v_t10->c2 = xcn_g->getVal(tri_1n[1],2)-xcn_g->getVal(tri_1n[0],2);

        v_t11->c0 = xcn_g->getVal(tri_1n[2],0)-xcn_g->getVal(tri_1n[0],0);
        v_t11->c1 = xcn_g->getVal(tri_1n[2],1)-xcn_g->getVal(tri_1n[0],1);
        v_t11->c2 = xcn_g->getVal(tri_1n[2],2)-xcn_g->getVal(tri_1n[0],2);
        Vec3D* n_t10_v1 = ComputeSurfaceNormal(v_t10,v_t11);
        double orient_t1_check = DotVec3D(r0,n_t10_v1);
        
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
        orient_t0_check = DotVec3D(r0,n_t0_v2);
//
        //n_t0 = n_t0_v2;
        v_t10->c0 = xcn_g->getVal(prism1[1],0)-xcn_g->getVal(prism1[0],0);
        v_t10->c1 = xcn_g->getVal(prism1[1],1)-xcn_g->getVal(prism1[0],1);
        v_t10->c2 = xcn_g->getVal(prism1[1],2)-xcn_g->getVal(prism1[0],2);

        v_t11->c0 = xcn_g->getVal(prism1[2],0)-xcn_g->getVal(prism1[0],0);
        v_t11->c1 = xcn_g->getVal(prism1[2],1)-xcn_g->getVal(prism1[0],1);
        v_t11->c2 = xcn_g->getVal(prism1[2],2)-xcn_g->getVal(prism1[0],2);
        Vec3D* n_t10_v2 = ComputeSurfaceNormal(v_t10,v_t11);
        orient_t1_check = DotVec3D(r0,n_t10_v2);
        
        //n_t10 = n_t10_v2;
//
        if(orient_t0_check > 0 && orient_t1_check > 0)
        {
            //std::cout << "Correct = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
            fwrong++;
        }
        if(orient_t0_check < 0 && orient_t1_check < 0)
        {
            //std::cout << "Incorrect = " << orient_t0_check  << " " << orient_t1_check  << std::endl;
            fright++;
        }
        
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
    
    int debug = 1;
//    if(world_rank == 0)
//    {
//        UnitTestEigenDecomp();
//    }
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
        int varia = 4;
        US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);

        int Nel_part = us3d->ien->getNrow();
        
        // ParallelState is an object that allows the user to get Offsets and Nlocs for that input array.
        
        ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
        ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
        ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);
        ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
        
        Array<double>* Uivar = new Array<double>(Nel_part,1);
        
        for(int i=0;i<Nel_part;i++)
        {
            Uivar->setVal(i,0,us3d->interior->getVal(i,varia));
        }
        
        clock_t t,t1;
        double tmax = 0.0;
        double tn = 0.0;
        t = clock();
        
        // ien -> element2node    map coming from parallel reading.
        // iee -> element2element map coming from parallel reading.
        // ief -> element2face    map coming from parallel reading.
        // ifn -> face2node       map coming from parallel reading.
        // ife -> face2element    map coming from parallel reading.
        
        Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief,
                                     us3d->ifn, us3d->ife, us3d->if_ref,
                                     parmetis_pstate, ien_pstate, ife_pstate,
                                     us3d->xcn, xcn_pstate, Uivar, comm);
        
        double duration = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
        double Ptime = 0.0;
        MPI_Allreduce(&duration, &Ptime, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
        std::cout << "Timing partitioning: " << duration << std::endl;
        }
        
        std::vector<int> LocElem = P->getLocElem();
        std::vector<double> Uvaria  = P->getLocElemVaria();
        
        std::map<int,double> Uvaria_map;
        double UvariaV = 0.0;
        for(int i=0;i<LocElem.size();i++)
        {
            int gid = LocElem[i];
            UvariaV   = Uvaria[i];
            Uvaria_map[gid] = UvariaV;
            
            std::cout << std::setprecision (15) << gid << " " << Uvaria_map[gid] << " " << UvariaV << " " << Uvaria[i] << std::endl;

        }
        
        std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uvaria_map,comm);

        int* bnd_map;
        int nBnd = 4;

        if(world_rank == 0)
        {
            std::cout << "Started creating mesh topology object... " << std::endl;
        }

        Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
        
        if(world_rank == 0)
        {
            std::cout << "Finished creating mesh topology object... "  << std::endl;
        }
        
        if(world_rank == 0)
        {
            std::cout << "Setting the ghost element data... "  << std::endl;
        }
        Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
        for(int i=0;i<us3d->ghost->getNrow();i++)
        {
            gB->setVal(i,0,us3d->ghost->getVal(i,varia));
        }
        t = clock();
        if(world_rank == 0)
        {
            std::cout << "Started reconstructing the gradient... " << std::endl;
        }
        
        std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D_v3(P,UauxNew,meshTopo,gB,comm);
        
        double Gtiming = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
        double Gmax_time = 0.0;
        MPI_Allreduce(&Gtiming, &Gmax_time, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            std::cout << "Finished reconstructing the gradient... " << std::endl;
            std::cout << "Timing gradient reconstruction... " << Gmax_time << std::endl;
        }
        
        Array<int>* lE2gE     = new Array<int>(dUdXi.size(),1);
        Array<double>* dUidxi = new Array<double>(dUdXi.size(),1);
        Array<double>* dUidyi = new Array<double>(dUdXi.size(),1);
        Array<double>* dUidzi = new Array<double>(dUdXi.size(),1);
        
        std::map<int,double> dUidxi_map;
        std::map<int,double> dUidyi_map;
        std::map<int,double> dUidzi_map;
        
        std::map<int,Array<double>* >::iterator grit;
        int i=0;
        for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
        {
            lE2gE->setVal(i,0,grit->first);
            dUidxi->setVal(i,0,grit->second->getVal(0,0));
            dUidyi->setVal(i,0,grit->second->getVal(1,0));
            dUidzi->setVal(i,0,grit->second->getVal(2,0));
            dUidxi_map[grit->first]=grit->second->getVal(0,0);
            dUidyi_map[grit->first]=grit->second->getVal(1,0);
            dUidzi_map[grit->first]=grit->second->getVal(2,0);
            
            i++;
        }
        
        std::map<int,double > dUdxauxNew  = P->CommunicateAdjacentDataUS3D(dUidxi_map,comm);
//        std::map<int,double > dUdyauxNew  = P->CommunicateAdjacentDataUS3D(dUidyi,comm);
//        std::map<int,double > dUdzauxNew  = P->CommunicateAdjacentDataUS3D(dUidzi,comm);
        
        //std::cout << "making sure " << dUdxauxNew.size() << " " << dUdXi.size() << std::endl;
        
        int nlElem = us3d->ien->getNrow();
        int nElem  = us3d->ien->getNglob();
        int nvg    = us3d->xcn->getNglob();
        
        std::vector<double> GuX_loc;
        std::vector<int> vids;
        std::vector<double> Uids;
        std::vector<double> Hids;
        int gid,lid;
        int nval = 6;
        std::set<int> vdone;
        
//
        Array<int>*  lE2gE_g;
        Array<double>*  GuX_g;
        Array<double>*  GuX_gr;
        Array<double>*  GuY_g;
        Array<double>*  GuZ_g;

        int nEl_glob = us3d->ien->getNglob();

        if(world_rank == 0)
        {
            lE2gE_g     = new Array<int>(nEl_glob,1);
            GuX_g       = new Array<double>(nEl_glob,1);
            GuX_gr      = new Array<double>(nEl_glob,1);
            GuY_g       = new Array<double>(nEl_glob,1);
            GuZ_g       = new Array<double>(nEl_glob,1);
        }
        else
        {
            lE2gE_g     = new Array<int>(1,1);
            GuX_g       = new Array<double>(1,1);
            GuX_gr      = new Array<double>(1,1);
            GuY_g       = new Array<double>(1,1);
            GuZ_g       = new Array<double>(1,1);
        }

        int* G_nlocs      = new int[world_size];
        int* red_G_nlocs  = new int[world_size];
        int* G_offsets    = new int[world_size];

        for(i=0;i<world_size;i++)
        {
            G_nlocs[i] = 0;
            
            if(i==world_rank)
            {
                G_nlocs[i] = dUidxi->getNrow();
            }
            else
            {
                G_nlocs[i] = 0;
            }
        }

        MPI_Allreduce(G_nlocs, red_G_nlocs, world_size, MPI_INT, MPI_SUM, comm);
        
        int offset = 0;
        for(i=0;i<world_size;i++)
        {
            G_offsets[i] = offset;
            offset = offset+red_G_nlocs[i];
        }

        MPI_Gatherv(&lE2gE->data[0],
                    lE2gE->getNrow(),
                    MPI_INT,
                    &lE2gE_g->data[0],
                    red_G_nlocs,
                    G_offsets,
                    MPI_INT, 0, comm);
        
        MPI_Gatherv(&dUidxi->data[0],
                    dUidxi->getNrow(),
                    MPI_DOUBLE,
                    &GuX_g->data[0],
                    red_G_nlocs,
                    G_offsets,
                    MPI_DOUBLE, 0, comm);
        
        
        if(world_rank == 0)
        {
            string filename = "GuX.dat";
            ofstream myfile;
            myfile.open(filename);

            for(int i=0;i<nEl_glob;i++)
            {
                int gid = lE2gE_g->getVal(i,0);
                GuX_gr->setVal(gid,0,GuX_g->getVal(i,0));
            }
            
            for(int i=0;i<nEl_glob;i++)
            {
                myfile << std::setprecision(16) <<  GuX_gr->getVal(i,0) << std::endl;
            }

            myfile.close();
        }
        
        
        /**/
        
        //delete dUdXi;
        /*
        Array<double>* dU2dXi2 = ComputedUdx_LSQ_US3D_v3(P,dUdxauxNew,meshTopo,gB,comm);
        Array<double>* dU2dYi2 = ComputedUdx_LSQ_US3D_v3(P,dUdyauxNew,meshTopo,gB,comm);
        Array<double>* dU2dZi2 = ComputedUdx_LSQ_US3D_v3(P,dUdzauxNew,meshTopo,gB,comm);

//      Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,gB,comm);
                
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
                    
        std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
        double max_v = *std::max_element(d2udx2_v.begin(), d2udx2_v.end());
        int dim = 3;
        

        Array<double>* metric = ComputeMetric(Verts,grad,hessian,max_v,loc_elem2verts_loc,us3d->ien->getNrow(),comm,dim);
        
        
        if(world_rank==0)
        {
            std::cout << "Started gathering metric data on rank 0..." <<std::endl;
        }
        MMG_Mesh* mmg = GetOptimizedMMG3DMeshOnRoot(P, us3d, hessian, metric, comm);
        if(world_rank==0)
        {
            std::cout << "Finished gathering metric data on rank 0..." <<std::endl;
        }

        //=================================================================
        //==================Output the data in Tecplot format==============
        //=================================================================
        
        if(debug == 1)
        {
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
        }
        
        
        
        delete us3d->interior;
        
        Array<double>*  xcn_g;
        Array<int>*     ief_g;
        Array<int>*     ien_g;
        Array<int>*     ifn_g;
        Array<int>*     if_ref_g;
        Array<int>*     ife_g;
        
        if(world_rank == 0)
        {
            xcn_g       = new Array<double>(us3d->xcn->getNglob(),3);
            ief_g       = new Array<int>(us3d->ief->getNglob(),6);
            ien_g       = new Array<int>(us3d->ien->getNglob(),8);
            if_ref_g    = new Array<int>(us3d->ifn->getNglob(),1);
            ifn_g       = new Array<int>(us3d->ifn->getNglob(),4);
            ife_g       = new Array<int>(us3d->ifn->getNglob(),2);
        }
        else
        {
            xcn_g    = new Array<double>(1,1);
            ief_g    = new Array<int>(1,1);
            ien_g    = new Array<int>(1,1);
            if_ref_g = new Array<int>(1,1);
            ifn_g    = new Array<int>(1,1);
            ife_g    = new Array<int>(1,1);
        }

        int* ien_nlocs      = new int[world_size];
        int* ien_offsets    = new int[world_size];
        int* ief_nlocs      = new int[world_size];
        int* ief_offsets    = new int[world_size];
        int* xcn_nlocs      = new int[world_size];
        int* xcn_offsets    = new int[world_size];
        int* ifn_nlocs      = new int[world_size];
        int* ifn_offsets    = new int[world_size];
        int* if_ref_nlocs   = new int[world_size];
        int* if_ref_offsets = new int[world_size];
        int* ife_nlocs      = new int[world_size];
        int* ife_offsets    = new int[world_size];
        
        for(int i=0;i<world_size;i++)
        {
            xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]  *3;
            xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;

            ief_nlocs[i]   = ien_pstate->getNlocs()[i]  *6;
            ief_offsets[i] = ien_pstate->getOffsets()[i]*6;

            ien_nlocs[i]   = ien_pstate->getNlocs()[i]  *8;
            ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
            
            ifn_nlocs[i]   = ife_pstate->getNlocs()[i]  *4;
            ifn_offsets[i] = ife_pstate->getOffsets()[i]*4;
            
            if_ref_nlocs[i]   = ife_pstate->getNlocs()[i]  *1;
            if_ref_offsets[i] = ife_pstate->getOffsets()[i]*1;
            
            ife_nlocs[i]   = ife_pstate->getNlocs()[i]  *2;
            ife_offsets[i] = ife_pstate->getOffsets()[i]*2;
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

        MPI_Gatherv(&us3d->ifn->data[0],
                    us3d->ifn->getNrow()*4,
                    MPI_INT,
                    &ifn_g->data[0],
                    ifn_nlocs,
                    ifn_offsets,
                    MPI_INT, 0, comm);
        
        MPI_Gatherv(&us3d->if_ref->data[0],
                    us3d->if_ref->getNrow()*1,
                    MPI_INT,
                    &if_ref_g->data[0],
                    if_ref_nlocs,
                    if_ref_offsets,
                    MPI_INT, 0, comm);

        MPI_Gatherv(&us3d->ife->data[0],
                    us3d->ife->getNrow()*2,
                    MPI_INT,
                    &ife_g->data[0],
                    ife_nlocs,
                    ife_offsets,
                    MPI_INT, 0, comm);
        
        delete P;
        delete d2udx2;
        delete d2udxy;
        delete d2udxz;
        delete d2udyx;
        delete d2udy2;
        delete d2udyz;
        delete d2udzx;
        delete d2udzy;
        delete d2udz2;
        
        if(world_rank == 0)
        {
            std::map<std::set<int>,int> tria_ref_map;
            std::map<std::set<int>,int> quad_ref_map;
            std::map<int,int> vert_ref_map;
            std::set<int> vert_ref_set;
            
            std::set<int> tria0;
            std::set<int> tria00;
            std::set<int> tria1;
            std::set<int> tria11;
            
            std::set<int> quad;
            int faceid;
            int nodeid;
            int nrow_ifn = ifn_g->getNrow();
            //Array<int>* ifn_ref  = new Array<int>(nrow_ifn,1);
            int ref;
            
            std::map<int,std::vector<int> > bnd_face_map;
            std::map<int,std::vector<int> >::iterator bmit;
            
            int r2=0;int r10=0;int r36=0;int r3=0;
            for(int i=0;i<nrow_ifn;i++)
            {
                ref = if_ref_g->getVal(i,0);
                //ifn_ref->setVal(i,0,ref);
                
                faceid = i;
                if(ref != 2)
                {
                    bnd_face_map[ref].push_back(faceid);
                }
                
                for(j=0;j<4;j++)
                {
                    nodeid = ifn_g->getVal(i,j); // This is actually node ID!!!!
                    //ifn_copy->setVal(i,j,ifn_g->getVal(i,j+1)-1);
                    
                    if(ref!=2)
                    {
                        if(vert_ref_set.find(nodeid)==vert_ref_set.end())
                        {
                            vert_ref_set.insert(nodeid);
                            vert_ref_map[nodeid] = if_ref_g->getVal(i,0);
                        }
                    }
                }
                
                tria0.insert(ifn_g->getVal(i,0));
                tria0.insert(ifn_g->getVal(i,1));
                tria0.insert(ifn_g->getVal(i,2));
                
                tria00.insert(ifn_g->getVal(i,0));
                tria00.insert(ifn_g->getVal(i,2));
                tria00.insert(ifn_g->getVal(i,3));
                
                tria1.insert(ifn_g->getVal(i,0));
                tria1.insert(ifn_g->getVal(i,1));
                tria1.insert(ifn_g->getVal(i,3));
                
                tria11.insert(ifn_g->getVal(i,1));
                tria11.insert(ifn_g->getVal(i,2));
                tria11.insert(ifn_g->getVal(i,3));
                
                quad.insert(ifn_g->getVal(i,0));
                quad.insert(ifn_g->getVal(i,1));
                quad.insert(ifn_g->getVal(i,2));
                quad.insert(ifn_g->getVal(i,3));
                
                if(tria_ref_map.find(tria0)==tria_ref_map.end() && ref!=2)
                {
                    tria_ref_map[tria0] = ref;
                    tria_ref_map[tria00] = ref;
                }
                if(tria_ref_map.find(tria1)==tria_ref_map.end() && ref!=2)
                {
                    tria_ref_map[tria1] = ref;
                    tria_ref_map[tria11] = ref;
                }
                
                if(quad_ref_map.find(quad)==quad_ref_map.end() && ref!=2)
                {
                    quad_ref_map[quad] = ref;
                }
                
                tria0.clear();
                tria00.clear();
                tria1.clear();
                tria11.clear();
                quad.clear();
            }
//            for(bmit=bnd_face_map.begin();bmit!=bnd_face_map.end();bmit++)
//            {
//                std::cout << " bmit " << bmit->first << " " << bmit->second.size() << std::endl;
//            }
         
            int wall_id = 3;
            int nLayer  = 230;
            
            if(nLayer>0)
            {
                int counter = 0;
                Mdata* Md = ReadMetricData();
                BLShellInfo* BLshell = FindOuterShellBoundaryLayerMesh(wall_id, nLayer, us3d,xcn_g,ien_g,ief_g,ife_g,ifn_g,xcn_pstate,ien_pstate,bnd_face_map,vert_ref_map,comm);
                        
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
                
                int nbHex       =  Md->arrHex.size();
                int nbPrisms    =  bnd_face_map[wall_id].size()*(nLayer)*2;
                int nbTets      = (nbHex-bnd_face_map[wall_id].size()*(nLayer))*6;
                int nbHexsNew   = (nbHex-bnd_face_map[wall_id].size()*(nLayer));
                
                
			
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
   		
                MMG5_pMesh mmgMesh_TET = NULL;
                MMG5_pSol mmgSol_TET   = NULL;
                
                MMG3D_Init_mesh(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppMet,&mmgSol_TET,
                MMG5_ARG_end);
                
                int nbVerts_TET=locTet_verts.size();

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
                        cshell++;
                    }
                    if(BLshell->ShellRef->getVal(gv,0)==-3)
                    {
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
                std::cout << "Check the orientation and modify when necessary so that we can cut each hex up into 6 tetrahedra..."<<std::endl;
                int num = H2T_chkorient(mmgMesh_TET,hexTabNew,nbHexsNew);
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
		
                int nTets_before_split  = mmgMesh_TET->ne;
                
                std::cout << "Cut each hexahedral up into 6 tetrahedra..."<<std::endl;
                //we need to do this in order to determine the orientation of the triangles at the shell interface and trace back how we need to tesselate the wall boundary into triangles so that they match eachothers orientation.
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
                
                std::map<std::set<int>,int > shelltri2fid=BLshell->ShellTri2FaceID;
                std::map<int,int> shellFace2bFace=BLshell->ShellFace2BFace;
                std::map<int,int> bFace2shellFace=BLshell->BFace2ShellFace;
                std::set<std::set<int> >::iterator itset;
                std::map<int,std::vector<int> > TriID2ShellFaceID;
                std::vector<std::vector<int> > u_tris;
                std::vector<std::vector<int> > u_tris_loc;
                int teller=0;
                std::set<int> unique_shell_vert;
                std::vector<int> U_shell_vert;
                std::map<int,int> glob2loc_shell_vert;
                std::map<int,int> loc2glob_shell_vert;
                int loc_s_v=0;
                for(itset=unique_shell_tris.begin();itset!=unique_shell_tris.end();itset++)
                {
                    std::set<int> shell_tri = *itset;
                    std::set<int>::iterator itsh;
                    std::vector<int> u_tri_vec(3);
                    std::vector<int> u_tri_vec_loc(3);
                    int bb = 0;
                    for(itsh=shell_tri.begin();itsh!=shell_tri.end();itsh++)
                    {
                        u_tri_vec[bb]=*itsh;
                        
                        if(unique_shell_vert.find(*itsh)==unique_shell_vert.end())
                        {
                            //u_tri_vec_loc.push_back(loc_s_v);
                            unique_shell_vert.insert(*itsh);
                            U_shell_vert.push_back(*itsh);
                            loc2glob_shell_vert[loc_s_v]=*itsh;
                            glob2loc_shell_vert[*itsh]=loc_s_v;
                            u_tri_vec_loc[bb]= loc_s_v;
                            loc_s_v++;
                        }
                        else
                        {
                            int local_s_vId = glob2loc_shell_vert[*itsh];
                            u_tri_vec_loc[bb]=local_s_vId;
                        }
                        bb++;
                    }
                    
                    u_tris.push_back(u_tri_vec);
                    u_tris_loc.push_back(u_tri_vec_loc);
                    
                    int shell_faceid        = shelltri2fid[shell_tri];
                    int bfaceID             = shellFace2bFace[shell_faceid];
                    int shell_faceid2       = bFace2shellFace[bfaceID];
                    std::map<int,int> v2v   = BLshell->ShellFace2ShellVert2OppositeBoundaryVerts[shell_faceid];
                    TriID2ShellFaceID[teller].push_back(shell_faceid);
                    BLshell->ShellFaceID2TriID[shell_faceid].push_back(teller);
                    teller++;
                }
                
                std::cout << "Extracting the prismatic boundary layer mesh..." <<std::endl;
                Mesh_Topology_BL* mesh_topo_bl2 =  ExtractBoundaryLayerMeshFromShell(u_tris, BLshell, wall_id, nLayer, us3d, xcn_g, ien_g, ief_g, ife_g, ifn_g, xcn_pstate, ien_pstate, bnd_face_map, tria_ref_map, quad_ref_map, comm);
                std::cout << "Outputting the prismatic boundary layer mesh in ---> BoundaryLayerMesh_0.dat" <<std::endl;
                 OutputBoundaryLayerPrisms(xcn_g, mesh_topo_bl2, comm, "BoundaryLayerMesh_");

                if(debug == 1)
                {
                    std::ofstream myfile_shell;
                    myfile_shell.open("shell.dat");
                    myfile_shell << "TITLE=\"new_volume.tec\"" << std::endl;
                    myfile_shell <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
                    myfile_shell <<"ZONE N = " << U_shell_vert.size() << ", E = " << teller << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;

                    for(int i=0;i<U_shell_vert.size();i++)
                    {
                        myfile_shell << xcn_g->getVal(loc2glob_shell_vert[i],0)
                        << " " << xcn_g->getVal(loc2glob_shell_vert[i],1)
                        << " " << xcn_g->getVal(loc2glob_shell_vert[i],2) << std::endl;
                    }
                    for(int i=0;i<teller;i++)
                    {
                        myfile_shell << u_tris_loc[i][0]+1
                        << " " << u_tris_loc[i][1]+1
                        << " " << u_tris_loc[i][2]+1 << std::endl;
                    }
                    myfile_shell.close();
                }
                
                
                // counting the number of boundary triangles and quads in the BL mesh.
                int nTriangles_BL  = 0;
                int nQuads_BL      = 0;
                std::map<int,std::vector<std::vector<int> > >::iterator itrr;
                for(itrr=mesh_topo_bl2->bcTria.begin();itrr!=mesh_topo_bl2->bcTria.end();itrr++)
                {
                    nTriangles_BL  = nTriangles_BL+itrr->second.size();

                }
                for(itrr=mesh_topo_bl2->bcQuad.begin();itrr!=mesh_topo_bl2->bcQuad.end();itrr++)
                {
                    nQuads_BL  = nQuads_BL+itrr->second.size();
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
                int yte = 0;
                int missing = 0;
                std::set<int> setones;

                double vxc = 0;
                double vyc = 0;
                double vzc = 0;
                
                double vxf = 0;
                double vyf = 0;
                double vzf = 0;
                int negf0=0;
                int negf1=0;
                int negf2=0;
                int negf3=0;
                int neg_bc=0;
                int pos_bc=0;
                std::vector<double> oris;
                
                std::cout << "Defining the tetrahedra mesh..."<<std::endl;
                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref != 0)
                    {
                        // Determine the vertices that are newly introduced by the tesselation of the hexes into tets.
                        int weight = 0;
                        std::vector<int> which;
                        std::vector<int> indexes;
                        vxc=0;vyc=0;vzc=0;
                        for(int s=0;s<4;s++)
                        {
                            vxc = vxc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[0];
                            vyc = vyc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[1];
                            vzc = vzc+mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[2];
                            
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
                        
                        vxc = vxc/4;
                        vyc = vyc/4;
                        vzc = vzc/4;
                        
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
                        
                        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria0];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v2+1;
                            tria[2] = v1+1;
                            
                            vxf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0])/3;
                            
                            vyf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1])/3;
                            
                            vzf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2])/3;
                            
                            Vec3D* r0 = new Vec3D;
                            double r0L = sqrt( (vxf-vxc)*(vxf-vxc)
                                              +(vyf-vyc)*(vyf-vyc)
                                              +(vzf-vzc)*(vzf-vzc));
                            r0->c0 = (vxf-vxc)/r0L;
                            r0->c1 = (vyf-vyc)/r0L;
                            r0->c2 = (vzf-vzc)/r0L;
                            Vec3D* n_f0 = new Vec3D;
                            n_f0->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f0->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f0->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f1 = new Vec3D;
                            n_f1->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f1->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f1->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f = ComputeSurfaceNormal(n_f0,n_f1);
                            double orient0 = DotVec3D(r0,n_f);
                            oris.push_back(orient0);
                            if(orient0<0)
                            {
                                negf0++;
                                neg_bc++;
                            }
                            else{
                                pos_bc++;
                            }
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        std::set<int> tria1;
                            
                        tria1.insert(v1);
                        tria1.insert(v2);
                        tria1.insert(v3);
                            
                        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria1];
                            int* tria = new int[3];
                            tria[0] = v1+1;
                            tria[1] = v2+1;
                            tria[2] = v3+1;
                            
                            vxf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0])/3;

                            vyf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1])/3;

                            vzf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2])/3;

                            Vec3D* r0 = new Vec3D;
                            double r0L = sqrt( (vxf-vxc)*(vxf-vxc)
                                              +(vyf-vyc)*(vyf-vyc)
                                              +(vzf-vzc)*(vzf-vzc));
                            r0->c0 = (vxf-vxc)/r0L;
                            r0->c1 = (vyf-vyc)/r0L;
                            r0->c2 = (vzf-vzc)/r0L;
                            Vec3D* n_f0 = new Vec3D;
                            n_f0->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0];
                            n_f0->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1];
                            n_f0->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2];
                            Vec3D* n_f1 = new Vec3D;
                            n_f1->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0];
                            n_f1->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1];
                            n_f1->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2];
                            Vec3D* n_f = ComputeSurfaceNormal(n_f0,n_f1);
                            double orient0 = DotVec3D(r0,n_f);
                            oris.push_back(orient0);
                            if(orient0<0)
                            {
                                negf1++;
                                neg_bc++;
                            }
                            else{
                                pos_bc++;
                            }
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra]     = tria;
                            bndtrisVolRef[tra]  = refer;
                            tra++;
                        }
                        
                        std::set<int> tria2;
                        
                        tria2.insert(v0);
                        tria2.insert(v3);
                        tria2.insert(v2);
                        
                        if(tria_ref_map.find(tria2)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria2];
                            int* tria = new int[3];
                            tria[0] = v0+1;
                            tria[1] = v3+1;
                            tria[2] = v2+1;
                            
                            vxf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0])/3;

                            vyf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1])/3;

                            vzf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2])/3;

                            Vec3D* r0 = new Vec3D;
                            double r0L = sqrt( (vxf-vxc)*(vxf-vxc)
                                              +(vyf-vyc)*(vyf-vyc)
                                              +(vzf-vzc)*(vzf-vzc));
                            r0->c0 = (vxf-vxc)/r0L;
                            r0->c1 = (vyf-vyc)/r0L;
                            r0->c2 = (vzf-vzc)/r0L;
                            Vec3D* n_f0 = new Vec3D;
                            n_f0->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f0->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f0->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f1 = new Vec3D;
                            n_f1->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f1->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f1->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f = ComputeSurfaceNormal(n_f0,n_f1);
                            double orient0 = DotVec3D(r0,n_f);
                            oris.push_back(orient0);
                            if(orient0<0)
                            {
                                negf2++;
                                neg_bc++;
                            }
                            else{
                                pos_bc++;
                            }
                            bound_tet[refer].push_back(tria);
                            bndtrisVol[tra] = tria;
                            bndtrisVolRef[tra] = refer;
                            tra++;
                        }
                        
                        std::set<int> tria3;
                                                   
                        tria3.insert(v0);
                        tria3.insert(v1);
                        tria3.insert(v3);
                       
                        if(tria_ref_map.find(tria3)!=tria_ref_map.end())
                        {
                           refer = tria_ref_map[tria3];
                           int* tria = new int[3];
                           tria[0] = v0+1;
                           tria[1] = v1+1;
                           tria[2] = v3+1;
                            
                            vxf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0])/3;

                            vyf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1])/3;

                            vzf=(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2]
                            +mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2])/3;

                            Vec3D* r0 = new Vec3D;
                            double r0L = sqrt( (vxf-vxc)*(vxf-vxc)
                                              +(vyf-vyc)*(vyf-vyc)
                                              +(vzf-vzc)*(vzf-vzc));
                            r0->c0 = (vxf-vxc)/r0L;
                            r0->c1 = (vyf-vyc)/r0L;
                            r0->c2 = (vzf-vzc)/r0L;
                            Vec3D* n_f0 = new Vec3D;
                            n_f0->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f0->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f0->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f1 = new Vec3D;
                            n_f1->c0 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[0]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[0];
                            n_f1->c1 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[1]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[1];
                            n_f1->c2 = mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].c[2]-mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].c[2];
                            Vec3D* n_f = ComputeSurfaceNormal(n_f0,n_f1);
                            double orient0 = DotVec3D(r0,n_f);
                            oris.push_back(orient0);
                            if(orient0<0)
                            {
                                negf3++;
                                neg_bc++;
                            }
                            else{
                                pos_bc++;
                            }
                           bound_tet[refer].push_back(tria);
                           bndtrisVol[tra]      = tria;
                           bndtrisVolRef[tra]   = refer;
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
                
                double min_oris = *std::min_element(oris.begin(),oris.end());
                double max_oris = *std::max_element(oris.begin(),oris.end());

                int nTriangles_Vol = 0;
                std::map<int,std::vector< int* > >::iterator itrrb;
                for(itrrb=bound_tet.begin();itrrb!=bound_tet.end();itrrb++)
                {
                    nTriangles_Vol  = nTriangles_Vol+itrrb->second.size();
                }
                
                
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
                    
                    m11n = m11n/nve->second.size();
                    m12n = m12n/nve->second.size();
                    m13n = m13n/nve->second.size();
                    m22n = m22n/nve->second.size();
                    m23n = m23n/nve->second.size();
                    m33n = m33n/nve->second.size();
                    
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
                
                if ( MMG3D_Set_meshSize(mmgMesh_hyb,nVertices_New,nbTets_New,nbPrisms,nTriangles_BL+nTriangles_Vol,nQuads_BL,0) != 1 )  exit(EXIT_FAILURE);
                
                if ( MMG3D_Set_solSize(mmgMesh_hyb,mmgSol_hyb,MMG5_Vertex,mmgMesh_hyb->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);

                int i     = 0;
                int tt    = 1;
                int qt    = 1;
                int qt0   = 0;
                int fnt   = 0;
                int qfid  = 0;
                int qfidL = 0;
                int qfidR = 0;
                int tfid  = 0;
                
                std::map<std::set<int>, int> tfacesmap;
                std::set<std::set<int> >tfaces;
                std::map<int,int> tlh;
                std::map<int,int> trh;
                
                std::map<std::set<int>, int> qfacesmap;
                std::set<std::set<int> >qfaces;
                std::map<int,int> qlh;
                std::map<int,int> qrh;
                std::cout << "Set the prisms in the mmgMesh..."<<std::endl;
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
                        if(tria_ref_map.find(tria0)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria0];
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
                        if(tria_ref_map.find(tria1)!=tria_ref_map.end())
                        {
                            refer = tria_ref_map[tria1];
                            mmgMesh_hyb->tria[tt].v[0] = prism[3]+1;
                            mmgMesh_hyb->tria[tt].v[1] = prism[4]+1;
                            mmgMesh_hyb->tria[tt].v[2] = prism[5]+1;
                            mmgMesh_hyb->tria[tt].ref  = refer;
                            tt++;
                        }
                        
                        
                        
                        quad0.insert(prism[0]);
                        quad0.insert(prism[2]);
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
                        if(quad_ref_map.find(quad0)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad0];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                            qt0++;
                        }
                        quad1.insert(prism[1]);
                        quad1.insert(prism[5]);
                        quad1.insert(prism[4]);
                        quad1.insert(prism[2]);
                        
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
                        if(quad_ref_map.find(quad1)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad1];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[1]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[4]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[2]+1;
                            mmgMesh_hyb->quadra[qt].ref  = refer;
                            qt++;
                            qt0++;
                        }
                        quad2.insert(prism[0]);
                        quad2.insert(prism[3]);
                        quad2.insert(prism[5]);
                        quad2.insert(prism[1]);
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
                        if(quad_ref_map.find(quad2)!=quad_ref_map.end())
                        {
                            refer = quad_ref_map[quad2];
                            mmgMesh_hyb->quadra[qt].v[0] = prism[0]+1;
                            mmgMesh_hyb->quadra[qt].v[1] = prism[3]+1;
                            mmgMesh_hyb->quadra[qt].v[2] = prism[5]+1;
                            mmgMesh_hyb->quadra[qt].v[3] = prism[1]+1;
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
                int sw_v = 0;
                std::vector<int*> added_tets;
                std::set<int> used_vert;
                std::vector<double> jtet;
                std::vector<double> vtet;
                Array<double>* Jv = new Array<double>(mmgMesh_TET->np,1);
                Array<double>* Jvc = new Array<double>(mmgMesh_TET->np,1);
                for(int i=0;i<mmgMesh_TET->np;i++)
                {
                    Jv->setVal(i,0,0.0);
                    Jvc->setVal(i,0,0);
                }
                
                std::cout << "Set the tetrahedra in the mmgMesh..."<<std::endl;

                for(int i=1;i<=nel_tets;i++)
                {
                    if(mmgMesh_TET->tetra[i].ref == 20)
                    {
                        tet++;
                        int* tetra = new int[4];
                        int* tetra2 = new int[4];
                        std::set<int> face00;
                        std::set<int> face11;
                        std::set<int> face22;
                        std::set<int> face33;
                        sw_v = 0;
                        double* Points = new double[4*3];
                        
                        for(int s=0;s<4;s++)
                        {
                            Points[s*3+0] =mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[0];
                            Points[s*3+1] =mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[1];
                            Points[s*3+2] =mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].c[2];
                            
                            if(mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[s]].ref == 0)
                            {
                                m11 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][0];
                                m12 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][1];
                                m13 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][2];
                                m22 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][3];
                                m23 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][4];
                                m33 = newvert2metric[mmgMesh_TET->tetra[i].v[s]][5];
                                if(used_vert.find(mmgMesh_TET->tetra[i].v[s])==used_vert.end())
                                {
                                    used_vert.insert(mmgMesh_TET->tetra[i].v[s]);
                                    sw_v = 1;
                                }
                                int locNew = globNew2locNew[mmgMesh_TET->tetra[i].v[s]];
                                tetra[s]   = locNew;
                                tetra2[s]  = mmgMesh_TET->tetra[i].v[s];
                                mmgMesh_hyb->tetra[tet].v[s] = locNew;

                                if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,locNew) != 1 ) exit(EXIT_FAILURE);
                                
                                here++;
                            }
                            else
                            {
                                int vg = lv2gv_tet_mesh[mmgMesh_TET->tetra[i].v[s]-1];
                                mmgMesh_hyb->tetra[tet].v[s] = vg+1;
                                tetra[s]   = vg;
                                tetra2[s]  = mmgMesh_TET->tetra[i].v[s];
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
                        
                        jtet.push_back(ComputeDeterminantJ_tet(Points));
                        vtet.push_back(ComputeTetVolume(Points));
                        Jv->setVal(tetra2[0]-1,0,Jv->getVal(tetra2[0]-1,0)+ComputeDeterminantJ_tet(Points));
                        Jv->setVal(tetra2[1]-1,0,Jv->getVal(tetra2[1]-1,0)+ComputeDeterminantJ_tet(Points));
                        Jv->setVal(tetra2[2]-1,0,Jv->getVal(tetra2[2]-1,0)+ComputeDeterminantJ_tet(Points));
                        Jv->setVal(tetra2[3]-1,0,Jv->getVal(tetra2[3]-1,0)+ComputeDeterminantJ_tet(Points));
                    
                        Jvc->setVal(tetra2[0]-1,0,Jvc->getVal(tetra2[0]-1,0)+1);
                        Jvc->setVal(tetra2[1]-1,0,Jvc->getVal(tetra2[1]-1,0)+1);
                        Jvc->setVal(tetra2[2]-1,0,Jvc->getVal(tetra2[2]-1,0)+1);
                        Jvc->setVal(tetra2[3]-1,0,Jvc->getVal(tetra2[3]-1,0)+1);
                        
                        if(sw_v==1)
                        {
                            added_tets.push_back(tetra2);
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
                
                std::set<int> un_add_vert;
                std::map<int,int> lv2gv_add;
                std::map<int,int> gv2lv_add;
                std::vector<int> l2g_add;
                Array<int>* added_elements = new Array<int>(added_tets.size(),4);
                int uadv = 1;
                for(int el=0;el<added_tets.size();el++)
                {
                    for(int ve=0;ve<4;ve++)
                    {
                        if(un_add_vert.find(added_tets[el][ve])==un_add_vert.end())
                        {
                            un_add_vert.insert(added_tets[el][ve]);
                            lv2gv_add[uadv]=added_tets[el][ve];
                            l2g_add.push_back(added_tets[el][ve]);
                            gv2lv_add[added_tets[el][ve]]=uadv;
                            added_elements->setVal(el,ve,uadv);
                            uadv++;
                        }
                        else
                        {
                            added_elements->setVal(el,ve,gv2lv_add[added_tets[el][ve]]);
                        }
                    }
                }
                
                if(debug == 1)
                {
                    std::ofstream myfile;
                    myfile.open("added_tets.dat");
                    myfile << "TITLE=\"new_volume.tec\"" << std::endl;
                    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
                    myfile <<"ZONE N = " << un_add_vert.size() << ", E = " << added_tets.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
                    std::map<int,int>::iterator iter_map;

                    for(int i=0;i<l2g_add.size();i++)
                    {
                        myfile << mmgMesh_TET->point[l2g_add[i]].c[0] << " " <<   mmgMesh_TET->point[l2g_add[i]].c[1] << " " << mmgMesh_TET->point[l2g_add[i]].c[2] <<  std::endl;
                    }
                    for(int i=0;i<added_tets.size();i++)
                    {
                        myfile << added_elements->getVal(i,0) << " " << added_elements->getVal(i,1) << " " << added_elements->getVal(i,2) << " " << added_elements->getVal(i,3) << std::endl;
                    }
                    myfile.close();
                    
                    std::ofstream myfile2;
                    myfile2.open("Jac_elements.dat");
                    myfile2 << "TITLE=\"new_volume.tec\"" << std::endl;
                    myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\", \"J\"" << std::endl;
                    myfile2 <<"ZONE N = " << mmgMesh_TET->np << ", E = " << tellert << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

                    for(int i=1;i<=mmgMesh_TET->np;i++)
                    {
                        myfile2 << mmgMesh_TET->point[i].c[0] << " " << mmgMesh_TET->point[i].c[1] << " " << mmgMesh_TET->point[i].c[2] << " " << Jv->getVal(i-1,0)/Jvc->getVal(i-1,0) << std::endl;
                    }
                    for(int i=1;i<=mmgMesh_TET->ne;i++)
                    {
                        if(mmgMesh_TET->tetra[i].ref==20)
                        {
                            myfile2 << mmgMesh_TET->tetra[i].v[0] << " " << mmgMesh_TET->tetra[i].v[1] << " " << mmgMesh_TET->tetra[i].v[2] << " " << mmgMesh_TET->tetra[i].v[3] << std::endl;
                        }
                    }
                    myfile2.close();
                }
                
             
                //MMG3D_Set_handGivenMesh(mmgMesh_hyb);
                if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, 2.0) != 1 )    exit(EXIT_FAILURE);

                //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
                MMG3D_Set_dparameter( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_DPARAM_hgradreq , -1 );
                std::cout<<"Start the adaptation of the tetrahedra..."<<std::endl;
                int ier = MMG3D_mmg3dlib(mmgMesh_hyb,mmgSol_hyb);
                std::cout<<"Finished the adaptation of the tetrahedra..."<<std::endl;

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
                
                
                std::cout<<"Started writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;
                OutputMesh_MMG(mmgMesh_TETCOPY,0,mmgMesh_TETCOPY->ne,"OuterVolume.dat");
                std::cout<<"Finished writing the adapted tetrahedra mesh in ---> OuterVolume.dat"<<std::endl;

                
                std::cout<<"Started writing the adapted hybrid mesh in US3D format..."<<std::endl;
                WriteUS3DGridFromMMG(mmgMesh_hyb, us3d, bnd_face_map, unique_shell_tris);
                std::cout<<"Finished writing the adapted hybrid mesh in US3D format..."<<std::endl;
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
                int nbTriangles = tria_ref_map.size();
                
                
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
                    if(tria_ref_map.find(tria0)!=tria_ref_map.end())
                    {
                        ref0 = tria_ref_map[tria0];
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
                    if(tria_ref_map.find(tria1)!=tria_ref_map.end())
                    {
                        ref1 = tria_ref_map[tria1];
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
                    if(tria_ref_map.find(tria2)!=tria_ref_map.end())
                    {
                        ref2 = tria_ref_map[tria2];
                        mmgMesh->tria[t].v[0] = mmgMesh->tetra[offset_NE+i].v[2];
                        mmgMesh->tria[t].v[1] = mmgMesh->tetra[offset_NE+i].v[3];
                        mmgMesh->tria[t].v[2] = mmgMesh->tetra[offset_NE+i].v[0];
                        mmgMesh->tria[t].ref  = ref2;
                        t++;
                    }

                    tria3.insert(mmgMesh->tetra[offset_NE+i].v[3]-1);
                    tria3.insert(mmgMesh->tetra[offset_NE+i].v[0]-1);
                    tria3.insert(mmgMesh->tetra[offset_NE+i].v[1]-1);
                    if(tria_ref_map.find(tria3)!=tria_ref_map.end())
                    {
                        ref3 = tria_ref_map[tria3];
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
                
               // WriteUS3DGridFromMMG(mmgMesh, us3d);
            }
        }
         */
        MPI_Finalize();
        
    }
     
    return 0;
}
