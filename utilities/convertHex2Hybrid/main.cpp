#include "../../src/adapt_io.h"
#include "../../src/adapt_boundary.h"
#include "../../src/hex2tet.h"
#include "../../src/adapt_bltopology.h"


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
    int i,j,k;
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;
    
    if(world_size > 1)
    {
        if(world_rank == 0)
        {
            std::cout << "This application only runs in serial so try to rerun using: mpirun -np 1 ./ConvertHex2HybridMesh " << std::endl;
        }
        
        MPI_Finalize();
        
        return -1;
    }
    
    const char* fn_grid   = "inputs/grid.h5";
    const char* fn_conn   = "inputs/conn.h5";
    int ReadFromStats = 0;
    US3D* us3d    = ReadUS3DGrid(fn_conn,fn_grid,ReadFromStats,comm,info);

    //US3D* us3d  = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve       = us3d->xcn->getNglob();
    
    int Nel_part  = us3d->ien->getNrow();
    int Nel_glob  = us3d->ien->getNglob();
    
    std::cout << "Starting to extract the boundary map... " << std::endl;
    BoundaryMap* bmap = new BoundaryMap(us3d->ifn, us3d->if_ref);
    std::cout << "Finished extracting the boundary map... " << std::endl;

//    if(debug == 1)
//    {
//        std::ofstream myfile;
//        myfile.open("cylinderOutput.dat");
//        myfile << "TITLE=\"new_volume.tec\"" << std::endl;
//        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//        myfile <<"ZONE N = " << us3d->xcn->getNrow() << ", E = " << us3d->ien->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
//
//        for(int i=0;i<us3d->xcn->getNrow();i++)
//        {
//            myfile << us3d->xcn->getVal(i,0) << " " << us3d->xcn->getVal(i,1)  << " " << us3d->xcn->getVal(i,2) <<  std::endl;
//        }
//        for(int i=0;i<us3d->ien->getNrow();i++)
//        {
//            myfile << us3d->ien->getVal(i,0)+1 << " "
//            << us3d->ien->getVal(i,1)+1 << " "
//            << us3d->ien->getVal(i,2)+1 << " "
//            << us3d->ien->getVal(i,3)+1 << " "
//            << us3d->ien->getVal(i,4)+1 << " "
//            << us3d->ien->getVal(i,5)+1 << " "
//            << us3d->ien->getVal(i,6)+1 << " "
//            << us3d->ien->getVal(i,7)+1 << std::endl;
//        }
//        myfile.close();
//    }
    
    
    int nPrmsNormal = 0;
    
    int wall_id     = 3;
    int wall_id2    = 14;
    
    if (argc < 2)
    {
        std::cout << "ERROR: Please provide the desired number of prisms in normal direction with respect to the wall." << std::endl;
        return 1;
    }
    else
    {
        nPrmsNormal = atoi(argv[1]);
    }

    if(nPrmsNormal>0)
    {
        int counter = 0;
        
        std::map<int,std::vector<int> > bnd_face_map = bmap->getBfaceMap();
        std::map<std::set<int>,int> tria_ref_map     = bmap->getTriaRefMap();
        std::map<std::set<int>,int> quad_ref_map     = bmap->getQuadRefMap();
        
        std::map<int,int> vert_ref_map = bmap->getNodeRefMap();
        BLShellInfo* BLshell = new BLShellInfo;
        
        FindOuterShellBoundaryLayerMesh_V2(BLshell,wall_id, nPrmsNormal,us3d->xcn,us3d->ien,us3d->iee,us3d->ief,us3d->ife,us3d->ifn, bnd_face_map,vert_ref_map,comm);
        
        FindOuterShellBoundaryLayerMesh_V2(BLshell,wall_id2, nPrmsNormal,us3d->xcn,us3d->ien,us3d->iee,us3d->ief,us3d->ife,us3d->ifn, bnd_face_map,vert_ref_map,comm);
        
        if(debug == 1)
        {
            std::ofstream myfile1;
            myfile1.open("BLOutput.dat");
            myfile1 << "TITLE=\"new_volume.tec\"" << std::endl;
            myfile1 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
            myfile1 <<"ZONE N = " << us3d->xcn->getNrow() << ", E = " << BLshell->elements_set.size() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;

            for(int i=0;i<us3d->xcn->getNrow();i++)
            {
                myfile1 << us3d->xcn->getVal(i,0) << " " << us3d->xcn->getVal(i,1)  << " " << us3d->xcn->getVal(i,2) <<  std::endl;
            }
            
            std::set<int>::iterator its;
            for(its=BLshell->elements_set.begin();its!=BLshell->elements_set.end();its++)
            {
                int gid = *its;
//                std::cout << gid << " " << us3d->ien->getVal(gid,0)+1 << " "
//                << us3d->ien->getVal(gid,1)+1 << " "
//                << us3d->ien->getVal(gid,2)+1 << " "
//                << us3d->ien->getVal(gid,3)+1 << " "
//                << us3d->ien->getVal(gid,4)+1 << " "
//                << us3d->ien->getVal(gid,5)+1 << " "
//                << us3d->ien->getVal(gid,6)+1 << " "
//                << us3d->ien->getVal(gid,7)+1 << " " << std::endl;
                
                myfile1 << us3d->ien->getVal(gid,0)+1 << " "
                << us3d->ien->getVal(gid,1)+1 << " "
                << us3d->ien->getVal(gid,2)+1 << " "
                << us3d->ien->getVal(gid,3)+1 << " "
                << us3d->ien->getVal(gid,4)+1 << " "
                << us3d->ien->getVal(gid,5)+1 << " "
                << us3d->ien->getVal(gid,6)+1 << " "
                << us3d->ien->getVal(gid,7)+1 << std::endl;
            }
            
            myfile1.close();
        }
        
        
        
        
        
        
        int nbHex       =  us3d->ien->getNrow();
        int nbPrisms    =  (bnd_face_map[wall_id2].size()+bnd_face_map[wall_id].size())*(nPrmsNormal)*2;
        int nbHexsNew   = (nbHex-(bnd_face_map[wall_id2].size()+bnd_face_map[wall_id].size())*(nPrmsNormal));
        
        int ith = 0;
        std::set<int> u_tet_vert;
        std::map<int,int> lv2gv_tet_mesh;
        std::map<int,int> gv2lv_tet_mesh;
        std::map<int,double*> metric_hex2tet;
        Array<int>* ien_hex2tet = new Array<int>(nbHexsNew,8);
        int sv = 0;
        std::vector<int> locTet_verts;
        
        std::cout << nbHexsNew << " " << nbHex << " " << nbPrisms << " " << bnd_face_map[wall_id].size() << " " << BLshell->elements_set.size() << std::endl;
        
        for(int i=0;i<us3d->ien->getNrow();i++)
        {
            if(BLshell->elements_set.find(i)==BLshell->elements_set.end())
            {
                for(int j=0;j<8;j++)
                {
                    int val = us3d->ien->getVal(i,j);
                    
                    if(u_tet_vert.find(val)==u_tet_vert.end())
                    {
                        u_tet_vert.insert(val);
                        locTet_verts.push_back(val);
                        gv2lv_tet_mesh[val]=sv;
                        lv2gv_tet_mesh[sv]=val;
                        ien_hex2tet->setVal(ith,j,sv);
                        
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
        int VertRef= 0;
        for(int i=0;i<nbVerts_TET;i++)
        {
            int gv=lv2gv_tet_mesh[i];
            
            if(BLshell->shellVrts.find(gv)!=BLshell->shellVrts.end())
            {
                VertRef = 333;
            }
            else
            {
                VertRef = 1;
                
            }
            
            if(BLshell->ShellRef->getVal(gv,0)==999)
            {
                cshell++;
            }
            if(BLshell->ShellRef->getVal(gv,0)==-3)
            {
                nshell++;
            }
            
            mmgMesh_TET->point[i+1].c[0] = us3d->xcn->getVal(gv,0);
            mmgMesh_TET->point[i+1].c[1] = us3d->xcn->getVal(gv,1);
            mmgMesh_TET->point[i+1].c[2] = us3d->xcn->getVal(gv,2);
            mmgMesh_TET->point[i+1].ref  = VertRef;
            

            
            double m11 = 1.0;//mv_g->getVal(gv,0);
            double m12 = 0.0;//mv_g->getVal(gv,1);
            double m13 = 0.0;//mv_g->getVal(gv,2);
            double m22 = 1.0;//mv_g->getVal(gv,3);
            double m23 = 0.0;//mv_g->getVal(gv,4);
            double m33 = 1.0;//mv_g->getVal(gv,5);
            
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
        
        std::cout << "Cut each hexahedral up into 6 tetrahedra..."<<std::endl;
        //we need to do this in order to determine the orientation of the triangles at the shell interface and trace back how we need to tesselate the wall boundary into triangles so that they match eachothers orientation.
        
        int ret = H2T_cuthex(mmgMesh_TET, &hed22, hexTabNew, adjahexNew, nbHexsNew);
        
        if(debug == 1)
        {
            std::ofstream myfile3;
            myfile3.open("TetVolOutput.dat");
            myfile3 << "TITLE=\"new_volume.tec\"" << std::endl;
            myfile3 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
            myfile3 <<"ZONE N = " << mmgMesh_TET->np << ", E = " << (int)mmgMesh_TET->ne/2 << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

            for(int i=0;i<mmgMesh_TET->np;i++)
            {
                myfile3 << mmgMesh_TET->point[i+1].c[0]  << " " << mmgMesh_TET->point[i+1].c[1]   << " " << mmgMesh_TET->point[i+1].c[2] <<  std::endl;
                // std::cout << mmgMesh_TET->point[i+1].ref << std::endl;
            }
            
            int start = (int) mmgMesh_TET->ne/2.0;
            for(int i=start;i<mmgMesh_TET->ne;i++)
            {
                myfile3 << mmgMesh_TET->tetra[i+1].v[0]  << " " << mmgMesh_TET->tetra[i+1].v[1] << " " << mmgMesh_TET->tetra[i+1].v[2] << " " << mmgMesh_TET->tetra[i+1].v[3]  <<  std::endl;
            }
            
            myfile3.close();
        }
        
        
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
                if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==333 )
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
                if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==333 )
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
                if(   mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[2]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==333 )
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
                if( mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[3]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[0]].ref==333
                   && mmgMesh_TET->point[mmgMesh_TET->tetra[i].v[1]].ref==333 )
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
            //std::cout << "shell_faceid " << shell_faceid << std::endl;
            teller++;
        }
        
        std::cout << "Extracting the prismatic boundary layer mesh... " << unique_shell_tris.size() <<std::endl;
        Mesh_Topology_BL* mesh_topo_bl2 = new Mesh_Topology_BL;
        mesh_topo_bl2->Nprisms = 0;

        
        ExtractBoundaryLayerMeshFromShell(mesh_topo_bl2,u_tris, BLshell, wall_id, nPrmsNormal, us3d->xcn, us3d->ien, us3d->ief, us3d->ife, us3d->ifn, bnd_face_map, tria_ref_map, quad_ref_map, comm);
        
        ExtractBoundaryLayerMeshFromShell(mesh_topo_bl2,u_tris, BLshell, wall_id2, nPrmsNormal, us3d->xcn, us3d->ien, us3d->ief, us3d->ife, us3d->ifn, bnd_face_map, tria_ref_map, quad_ref_map, comm);
        
        std::cout << "mesh_topo_bl2->Nprisms " << mesh_topo_bl2->Nprisms << std::endl;
        if(debug == 1)
        {
            std::cout << "Outputting the prismatic boundary layer mesh in ---> BoundaryLayerMesh_0.dat" <<std::endl;
            OutputBoundaryLayerPrisms(us3d->xcn, mesh_topo_bl2, comm, "BoundaryLayerMesh_");
            std::ofstream myfile_shell;
            myfile_shell.open("shell.dat");
            myfile_shell << "TITLE=\"new_volume.tec\"" << std::endl;
            myfile_shell <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
            myfile_shell <<"ZONE N = " << U_shell_vert.size() << ", E = " << teller << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;

            for(int i=0;i<U_shell_vert.size();i++)
            {
                myfile_shell << us3d->xcn->getVal(loc2glob_shell_vert[i],0)
                << " " << us3d->xcn->getVal(loc2glob_shell_vert[i],1)
                << " " << us3d->xcn->getVal(loc2glob_shell_vert[i],2) << std::endl;
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
        std::map<int,std::vector<int> > newvert2vert;
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
                                newvert2vert[mmgMesh_TET->tetra[i].v[s]].push_back(mmgMesh_TET->tetra[i].v[t]);
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
        std::cout << "unique_new_verts.size = > " << unique_new_verts.size() << std::endl;

        double min_oris = *std::min_element(oris.begin(),oris.end());
        double max_oris = *std::max_element(oris.begin(),oris.end());

        int nTriangles_Vol = 0;
        std::map<int,std::vector< int* > >::iterator itrrb;
        for(itrrb=bound_tet.begin();itrrb!=bound_tet.end();itrrb++)
        {
            nTriangles_Vol  = nTriangles_Vol+itrrb->second.size();
        }
        
        
        MMG5_pMesh mmgMesh_hyb = NULL;
        MMG5_pSol mmgSol_hyb   = NULL;
        
        MMG3D_Init_mesh(MMG5_ARG_start,
        MMG5_ARG_ppMesh,&mmgMesh_hyb,MMG5_ARG_ppMet,&mmgSol_hyb,
        MMG5_ARG_end);
        
        int nVertices_New  = mmgMesh_TET->np+unique_prism_verts.size()-cshell;
        int nbTets_New     = tellert;
        
        std::cout << "check cshell " << nVertices_New << " " << cshell << " " << unique_prism_verts.size() << std::endl;
        
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
                mmgMesh_hyb->prism[i+1].v[1] = prism[1]+1;
                mmgMesh_hyb->prism[i+1].v[2] = prism[2]+1;
                mmgMesh_hyb->prism[i+1].v[3] = prism[3]+1;
                mmgMesh_hyb->prism[i+1].v[4] = prism[4]+1;
                mmgMesh_hyb->prism[i+1].v[5] = prism[5]+1;
                
                mmgMesh_hyb->prism[i+1].ref  = 0;
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
        int nbVertices = us3d->xcn->getNrow();
        for(int i=0;i<nbVertices;i++)
        {
            mmgMesh_hyb->point[i+1].c[0] = us3d->xcn->getVal(i,0);
            mmgMesh_hyb->point[i+1].c[1] = us3d->xcn->getVal(i,1);
            mmgMesh_hyb->point[i+1].c[2] = us3d->xcn->getVal(i,2);
        }
        std::cout << "mmgMesh_hyb->np - " << mmgMesh_hyb->np << std::endl;
        std::map<int,std::vector<int> >::iterator v2m;
        std::map<int,int> locNew2globNew;
        std::map<int,int> globNew2locNew;
        int p = 0;
        
        for(v2m=newvert2vert.begin();v2m!=newvert2vert.end();v2m++)
        {
            mmgMesh_hyb->point[nbVertices+p+1].c[0] = mmgMesh_TET->point[v2m->first].c[0];
            mmgMesh_hyb->point[nbVertices+p+1].c[1] = mmgMesh_TET->point[v2m->first].c[1];
            mmgMesh_hyb->point[nbVertices+p+1].c[2] = mmgMesh_TET->point[v2m->first].c[2];

            locNew2globNew[nbVertices+p+1] = v2m->first;
            globNew2locNew[v2m->first]     = nbVertices+p+1;
            //if ( MMG3D_Set_tensorSol(mmgSol_hyb, m11,m12,m13,m22,m23,m33,nbVertices+p+1) != 1 ) exit(EXIT_FAILURE);

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
        if ( MMG3D_Set_dparameter(mmgMesh_hyb,mmgSol_hyb,MMG3D_DPARAM_hgrad, 3.0) != 1 )    exit(EXIT_FAILURE);

        //MMG3D_Set_iparameter ( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_IPARAM_nosizreq , 1 );
        MMG3D_Set_dparameter( mmgMesh_hyb,  mmgSol_hyb,  MMG3D_DPARAM_hgradreq , -1 );
//        std::cout<<"Start the adaptation of the tetrahedra..."<<std::endl;
//        int ier = MMG3D_mmg3dlib(mmgMesh_hyb,mmgSol_hyb);
//        std::cout<<"Finished the adaptation of the tetrahedra..."<<std::endl;

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

        MMG3D_Free_all(MMG5_ARG_start,
                       MMG5_ARG_ppMesh,&mmgMesh_TETCOPY,MMG5_ARG_ppSols,&mmgSol_TETCOPY,
                       MMG5_ARG_end);
        
        MMG3D_Free_all(MMG5_ARG_start,
                       MMG5_ARG_ppMesh,&mmgMesh_TET,MMG5_ARG_ppSols,&mmgSol_TET,
                       MMG5_ARG_end);
        
        std::cout<<"Started writing the adapted hybrid mesh in US3D format..."<<std::endl;
//        WriteUS3DGridFromMMG_it0_NEW(mmgMesh_hyb, mmgSol_hyb, us3d);
        WriteUS3DGridFromMMG_it0(mmgMesh_hyb, mmgSol_hyb, us3d);
        std::cout<<"Finished writing the adapted hybrid mesh in US3D format..."<<std::endl;
        //
    
        
        
  
    }
    else
    {
        std::cout << "Please specify the number of prisms in normal direction from the wall." << std::endl;
    }
     
    
    
    MPI_Finalize();
    
}
