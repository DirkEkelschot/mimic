
#include <vector>
#include <math.h>
#include <unordered_set>
#include <set>
#include <map>
#include <string>
#include <utility>
#include <math.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include <sstream>
#include <iomanip>

#include "adapt_io.h"
#include "adapt_compute.h"
#include "adapt_operations.h"
#include "adapt_math.h"
#include "adapt_output.h"
#include "adapt_partition.h"
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"
#include "parmmg/libparmmg.h"
#include "hex2tet.h"



int main(int argc, char** argv) {
        

    MPI_Init(NULL, NULL);

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::clock_t startr;
    startr = std::clock();
    //  GetXadjandAdjcyArrays(iee,ien,comm);
    //  Example3DPartitioningWithParVarParMetis();
    //  ExampleUS3DPartitioningWithParVarParMetis();
    //Example3DPartitioningWithParVarParMetis();
//============================================================

    //const char* fn_conn="grids/piston/conn.h5";
    const char* fn_conn="../adapt_step1/conn_66k.h5";
    const char* fn_grid="../adapt_step1/grid_66k.h5";

    // Read the hex mesh
    Array<double>* xcn    = ReadDataSetFromFile<double>(fn_grid,"xcn");
    Array<int>* ien       = ReadDataSetFromFile<int>(fn_conn,"ien");

    std::ifstream fin;
    fin.open("../adapt_step1/met.dat");

    // Read the file row by row
    std::vector<double> row(12);
    std::vector<std::vector<double> > arr;
    int t=0;
    while(fin >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5] >> row[6] >> row[7] >> row[8] >> row[9] >> row[10] >> row[11])
    {
       arr.push_back(row);
       t++;
    }
    int L = arr.size();
    std::ifstream finhex;
    finhex.open("../adapt_step1/pre_elem_rank_0.dat");

    // Read the file row by row
    std::vector<int> rowHex(8);
    std::vector<std::vector<int> > arrHex;

    while(finhex >> rowHex[0] >> rowHex[1] >> rowHex[2] >> rowHex[3] >> rowHex[4] >> rowHex[5] >> rowHex[6] >> rowHex[7])
    {
       arrHex.push_back(rowHex);
    }
    
    int nbHex      = arrHex.size();
    int nbVertices = arr.size();
    
/** 2) Build mesh in MMG5 format */
    /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
        file formatted or manually set your mesh using the MMG3D_Set* functions */
    int k =0;
    int ier;
    /** Manually set of the mesh */
    /** a) give the size of the mesh: 12 vertices, 12 tetra,0 prisms, 20
     * triangles, 0 quads, 0 edges */

    MMG5_pMesh mmgMesh = NULL;
    MMG5_pSol mmgSol   = NULL;
    
    int ref = 0;
    std::cout << "nbVertices =" << " " << arr.size() << std::endl;
    MMG3D_Init_mesh(MMG5_ARG_start,
    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
    MMG5_ARG_end);

    if ( MMG3D_Set_meshSize(mmgMesh,nbVertices,nbHex*6,0,0,0,0) != 1 )  exit(EXIT_FAILURE);

//    mmgMesh->nenil = 1;
//    for ( k=mmgMesh->nenil; k<mmgMesh->nemax-1; k++)
//      mmgMesh->tetra[k].v[3] = k+1;

    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);


//    double hmax = 0.001;
//    double hmin = 0.001;
    //double hmax = 0.01;
    for(int i=0;i<arr.size();i++)
    {
//        mmgMesh->point[i+1].c[0] = xcn->getVal(i,0);
//        mmgMesh->point[i+1].c[1] = xcn->getVal(i,1);
//        mmgMesh->point[i+1].c[2] = xcn->getVal(i,2);
        mmgMesh->point[i+1].c[0] = arr[i][0];
        mmgMesh->point[i+1].c[1] = arr[i][1];
        mmgMesh->point[i+1].c[2] = arr[i][2];
        
        mmgMesh->point[i+1].ref = 1;
        //std::cout << mmgMesh->point[i+1].c[0] << " " << mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] << std::endl;
//        double x = xcn->getVal(i,0);
//        double y = xcn->getVal(i,1);
//        double z = xcn->getVal(i,2);
//        double a = 60*fabs(0.75-sqrt(x*x+y*y));
//        double hx = min(0.002*pow(5,a),hmax);
//        double hy = min(0.05*pow(2,a),hmax);
//        double hz = hmax;
//        double rat;
//        if(x==0)
//        {
//            rat = 1.0;
//        }
//        else
//        {
//            rat = y/x;
//        }
//
//        double theta = atan(rat);
        //std::cout << hx << " " << hy << " " << hz << " " << theta << " " << rat << std::endl;
//        double m11 = (1.0/(hx*hx))*cos(theta)*cos(theta)+(1.0/(hy*hy))*sin(theta)*sin(theta);
//        double m12 = (1.0/(hx*hx)-1.0/(hy*hy))*cos(theta)*sin(theta);
//        double m13 = 0.0;
//        double m22 = (1.0/(hx*hx))*sin(theta)*sin(theta)+(1.0/(hy*hy))*cos(theta)*cos(theta);
//        double m23 = 0.0;
//        double m33 = 1.0/(hmax*hmax);

//        double m11 = 1.0/(hmax*hmax);
//        double m12 = 0.0;
//        double m13 = 0.0;
//        double m22 = 1.0/(hmax*hmax);
//        double m23 = 0.0;
//        double m33 = 1.0/(hmin*hmin);
        
        
        
        double m11 = arr[i][3];
        double m12 = arr[i][4];
        double m13 = arr[i][5];
        double m22 = arr[i][7];
        double m23 = arr[i][8];
        double m33 = arr[i][11];
        if((i+1)==8405 || i == 8405)
        {
            std::cout << m11 << " " << m12 << " " << m13 << " " <<m22<< " "<< m23 << " " << m33 << std::endl;
        }
        //
        
        if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,i+1) != 1 ) exit(EXIT_FAILURE);
    }




    Array<int>* ien_copy = new Array<int>(ien->getNrow(), ien->getNcol()-1);
    int* hexTab = new int[9*(nbHex+1)];
    for(int i=0;i<nbHex;i++)
    {
        int hexTabPosition = 9*(i+1);
        for(int j=0;j<8;j++)
        {
            //int val = ien->getVal(i,j+1);
            int val = arrHex[i][j];
            ien_copy->setVal(i,j,val);
            hexTab[hexTabPosition+j] = val;
        }
        hexTab[hexTabPosition+8] = ref;

    }
    
//        for(int i=0;i<(nbHex+1);i++)
//        {
//            for(int j=0;j<9;j++)
//            {
//                std::cout << hexTab[i*9+j] << " ";
//            }
//            std::cout << std::endl;
//        }

    // Make sure the orientation of the hex is correct in order to split up it up in tetrahedra.
    int num = H2T_chkorient(mmgMesh,hexTab,nbHex);

    /*
    * Hexahedron
    *
    *
    *   4----------5          .----------.           .----------.
    *   |\         |\         |\         |\          |\      _/ |\
    *   | \        | \        | \     3  | \         | \   _/   | \
    *   |  \       |  \       |  \  5    |  \        |  \ /     |  \
    *   |   7------+---6      |   .------+---.       |  /.------+---.
    *   |   |      |   |      | 1 |      | 4 |       | / |\_    |   |
    *   0---+------1   |      .---+------.   |       ./--+--\_--.   |
    *    \  |       \  |       \  |     2 \  |        \ \|__  \_ \  |
    *     \ |        \ |        \ |  0     \ |         \ |  \__ \_\ |
    *      \|         \|         \|         \|          \|     \__ \|
    *       3----------2          .----------.           .----------.
    *
    */

    int* adjahex = NULL;
    adjahex = (int*)calloc(6*nbHex+7,sizeof(int));
    assert(adjahex);
//
    if(!H2T_hashHexa(hexTab,adjahex,nbHex)) return H2T_STRONGFAILURE;
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
     for (k=6*nbHex; k<hed2.nhmax; k++)
       hed2.item[k].nxt = k+1;
//
    int ret = H2T_cuthex(mmgMesh, &hed2, hexTab, adjahex, nbHex);


//    std::vector<int> T1(4),T2(4),T3(4),T4(4),T5(4),T6(4);
//    std::vector<std::vector<int> > tets;
//    std::vector<int> H;
//    for(int i=0;i<nbHex;i++)
//    {
//        for(int j=0;j<8;j++)
//        {
//            H.push_back(hexTab[i*9+j]);
//        }
//        //std::cout << nbHex << " " << i << " " << H.size() << std::endl;
//        //======================================================================
//        T1[0] = H[0]; T1[1] = H[1]; T1[2] = H[3]; T1[3] = H[7];
//        T2[0] = H[7]; T2[1] = H[2]; T2[2] = H[6]; T2[3] = H[1];
//        T3[0] = H[1]; T3[1] = H[4]; T3[2] = H[5]; T3[3] = H[7];
//        T4[0] = H[7]; T4[1] = H[4]; T4[2] = H[0]; T4[3] = H[1];
//        T5[0] = H[1]; T5[1] = H[6]; T5[2] = H[7]; T5[3] = H[5];
//        T6[0] = H[1]; T6[1] = H[3]; T6[2] = H[7]; T6[3] = H[2];
//
//        tets.push_back(T1);
//        tets.push_back(T2);
//        tets.push_back(T3);
//        tets.push_back(T4);
//        tets.push_back(T5);
//        tets.push_back(T6);
//
//        H.clear();
//    }
////
//    std::cout << tets.size() << " " << mmgMesh->ne << std::endl;
//    for(int i=0;i<tets.size();i++)
//    {
//        mmgMesh->tetra[i+1].v[0] = tets[i][0];
//        mmgMesh->tetra[i+1].v[1] = tets[i][1];
//        mmgMesh->tetra[i+1].v[2] = tets[i][2];
//        mmgMesh->tetra[i+1].v[3] = tets[i][3];
//        mmgMesh->tetra[i+1].ref = 1;
//    }
    MMG3D_Set_handGivenMesh(mmgMesh);
    
    //if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hausd, 0.00001) != 1 )
    //    exit(EXIT_FAILURE);
    if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 1.3) != 1 )
        exit(EXIT_FAILURE);
    
//
    ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

    //if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

    std::ofstream myfile2;
    myfile2.open("mmgMesh_tecplot2.dat");
    myfile2 << "TITLE=\"volume_part_"  + std::to_string(0) +  ".tec\"" << std::endl;
    //myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\", \"m11\", \"m12\", \"m13\", \"m22\", \"m23\", \"m33\"" << std::endl;
    myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    //std::cout << " verts check " << LVerts.size() << " " << hx.size() << std::endl;
    myfile2 <<"ZONE N = " << mmgMesh->np << ", E = " << mmgMesh->ne << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
    for(int i=0;i<mmgMesh->np;i++)
    {
        myfile2 << mmgMesh->point[i+1].c[0] << " " <<mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] <<  std::endl;

        //myfile2 << mmgMesh->point[i+1].c[0] << " " <<mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] << " " << metric[i][0]<< " " << metric[i][1]<< " " << metric[i][2]<< " " << metric[i][3]<< " " << metric[i][4] << " " << metric[i][5] <<  std::endl;
    }

    for(int i=1;i<=mmgMesh->ne;i++)
    {
        myfile2 << mmgMesh->tetra[i].v[0] << " " << mmgMesh->tetra[i].v[1] << " " << mmgMesh->tetra[i].v[2] << " " << mmgMesh->tetra[i].v[3] << std::endl;
    }

    myfile2.close();

    MPI_Finalize();

    return 0;
     
}
