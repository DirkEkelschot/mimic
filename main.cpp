#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
//#include "petsc.h"
#include <mpi.h>

//#include <hdf5.h>
#include <math.h>
#include "us3d_io.h"
#include "us3d_compute.h"
#include "us3d_math.h"
#include <fstream>
//#include "us3d_print.h"
//#include "us3d_datastruct.h"
#include <utility>
#include <string>
#include <vector>
#include <math.h>
#include <map>
#include <unordered_set>
#include <set>
#include "parmetis.h"

int mpi_size, mpi_rank;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


extern "C" {
  void dgeev_(char const * __restrict JOBVL, char const * __restrict JOBVR, int const * __restrict n, double * __restrict A, int const * lda, double * __restrict WR, double * __restrict WI, double * __restrict VL, int const * __restrict ldvl, double * __restrict VR, int const * __restrict ldvr, double * __restrict Work, int const * __restrict lwork, int       * __restrict info );
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  void dgemm_(char * transA,  char * transB, int * m, int * n, int * k, double * alpha, double * A, int * lda, double * B, int * ldb, double * beta, double * C, int * ldc );
  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

}


using namespace std;


struct HalfEdge
{
   HalfEdge * oppositeHalfEdge;
   HalfEdge * prevHalfEdge;
   HalfEdge * nextHalfEdge;
   int vertex;
   int opposite_vertex;
   int face_g;
   int face_l;
};



void EigenDecomp(int n, double * A,  double * WR, double * WI, double * V, double * iV )
{
  char JOBVL = 'V';
  char JOBVR = 'N';
  int size = 10*n;
  double WORK [size];
  int info;
  int i,j;
  int Pivot[n];
  
  // Copy A into V 
  memcpy( V, A, n*n*sizeof(double) );
  
  // Factor A, right eigenvectors are in iV though column major
  dgeev_( &JOBVL, &JOBVR, &n, V, &n, WR, WI, iV, &n, NULL, &n, WORK, &size, &info );
  
  // Copy right eigenvectors into V (with transpose)
  for ( i = 0; i < n; i++)
    for ( j = 0; j < n; j++)
      V[i*n+j] = iV[j*n+i];
  
  // Compute inverse of V1 
  memcpy( iV, V, n*n*sizeof(double) );
  dgetrf_(&n, &n, iV, &n, Pivot, &info);
  dgetri_(&n, iV, &n, Pivot, WORK, &size, &info);
}




//typedef map<pair<int, int>, HalfEdge *> MapType;

map< pair<int,int>, HalfEdge* > GetHalfEdges(int* element_verts, int* M, int nloc, int offset)
{
    map< pair<int, int>, HalfEdge* > Edges;
    int cnt=0;
    int dupl=0;
    int u,v,un,vn,up,vp;
    for(int i=0;i<nloc;i++)
    {
        int id_b = M[i];
        int id_e = M[i+1];
        
        int type = id_e-id_b;
        //std::cout << "type - " << type <<std::endl;
        if (type == 3)
        {
            /*
            std::cout << "T :: ";
            for(int j=id_b;j<id_e;j++)
            {
                std::cout << element_verts[j] << " ";
            }
            */
             
            for(int j=id_b;j<id_e;j++)
            {
            
                if(j == id_b)
                {
                    u  = element_verts[j];
                    v  = element_verts[j+1];
                
                    up = element_verts[id_e-1];
                    vp = element_verts[j];
                
                    un = element_verts[j+1];
                    vn = element_verts[j+2];
                }
                else if(j == id_e-1)
                {
                    u  = element_verts[j];
                    v  = element_verts[id_b];
                    
                    up = element_verts[j-1];
                    vp = element_verts[j];
                    
                    un = element_verts[id_b];
                    vn = element_verts[id_b+1];
                }
                else
                {
                    u  = element_verts[j];
                    v  = element_verts[id_e-1];
                    
                    up = element_verts[id_b];
                    vp = element_verts[j];
                    
                    un = element_verts[id_e-1];
                    vn = element_verts[id_b];
                }
                
                pair<int,int> edge = make_pair(u,v);
                HalfEdge* HE = new HalfEdge();
                HalfEdge* HE_opp = new HalfEdge();
                HE->oppositeHalfEdge = HE_opp;
                HE->vertex = u;
                HE->opposite_vertex = v;
                HE->face_g = offset+i;
                HE->face_l = i;
                Edges[edge]=HE;
                
                pair<int,int> prev_edge = make_pair(up,vp);
                pair<int,int> next_edge = make_pair(un,vn);

                Edges[ edge ]->prevHalfEdge = Edges[prev_edge];
                Edges[ edge ]->nextHalfEdge = Edges[next_edge];
                
                pair<int,int> rev_edge;
                rev_edge = make_pair(v,u);
                
                if ( Edges.find( rev_edge ) != Edges.end() )
                {
                    Edges[edge]->oppositeHalfEdge = Edges[rev_edge];
                    Edges[rev_edge]->oppositeHalfEdge = Edges[edge];
                }
            }
        }
        
        if (type == 4)
        {
            /*
            std::cout << "Q :: ";
            for(int j=id_b;j<id_e;j++)
            {
                std::cout << element_verts[j] << " ";
            }
            
            std::cout << std::endl;
            */
            
            for(int j=id_b;j<id_e;j++)
            {
                if(j == id_b)
                {
                    u  = element_verts[j];
                    v  = element_verts[j+1];
                
                    up = element_verts[id_e-1];
                    vp = element_verts[j];
                
                    un = element_verts[j+1];
                    vn = element_verts[j+2];
                }
                else if(j == id_e-1)
                {
                    u  = element_verts[j];
                    v  = element_verts[id_b];
                    
                    up = element_verts[j-1];
                    vp = element_verts[j];
                    
                    un = element_verts[id_b];
                    vn = element_verts[id_b+1];
                }
                else if(j == id_b+1)
                {
                    u  = element_verts[j];
                    v  = element_verts[j+1];
                    
                    up = element_verts[id_b];
                    vp = element_verts[j];
                    
                    un = element_verts[j+1];
                    vn = element_verts[j+2];
                }
                else if(j == id_b+2)
                {
                    u  = element_verts[j];
                    v  = element_verts[j+1];
                    
                    up = element_verts[j-1];
                    vp = element_verts[j];
                    
                    un = element_verts[id_e-1];
                    vn = element_verts[id_b];
                }
                
                pair<int,int> edge;
                edge = make_pair(u,v);
                HalfEdge* HE = new HalfEdge();
                HalfEdge* HE_opp = new HalfEdge();

                HE->oppositeHalfEdge = HE_opp ;
                HE->vertex = u;
                HE->opposite_vertex = v;
                HE->face_g = offset+i;
                HE->face_l = i;
                
                Edges[edge]=HE;
                
                pair<int,int> prev_edge = make_pair(up,vp);
                pair<int,int> next_edge = make_pair(un,vn);

                Edges[ edge ]->prevHalfEdge = Edges[prev_edge];
                Edges[ edge ]->nextHalfEdge = Edges[next_edge];
                
                pair<int,int> rev_edge;
                rev_edge = make_pair(v,u);
                
                if ( Edges.find( rev_edge ) != Edges.end())
                {
                    Edges[edge]->oppositeHalfEdge = Edges[rev_edge];
                    Edges[rev_edge]->oppositeHalfEdge = Edges[edge];
                }
            }
        }
    }
    
    //std::cout << "data  " << cnt << " " << dupl << std::endl;
    
    return Edges;
}



void UnitTestEigenDecomp()
{
    double *M = new double[3*3];
    M[0] = 0.25;M[1]=0.0;M[2]=0.0;
    M[3] = 0.0;M[4]=0.25;M[5]=0.0;
    M[6] = 0.0;M[7]=0.0;M[8]=0.25;
    
    double * WR = new double[3];
    double * WI = new double[3];
    double * V = new double[3*3];
    double * iV = new double[3*3];
    EigenDecomp(3, M,  WR,  WI, V, iV );
    
    for(int i=0;i<3;i++)
    {
        std::cout << "eigenvalues = " << WR[i] << " + " << WI[i] << "i" << std::endl;
    }
    
    
    
    std::cout << "V ==> "  << std::endl;;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << V[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
    std::cout << "Vinv ==> " << std::endl;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << iV[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void UnitTestJacobian()
{
    double* Hex = new double[8*3];
    
    Hex[0*3+0] = 0;     Hex[0*3+1] = 0;     Hex[0*3+2] = 0;
    Hex[1*3+0] = 0.5;   Hex[1*3+1] = 0;     Hex[1*3+2] = 0;
    Hex[2*3+0] = 0.5;   Hex[2*3+1] = 0.5;   Hex[2*3+2] = 0;
    Hex[3*3+0] = 0;     Hex[3*3+1] = 0.5;   Hex[3*3+2] = 0;
    
    Hex[4*3+0] = 0;     Hex[4*3+1] = 0;     Hex[4*3+2] = 0.5;
    Hex[5*3+0] = 0.5;   Hex[5*3+1] = 0;     Hex[5*3+2] = 0.5;
    Hex[6*3+0] = 0.5;   Hex[6*3+1] = 0.5;   Hex[6*3+2] = 0.5;
    Hex[7*3+0] = 0;     Hex[7*3+1] = 0.5;   Hex[7*3+2] = 0.5;
    
    double* JP1 = ComputeJAtCenter(Hex,8);
    
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << JP1[i*3+j] << " ";
        }
        std::cout << " " << std::endl;
    }
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
                 -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
                 +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    std::cout << DetJ << std::endl;
}





void PartitionBigMesh(Array<double>* xcn, Array<int>* ien, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int eltype = 8;
    int nel    = ien->nglob;
    int npo    = 0;
    
    for(int i = 0;i<nel;i++)
    {
        npo += 8;
    }
    
    int nloc     = int(nel/world_size) + ( world_rank < nel%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(nel/world_size) + MIN(world_rank, nel%world_size);
    int* elmdist = new int[world_size+1];

    int npo_loc=0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += 8;
    }
    
    int* locs        = new int[world_size];
    int* npo_locs    = new int[world_size];
    int* npo_offset  = new int[world_size+1];
    npo_offset[0]=0;
    
    
    
    for(int i=0;i<world_size;i++)
    {
        if (i==world_rank)
        {
            locs[i]     = nloc;
            npo_locs[i] = npo_loc;
        }
        else
        {
            locs[i]     = 0;
            npo_locs[i] = 0;
        }
    }
    
    
     
     
    int* offsets = new int[world_size+1];
    for(int i=0;i<world_size+1;i++)
    {
        if (i==world_rank)
        {
            elmdist[i]    = offset;
        }
        else
        {
            elmdist[i]    = 0;
        }
    }
    
    
    
    
    int* red_locs       = new int[world_size];
    int* red_npo_locs   = new int[world_size];
    int* red_elmdist    = new int[world_size+1];
    std::cout << "world_size:: " << world_size << std::endl;
    for(int i=0;i<world_size;i++)
    {
        red_locs[i]    = 0;
        red_elmdist[i] = 0;
    }
    
    MPI_Allreduce(locs,     red_locs,     world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(npo_locs, red_npo_locs, world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(elmdist,  red_elmdist,  world_size+1, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    
    for(int i=0;i<world_size;i++)
    {
        npo_offset[i+1] = npo_offset[i]+red_npo_locs[i];
    }
    
    red_elmdist[world_size] = nel;
    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    
    eptr[0]  = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1] = eptr[i]+8;
        
        for(int j=0;j<8;j++)
        {
            eind[i*8+j] = ien->getVal(offset+i,j+1)-1;
        }
    }
    
    //std::cout << "extra check " << eptr[nloc] << " " << npo_loc << std::endl;
    //std::cout << "ntot " <<  nloc*8 << " " << npo_loc << std::endl;
    
    /*
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
    idx_t *adjwgt    = NULL;
    idx_t *vsize     = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt;
    
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0;i<np*ncon[0];i++)
    {
        tpwgts[i]    = 1.0/np;
    }
    
    idx_t nparts_[]  = {np};
    idx_t *nparts    = nparts_;
    
    idx_t part_[]    = {nloc};
    idx_t *part      = part_;
    */
    /*
    if(world_rank == 0)
    {
        for(int i=0;i<world_size+1;i++)
        {
            std::cout << "efwf " << world_rank << " " << red_elmdist[i] << std::endl;
        }
    }
    if(world_rank == 1)
    {
        for(int i=0;i<nloc;i++)
        {
            //std::cout << eptr[i] << std::endl;
            //std::cout << "element - " << i << " :: ";
            for(int j=0;j<8;j++)
            {
                std::cout << eind[i*8+j] << " ";
            }
            
            std::cout << std::endl;
            
        }
    }
    */
    //ParMETIS_V3_PartMeshKway(red_elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
    //ParMETIS_V3_Mesh2Dual(red_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    //idx_t *nparts2 = nparts_;
    
    //ParMETIS_V3_AdaptiveRepart(red_elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
    
    //std::cout << "part = " << part[0] << std::endl;
    /*
    if(world_rank == 1)
    {
        std::cout << std::endl;
        
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
   */ //===================================================================================================================================
    /*
    int rank = 0;
    if(world_rank == rank)
    {

        int cnt = 0;
        
        int nshow = 20;
        
        for(int i=0;i<nloc;i++)
        {
            if(cnt < nshow)
            {
                //std::cout << "element " << i+offset << " :: -> ";
            }
            for(int j=xadj[i];j<xadj[i+1];j++)
            {
                if(cnt < nshow)
                {
                    std::cout << adjncy[j] << ", ";
                }
            }
            if(cnt < nshow)
            {
                std::cout << " " << std::endl;
            }
            cnt++;
        }
    }
    */
}

void Example3DPartitioning()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    int nel = 8;
    int * eltype = new int[nel];
    
    eltype[0] = 8;
    eltype[1] = 8;
    eltype[2] = 8;
    eltype[3] = 8;
    eltype[4] = 8;
    eltype[5] = 8;
    eltype[6] = 8;
    eltype[7] = 8;
    
    int npo   = 0;
    
    for(int i = 0;i < nel; i++)
    {
        npo += eltype[i];
    }
    
    int * test = new int[npo];
    
    test[0] = 0;test[1] = 1;test[2] = 6;test[3] = 5;        test[4]  = 14+0;test[5] = 14+1;test[6] = 14+6;test[7] = 14+5;
    test[8] = 1;test[9] = 2;test[10] = 7;test[11] = 6;      test[12] = 14+1;test[13] = 14+2;test[14] = 14+7;test[15] = 14+6;
    test[16] = 2;test[17] = 3;test[18] = 8;test[19] = 7;    test[20] = 14+2;test[21] = 14+3;test[22] = 14+8;test[23] = 14+7;
    test[24] = 3;test[25] = 4;test[26] = 9;test[27] = 8;    test[28] = 14+3;test[29] = 14+4;test[30] = 14+9;test[31] = 14+8;
    test[32] = 5;test[33] = 6;test[34] = 11;test[35] = 10;  test[36] = 14+5;test[37] = 14+6;test[38] = 14+11;test[39] = 14+10;
    test[40] = 6;test[41] = 7;test[42] = 12;test[43] = 11;  test[44] = 14+6;test[45] = 14+7;test[46] = 14+12;test[47] = 14+11;
    test[48] = 7;test[49] = 8;test[50] = 13;test[51] = 12;  test[52] = 14+7;test[53] = 14+8;test[54] = 14+13;test[55] = 14+12;
    test[56] = 8;test[57] = 9;test[58] = 14;test[59] = 13;  test[60] = 14+8;test[61] = 14+9;test[62] = 14+14;test[63] = 14+13;

    
    int nloc     = int(nel/world_size) + ( world_rank < nel%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(nel/world_size) + MIN(world_rank, nel%world_size);
    int* elmdist = new int[world_size];

    int npo_loc=0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += eltype[offset+i];
    }
    
    int* locs        = new int[world_size];
    int* npo_locs    = new int[world_size];
    int* npo_offset  = new int[world_size+1];
    npo_offset[0]=0;
    
    for(int i=0;i<world_size;i++)
    {
        if (i==world_rank)
        {
            locs[i]     = nloc;
            npo_locs[i] = npo_loc;
        }
        else
        {
            locs[i]     = 0;
            npo_locs[i] = 0;
        }
    }
    
    int* offsets = new int[world_size+1];
    for(int i=0;i<world_size+1;i++)
    {
        if (i==world_rank)
        {
            elmdist[i]    = offset;
        }
        else
        {
            elmdist[i]    = 0;
        }
    }
    
    int* red_locs       = new int[world_size];
    int* red_npo_locs   = new int[world_size];
    int* red_elmdist    = new int[world_size+1];
    
    for(int i=0;i<world_size;i++)
    {
        red_locs[i]    = 0;
        red_elmdist[i] = 0;
    }
    
    MPI_Allreduce(locs,     red_locs,     world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(npo_locs, red_npo_locs, world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(elmdist,  red_elmdist,  world_size+1, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    
    for(int i=0;i<world_size;i++)
    {
        npo_offset[i+1] = npo_offset[i]+red_npo_locs[i];
    }
    
    red_elmdist[world_size] = nel;
    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    std::cout << "nloc = " << nloc << std::endl;
    eptr[0]  = 0;
    
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+eltype[offset+i];
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = test[npo_offset[world_rank]+j];
            //std::cout << eind[j] << " ";
        }
        //std::cout << std::endl;
    }
    
    /*
    //map< pair<int, int>, HalfEdge* > HE = GetHalfEdges(test,eptr,nloc,offset);
    //map< pair<int, int>, HalfEdge* >::iterator it;
    

    if (world_rank == 0)
    {
        int cnt = 0;
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        for(it = HE.begin(); it != HE.end(); it++)
        {
            //std::cout << nloc << " " << cnt << std::endl;
            std::cout << cnt << " " << it->first.first << " " << it->first.second << " (" << it->second->oppositeHalfEdge->vertex << " " << it->second->oppositeHalfEdge->opposite_vertex <<")" << std::endl;
            cnt++;
            
        }
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    */
    /*
    if (world_rank == 2)
    {
        std::cout << "==========================================="<<std::endl;
        for(int i=0;i<nloc;i++)
        {
            for(int j=eptr[i];j<eptr[i+1];j++)
            {
                int vid = eind[j];
                std::cout << vid << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "==========================================="<<std::endl;
    }
    */
    /*
    if (world_rank == 3)
    {
        for(int i=0;i<nloc;i++)
        {
            for(int j=eptr[i];j<eptr[i+1];j++)
            {
                std::cout << eind[j] << " ";
            }
            
            std::cout << std::endl;
        }
    }
    */
    //===================================================================================================================================
    
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
    idx_t *adjwgt    = NULL;
    idx_t *vsize     = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt;
    
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }
    
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    
    idx_t part_[]    = {nloc};
    idx_t *part      = part_;
    
    
    ParMETIS_V3_PartMeshKway(red_elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
    if(world_rank == 1
       )
    {
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    
    
    ParMETIS_V3_Mesh2Dual(red_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
    ParMETIS_V3_AdaptiveRepart(red_elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
    int rank = 0;
    if(world_rank == 1)
    {
        std::cout << std::endl;
        
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    //===================================================================================================================================
    
    if(world_rank == rank)
    {
    
        //std::cout << "rank :: " << world_rank << std::endl;
        for(int i=0;i<nloc+1;i++)
        {
            std::cout << xadj[i] << " --> ";
        }
        
        for(int j=0;j<xadj[nloc];j++)
        {
            std::cout << adjncy[j] << " ";
        }
    }
}



void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs,MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    if (world_rank == 0)
    {
        std::cout << "printing boundary data..." << std::endl;
        for(int i=0;i<zdefs->nloc;i++)
        {
            for(int j=0;j<znames->ncol;j++)
            {
                std::cout << znames->getVal(i,j) << "";
            }
            std::cout << " :: ";
            for(int j=0;j<zdefs->nloc;j++)
            {
                std::cout << zdefs->getVal(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
    
}

void WriteBoundaryDataInParallel(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, Array<double>* Variables, MPI_Comm comm)
{
    /*
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    for(int bc=3;bc<zdefs->nrow;bc++)
    {
        int b_offset = zdefs->getVal(2,4);
        map< int, int > Loc2GlobBound;
        map< int, Vert > BC_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        
        int nloc     = int(n_bc_faces/world_size) + ( world_rank < n_bc_faces%world_size );
        //  compute offset of rows for each proc;
        int offset   = world_rank*int(n_bc_faces/world_size) + MIN(world_rank, n_bc_faces%world_size);
        
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    cnt++;
                }
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + "_" + std::to_string(world_rank) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\" \"dJ\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
    }
    
    */
}



void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts)
{
    for(int bc=3;bc<zdefs->nglob;bc++)
    {
        int b_offset = zdefs->getVal(2,4);
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        int* pltJ = new int[n_bc_faces*4];
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    J_verts[cnt]  = detJ_verts[val-1];
                    V_verts[cnt]  = vol_verts[val-1];
                    Jno_verts[cnt]  = Jnorm_verts[val-1];
                    cnt++;
                }
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dJ\", \"Vol\", \"dJ/Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << J_verts[(i)] << "   " << V_verts[(i)] << "   " << Jno_verts[(i)] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
    }
}





int largest(int arr[], int n)
{
    int i;
      
    // Initialize maximum element
    int max = arr[0];
  
    // Traverse array elements
    // from second and compare
    // every element with current max
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
  
    return max;
}


void WriteBoundaryDataInSerial3(Array<double>* xcn)
{
    string filename = "boundary_nodes.dat";
    ofstream myfile;
    myfile.open(filename);
    
    
    for(int i=0;i<3097156;i++)
    {
        myfile << xcn->getVal(i,0) << " " << xcn->getVal(i,1) << " " << xcn->getVal(i,2) << std::endl;
    }
    myfile.close();
}



void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts)
{
    for(int bc=3;bc<zdefs->nloc;bc++)
    {
        int b_offset = zdefs->getVal(2,4);
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        int* pltJ = new int[n_bc_faces*4];
        int teller = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                pltJ[teller] = val;
                //std::cout << val << " ";
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    V_verts[cnt]   = vol_verts[val-1];
                    cnt++;
                }
                
                teller=teller+1;
            }
            //std::cout << std::endl;
            tel++;
        }
        cout << "\nlargest id = " << largest(pltJ,n_bc_faces*4) << std::endl;
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        std::cout << "number of boundary nodes -> " << BC_verts.size() << std::endl;
        
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << V_verts[i] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
    }
}









ParVar* ComputeParallelStateArray(int nel,MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int nloc     = int(nel/world_size) + ( world_rank < nel%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(nel/world_size) + MIN(world_rank, nel%world_size);
    
    int* locs        = new int[world_size];
    
    for(int i=0;i<world_size;i++)
    {
        if (i==world_rank)
        {
            locs[i]     = nloc;
        }
        else
        {
            locs[i]     = 0;
        }
    }
    
    int* offsets = new int[world_size+1];
    for(int i=0;i<world_size;i++)
    {
        if (i==world_rank)
        {
            offsets[i]    = offset;
        }
        else
        {
            offsets[i]    = 0;
        }
    }
    
    int* red_locs       = new int[world_size];
    int* red_offsets    = new int[world_size+1];
    
    for(int i=0;i<world_size;i++)
    {
        red_locs[i]    = 0;
        red_offsets[i] = 0;
    }
    
    MPI_Allreduce(locs,     red_locs,     world_size,   MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(offsets,  red_offsets,  world_size, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
    
    red_offsets[world_size] = nel;
    ParVar* pv = new ParVar;
    
    pv->size    = world_size;
    pv->nlocs   = red_locs;
    pv->offsets = red_offsets;
    
    return pv;
    
}


int FindRank(int* arr, int size, int val)
{
    int start = 0;
    int last  = size;
    
    int mid = (start+last)/2;
    
    int found = 1;
    while (start<=last)
    {
        if (arr[mid]<=val)
        {
            start = mid + 1;
        }
        else
        {
            last = mid - 1;
        }
        mid = (start+last)/2;
    }
    return mid;
}



map< int, std::set<int> > GetRequestedVertices(Array<int>* ien, Array<double>* xcn, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    map< int, std::set<int> > Request;
    int rank_f = 0;
    std::vector<std::vector<int> > req;
    
    ParVar* pv_xcn = ComputeParallelStateArray(xcn->nglob, comm);
    int val;
    
    std::clock_t start;
    double duration;
    start = std::clock();
    
    if(world_rank == 0)
    {
        for(int i=0;i<world_size+1;i++)
        {
            std::cout << "offsets ien " << pv_xcn->offsets[i] << std::endl;
        }
    }
    
    //std::cout << "size decreases right? " << world_rank << " " << ien->nloc << std::endl;
    //ofstream myfile,myfile2;
    //myfile.open("val_"+std::to_string(world_rank)+".dat");

    //myfile2.open("rank_"+std::to_string(world_rank)+".dat");
    
    for(int i=0;i<ien->nloc;i++)
    {
        //std::cout << pv_xcn->offsets[world_rank] << " " << pv_xcn->offsets[world_rank+1] << "   ::  ";
        for(int j=1;j<ien->ncol;j++)
        {
            val = ien->getVal(i,j);
            //rank_f = FindRank(pv_xcn->offsets,world_size+1,val);
            //myfile <<val<< " ";

            //myfile2 <<rank_f << " ";
            if(val>pv_xcn->offsets[world_rank+1] || val < pv_xcn->offsets[world_rank])
            {
                rank_f = FindRank(pv_xcn->offsets,world_size+1,val);
                
                if ( Request[rank_f].find( val ) == Request[rank_f].end() )
                {
                    Request[rank_f].insert(val);
                }
            }
        }
        //myfile << std::endl;
        //myfile2 << std::endl;
    }
    
    
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    //std::cout << world_rank << " loopie loop = " << duration << std::endl;
    
    return Request;
}


TmpStruct* GetSchedule(Array<int>* ien, Array<double>* xcn, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    

    
    
    map< int, std::set<int> > Request = GetRequestedVertices(ien,xcn,comm);
    map<int,  std::set<int> >::iterator it;
    int num_req_proc = Request.size();
    
    /*
    for(it=Request.begin();it!=Request.end();it++)
    {
        std::cout << "sizing :: " << " " << world_rank << " " << it->second.size() << std::endl;
    }
    */
    std::clock_t start;
    double duration;
    start = std::clock();
    int* collect = new int[num_req_proc+1];
    int* sizing  = new int[num_req_proc];
    int tot_req_proc = 0;
    
    
    int* num_req_procs = new int[world_size];
    int* red_num_req_procs = new int[world_size];
    int* num_req_sizing = new int[world_size];
    int* red_num_req_sizing = new int[world_size];
    int* proc_offset = new int[world_size];
    int* proc_offset_sizing = new int[world_size];
    for(int i=0;i<world_size;i++)
    {
        if(i==world_rank)
        {
            num_req_procs[i] = num_req_proc+1;
            num_req_sizing[i] = num_req_proc;
        }
        else
        {
            num_req_procs[i]=0;
            num_req_sizing[i] = 0;
        }
    }
    
    MPI_Allreduce(num_req_procs, red_num_req_procs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(num_req_sizing, red_num_req_sizing, world_size, MPI_INT, MPI_SUM, comm);
    /*
    for(int i=0;i<world_size;i++)
    {
        std::cout << world_rank << " " << red_num_req_sizing[i] << std::endl;
    }
    
    std::cout << "=====" << std::endl;
     */
    int offset = 0;
    int offset_sizing = 0;
    for(int i=0;i<world_size;i++)
    {
        proc_offset[i] = offset;
        proc_offset_sizing[i] = offset_sizing;
        
        offset = offset + red_num_req_procs[i];
        offset_sizing = offset_sizing+red_num_req_sizing[i];
    }
    
    
    
    MPI_Allreduce(&num_req_proc, &tot_req_proc, 1, MPI_INT, MPI_SUM, comm);
    int* reduce_req_procs = new int[world_size+tot_req_proc];
    int* reduce_req_sizing = new int[world_size+tot_req_proc];
    collect[0] = world_rank;
    sizing[0] = -1;
    int t = 1;
    for(it = Request.begin(); it != Request.end(); it++)
    {
        collect[t] = it->first;
        sizing[t] = it->second.size();
        t++;
    }
    
    //std::cout << world_rank << " " << num_req_proc+1 << std::endl;
    MPI_Gatherv(&collect[0], (num_req_proc+1), MPI_INT, &reduce_req_procs[0], red_num_req_procs, proc_offset, MPI_INT, 0, comm);
    MPI_Gatherv(&sizing[0], (num_req_proc+1), MPI_INT, &reduce_req_sizing[0], red_num_req_procs, proc_offset, MPI_INT, 0, comm);
    MPI_Bcast(reduce_req_procs, world_size+tot_req_proc, MPI_INT, 0, comm);
    MPI_Bcast(reduce_req_sizing, world_size+tot_req_proc, MPI_INT, 0, comm);
    
    Array<int>* schedule = new Array<int>(world_size+tot_req_proc,1);
    //JaggedArray<int> schedule_jag = new JaggedArray<int>(nrow);
    schedule->data = reduce_req_procs;
    //schedule->sizing = reduce_req_sizing;
    
    if (world_rank == 0)
    {
        for(int i=0;i<tot_req_proc;i++)
        {
            std::cout<< reduce_req_sizing[i] << std::endl;
        }
        
        std::cout << "++++++++++++++++++++++++" << std::endl;
    }
    
    
    
    TmpStruct* t_struct = new TmpStruct;
    t_struct->data = reduce_req_procs;
    t_struct->sizing = reduce_req_sizing;
    t_struct->offsets = proc_offset;
    t_struct->nlocs = red_num_req_procs;
    t_struct->offsets_sizing = proc_offset_sizing;
    t_struct->nlocs_sizing = red_num_req_sizing;
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout << world_rank << " scheduling inside = " << duration << std::endl;
    
    return t_struct;
}


void TestReadInParallelToRoot(MPI_Comm comm, MPI_Info info)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::clock_t start;
    double duration;

    start = std::clock();
    Array<int>*   iee    = ReadDataSetFromFile<int>("grids/conn.h5","iee");
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " reading_serial = " << duration << std::endl;
    
    start = std::clock();
    Array<int>*   iee_r  = ReadDataSetFromFileInParallelToRoot<int>("grids/conn.h5","iee",comm,info);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " reading_par= " << duration << std::endl;
    
    //little test testing the reading
    if(world_rank == 0)
    {
        for(int i=0;i<iee_r->nloc;i++)
        {
            for(int j=0;j<iee_r->ncol;j++)
            {
                if((iee->getVal(i,j)-iee_r->getVal(i,j))!=0)
                {
                    std::cout << "not the same" << std::endl;
                }
            }
        }
    }
}


void TestBrutePartioningUS3D()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    const char* fn_conn="grids/conn.h5";
    const char* fn_grid="grids/grid.h5";
    
    Array<double>*   xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
    Array<int>*      ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    
    map< int, std::set<int> > Request = GetRequestedVertices(ien,xcn,comm);
    map<int,  std::set<int> >::iterator it;
    
    
    TmpStruct* schedule = GetSchedule(ien,xcn,comm);

    std::map< int, set<int> > s_send;
    std::map< int, set<int> > s_recv;
    std::map< int, std::vector<int> > s_send_size;
    std::map< int, std::vector<int> > s_recv_size;

    std::map< int, std::map< int, int> > s_alloc;
    std::map< int, std::map< int, int> > s_recv_alloc;
    //std::cout << "we are jagged " << std::endl;
    int nval_tot = 0;
    for(int i=0;i<world_size;i++)
    {

       int nloc   = schedule->nlocs[i];
       int offset = schedule->offsets[i];

       int nloc_sizing   = schedule->nlocs_sizing[i];
       int offset_sizing = schedule->offsets_sizing[i];
        
       for(int j=offset+1;j<(offset+nloc);j++)
       {
           s_send[i].insert(schedule->data[j]);
           s_recv[schedule->data[j]].insert(i);

           s_send_size[i].push_back(schedule->sizing[j]);
           s_recv_size[schedule->data[j]].push_back(schedule->sizing[j]);

           s_alloc[i][schedule->data[j]] = schedule->sizing[j];
           s_recv_alloc[schedule->data[j]][i] = schedule->sizing[j];
       }
    }



    std::map< int, int> loc_alloc = s_recv_alloc[world_rank];


    //std::cout << " loc_alloc.size()  = " << loc_alloc.size() << std::endl;
    int* recv_offset = new int[loc_alloc.size()+1];
    recv_offset[0]   = 0;
    int* recv_loc    = new int[loc_alloc.size()];;
    int recv_size    = 0;
    int i = 0;

    std::map< int, int> offset_map;
    std::map< int, int> loc_map;
    std::map< int, int>::iterator it_loc;
    
    for(it_loc=loc_alloc.begin();it_loc!=loc_alloc.end();it_loc++)
    {
        recv_loc[i]           = it_loc->second;
        recv_offset[i+1]      = recv_offset[i]+recv_loc[i];
        recv_size = recv_size+it_loc->second;

        loc_map[it_loc->first]=recv_loc[i];
        offset_map[it_loc->first]=recv_offset[i];

        i++;
    }

    int* recv_collector = new int[recv_size];

    int n_req_recv;

    for(int q=0;q<world_size;q++)
    {
        if (world_rank == q)
        {
            int tel = 0;
            int i = 0;
            for (it = Request.begin(); it != Request.end(); it++)
            {
                int n_req = it->second.size();
                int* req_arr_send = new int[n_req];

                int dest   = it->first;
                std::set<int>::iterator it_set;

                for(it_set=it->second.begin();it_set != it->second.end();++it_set)
                {
                    req_arr_send[i] = *it_set;
                }
                std::cout << "wr = " << q << " dest " << dest << " " << n_req << " ";

                MPI_Send(&n_req, 1, MPI_INT, dest, dest, comm);
                MPI_Send(&req_arr_send[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                tel = tel + 1;
                i++;
            }
        }
        else if (s_send[q].find( world_rank ) != s_send[q].end())
        {
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, world_rank, comm, MPI_STATUS_IGNORE);

            MPI_Recv(&recv_collector[offset_map[q]], n_req_recv, MPI_INT, q, 100+world_rank*2, comm, MPI_STATUS_IGNORE);

        }
    }
}

void TestFindRank()
{
    int* arr = new int[10];
    int* arr_recv = new int[10];
    
    
    arr[0] = 0;
    arr[1] = 10;
    arr[2] = 12;
    arr[3] = 24;
    arr[4] = 42;
    arr[5] = 55;
    arr[6] = 65;
    arr[7] = 66;
    arr[8] = 78;
    arr[9] = 81;
    
    
    int value = 25;
    
    int res = FindRank(arr,10,value);
    if(res == 3)
    {
        std::cout << "TestFindRank() has passed!" << std::endl;
    }
    else
    {
        std::cout << "TestFindRank() has failed!" << std::endl;

    }
    
}






void Example3DPartitioningWithParVarParMetis()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    int nel = 8;
    int * eltype = new int[nel];
    
    eltype[0] = 8;
    eltype[1] = 8;
    eltype[2] = 8;
    eltype[3] = 8;
    eltype[4] = 8;
    eltype[5] = 8;
    eltype[6] = 8;
    eltype[7] = 8;
    
    int npo   = 0;
    
    for(int i = 0;i < nel; i++)
    {
        npo += eltype[i];
    }
    
    int * test = new int[npo];
    
    Array<int>* e2n = new Array<int>(8,8);
    
    e2n->setVal(0,0,0);e2n->setVal(0,1,1);e2n->setVal(0,2,6);e2n->setVal(0,3,5);    e2n->setVal(0,4,14+0);e2n->setVal(0,5,14+1);e2n->setVal(0,6,14+6);e2n->setVal(0,7,14+5);
    e2n->setVal(1,0,1);e2n->setVal(1,1,2);e2n->setVal(1,2,7);e2n->setVal(1,3,6);    e2n->setVal(1,4,14+1);e2n->setVal(1,5,14+2);e2n->setVal(1,6,14+7);e2n->setVal(1,7,14+6);
    e2n->setVal(2,0,2);e2n->setVal(2,1,3);e2n->setVal(2,2,8);e2n->setVal(2,3,7);    e2n->setVal(2,4,14+2);e2n->setVal(2,5,14+3);e2n->setVal(2,6,14+8);e2n->setVal(2,7,14+7);
    e2n->setVal(3,0,3);e2n->setVal(3,1,4);e2n->setVal(3,2,9);e2n->setVal(3,3,8);    e2n->setVal(3,4,14+3);e2n->setVal(3,5,14+4);e2n->setVal(3,6,14+9);e2n->setVal(3,7,14+8);
    e2n->setVal(4,0,5);e2n->setVal(4,1,6);e2n->setVal(4,2,11);e2n->setVal(4,3,10);  e2n->setVal(4,4,14+5);e2n->setVal(4,5,14+6);e2n->setVal(4,6,14+11);e2n->setVal(4,7,14+10);
    e2n->setVal(5,0,6);e2n->setVal(5,1,7);e2n->setVal(5,2,12);e2n->setVal(5,3,11);  e2n->setVal(5,4,14+6);e2n->setVal(5,5,14+7);e2n->setVal(5,6,14+12);e2n->setVal(5,7,14+11);
    e2n->setVal(6,0,7);e2n->setVal(6,1,8);e2n->setVal(6,2,13);e2n->setVal(6,3,12);  e2n->setVal(6,4,14+7);e2n->setVal(6,5,14+8);e2n->setVal(6,6,14+13);e2n->setVal(6,7,14+12);
    e2n->setVal(7,0,8);e2n->setVal(7,1,9);e2n->setVal(7,2,14);e2n->setVal(7,3,13);  e2n->setVal(7,4,14+8);e2n->setVal(7,5,14+9);e2n->setVal(7,6,14+14);e2n->setVal(7,7,14+14);

    ParVar_ParMetis* pv_parmetis = CreateParallelDataParmetis(e2n,comm,8);
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;
    idx_t *adjwgt    = NULL;
    idx_t *vsize     = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt;
    
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }
    
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    
    idx_t part_[]    = {pv_parmetis->nlocs[world_rank]};
    idx_t *part      = part_;
    
    
    ParMETIS_V3_PartMeshKway(pv_parmetis->elmdist, pv_parmetis->eptr, pv_parmetis->eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
//    if(world_rank == 1)
//    {
//        for(int i = 0; i < pv_parmetis->nlocs[world_rank]; i++)
//        {
//            std::cout << part[i] << std::endl;
//        }
//    }
    
    ParMETIS_V3_Mesh2Dual(pv_parmetis->elmdist,pv_parmetis->eptr,pv_parmetis->eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
    ParMETIS_V3_AdaptiveRepart(pv_parmetis->elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
    
//    if(world_rank == 1)
//    {
//        for(int i = 0; i < pv_parmetis->nlocs[world_rank]; i++)
//        {
//            std::cout << part[i] << std::endl;
//        }
//    }
    
    int rank = 0;
    
    int* arr_res = new int[9];
    arr_res[0] = 1;arr_res[1] = 4;arr_res[2] = 0;arr_res[3] = 2;
    arr_res[4] = 5;arr_res[5] = 1;arr_res[6] = 3;arr_res[7] = 6;arr_res[8] = 2;
    int fail = 0;
    
    if(world_rank == rank)
    {
    
//        for(int i=0;i<pv_parmetis->nlocs[world_rank]+1;i++)
//        {
//            std::cout << xadj[i] << " --> ";
//        }
        
        for(int j=0;j<xadj[pv_parmetis->nlocs[world_rank]];j++)
        {
//            std::cout << adjncy[j] << " ";
            
            if(adjncy[j]-arr_res[j]!=0)
            {
                fail = 1;
            }
        }
        
        if (fail == 0)
        {
            std::cout << "The ParVar_ParMetis test has passed!" << std::endl;
        }
        else
        {
            std::cout << "The ParVar_ParMetis test has failed!" << std::endl;
        }
    }
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
    
    //GetXadjandAdjcyArrays(iee,ien,comm);
    
//============================================================
    
    const char* fn_conn="grids/conn.h5";
    const char* fn_grid="grids/grid.h5";
    
    //Array<int>*    zdefs = ReadDataSetFromGroupFromFile<int>(fn_conn,"zones","zdefs");
    //Array<char>*  znames = ReadDataSetFromGroupFromFile<char>(fn_conn,"zones","znames");
    //PlotBoundaryData(znames,zdefs,comm);
    //TestBrutePartioningUS3D();
    
    if (world_rank == 0)
    {
        TestFindRank();
    }
    Example3DPartitioningWithParVarParMetis();

    
    Array<double>*   xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
    Array<int>*      ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    
    //std::cout << ien->nglob << std::endl;
    int N = ien->nglob;
    //ParVar_ParMetis* pv_parmetis = CreateParallelDataParmetis(N,comm,8);
    
//    if (world_rank == 0)
//    {
//        for(int i=0;i<world_size+1;i++)
//        {
//            if(i < world_size)
//            {
//                std::cout << pv_parmetis->npo_locs[i] << " " << pv_parmetis->nlocs[i] << " ";
//            }
//            std::cout << pv_parmetis->elmdist[i] << " ";
//            std::cout << std::endl;
//        }
//    }
    
    //Example3DPartitioning();
    //Example3DPartitioning();
    //WriteBoundaryDataInSerial2(xcn,zdefs,ifn,Vollie_verts);

    /*
    Array<double>*   var_el       = new Array<double>(ien->nrow,0);
    Array<double>*   var_loc_el   = new Array<double>(nloc,0);

    // In order to isolate a single boundary surface to plot in tecplot.
    // !!!!  This routine relies on the fact that we only have quads on the suface.  !!!!
    
    Array<double>*   var_vert = new Array<double>(xcn->nrow,0);
    Array<int>*      av_vert  = new Array<int>(xcn->nrow,0);
    Array<double>*   Jvert = new Array<double>(xcn->nrow,0);
    
    for(int i=0;i<xcn->nrow;i++)
    {
        av_vert->setVal(i,0,0);
        var_vert->setVal(i,0,0.0);
    }
    int Nnodes = xcn->nrow;
    int Nelement = ien->nrow;
    
    
    std::vector<std::vector<int>> N2E(Nnodes);
    std::vector<std::vector<int>> E2N(Nelement);
    std::vector<int> ElType(Nelement);
    int tel =0;
    std::clock_t start;
    double duration;
    start = std::clock();
    for(int i=0;i<Nelement;i++)
    {
        ElType[i]=ien->getVal(i,0);
        
        for(int j=0;j<8;j++)
        {
            //E2N[i].push_back(ien->getVal(i,j+1)-1);
            
            if(ien->getVal(i,j+1)-1 < Nnodes)
            {
                N2E[ien->getVal(i,j+1)-1].push_back(i);
            }
        }
    }
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer = " << duration << std::endl;
    
    
    for(int i=0;i<Nelement;i++)
    {
        ElType[i]=ien->getVal(i,0);
        
        for(int j=0;j<8;j++)
        {
            //E2N[i].push_back(ien->getVal(i,j+1)-1);
            
            if(ien->getVal(i,j+1)-1 < Nnodes)
            {
                N2E[ien->getVal(i,j+1)-1].push_back(i);
            }
        }
    }
    

    start = std::clock();
    double* Jaccie = ComputeDeterminantofJacobian(xcn,ien,nloc,offset,var_el);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer jaccie = " << duration << std::endl;
    
    start = std::clock();
    double* Vollie = ComputeVolumeCells(xcn,ien);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer vollie = " << duration << std::endl;
    int id;
    double sum = 0.0;double vu = 0.0;
    double sumv = 0.0;double vuv = 0.0;
    double* Jvert2 = new double[Nnodes];
    double* Vvert2 = new double[Nnodes];
    double* Jnorm  = new double[Nnodes];
    for(int i=0;i<N2E.size();i++)
    {
        sum = 0.0;
        sumv = 0.0;
        for(int j=0;j<N2E[i].size();j++)
        {
            id = N2E[i][j];
            sum = sum+Jaccie[id];
            sumv = sumv+Vollie[id];
            
            //std::cout << "id = " << id << " val = ";
            
            //std::cout << Jaccie[id] << " ";
        }
        //std::cout << std::endl;
        vu = sum/N2E[i].size();
        vuv = sumv/N2E[i].size();
        Jvert2[i] = vu;
        Vvert2[i] = vuv;
        Jnorm[i]  = vu/vuv;
    }
        
    WriteBoundaryDataInSerial(xcn,zdefs,ifn,Jvert2,Vvert2,Jnorm);
    */
    /*
    double *rbuf = new double[nloc];
    double *recv = new double[ien->nrow];
    
    for(int i=0;i<nloc;i++)
    {
        rbuf[i] = var->getVal(i,0);
    }
    start = std::clock();
    MPI_Gather(rbuf, nloc, MPI_DOUBLE, recv, nloc, MPI_DOUBLE, 0, comm);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer gather = " << duration << std::endl;
    
    Array<double>detJ_verts(xcn->nrow,0);
    std::vector< std::vector<double> > detJverts(xcn->nrow);
    start = std::clock();
    int Vid;
    double value;
    for(int i=0;i<ien->nrow;i++)
    {
        value = recv[i];
        for(int j=0;j<ien->ncol;j++)
        {
            Vid = ien->getVal(i,j);
            detJverts[Vid].push_back(value);
        }
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer reducing to verts = " << duration << std::endl;
     */
    //for(int i=0;i<ien->nrow;i++)
    //{
        //std::cout << i << " " << recv[i] << std::endl;
    //}
    
    // MPI_Gather(var_loc->data, nloc, MPI_DOUBLE, &var->data[offset], nloc, MPI_DOUBLE, 0, comm);
    
    
    //WriteBoundaryDataInSerial(xcn,zdefs,ifn,detJ_verts);
    
    /*
    for(int i=b_start;i<b_end;i++)
    {
        
        std::cout << "I = " << i << " :: ";
        for(int j=0;j<4;j++)
        {
            std::cout << Loc[i*4+j] << " ";
        }
        
        std::cout << " <--::--> ";
        
        for(int j=0;j<4;j++)
        {
            std::cout << ifn->getVal(i,j+1) << " ";
        }
        std::cout << std::endl;
    }
    */
    
    /*
    ofstream myfile;
    myfile.open("boundary.ascii");
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"J\" ";
    for(int i=2;i<3;i++)
    {
        int b_start = zdefs->getVal(i,3)-b_offset-1;
        int b_end   = zdefs->getVal(i,4)-b_offset-1;
        
        std::cout << "b_start " << b_start << " " << b_end << std::endl;
        
        
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=1;k<5;k++)
            {
                std::cout << ifn->getVal(j+b_offset,k) << " ";
                myfile << xcn->getVal(ifn->getVal(j+b_offset,k),0) << " " << xcn->getVal(ifn->getVal(j+b_offset,k),1) << " " << xcn->getVal(ifn->getVal(j+b_offset,k),2) << " ";
            }
            myfile << std::endl;
            std::cout << std::endl;
            //std::cout << bound->getVal(j,0) << " " << bound->getVal(j,1) << " " << bound->getVal(j,2) << std::endl;
        }
    }
    myfile.close();
     */
    //std::cout << bound->nrow << std::endl;
    //PlotBoundaryData(znames,zdefs);
    
    //                                                                                                            map< pair<int, int>, HalfEdge* > HE = GetHalfEdges(test,eptr,nloc,offset);
    //map< pair<int, int>, HalfEdge* >::iterator it;
    
    //UnitTestJacobian();
    //UnitTestEigenDecomp();
    MPI_Finalize();
    return 0;
     
}
