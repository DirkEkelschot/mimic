#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
//#include "petsc.h"
#include "array.h"
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
    int nel    = ien->nrow;
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
    
    ParMETIS_V3_Mesh2Dual(red_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
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
    
    int rank = 0;
    if(world_rank == rank)
    {
    
        //std::cout << "rank :: " << world_rank << std::endl;
        /*for(int i=0;i<nloc+1;i++)
        {
            std::cout << xadj[i] << " --> ";
        }
         */
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
}

void Example3DPartitioning(MPI_Comm comm)
{
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
    std::cout << nloc << std::endl;
    eptr[0]  = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+eltype[offset+i];
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = test[npo_offset[world_rank]+j];
            std::cout << eind[j] << " ";
        }
        std::cout << std::endl;
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
    
    /*
    ParMETIS_V3_PartMeshKway(red_elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    
    if(world_rank == 2)
    {
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    */
    
    ParMETIS_V3_Mesh2Dual(red_elmdist,eptr,eind,numflag,ncommonnodes,&xadj,&adjncy,&comm);
    
    idx_t *nparts2 = nparts_;
    
    //ParMETIS_V3_AdaptiveRepart(red_elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts2, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    
    if(world_rank == 1)
    {
        std::cout << std::endl;
        
        for(int i = 0; i < nloc; i++)
        {
            std::cout << part[i] << std::endl;
        }
    }
    //===================================================================================================================================
    int rank = 0;
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



void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs)
{
    for(int i=0;i<zdefs->nrow;i++)
    {
        for(int j=0;j<znames->ncol;j++)
        {
            std::cout << znames->getVal(i,j) << "";
        }
        std::cout << " :: ";
        for(int j=0;j<zdefs->ncol;j++)
        {
            std::cout << zdefs->getVal(i,j) << " ";
        }
        std::cout << std::endl;
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
    for(int bc=3;bc<zdefs->nrow;bc++)
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











void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts)
{
    for(int bc=3;bc<zdefs->nrow;bc++)
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
                    
                    V_verts[cnt]   = vol_verts[val-1];
                    cnt++;
                }
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
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











int main(int argc, char** argv) {
    
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::clock_t start;
    double duration;
    
    
    //Array<int> iee   = ReadDataSetFromFileInParallel<int>("../data_files/conn.h5","ife",comm,info);
    //Array<int> ife   = ReadDataSetFromFileInParallel<int>("../data_files/conn.h5","ife",comm,info);
    //Array<double> bound_par = ReadDataSetFromRunInFileInParallel<double>("../data_files/data_r.h5","run_1","boundaries",comm,info);
    //Array<double>* bound = ReadDataSetFromRunInFile<double>("../data_files/adept_geom/data_r.h5","run_1","boundaries");
    
    
    
    // Read in adept geometry.
    
    /*
    Array<int>*    zdefs = ReadDataSetFromGroupFromFile<int>("../data_files/adept_geom/conn.h5","zones","zdefs");
    Array<char>*  znames = ReadDataSetFromGroupFromFile<char>("../data_files/adept_geom/conn.h5","zones","znames");
    Array<int>*      ifn = ReadDataSetFromFile<int>("../data_files/adept_geom/grid.h5","ifn");
    Array<double>*   xcn = ReadDataSetFromFile<double>("../data_files/adept_geom/grid.h5","xcn");
    //Array<int>*      ien = ReadDataSetFromFileInParallel<int>("../data_files/adept_geom/conn.h5","ien",comm,info);
    Array<int>*      ien = ReadDataSetFromFile<int>("../data_files/adept_geom/conn.h5","ien");
     */
     
    Array<int>*    zdefs = ReadDataSetFromGroupFromFile<int>("../data_files/conn_kl4.h5","zones","zdefs");
    Array<char>*  znames = ReadDataSetFromGroupFromFile<char>("../data_files/conn_kl4.h5","zones","znames");
    Array<int>*      ifn = ReadDataSetFromFile<int>("../data_files/grid_kl4.h5","ifn");
    
    //Array<double>*   xcn = ReadDataSetFromFile<double>("../data_files/grid_pstn.h5","xcn");
    //Array<int>*      ien = ReadDataSetFromFile<int>("../data_files/conn_pstn.h5","ien");
    
    
    Array<double>*   xcn = ReadDataSetFromFile<double>("../data_files/adept_geom/grid.h5","xcn");
    Array<int>*      ien = ReadDataSetFromFile<int>("../data_files/adept_geom/conn.h5","ien");
    
    int Nnodes = xcn->nrow;
    int Nelement = ien->nrow;
    //std::cout << "Nelement = " << Nelement << " " << Nnodes << std::endl;
    int nloc     = int(ien->nrow/world_size) + ( world_rank < ien->nrow%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(ien->nrow/world_size) + MIN(world_rank, ien->nrow%world_size);
    /*
    // returns array with volumes for elements with length of local number of elements;
    start = std::clock();
    double* Vollie = ComputeVolumeCells(xcn,ien,comm);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer vollie cell = " << duration << std::endl;
    
    start = std::clock();
    double* Vollie_verts = ComputeVolumeCellsReducedToVerts(xcn,ien);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer vollie vert = " << duration << std::endl;
    */
    /*
    double *recv = new double[ien->nrow];
    
    Array<double>* Volumes = new Array<double>(ien->nrow,0);
    Array<double>* Jvert2  = new Array<double>(ien->nrow,0);
    Array<double>* Vvert2  = new Array<double>(ien->nrow,0);
    
    MPI_Gather(Vollie, nloc, MPI_DOUBLE, recv, nloc, MPI_DOUBLE, 0, comm);
    */
    
    //Example3DPartitioning(comm);
    
    
    start = std::clock();
    PartitionBigMesh(xcn,ien,comm);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer partitioning = " << duration << std::endl;
    
    start = std::clock();
    
    double* xadjn = new double[nloc*8];
    for(int i=0;i<nloc;i++)
    {
        for(int j=0;j<8;j++)
        {
            xadjn[i*8+j] = ien->getVal(offset+i,j+1)-1;
        }
    }
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << world_rank << " timer partitioning myself = " << duration << std::endl;
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
    
    //Example3DPartitioning(comm);
    //UnitTestJacobian();
    //UnitTestEigenDecomp();
    MPI_Finalize();
    return 0;
     
}
