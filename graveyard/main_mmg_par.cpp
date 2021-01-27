#include "adapt_io.h"
#include "adapt_compute.h"
#include "adapt_operations.h"
#include "adapt_math.h"
#include <math.h>
#include "adapt_output.h"
#include "adapt_partition.h"
//#include "adapt.h"
#include <iomanip>
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"
#include "parmmg/libparmmg.h"


int mpi_size, mpi_rank;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

namespace lapack {
    extern "C" {
      void dgeev_(char const * __restrict JOBVL, char const * __restrict JOBVR, int const * __restrict n, double * __restrict A, int const * lda, double * __restrict WR, double * __restrict WI, double * __restrict VL, int const * __restrict ldvl, double * __restrict VR, int const * __restrict ldvr, double * __restrict Work, int const * __restrict lwork, int       * __restrict info );
      // LU decomoposition of a general matrix
      void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

      void dgemm_(char * transA,  char * transB, int * m, int * n, int * k, double * alpha, double * A, int * lda, double * B, int * ldb, double * beta, double * C, int * ldc );
      // generate inverse of a matrix given its LU decomposition
      void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);


        void dgeqrf_(int* M, int* N,
                     double* A, int* LDA, double* TAU,
                     double* WORK, int* LWORK, int* INFO );

        void dormqr_(char *side, char *trans, int* m, int *n,
                     int *k, double* A, int* lda, double* tau, double* C,
                     int* ldc, double *WORK, int* lwork, int* info);

        void dtrtrs_(char *UPLO, char *TRANS, char *DIAG, int* N, int *NRHS, double* A, int* lda, double* B, int* ldb, int* info);
    }

    int geqrf(int m, int n,
              double* A, int lda, double *tau)
    {
        int info=0;
        int lwork=-1;
        double iwork;
        dgeqrf_(&m, &n, A, &lda, tau,
                        &iwork, &lwork, &info);
        lwork = (int)iwork;
        double* work = new double[lwork];
        dgeqrf_(&m, &n, A, &lda, tau,
                        work, &lwork, &info);
        delete[] work;
        return info;
    }

    int ormqr(char side, char trans, int m, int n, int k,
              double *A, int lda, double *tau, double* C, int ldc)
    {
        int info=0;
        int lwork=-1;
        double iwork;
        dormqr_(&side, &trans, &m, &n, &k,
                A, &lda, tau, C, &ldc, &iwork, &lwork, &info);
        lwork = (int)iwork;
        double* work = new double[lwork];
        dormqr_(&side, &trans, &m, &n, &k,
                A, &lda, tau, C, &ldc, work, &lwork, &info);
        delete[] work;
        return info;
    }

    int trtrs(char uplo, char trans, char diag,
              int n, int nrhs,
              double* A, int lda, double* B, int ldb)
    {
        int info = 0;
        dtrtrs_(&uplo, &trans, &diag, &n, &nrhs,
                A, &lda, B, &ldb, &info);
        return info;
    }
}


using namespace std;


bool isDiagonalMatrix(Array<double>* Msq)
{
    
    for (int i = 0; i < Msq->getNrow(); i++)
        for (int j = 0; j < Msq->getNrow(); j++)
            // condition to check other elements
            // except main diagonal are zero or not.
            if ((i != j) && (Msq->getVal(i,j) != 0))
                return false;
    return true;
}

Array<double>* MatInv(Array<double>* A)
{
    int n = A->getNrow();
    int size = n*n;
    double WORK [size];
    int info;
    int Pivot[n];
    Array<double>* R = new Array<double>(n,n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            R->setVal(i,j,A->getVal(i,j));
        }
    }
    
    lapack::dgetrf_(&n, &n, R->data, &n, Pivot, &info);
    lapack::dgetri_(&n, R->data, &n, Pivot, WORK, &size, &info);
    
    return R;
}

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

int largest(int* arr, int n)
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


int binarySearch(int* arr, int low, int high, int key)
{
    if (high < low)
        return -1;
    int mid = (low + high) / 2; /*low + (high - low)/2;*/
    if (key == arr[mid])
        return mid;
    if (key > arr[mid])
        return binarySearch(arr, (mid + 1), high, key);
    
    return binarySearch(arr, low, (mid - 1), key);
}

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
  lapack::dgeev_( &JOBVL, &JOBVR, &n, V, &n, WR, WI, iV, &n, NULL, &n, WORK, &size, &info );
  
  // Copy right eigenvectors into V (with transpose)
  for ( i = 0; i < n; i++)
    for ( j = 0; j < n; j++)
      V[i*n+j] = iV[j*n+i];
  
  // Compute inverse of V1 
  memcpy( iV, V, n*n*sizeof(double) );
  lapack::dgetrf_(&n, &n, iV, &n, Pivot, &info);
  lapack::dgetri_(&n, iV, &n, Pivot, WORK, &size, &info);
}






//typedef map<pair<int, int>, HalfEdge *> MapType;

map< pair<int,int>, HalfEdge* > GetHalfEdges(int* element_verts, int* M, int nloc, int offset)
{
    map< pair<int, int>, HalfEdge* > Edges;
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
    M[0] = 0.25;M[1]=0.1;M[2]=0.14;
    M[3] = 0.2;M[4]=0.25;M[5]=0.0;
    M[6] = 0.4;M[7]=0.3;M[8]=0.25;
    
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








void PlotBoundaryData(Array<char>* znames, Array<int>* zdefs,MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int nrow = zdefs->getNrow();
    int ncol = znames->getNcol();
    if (world_rank == 0)
    {
        std::cout << "printing boundary data..." << nrow << " " << zdefs->getNcol() << std::endl;
        for(int i=0;i<nrow;i++)
        {
            for(int j=0;j<ncol;j++)
            {
                std::cout << znames->getVal(i,j) << "";
            }
            std::cout << " :: ";
            for(int j=0;j<zdefs->getNcol();j++)
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





map< int, std::set<int> > GetRequestedVertices(ParArray<int>* ien, ParArray<double>* xcn, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    map< int, std::set<int> > Request;
    int rank_f = 0;
    std::vector<std::vector<int> > req;
    
    ParVar* pv_xcn = ComputeParallelStateArray(xcn->getNglob(), comm);
    int val;
    
    for(int i=0;i<ien->getNrow();i++)
    {
        for(int j=1;j<ien->getNcol();j++)
        {
            val = ien->getVal(i,j);
            
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
    delete pv_xcn;
    return Request;
}





TmpStruct* GetSchedule(ParArray<int>* ien, ParArray<double>* xcn, MPI_Comm comm)
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

    int* collect = new int[num_req_proc+1];
    int* sizing  = new int[num_req_proc+1];
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
            num_req_procs[i]  = num_req_proc+1;
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
        sizing[t]  = it->second.size();
        t++;
    }
    
    //std::cout << world_rank << " " << num_req_proc+1 << std::endl;
    MPI_Gatherv(&collect[0], (num_req_proc+1), MPI_INT, &reduce_req_procs[0], red_num_req_procs, proc_offset, MPI_INT, 0, comm);
    MPI_Gatherv(&sizing[0], (num_req_proc+1), MPI_INT, &reduce_req_sizing[0], red_num_req_procs, proc_offset, MPI_INT, 0, comm);
    MPI_Bcast(reduce_req_procs, world_size+tot_req_proc, MPI_INT, 0, comm);
    MPI_Bcast(reduce_req_sizing, world_size+tot_req_proc, MPI_INT, 0, comm);
    
    if(world_rank == 1)
    {
        for(int i=0;i<(num_req_proc+1);i++)
        {
            std::cout << "reduce_req_procs[ " << collect[i] << std::endl;
        }
    }
    
    Array<int>* schedule = new Array<int>(world_size+tot_req_proc,1);
    //JaggedArray<int> schedule_jag = new JaggedArray<int>(nrow);
    schedule->data = reduce_req_procs;
    //schedule->sizing = reduce_req_sizing;
    
    if (world_rank == 0)
    {
        for(int i=0;i<tot_req_proc+world_size;i++)
        {
            std::cout << reduce_req_procs[i] << " " << reduce_req_sizing[i] << std::endl;
        }

        std::cout << "++++++++++++++++++++++++" << std::endl;
        
        for(int i=0;i<world_size;i++)
        {
            std::cout << red_num_req_procs[i] << " "  << proc_offset_sizing[i] << std::endl;
        }
    }

    
    delete schedule;
    
    TmpStruct* t_struct = new TmpStruct;
    t_struct->data = reduce_req_procs;
    t_struct->sizing = reduce_req_sizing;
    t_struct->offsets = proc_offset;
    t_struct->nlocs = red_num_req_procs;
    t_struct->offsets_sizing = proc_offset_sizing;
    t_struct->nlocs_sizing = red_num_req_sizing;
    
    
    return t_struct;
}


void TestReadInParallelToRoot(MPI_Comm comm, MPI_Info info)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    

    //start = std::clock();
    Array<int>*   iee    = ReadDataSetFromFile<int>("grids/piston/conn.h5","iee");
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout << world_rank << " reading_serial = " << duration << std::endl;
    //start = std::clock();
    Array<int>*   iee_r  = ReadDataSetFromFileInParallelToRoot<int>("grids/piston/conn.h5","iee",comm,info);
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout << world_rank << " reading_par= " << duration << std::endl;
    
    //little test testing the reading
    
    if(world_rank == 0)
    {
        int tel = 0;
        for(int i=0;i<iee_r->getNrow();i++)
        {
            for(int j=0;j<iee_r->getNcol();j++)
            {
                if((iee->getVal(i,j)-iee_r->getVal(i,j))!=0)
                {
                    
                    std::cout << "not the same" << std::endl;
                    tel = tel + 1;
                }
            }
        }
        if (tel == 0)
        {
            std::cout << "Parallel reading test to root has passed!!!" << std::endl;
        }
    }
    
    delete iee_r;
    delete iee;
}


int* TestBrutePartioningUS3D()
{
    std::clock_t start;
    double duration;
    start = std::clock();
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    const char* fn_conn="grids/piston/conn.h5";
    const char* fn_grid="grids/piston/grid.h5";
    
    ParArray<double>*   xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
    ParArray<int>*      ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    //Array<double>*   ief = ReadDataSetFromFileInParallel<double>(fn_grid,"ief",comm,info);
    //ParArray<double>*   ifn = ReadDataSetFromFileInParallel<double>(fn_grid,"ifn",comm,info);

    
    map< int, std::set<int> > Request = GetRequestedVertices(ien,xcn,comm);
    map<int,  std::set<int> >::iterator it;
    
    TmpStruct* schedule = GetSchedule(ien,xcn,comm);

    std::map< int, set<int> > s_send;
    std::map< int, set<int> > s_recv;
    std::map< int, std::map< int, int> > s_recv_alloc;
    //std::cout << "we are jagged " << std::endl;
    for(int i=0;i<world_size;i++)
    {
       int nloc   = schedule->nlocs[i];
       int offset = schedule->offsets[i];
        
       for(int j=offset+1;j<(offset+nloc);j++)
       {
           s_send[i].insert(schedule->data[j]);
           s_recv[schedule->data[j]].insert(i);
           
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
        recv_loc[i]      = it_loc->second;
        recv_offset[i+1] = recv_offset[i]+recv_loc[i];
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
                MPI_Send(&n_req, 1, MPI_INT, dest, dest, comm);
                MPI_Send(&req_arr_send[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                tel = tel + 1;
                i++;
                delete[] req_arr_send;
            }
        }
        else if (s_send[q].find( world_rank ) != s_send[q].end())
        {
            
            MPI_Recv(&n_req_recv, 1, MPI_INT, q, world_rank, comm, MPI_STATUS_IGNORE);

            MPI_Recv(&recv_collector[offset_map[q]], n_req_recv, MPI_INT, q, 100+world_rank*2, comm, MPI_STATUS_IGNORE);

        }
    }
    
    delete schedule;
    delete xcn;
    delete ien;
    delete[] recv_loc;
    delete[] recv_offset;
    //delete ief;
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
        
    std::cout << "duration " << duration  << std::endl;
    return recv_collector;
    

}

void TestFindRank(MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int* arr = new int[10];
    
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
    
    int value = 82;
    
    int res = FindRank(arr,10,value);
    
    if(res == 3)
    {
        std::cout << "TestFindRank() has passed! " << res << std::endl;
    }
    else
    {
        
        std::cout << rank << "  TestFindRank() has failed! " << res << std::endl;
        std::cout << std::endl;
    }
    
}



void Example3DPartitioningWithParVarParMetis()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    ParArray<int>*   ien = ReadDataSetFromFileInParallel<int>("tools/test_grid.h5","ien",comm,info);
    
    int nel = ien->getNglob();
    int * eltype = new int[nel];

    int npo   = 0;

    for(int i = 0;i < nel; i++)
    {
        eltype[i] = 8;
        npo += eltype[i];
    }

    ParVar_ParMetis* pv_parmetis = CreateParallelDataParmetis(ien,comm,8);

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

    idx_t *elmwgt = NULL;

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
            std::cout << adjncy[j] << " ";

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
            std::cout << "The ParVar_ParMetis test has failed" << std::endl;
            std::cout << "=============Suggetions============" << std::endl;
            std::cout << " 1) Make sure you run this test using 4 processors!!!" << std::endl;
            std::cout << " 2) Make sure you run make_test_array.py in the tools folder which generates the test mesh." << std::endl;
            std::cout << "===================================" << std::endl;

        }
    }
    
    delete[] eltype;
    delete[] arr_res;
    //delete nparts2;
    //delete nparts;
    delete[] tpwgts;
    //delete wgtflag;
    
    //delete ncon;
    //delete itr;
    delete[] xadj;
    delete[] adjncy;
    
}

//
//This function generates the adjacency map based on the element2face and face2element maps which are
//For now this map is read in from us3d grid data files but could potentially be used in general.
//

/*

void GetAdjacencyForUS3D(Array<int>* ife, ParArray<int>* ief, int nelem, MPI_Comm comm)
{
    

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    
    int fid = 0;
    int eid = 0;
    ParallelState* pv = ief->getParallelState();
    std::map<int, set<int> > e2e;
    int nrow = ief->getNrow();
    int ncol = ief->getNcol();
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            fid = fabs(ief->getVal(i,j+1))-1;
            
            for(int k=0;k<ncol;k++)
            {
                eid = ife->getVal(fid,k)-1;
                //std::cout << eid << " ";
                if(e2e[i].find(eid) == e2e[i].end() && (eid<nelem) && eid!=i+pv->getOffset(world_rank))
                {
                    e2e[i].insert(eid);
                }
            }
        }
    }
    
    set<int>::iterator it;
    for(int i=0;i<nrow;i++)
    {
        std::cout << world_rank << " " << i << " :: ";
        for (it=e2e[i].begin(); it != e2e[i].end(); ++it)
        {
            cout << *it << " " ;
        }

        std::cout << std::endl;
    }
}
*/
//
//This function generates the adjacency map based on the element2face map.adjacency
//For now this map is read in from us3d grid data files but could potentially be used in general.
//

void GetAdjacencyForUS3D_V2(ParArray<int>* ief, int nelem, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int fid = 0;
    
    std::map<int, set<int> > f2e;
    std::map<int, std::vector<int> > f2e_vec;
    std::map<int, set<int> > e2f;
    std::map<int, std::vector<int> > e2f_vec;

    int nrow = ief->getNrow();
    int ncol = ief->getNcol();
    int offset = ief->getOffset(world_rank);
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            fid = fabs(ief->getVal(i,j+1))-1;
            
            if(f2e[fid].find(i) == f2e[fid].end())
            {
                f2e[fid].insert(offset+i);
                f2e_vec[fid].push_back(offset+i);
            }
            if(e2f[i].find(fid) == e2f[i].end())
            {
                e2f[offset+i].insert(fid);
                e2f_vec[offset+i].push_back(fid);
            }
        }
    }
    
    set<int>::iterator it;
    std::map<int, std::vector<int> > e2e;
    for(int i=0;i<f2e_vec.size();i++)
    {
        if(f2e_vec[i].size()==2)
        {
            e2e[f2e_vec[i][0]].push_back(f2e_vec[i][1]);
            e2e[f2e_vec[i][1]].push_back(f2e_vec[i][0]);
        }
    }
    
    std::cout << "e2e.size() " << e2e.size() << std::endl;
    for(int i=0;i<e2e.size();i++)
    {
        std::cout << world_rank << " element = " << i << " :: ";
        for(int j=0;j<e2e[i].size();j++)
        {
            std::cout << e2e[i][j] << " ";
        }
        std::cout << std::endl;
    }
}







int ExampleUS3DPartitioningWithParVarParMetis()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    const char* fn_conn="grids/piston/conn.h5";
    ParArray<int>*   ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    int N = ien->getNglob();
    
    ParArray<int>* ien_copy = new ParArray<int>(N, ncol, comm);
    int val = 0.0;
//
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            val = ien->getVal(i,j+1)-1;
            ien_copy->setVal(i,j,val);
        }
    }

    ParVar_ParMetis* pv_parmetis = CreateParallelDataParmetis(ien_copy,comm,8);
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj      = NULL;
    idx_t *adjncy    = NULL;

    idx_t *elmwgt = new idx_t[pv_parmetis->nlocs[world_rank]];
    idx_t *adjwgt = new idx_t[pv_parmetis->nlocs[world_rank]];

    for(int i=0;i<pv_parmetis->nlocs[world_rank];i++)
    {
        elmwgt[i] = 1;
        adjwgt[i] = 1;
    }
    
    
    

//    int* part = new int[nloc];
//    
//    for(int i=0;i<nloc;i++)
//    {
//        part[i] = 0.0;
//    }
    
    int status = ParMETIS_V3_Mesh2Dual(pv_parmetis->elmdist, pv_parmetis->eptr, pv_parmetis->eind, numflag, ncommonnodes, &xadj, &adjncy, &comm);
    
//    for(int i=0;i<nloc;i++)
//    {
//        std::cout << "rank = " << world_rank << "  el = " << i+pv_parmetis->elmdist[world_rank] << " ---> ";
//        for(int j=xadj[i];j<xadj[i+1];j++)
//        {
//            std::cout << adjncy[j] << " ";
//        }
//        std::cout << std::endl;
//
//    }
    //status = ParMETIS_V3_AdaptiveRepart(pv_parmetis->elmdist, xadj, adjncy, elmwgt, adjwgt, vsize, wgtflag, numflag, ncon, nparts, tpwgts, ubvec, itr, options, &edgecut, part, &comm);
    delete[] elmwgt;
    delete[] adjwgt;
    delete ien;
    delete pv_parmetis;
    delete[] adjncy;
    delete[] xadj;
    return status;
    
}


//
LocalPartitionData* GetPartitionData(Array<double>* xcn, ParArray<int>* ien, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    
    std::map<int, int> Loc2GlobElement;
    std::map<int, int> Glob2LocElement;
    
    std::map<int, int> Loc2GlobNodes;
    std::map<int, int> Glob2LocNodes;
    
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    Array<int>* ien_loc = new Array<int>(nrow,ncol-1);
    
    int cnt = 0;
    int val = 0;
    for(int i=0;i<nrow;i++)
    {
        Loc2GlobElement[i] = ien->getOffset(rank)+i;
        Glob2LocElement[ien->getOffset(rank)+i] = i;
        
        for(int j=0;j<(ncol-1);j++)
        {
            val = ien->getVal(i,j+1)-1;

            if ( Loc2GlobNodes.find( val ) == Loc2GlobNodes.end())
            {
                Loc2GlobNodes[cnt] = val;
                Glob2LocNodes[val] = cnt;
                ien_loc->setVal(i,j,cnt);
                cnt++;
            }
            else{
                ien_loc->setVal(i,j,Glob2LocNodes[val]);
            }
        }
    }
    
    LocalPartitionData* lpd = new LocalPartitionData;
    lpd->loc2glob_el        = Loc2GlobElement;
    lpd->glob2loc_el        = Glob2LocElement;
    lpd->ien_loc            = ien_loc;
    lpd->loc2glob_vrt       = Loc2GlobNodes;
    lpd->glob2loc_vrt       = Glob2LocNodes;
    
    return lpd;
}






void MergeTest()
{
    std::vector<int> vec;
    vec.push_back(1);
    vec.push_back(3);
    vec.push_back(2);
    vec.push_back(14);
    vec.push_back(8);
    vec.push_back(2);
    vec.push_back(3);
    vec.push_back(14);
    vec.push_back(9);
    vec.push_back(15);
    vec.push_back(1);
    vec.push_back(14);
    sort(vec.begin(),vec.end());
    std::cout << "==vec1==" << std::endl;
    for(int i = 0;i<vec.size();i++)
    {
        std::cout << vec[i] << std::endl;
    }
    std::cout << "====" << std::endl;
    std::vector<int> vec2;
    vec2.push_back(2);
    vec2.push_back(6);
    vec2.push_back(5);
    vec2.push_back(7);
    vec2.push_back(9);
    vec2.push_back(4);
    vec2.push_back(31);
    vec2.push_back(4);
    sort(vec2.begin(),vec2.end());
    std::cout << "==vec2==" << std::endl;
    for(int i = 0;i<vec2.size();i++)
    {
        std::cout << vec2[i] << std::endl;
    }
    std::cout << "====" << std::endl;
    
    std::vector<int> res = merge_vec(vec,vec2);
    
    std::cout << "====" << std::endl;
    for(int i = 0;i<res.size();i++)
    {
        std::cout << res[i] << std::endl;
    }
    std::cout << "====" << std::endl;
}



void FindDuplicateTest()
{
    std::vector<int> vec;
    vec.push_back(1);
    vec.push_back(3);
    vec.push_back(2);
    vec.push_back(14);
    vec.push_back(8);
    vec.push_back(2);
    vec.push_back(3);
    vec.push_back(14);
    vec.push_back(9);
    vec.push_back(15);
    vec.push_back(1);
    vec.push_back(14);
    std::vector<int> res = FindDuplicates(vec);
    
    std::cout << "====" << std::endl;
    for(int i = 0;i<res.size();i++)
    {
        std::cout << res[i] << std::endl;
    }
    std::cout << "====" << std::endl;
}





void ParallelSortTest()
{
    int i=0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int levels = log2(size);
    
    const char* fn_conn="grids/piston/conn.h5";
    
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    
    int nrow = ief->getNrow();
    int ncol = ief->getNcol();
    
    std::vector<int> ief_copy(nrow*(ncol-1));
    int fid;
    int cnt=0;
    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            fid = fabs(ief->getVal(i,j+1))-1;
            //std::cout << cnt << " " << fid << std::endl;
            ief_copy[cnt] = fid;
            cnt++;
        }
    }
    
    int lsize = ief_copy.size();
    int* ief_arr = new int[lsize];
    
    for(i=0;i<lsize;i++)
    {
        ief_arr[i] = ief_copy[i];
    }
    
    int *glob_arr = NULL;
    int nglob = ief->getNglob();

    if (rank == 0)
    {
        glob_arr = new int[nglob*6];
    }
    
    int* sorted = mergeSort(levels, rank, ief_arr, lsize, comm, glob_arr);
    
    if(rank == 0)
    {
        std::cout << "hoi " << rank << " " << nglob*6 << std::endl;
       for(int i=0;i<nglob*6;i++)
        {
           std::cout << rank << " " << i << " " << sorted[i] << std::endl;
        }
    }
    
    
    delete ief;
    ief_copy.erase(ief_copy.begin(),ief_copy.end());
    
}


void ParallelSortTest_v2()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int levels = log2(size);
    
    const char* fn_conn="grids/piston/conn.h5";
    //const char* fn_grid="grids/piston/grid.h5";
    
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    
    int nrow = ief->getNrow();
    int ncol = ief->getNcol();
    
    int *data = new int[nrow*(ncol-1)];
    std::vector<int> ief_copy(nrow*(ncol-1));
    int fid;
    int cnt=0;
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            fid = fabs(ief->getVal(i,j+1))-1;
            ief_copy[cnt] = fid;
            cnt++;
        }
    }
    int lsize = ief_copy.size();
    int* ief_arr = new int[lsize];
    for(int i=0;i<lsize;i++)
    {
        ief_arr[i] = ief_copy[i];
    }
    //int *glob_arr;
//    if (rank == 0)
//    {
//        glob_arr = new int[ief->nglob*6];
//
//    }
    int nglob = ief->getNglob();
    std::vector<int> glob_arr(nglob*6);
    
    std::vector<int> sorted = mergeSort_vec(levels, rank, ief_copy, lsize, comm, glob_arr);
    
    if(rank == 0)
    {
        std::cout << "hoi " << rank << " " << nglob*6 << std::endl;
       for(int i=0;i<nglob*6;i++)
        {
           std::cout << "vec " << rank << " " << i << " " << sorted[i] << std::endl;
        }
    }
    
    delete[] ief_arr;
    delete[] data;
    delete ief;
    
    ief_copy.erase(ief_copy.begin(),ief_copy.end());
    
}

/*
std::map<int, std::map<int, int> > getElement2ElementTopology(Partition* P, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int* xadj                                       = P->getXadj();
    int* adjcny                                     = P->getAdjcny();
    int nElemLoc                                    = P->getNlocElem();
    std::map<int,std::vector<int> > elem2globface   = P->getElem2GlobalFace();
    
    std::set<int> elem_set = P->getElemSet();

    std::map<int,std::map<int,int> > Element2ElementTopology;
    int nloc = P->getPart()->getNloc(world_rank);
    int of = P->getPart()->getOffset(world_rank);
    //std::set<int> owned_faces;
    //std::cout << elem2globface.size() << std::endl;
    for(int i = 0;i<nloc;i++)
    {
       std::set<int> owned_faces;
//        if(world_rank == 0)
//        {
//           std::cout << " element " << i << " ";
//        }
        
       for(int n=0;n<6;n++)
       {
           if(world_rank == 0 && i<10)
           {
               std::cout << elem2globface[i+of][n] << " ";
           }
          owned_faces.insert(elem2globface[i+of][n]);
       }
        if(world_rank == 0 && i<10)
        {
            std::cout << std::endl;
            
        }
        
       int start = xadj[i];
       int end   = xadj[i+1];
       for(int j=start;j<end;j++)
       {
          int adjEl_id = adjcny[j];
	
          //if(elem_set.find(adjEl_id) != elem_set.end())
          //{
             for(int n=0;n<6;n++)
             {
                int f_id_ex = elem2globface[adjEl_id][n];
//                 if(world_rank == 0)
//                 {
//                     std::cout << n << " " << f_id_ex << std::endl;
//                 }
                if(owned_faces.find(f_id_ex)!=owned_faces.end())
                {
                    Element2ElementTopology[i+of][n]=adjEl_id;
                }
             }
          //}
           if(Element2ElementTopology.size()<=3)
           {
               std::cout << "rank = " << world_rank << " element " << i << " " << Element2ElementTopology.size() << std::endl;
           }
           
        }
        owned_faces.clear();
    }
    elem_set.erase(elem_set.begin(),elem_set.end());
    delete[] xadj;
    delete[] adjcny;
    return Element2ElementTopology;
}

Array<double>* SolveQR(double* A, int m, int n, Array<double>* b)
{
    
    int nrow = m; 
    int ncol = n; 
    
    
    int info = 0; 
    
    double tau[ncol];
    int lda = nrow;
    int lwork=1;
    double iwork;
    Array<double>* b_copy = new Array<double>(m,1);
    for(int i=0;i<m;i++)
    {
	b_copy->setVal(i,0,b->getVal(i,0));
    }
    // DGEQRF for Q*R=A, i.e., A and tau hold R and Householder reflectors

    lapack::geqrf(nrow, ncol, A, nrow, tau);
   
    lapack::ormqr('L', 'T', nrow, 1, ncol, A, nrow, tau, b_copy->data, nrow);
    
    lapack::trtrs('U', 'N', 'N', ncol, 1, A, nrow, b_copy->data, nrow);
    
    Array<double>* out = new Array<double>(ncol,1);
    for(int i = 0;i<ncol;i++)
    {    
        out->setVal(i,0,b_copy->getVal(i,0));
    }  
    delete b_copy; 
    return out; 
}

Array<double>* ComputeGradient(Partition* P, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    Array<double>* grad = new Array<double>(Nel,3);
    
    std::map<int,std::vector<int> > gE2lV       = P->getGlobElem2LocVerts();

    std::vector<Vert> locVerts                  = P->getLocalVerts();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    std::map<int,std::map<int,int> >::iterator itadj;

    double u_ip1jk = 0.0;
    double u_im1jk = 0.0;
    double u_ijp1k = 0.0;
    double u_ijm1k = 0.0;
    double u_ijkp1 = 0.0;
    double u_ijkm1 = 0.0;
    int offset = P->getPart()->getOffset(rank);
    int* xadj = P->getXadj();
    int* adjcny = P->getAdjcny();
    int nloc = P->getPart()->getNrow();

    
    for(int i = 0;i<nloc;i++)
    {
        int start = xadj[i];
        int end   = xadj[i+1];
        int nadj  = xadj[i+1]-xadj[i];
        
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        for(int q=0;q<nadj;q++)
        {
            
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        double u_ijk = P->getU0atGlobalElem(i+offset);
        std::vector<int> vijkIDs = gE2lV[i+offset];
        
        double* Pijk = new double[np*3];
        for(int k=0;k<vijkIDs.size();k++)
        {
            loc_vid     = vijkIDs[k];
            Pijk[k*3+0] = locVerts[loc_vid].x;
            Pijk[k*3+1] = locVerts[loc_vid].y;
            Pijk[k*3+2] = locVerts[loc_vid].z;
        }
        

        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        
        
        
        int tel   = 0;
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            
            double u_po = P->getU0atGlobalElem(adjEl_id);
            double* Po = new double[np*3];
            
            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                Po[k*3+0] = locVerts[loc_vid].x;
                Po[k*3+1] = locVerts[loc_vid].y;
                Po[k*3+2] = locVerts[loc_vid].z;
            }
            
            Vert* Vpo = ComputeCenterCoord(Po,8);
            
            Vrt->setVal(tel,0,Vpo->x-Vijk->x);
            Vrt->setVal(tel,1,Vpo->y-Vijk->y);
            Vrt->setVal(tel,2,Vpo->z-Vijk->z);
            
            b->setVal(tel,0,u_po-u_ijk);
            
            delete[] Po;
            delete Vpo;
            tel++;
        }
        
        double* A_cm = new double[nadj*3];
        double sum =0.0;
        for(int s=0;s<nadj;s++)
        {
            for(int j=0;j<3;j++)
            {
                A_cm[j*nadj+s] = Vrt->getVal(s,j);
            }
        }
        
        //Array<double>* x = SolveQR(A_cm,nadj,3,b);
        Array<double>* x = new Array<double>(3,1);
        x->setVal(0,0,0.0);
        x->setVal(1,0,0.0);
        x->setVal(2,0,0.0);
        grad->setVal(i,0,x->getVal(0,0));
        grad->setVal(i,1,x->getVal(1,0));
        grad->setVal(i,2,x->getVal(2,0));
        
        
        delete[] Pijk;
        vijkIDs.clear();
        delete Vrt_T;
        delete Vrt;
        delete b;
        delete[] A_cm;
        delete x;
        delete Vijk;
        
    }
   
    delete[] xadj;
    delete[] adjcny;
    gE2lV.clear();
    locVerts.clear();
    
    return grad;
}



Array<double>* ComputeHessian(Partition* P, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int nloc = P->getPart()->getNrow();
    Array<double>* hessian = new Array<double>(Nel,3);
    
    std::map<int,std::vector<int> > gE2lV       = P->getGlobElem2LocVerts();

    std::vector<Vert> locVerts                  = P->getLocalVerts();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    std::map<int,std::map<int,int> >::iterator itadj;

    int offset = P->getPart()->getOffset(rank);
    int* xadj = P->getXadj();
    int* adjcny = P->getAdjcny();
  
    for(int i = 0;i<nloc;i++)
    {
        
        int start = xadj[i];
        int end   = xadj[i+1];
        int nadj  = xadj[i+1]-xadj[i];
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
         
        for(int q=0;q<nadj;q++)
        {
            
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        double u_ijk = P->getUauxatGlobalElem(i+offset);
        
        std::vector<int> vijkIDs = gE2lV[i+offset];
        
        double* Pijk = new double[np*3];
        for(int k=0;k<vijkIDs.size();k++)
        {
            loc_vid     = vijkIDs[k];
            Pijk[k*3+0] = locVerts[loc_vid].x;
            Pijk[k*3+1] = locVerts[loc_vid].y;
            Pijk[k*3+2] = locVerts[loc_vid].z;
        }
        
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        
        delete[] Pijk;
        
        int tel   = 0;
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            
            double u_po = P->getUauxatGlobalElem(adjEl_id);
            double* Po = new double[np*3];
            
            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                Po[k*3+0] = locVerts[loc_vid].x;
                Po[k*3+1] = locVerts[loc_vid].y;
                Po[k*3+2] = locVerts[loc_vid].z;
            }
            
            Vert* Vpo = ComputeCenterCoord(Po,8);
            
            Vrt->setVal(tel,0,Vpo->x-Vijk->x);
            Vrt->setVal(tel,1,Vpo->y-Vijk->y);
            Vrt->setVal(tel,2,Vpo->z-Vijk->z);
            
       
            b->setVal(tel,0,u_po-u_ijk);
            
            delete[] Po;
            delete Vpo;
            tel++;
        }
        
        double* A_cm = new double[nadj*3];
        double sum =0.0;
        for(int s=0;s<nadj;s++)
        {
            for(int j=0;j<3;j++)
            {
                A_cm[j*nadj+s] = Vrt->getVal(s,j);
            }
        }
        
        Array<double>* x = SolveQR(A_cm,nadj,3,b);
        
        hessian->setVal(i,0,x->getVal(0,0));
        hessian->setVal(i,1,x->getVal(1,0));
        hessian->setVal(i,2,x->getVal(2,0));
        
    
        vijkIDs.clear();
        delete Vijk;
        delete[] A_cm;
        delete Vrt_T;
        delete Vrt;
        delete b;
    }
    
    //delete[] xadj;
    //delete[] adjcny;
        
    return hessian;
}

*/

/*
void QRdecomTest()
{
    Array<double>* Vrt = new Array<double>(4,3);
    int m = Vrt->getNrow();
    int n = Vrt->getNcol();
    Vrt->setVal(0,0,1.0);Vrt->setVal(0,1,1.0);Vrt->setVal(0,2,1.0);
    Vrt->setVal(1,0,1.0);Vrt->setVal(1,1,1.0);Vrt->setVal(1,2,0.0);
    Vrt->setVal(2,0,-1.0);Vrt->setVal(2,1,0.0);Vrt->setVal(2,2,1.0);
    Vrt->setVal(3,0,0.0);Vrt->setVal(3,1,0.0);Vrt->setVal(3,2,1.0);
    Array<double>* A_T = new Array<double>(3,4);

    Array<double>* A = new Array<double>(4,3);

    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            A_T->setVal(j,i,Vrt->getVal(i,j));
            A->setVal(i,j,Vrt->getVal(i,j));
        }
    }

    Array<double>* b = new Array<double>(4,1);
    b->setVal(0,0,1.0);
    b->setVal(1,0,2.0);
    b->setVal(2,0,3.0);
    b->setVal(3,0,4.0);

    Array<double>* b_copy = new Array<double>(4,1);
    b_copy->setVal(0,0,1.0);
    b_copy->setVal(1,0,2.0);
    b_copy->setVal(2,0,3.0);
    b_copy->setVal(3,0,4.0);

    double* A_cm = new double[n*m];
    double* A_rm = new double[n*m];

    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            A_rm[i*n+j] = A->getVal(i,j);
            A_cm[j*m+i] = A->getVal(i,j);
            std::cout << A->getVal(i,j) << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << m << " " << n << std::endl;
    for(int i=0;i<m*n;i++)
    {
        std::cout << A_cm[i] << " " << A->data[i] << " " << A_rm[i] << std::endl;
    }
    std::cout << "========================================" << std::endl;
    
    Array<double>* x = SolveQR(A_cm,m,n,b);
    double* ref = new double[3];
    ref[0] = -0.666667;
    ref[1] = 1.0;
    ref[2] = 7.0/3.0;
    
    for(int i=0;i<3;i++)
    {
        std::cout << "should be zero " << x->getVal(i,0)-ref[i] << std::endl;
    }
    
    Array<double>* approx = new Array<double>(4,1);
    double res;
    for(int i=0;i<4;i++)
    {
        res = 0.0;
        for(int j=0;j<3;j++)
        {
            res = res + A->getVal(i,j)*x->getVal(j,0);
        }
        approx->setVal(i,0,res);
    }
    double LQ = 0.0;
    for(int i=0;i<m;i++)
    {
        LQ = LQ+
        (approx->getVal(i,0)-b_copy->getVal(i,0))*(approx->getVal(i,0)-b_copy->getVal(i,0));
        
    }
    
    std::cout << "LQ = " << LQ << " sqrt(LQ) " << sqrt(LQ) << std::endl;
}



void QRdecomGradRecTest()
{
    Array<double>* Vrt = new Array<double>(4 ,3);
    int m = Vrt->getNrow();
    int n = Vrt->getNcol();
    
    Vrt->setVal(0,0,-1.0);  Vrt->setVal(0,1,0.0);    Vrt->setVal(0,2,0.000001);
    Vrt->setVal(1,0,0.0);   Vrt->setVal(1,1,-1.0);   Vrt->setVal(1,2,0.000001);
    Vrt->setVal(2,0,1.0);   Vrt->setVal(2,1,0.0);    Vrt->setVal(2,2,0.000001);
    Vrt->setVal(3,0,0.0);   Vrt->setVal(3,1,1.0);    Vrt->setVal(3,2,0.000001);
    
    Array<double>* b = new Array<double>(4,1);
    b->setVal(0,0,100-200.0);
    b->setVal(1,0,100-200.0);
    b->setVal(2,0,300-200.0);
    b->setVal(3,0,300-200.0);
    
    Array<double>* biff = new Array<double>(4,1);
    biff->setVal(0,0,100.0-200.0);
    biff->setVal(1,0,100.0-200.0);
    biff->setVal(2,0,300.0-200.0);
    biff->setVal(3,0,300.0-200.0);

    Array<double>* b_copy = new Array<double>(4,1);
    b_copy->setVal(0,0,100.0-200.0);
    b_copy->setVal(1,0,100.0-200.0);
    b_copy->setVal(2,0,300.0-200.0);
    b_copy->setVal(3,0,300.0-200.0);

    double* A_cm = new double[m*n];
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            A_cm[j*m+i] = Vrt->getVal(i,j);
        }
        
    }
    Array<double>* x = SolveQR(A_cm,m,n,b);
    Array<double>* r = MatMul(Vrt,x);
    double error;
     for(int i=0;i<m;i++)
       {
           error = r->getVal(i,0)-b_copy->getVal(i,0);
           std::cout << r->getVal(i,0) << " " << b_copy->getVal(i,0) << " " << r->getVal(i,0)-b_copy->getVal(i,0) << std::endl;
           if(error>0.001)
           {
               throw std::runtime_error("error :: /GradRecTest() has failed/");
           }
           
       }
}

*/

struct PartTmp {
    int* xadj;
    int* adjcny;
    Array<int>* lPart;
    Array<int>* gPart;
};

PartTmp* getPartionIDS(ParArray<int>* ien, ParallelState* pstate, MPI_Comm comm, int ElType)
{
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int Nel = ien->getNglob();
    
    PartTmp* pTmp = new PartTmp;
    pTmp->gPart = new Array<int>(Nel,1);
    
    ParallelState_Parmetis* pstate_parmetis = new ParallelState_Parmetis(ien,comm,ElType);
    
    
    int nrow = ien->getNrow();
    int nloc = nrow;
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {3};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.05};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt = NULL;

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
    int* part_arr = new int[nloc];
    real_t itr_[]    = {1.05};
    real_t *itr      = itr_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
//
    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist(),
                          pstate_parmetis->getEptr(),
                          pstate_parmetis->getEind(),
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);
//
    
//    for(int u=0;u<nloc;u++)
//    {
//        part_arr[u] = rank;
//    }

//
    ParMETIS_V3_AdaptiveRepart(pstate_parmetis->getElmdist(),
                           xadj_par, adjncy_par,
                           elmwgt, adjwgt,
                       vsize, wgtflag,
                   numflag, ncon, nparts,
                   tpwgts, ubvec, itr, options,
                   &edgecut, part_arr, &comm);
/*
    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist(),
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);
*/

    pTmp->lPart = new Array<int>(nloc,1);

    pTmp->lPart->data = part_arr;
    pTmp->xadj   = xadj_par;
    pTmp->adjcny = adjncy_par;

    MPI_Allgatherv(&pTmp->lPart->data[0],
                   nloc, MPI_INT,
                   &pTmp->gPart->data[0],
                   pstate->getNlocs(),
                   pstate->getOffsets(),
                   MPI_INT,comm);
    
    return pTmp;
}


void GradRecTest()
{
    Array<double>* Vrt = new Array<double>(4,3);

    int m = Vrt->getNrow();
    int n = Vrt->getNcol();
    
    Vrt->setVal(0,0,-1.0);  Vrt->setVal(0,1,0.0);    Vrt->setVal(0,2,1.0e-12);
    Vrt->setVal(1,0,0.0);   Vrt->setVal(1,1,-1.0);   Vrt->setVal(1,2,1.0e-12);
    Vrt->setVal(2,0,1.0);   Vrt->setVal(2,1,0.0);    Vrt->setVal(2,2,1.0e-12);
    Vrt->setVal(3,0,0.0);   Vrt->setVal(3,1,1.0);    Vrt->setVal(3,2,1.0e-12);

//    Vrt->setVal(0,0,-1.0);  Vrt->setVal(0,1,0.0);    Vrt->setVal(0,2,0.0);
//    Vrt->setVal(1,0,0.0);   Vrt->setVal(1,1,-1.0);   Vrt->setVal(1,2,0.0);
//    Vrt->setVal(2,0,1.0);   Vrt->setVal(2,1,0.0);    Vrt->setVal(2,2,0.0);
//    Vrt->setVal(3,0,0.0);   Vrt->setVal(3,1,1.0);    Vrt->setVal(3,2,0.0);
//    
    Array<double>* Vrt_T = new Array<double>(3,4);

    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            Vrt_T->setVal(j,i,Vrt->getVal(i,j));
        }
    }
//
    Array<double>* b = new Array<double>(4,1);
    b->setVal(0,0,100-200.0);
    b->setVal(1,0,100-200.0);
    b->setVal(2,0,300-200.0);
    b->setVal(3,0,300-200.0);
    
    Array<double>* R = MatMul(Vrt_T,Vrt);
    Array<double>*Rinv = MatInv(R);
    Array<double>* Rn = MatMul(Rinv,Vrt_T);
//
    Array<double>* x = MatMul(Rn,b);

    for(int i=0;i<x->getNrow();i++)
    {
        std::cout << x->getVal(i,0) << std::endl;
    }
    Array<double>* r = MatMul(Vrt,x);
    int pass = 0;
    double error;
    for(int i=0;i<m;i++)
    {
        error = r->getVal(i,0)-b->getVal(i,0);
        std::cout << r->getVal(i,0) << " " << b->getVal(i,0) << " " << r->getVal(i,0)-b->getVal(i,0) << std::endl;
        if(error>0.001)
        {
            throw std::runtime_error("error :: /GradRecTest() has failed/");
        }
        
    }
    
}




/* Function to get a local mesh from a global one, */
void get_local_mesh(int np, int ne, int nt, int *pmask, int *inv_pmask,
                    int *emask, int *tmask, int *inv_tmask,
                    double *pcoor, double *pcoor_all, int *pref, int *pref_all,
                    int *evert, int *evert_all, int *eref, int *eref_all,
                    int *tvert, int *tvert_all, int *tref, int *tref_all,
                    double *met, double *met_all,int ncomm,
                    int *ntifc, int **ifc_tria_loc, int **ifc_tria_glob,
                    int *npifc, int **ifc_nodes_loc, int **ifc_nodes_glob) {

  int k,d,icomm;

  for( k=0; k<np; k++ ) {
    inv_pmask[pmask[k]-1] = k+1;
    for( d=0; d<3; d++ )
      pcoor[3*k+d] = pcoor_all[3*(pmask[k]-1)+d];
    pref[k] = pref_all[pmask[k]-1];
  }
  for( k=0; k<ne; k++ ) {
    for( d=0; d<4; d++ )
      evert[4*k+d] = inv_pmask[evert_all[4*(emask[k]-1)+d]-1];
    eref[k] = eref_all[emask[k]-1];
  }
  for( k=0; k<nt; k++ ) {
    inv_tmask[tmask[k]-1] = k+1;
    for( d=0; d<3; d++ )
      tvert[3*k+d] = inv_pmask[tvert_all[3*(tmask[k]-1)+d]-1];
    tref[k] = tref_all[tmask[k]-1];
  }
  for( k=0; k<np; k++ ) {
    met[k] = met_all[pmask[k]-1];
  }
  for( icomm=0; icomm<ncomm; icomm++ ) {
    for( k=0; k<ntifc[icomm]; k++ ) {
      ifc_tria_loc[icomm][k] = inv_tmask[ifc_tria_glob[icomm][k]-1];
    }
    for( k=0; k<npifc[icomm]; k++ ) {
      ifc_nodes_loc[icomm][k] = inv_pmask[ifc_nodes_glob[icomm][k]-1];
    }
  }
}



struct MmgTestData{
    
    MMG5_pMesh mmgMesh;
    MMG5_pSol mmgSol;
    
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > E2N;
    std::map<int,std::vector<int> > F2N;
    
    double* vert_coor_all;
    int* vert_ref_all;
    int* tetra_vert_all;
    int* tetra_ref_all;
    int* tria_vert_all;
    int* tria_ref_all;
    double* met_all;
};


void UpdateConnectivityMmgMesh(MmgTestData* MmgTdata,MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<int> T;
    int fid = 0;
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > F2N;
    std::map<int,std::vector<int> > E2N;
    std::set<int> tri1;
    std::set<int> tri2;
    std::set<int> tri3;
    std::set<int> tri4;
    std::map<int,std::set<int> > Tr;
    std::map<std::set<int>, int> TrID;
    std::set<std::set<int> > faces;
    int e = 0;
    //std::cout << "Number of elements " << MmgTdata->mmgMesh->ne << std::endl;
    for(int i=1;i<=MmgTdata->mmgMesh->ne;i++)
    {
        e = i-1;
        
        int vid0 = MmgTdata->mmgMesh->tetra[i].v[0]-1;
        int vid1 = MmgTdata->mmgMesh->tetra[i].v[1]-1;
        int vid2 = MmgTdata->mmgMesh->tetra[i].v[2]-1;
        int vid3 = MmgTdata->mmgMesh->tetra[i].v[3]-1;
    
        T.push_back(vid0);
        T.push_back(vid1);
        T.push_back(vid2);
        T.push_back(vid3);
        E2N[e] =  T;
//        if(world_rank == 0)
//        {
//            std::cout << "E2N[e],.size " << E2N[e].size() << << std::endl;
//        }
        tri1.insert(T[0]);tri1.insert(T[1]);tri1.insert(T[2]);
        if( faces.count(tri1) != 1 )
        {
            F2N[fid].push_back(T[0]);
            F2N[fid].push_back(T[1]);
            F2N[fid].push_back(T[2]);
            faces.insert(tri1);
            Tr[fid]=tri1;
            TrID[tri1]=fid;
            E2F[e].push_back(fid);
            F2E[fid].push_back(e);
            fid++;
        }
        else
        {
            E2F[e].push_back(TrID[tri1]);
            F2E[TrID[tri1]].push_back(e);
        }
        
        
        tri2.insert(T[1]);tri2.insert(T[2]);tri2.insert(T[3]);
        if( faces.count(tri2) != 1 )
        {
            F2N[fid].push_back(T[1]);
            F2N[fid].push_back(T[2]);
            F2N[fid].push_back(T[3]);
            faces.insert(tri2);
            Tr[fid]=tri2;
            TrID[tri2]=fid;
            E2F[e].push_back(fid);
            F2E[fid].push_back(e);
            fid++;
        }
        else
        {
            E2F[e].push_back(TrID[tri2]);
            F2E[TrID[tri2]].push_back(e);
        }
        
        
        tri3.insert(T[2]);tri3.insert(T[3]);tri3.insert(T[0]);
        if( faces.count(tri3) != 1 )
        {
            F2N[fid].push_back(T[2]);
            F2N[fid].push_back(T[3]);
            F2N[fid].push_back(T[0]);
            faces.insert(tri3);
            Tr[fid]=tri3;
            TrID[tri3]=fid;
            E2F[e].push_back(fid);
            F2E[fid].push_back(e);
            fid++;
        }
        else
        {
            E2F[e].push_back(TrID[tri3]);
            F2E[TrID[tri3]].push_back(e);
        }
        
        
        
        tri4.insert(T[3]);tri4.insert(T[0]);tri4.insert(T[1]);
        if( faces.count(tri4) != 1 )
        {
            F2N[fid].push_back(T[3]);
            F2N[fid].push_back(T[0]);
            F2N[fid].push_back(T[1]);
            faces.insert(tri4);
            Tr[fid]=tri4;
            TrID[tri4]=fid;
            E2F[e].push_back(fid);
            F2E[fid].push_back(e);
            fid++;
        }
        else
        {
            E2F[e].push_back(TrID[tri4]);
            F2E[TrID[tri4]].push_back(e);
        }
        
        T.clear();
        tri1.clear();
        tri2.clear();
        tri3.clear();
        tri4.clear();
    }
    
    MmgTdata->E2N = E2N;
    MmgTdata->F2N = F2N;
    MmgTdata->E2F = E2F;
    MmgTdata->F2E = F2E;
    MmgTdata->E2N = E2N;
}




MmgTestData* GetStructuredBlockMmgMesh(int N,MMG5_pMesh& mmgMesh, MMG5_pSol& mmgSol)
{
    
    MmgTestData* MmgTdata = new MmgTestData;
    int nvertices = N*N*N;
    int Nel = (N-1)*(N-1)*(N-1);

    double dx = 1.0/(N-1);
    double dy = 1.0/(N-1);
    double dz = 1.0/(N-1);
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    std::vector<Vert> verts;
    int cnt = 0;
    
    MMG3D_Init_mesh(MMG5_ARG_start,
    MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
    MMG5_ARG_end);
    
    
//    double* vert_coor_all;
//    int vert_ref_all;
//    int* tetra_vert_all;
//    int* tetra_ref_all;
//    int* tria_vert_all;
//    int* tria_ref_all;
//
    MmgTdata->vert_coor_all = new double[nvertices*3];
    MmgTdata->vert_ref_all = new int[nvertices];

    for(int i=0;i<N;i++)
    {
        vy = 0;
        for(int j=0;j<N;j++)
        {
            vx = 0;
            for(int k=0;k<N;k++)
            {
                Vert V;
                V.x = vx;
                V.y = vy;
                V.z = vz;
                verts.push_back(V);
                
                MmgTdata->vert_coor_all[cnt*3+0] = V.x;
                MmgTdata->vert_coor_all[cnt*3+1] = V.y;
                MmgTdata->vert_coor_all[cnt*3+2] = V.z;
                MmgTdata->vert_ref_all[cnt] = 1;
                
                vx = vx+dx;
                cnt++;
            }
            vy = vy+dy;
        }
        vz = vz+dz;
    }

    std::vector<std::vector<int> > hexes;
    std::vector<std::vector<int> > tets;
    //Array<int>* E2F = new Array<int>(Nel*6,4);
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > F2N;
    std::map<int,std::vector<int> > E2N;
    std::set<std::set<int> > triangles;
    std::set<int> tri1;
    std::set<int> tri2;
    std::set<int> tri3;
    std::set<int> tri4;
    std::map<int,std::set<int> > Tr;
    std::map<std::set<int>, int> TrID;
    int fid = 0;
    int e   = 0;
    for(int i=0;i<N-1;i++)
    {
        for(int j=0;j<N-1;j++)
        {
            for(int k=0;k<N-1;k++)
            {
                std::vector<int> H(8);
                H[0] = N*i+j+N*N*(k);
                H[1] = N*i+j+1+N*N*(k);
                H[2] = N*(i+1)+j+N*N*(k);
                H[3] = N*(i+1)+j+1+N*N*(k);

                H[4] = N*i+j+N*N*(k+1);
                H[5] = N*i+j+1+N*N*(k+1);
                H[6] = N*(i+1)+j+N*N*(k+1);
                H[7] = N*(i+1)+j+1+N*N*(k+1);
                
//                for(int i=0;i<8;i++)
//                {
//                    std::cout << H[i] << " ";
//                }
//                std::cout << std::endl;
                
                std::vector<int> T1(4);
                std::vector<int> T2(4);
                std::vector<int> T3(4);
                std::vector<int> T4(4);
                std::vector<int> T5(4);
                std::vector<int> T6(4);
                
                
                //======================================================================
                T1[0] = H[0];
                T1[1] = H[4];
                T1[2] = H[1];
                T1[3] = H[6];
                E2N[e]=T1;
                //======================================================================
                tri1.insert(T1[0]);tri1.insert(T1[1]);tri1.insert(T1[2]);
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T1[0]);
                    F2N[fid].push_back(T1[1]);
                    F2N[fid].push_back(T1[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T1[1]);tri2.insert(T1[2]);tri2.insert(T1[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T1[1]);
                    F2N[fid].push_back(T1[2]);
                    F2N[fid].push_back(T1[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T1[2]);tri3.insert(T1[3]);tri3.insert(T1[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T1[2]);
                    F2N[fid].push_back(T1[3]);
                    F2N[fid].push_back(T1[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T1[3]);tri4.insert(T1[0]);tri4.insert(T1[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T1[3]);
                    F2N[fid].push_back(T1[0]);
                    F2N[fid].push_back(T1[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================
                T2[0] = H[1];
                T2[1] = H[7];
                T2[2] = H[5];
                T2[3] = H[6];
                E2N[e]=T2;
                //======================================================================
                tri1.insert(T2[0]);tri1.insert(T2[1]);tri1.insert(T2[2]);
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T2[0]);
                    F2N[fid].push_back(T2[1]);
                    F2N[fid].push_back(T2[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T2[1]);tri2.insert(T2[2]);tri2.insert(T2[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T2[1]);
                    F2N[fid].push_back(T2[2]);
                    F2N[fid].push_back(T2[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T2[2]);tri3.insert(T2[3]);tri3.insert(T2[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T2[2]);
                    F2N[fid].push_back(T2[3]);
                    F2N[fid].push_back(T2[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T2[3]);tri4.insert(T2[0]);tri4.insert(T2[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T2[3]);
                    F2N[fid].push_back(T2[0]);
                    F2N[fid].push_back(T2[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================


                T3[0] = H[1];
                T3[1] = H[2];
                T3[2] = H[3];
                T3[3] = H[6];
                E2N[e]=T3;
                //======================================================================
                tri1.insert(T3[0]);tri1.insert(T3[1]);tri1.insert(T3[2]);
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T3[0]);
                    F2N[fid].push_back(T3[1]);
                    F2N[fid].push_back(T3[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T3[1]);tri2.insert(T3[2]);tri2.insert(T3[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T3[1]);
                    F2N[fid].push_back(T3[2]);
                    F2N[fid].push_back(T3[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T3[2]);tri3.insert(T3[3]);tri3.insert(T3[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T3[2]);
                    F2N[fid].push_back(T3[3]);
                    F2N[fid].push_back(T3[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T3[3]);tri4.insert(T3[0]);tri4.insert(T3[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T3[3]);
                    F2N[fid].push_back(T3[0]);
                    F2N[fid].push_back(T3[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================

                T4[0] = H[0];
                T4[1] = H[2];
                T4[2] = H[1];
                T4[3] = H[6];
                E2N[e]=T4;
                //======================================================================
                tri1.insert(T4[0]);tri1.insert(T4[1]);tri1.insert(T4[2]);
                
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T4[0]);
                    F2N[fid].push_back(T4[1]);
                    F2N[fid].push_back(T4[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T4[1]);tri2.insert(T4[2]);tri2.insert(T4[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T4[1]);
                    F2N[fid].push_back(T4[2]);
                    F2N[fid].push_back(T4[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T4[2]);tri3.insert(T4[3]);tri3.insert(T4[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T4[4]);
                    F2N[fid].push_back(T4[3]);
                    F2N[fid].push_back(T4[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T4[3]);tri4.insert(T4[0]);tri4.insert(T4[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T4[3]);
                    F2N[fid].push_back(T4[0]);
                    F2N[fid].push_back(T4[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================

                T5[0] = H[3];
                T5[1] = H[7];
                T5[2] = H[1];
                T5[3] = H[6];
                E2N[e]=T5;
                //======================================================================
                tri1.insert(T5[0]);tri1.insert(T5[1]);tri1.insert(T5[2]);
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T5[0]);
                    F2N[fid].push_back(T5[1]);
                    F2N[fid].push_back(T5[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T5[1]);tri2.insert(T5[2]);tri2.insert(T5[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T5[1]);
                    F2N[fid].push_back(T5[2]);
                    F2N[fid].push_back(T5[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T5[2]);tri3.insert(T5[3]);tri3.insert(T5[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T5[2]);
                    F2N[fid].push_back(T5[3]);
                    F2N[fid].push_back(T5[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T5[3]);tri4.insert(T5[0]);tri4.insert(T5[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T5[3]);
                    F2N[fid].push_back(T5[0]);
                    F2N[fid].push_back(T5[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================

                T6[0] = H[5];
                T6[1] = H[4];
                T6[2] = H[1];
                T6[3] = H[6];
                E2N[e]=T6;
                //======================================================================
                tri1.insert(T6[0]);tri1.insert(T6[1]);tri1.insert(T6[2]);
                if( triangles.count(tri1) != 1 )
                {
                    F2N[fid].push_back(T6[0]);
                    F2N[fid].push_back(T6[1]);
                    F2N[fid].push_back(T6[2]);
                    triangles.insert(tri1);
                    Tr[fid]=tri1;
                    TrID[tri1]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri1]);
                    F2E[TrID[tri1]].push_back(e);
                }
                
                
                tri2.insert(T6[1]);tri2.insert(T6[2]);tri2.insert(T6[3]);
                if( triangles.count(tri2) != 1 )
                {
                    F2N[fid].push_back(T6[1]);
                    F2N[fid].push_back(T6[2]);
                    F2N[fid].push_back(T6[3]);
                    triangles.insert(tri2);
                    Tr[fid]=tri2;
                    TrID[tri2]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri2]);
                    F2E[TrID[tri2]].push_back(e);
                }
                
                
                tri3.insert(T6[2]);tri3.insert(T6[3]);tri3.insert(T6[0]);
                if( triangles.count(tri3) != 1 )
                {
                    F2N[fid].push_back(T6[2]);
                    F2N[fid].push_back(T6[3]);
                    F2N[fid].push_back(T6[0]);
                    triangles.insert(tri3);
                    Tr[fid]=tri3;
                    TrID[tri3]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri3]);
                    F2E[TrID[tri3]].push_back(e);
                }
                
                
                
                tri4.insert(T6[3]);tri4.insert(T6[0]);tri4.insert(T6[1]);
                if( triangles.count(tri4) != 1 )
                {
                    F2N[fid].push_back(T6[3]);
                    F2N[fid].push_back(T6[0]);
                    F2N[fid].push_back(T6[1]);
                    triangles.insert(tri4);
                    Tr[fid]=tri4;
                    TrID[tri4]=fid;
                    E2F[e].push_back(fid);
                    F2E[fid].push_back(e);
                    fid++;
                }
                else
                {
                    E2F[e].push_back(TrID[tri4]);
                    F2E[TrID[tri4]].push_back(e);
                }
                
                tri1.clear();
                tri2.clear();
                tri3.clear();
                tri4.clear();
                e++;
                //======================================================================

                tets.push_back(T1);
                tets.push_back(T2);
                tets.push_back(T3);
                tets.push_back(T4);
                tets.push_back(T5);
                tets.push_back(T6);
                
                
                
                hexes.push_back(H);
            }
        }
    }

    //std::cout << "sizng " << triangles.size() << " " << F2N.size() << std::endl;
    if ( MMG3D_Set_meshSize(mmgMesh,nvertices,tets.size(),0,F2N.size(),0,0) != 1 )  exit(EXIT_FAILURE);
    MmgTdata->met_all = new double[nvertices];
    for(int i=0;i<nvertices;i++)
    {
        mmgMesh->point[i+1].c[0]  = verts[i].x;  mmgMesh->point[i+1].c[1]  = verts[i].y; mmgMesh->point[i+1].c[2]  = verts[i].z; mmgMesh->point[i+1].ref  = 0;
        MmgTdata->met_all[i]=0.01;

    }
    
    std::cout << "Test mesh stats:" << std::endl;
    std::cout << " Number of tetrahedra = " << tets.size() << std::endl;
    std::cout << " Number of triangles = " << triangles.size() << std::endl;
    std::cout << " Number of triangles map = " << Tr.size() << std::endl;
    std::cout << " size F2E = " << F2E.size() << std::endl;
    std::cout << " Number of vertices = " << nvertices << std::endl;
    
    MmgTdata->tetra_vert_all = new int[tets.size()*4];
    MmgTdata->tetra_ref_all = new int[tets.size()];
    
    for(int i=1;i<=tets.size();i++)
    {
        
        mmgMesh->tetra[i].v[0] = tets[i-1][0]+1;
        mmgMesh->tetra[i].v[1] = tets[i-1][1]+1;
        mmgMesh->tetra[i].v[2] = tets[i-1][2]+1;
        mmgMesh->tetra[i].v[3] = tets[i-1][3]+1;
        
        MmgTdata->tetra_vert_all[(i-1)*4+0] = tets[i-1][0]+1;
        MmgTdata->tetra_vert_all[(i-1)*4+1] = tets[i-1][1]+1;
        MmgTdata->tetra_vert_all[(i-1)*4+2] = tets[i-1][2]+1;
        MmgTdata->tetra_vert_all[(i-1)*4+3] = tets[i-1][3]+1;
        MmgTdata->tetra_ref_all[(i-1)] = 1;
       // mmgMesh->tetra[i].ref  = 1;

    }
    
    MmgTdata->tria_vert_all = new int[F2N.size()*3];
    MmgTdata->tria_ref_all = new int[F2N.size()];

    std::map<int,std::vector<int> >::iterator f2nit;
    int q=0;
    for(f2nit=F2N.begin();f2nit!=F2N.end();f2nit++)
    {
        for(int j=0;j<3;j++)
        {
           mmgMesh->tria[q+1].v[j] = f2nit->second[j]+1;
           MmgTdata->tria_vert_all[q*3+j] = f2nit->second[j]+1;
        }
        
        MmgTdata->tria_ref_all[q]   = 0;
        mmgMesh->tria[q+1].ref      = 0;
        q++;

    }

    MMG3D_saveMesh(mmgMesh,"meshname.mesh");
    
    MmgTdata->mmgMesh   = mmgMesh;
    MmgTdata->E2N       = E2N;
    MmgTdata->F2N       = F2N;
    MmgTdata->E2F       = E2F;
    MmgTdata->F2E       = F2E;
    MmgTdata->E2N       = E2N;
    return MmgTdata;
}

struct TestMesh{
    
    int nv;
    int ne;
    int nt;
    
    Array<int>* ien_glob;
    ParArray<int>* ien_loc;
    ParallelState* pstate;
    
    double* vert_coor_all;
    int* vert_ref_all;
    int* tetra_vert_all;
    int* tetra_ref_all;
    int* tria_vert_all;
    int* tria_ref_all;
    double* met_all;
    
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > E2N;
    std::map<int,std::vector<int> > F2N;
    
};

TestMesh* GetStructureBlockMesh(int Nx, int Ny, int Nz, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int nvertices = Nx*Ny*Nz;
    int Nel = (Nx-1)*(Ny-1)*(Nz-1);

    double dx = 1.0/(Nx-1);
    double dy = 1.0/(Ny-1);
    double dz = 1.0/(Nz-1);
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    std::vector<Vert> verts;
    double* vert_coor_all   = new double[nvertices*3];
    int* vert_ref_all       = new int[nvertices];

    TestMesh* tmesh = new TestMesh;
    tmesh->nv = Nx*Ny*Nz;
    tmesh->ne = (Nx-1)*(Ny-1)*(Nz-1);
    
    int cnt = 0;
    for(int i=0;i<Nz;i++)
    {
       vy = 0;
       for(int j=0;j<Ny;j++)
       {
           vx = 0;
           for(int k=0;k<Nx;k++)
           {
               Vert V;
               V.x = vx;
               V.y = vy;
               V.z = vz;
               verts.push_back(V);
               
               vert_coor_all[cnt*3+0] = V.x;
               vert_coor_all[cnt*3+1] = V.y;
               vert_coor_all[cnt*3+2] = V.z;
               vert_ref_all[cnt] = 1;
               //std::cout << V.x << ", " << V.y << ", " << V.z << std::endl;
               
               vx = vx+dx;
               cnt++;
           }
           vy = vy+dy;
       }
       vz = vz+dz;
    }
        
    tmesh->vert_coor_all = vert_coor_all;
    tmesh->vert_ref_all  = vert_ref_all;
    
    
    int e   = 0;
    std::vector<std::vector<int> > tets;
    std::vector<std::vector<int> > hexes;
    for(int k=0;k<(Nz-1);k++)
    {
        for(int j=0;j<(Ny-1);j++)
        {
            for(int i=0;i<(Nx-1);i++)
            {
                std::vector<int> H(8);

                H[0] = j*(Nx)+i+(Nx)*(Ny)*(k);
                H[1] = j*(Nx)+i+1+(Nx)*(Ny)*(k);
                H[2] = (Nx)*(j+1)+i+(Nx)*(Ny)*(k);
                H[3] = (Nx)*(j+1)+i+1+(Nx)*(Ny)*(k);

                H[4] = (Nx)*j+i+(Nx)*(Ny)*(k+1);
                H[5] = (Nx)*j+i+1+(Nx)*(Ny)*(k+1);
                H[6] = (Nx)*(j+1)+i+(Nx)*(Ny)*(k+1);
                H[7] = (Nx)*(j+1)+i+1+(Nx)*(Ny)*(k+1);
                std::vector<int> T1(4);
                std::vector<int> T2(4);
                std::vector<int> T3(4);
                std::vector<int> T4(4);
                std::vector<int> T5(4);
                std::vector<int> T6(4);
                //std::cout << "Element " <<  H[0] << " " << H[1] << " "<< H[2] << " " << H[3] << " " << H[4] << " " << H[5] << " " << H[6] << " " << H[7] << std::endl;
                //========================================================
                //========================================================
                T1[0] = H[0]+1;
                T1[1] = H[4]+1;
                T1[2] = H[1]+1;
                T1[3] = H[6]+1;
                //========================================================
                //========================================================
                T2[0] = H[1]+1;
                T2[1] = H[7]+1;
                T2[2] = H[5]+1;
                T2[3] = H[6]+1;
                //========================================================
                //========================================================
                T3[0] = H[1]+1;
                T3[1] = H[2]+1;
                T3[2] = H[3]+1;
                T3[3] = H[6]+1;
                //========================================================
                //========================================================
                T4[0] = H[0]+1;
                T4[1] = H[2]+1;
                T4[2] = H[1]+1;
                T4[3] = H[6]+1;
                //========================================================
                //========================================================
                T5[0] = H[3]+1;
                T5[1] = H[7]+1;
                T5[2] = H[1]+1;
                T5[3] = H[6]+1;
                //========================================================
                //========================================================
                T6[0] = H[5]+1;
                T6[1] = H[4]+1;
                T6[2] = H[1]+1;
                T6[3] = H[6]+1;
                //========================================================
                //========================================================

                tets.push_back(T1);
                tets.push_back(T2);
                tets.push_back(T3);
                tets.push_back(T4);
                tets.push_back(T5);
                tets.push_back(T6);

                T1.clear();
                T2.clear();
                T3.clear();
                T4.clear();
                T5.clear();
                T6.clear();

                hexes.push_back(H);
                H.clear();
            }
        }
    }
//
    tmesh->ne = tets.size();
    int* tetra_vert_all = new int[tmesh->ne*4];
    int* tetra_ref_all = new int[tmesh->ne];

    std::set<int> tri;
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > E2N;
    std::map<int,std::vector<int> > F2N;
    Array<int>* ien_glob = new Array<int>(tmesh->ne,4);
    ParArray<int>* ien_loc = new ParArray<int>(tmesh->ne,4,comm);
    ParallelState* pstate  = new ParallelState(ien_loc->getNglob(),comm);
    //
    int offset = pstate->getOffset(world_rank);
    std::map<int, std::set<int> > fid2t;
    std::map<std::set<int>,int> t2fid;
//
    int fid = 1;
    for(int i=0;i<tmesh->ne;i++)
    {
        for(int j=0;j<4;j++)
        {
            ien_glob->setVal(i,j,tets[i][j]);
            tetra_vert_all[i*4+j] = tets[i][j];
            
            E2N[i+1].push_back(tets[i][j]);
        }
        
        tetra_ref_all[i] = 1;

        tri.insert(ien_glob->getVal(i,0));
        tri.insert(ien_glob->getVal(i,1));
        tri.insert(ien_glob->getVal(i,2));

        if(t2fid.find(tri)==t2fid.end())
        {
            t2fid[tri] = fid;
            fid2t[fid] = tri;
            fid++;
        }
        tri.clear();


        tri.insert(ien_glob->getVal(i,1));
        tri.insert(ien_glob->getVal(i,2));
        tri.insert(ien_glob->getVal(i,3));
        if(t2fid.find(tri)==t2fid.end())
        {
            t2fid[tri] = fid;
            fid2t[fid] = tri;
            fid++;
        }
        tri.clear();

        tri.insert(ien_glob->getVal(i,2));
        tri.insert(ien_glob->getVal(i,3));
        tri.insert(ien_glob->getVal(i,0));
        if(t2fid.find(tri)==t2fid.end())
        {
            t2fid[tri] = fid;
            fid2t[fid] = tri;
            fid++;
        }
        tri.clear();

        tri.insert(ien_glob->getVal(i,3));
        tri.insert(ien_glob->getVal(i,0));
        tri.insert(ien_glob->getVal(i,1));
        if(t2fid.find(tri)==t2fid.end())
        {
            t2fid[tri] = fid;
            fid2t[fid] = tri;
            fid++;
        }
        tri.clear();
    }

    std::map<int, std::set<int> >::iterator fit;
    std::map<int,int> loc2glob_face;
    int i=0;
    for(fit=fid2t.begin();fit!=fid2t.end();fit++)
    {
        std::set<int>::iterator fis;
        int j=0;
//        if(world_rank == 0)
//        {
//            std::cout << "face = " << fit->first << " -> ";
//
//        }
        for(fis=fit->second.begin();fis!=fit->second.end();fis++)
        {
            F2N[fit->first].push_back(*fis);
//            if(world_rank == 0)
//            {
//                std::cout << *fis << " ";
//            }
        }
//        if(world_rank == 0)
//        {
//            std::cout << std::endl;
//        }
        loc2glob_face[i]=fit->first;
        i++;
    }
    tmesh->ne = ien_glob->getNrow();
    tmesh->nt = F2N.size();
    //std::cout << "number of tris = " <<  tmesh->ne << " " << tmesh->nt << " " << tmesh->nv << std::endl;
    int* tria_vert_all = new int[tmesh->nt*3];
    int* tria_ref_all = new int[tmesh->nt];

//    for(int i=0;i<tmesh->nt;i++)
//    {
//        for(int j=0;j<3;j++)
//        {
//            tria_vert_all[i*3+j] = F2N[loc2glob_face[i]][j];
//        }
//
//        tria_ref_all[i] = 0;
//
//    }
    i=0;
    for(fit=fid2t.begin();fit!=fid2t.end();fit++)
    {
        std::set<int>::iterator fis;
        int j=0;
        for(fis=fit->second.begin();fis!=fit->second.end();fis++)
        {
            tria_vert_all[i*3+j] = *fis;
            j++;
        }
        i++;
    }
    std::map<std::set<int>,int>::iterator itts;

    for(int i=0;i<ien_glob->getNrow();i++)
    {
        tri.insert(ien_glob->getVal(i,0));
        tri.insert(ien_glob->getVal(i,1));
        tri.insert(ien_glob->getVal(i,2));
        E2F[i+1].push_back(t2fid[tri]);
        F2E[t2fid[tri]].push_back(i+1);
        tri.clear();

        tri.insert(ien_glob->getVal(i,1));
        tri.insert(ien_glob->getVal(i,2));
        tri.insert(ien_glob->getVal(i,3));
        E2F[i+1].push_back(t2fid[tri]);
        F2E[t2fid[tri]].push_back(i+1);
        tri.clear();


        tri.insert(ien_glob->getVal(i,2));
        tri.insert(ien_glob->getVal(i,3));
        tri.insert(ien_glob->getVal(i,0));
        E2F[i+1].push_back(t2fid[tri]);
        F2E[t2fid[tri]].push_back(i+1);
        tri.clear();

        tri.insert(ien_glob->getVal(i,3));
        tri.insert(ien_glob->getVal(i,0));
        tri.insert(ien_glob->getVal(i,1));
        E2F[i+1].push_back(t2fid[tri]);
        F2E[t2fid[tri]].push_back(i+1);
        tri.clear();
    }

    /** d) give solutions values and positions */

    for(int i=0;i<ien_loc->getNrow();i++)
    {
        for(int j=0;j<4;j++)
        {
            ien_loc->setVal(i,j,tetra_vert_all[(i+offset)*4+j]);
        }
    }

    double* met_all = new double[tmesh->nv];
    for(int i=0;i<tmesh->nv;i++)
    {
        met_all[i] = 1.0;
    }

    tmesh->ien_loc = ien_loc;
    tmesh->ien_glob = ien_glob;
    tmesh->pstate = pstate;

    tmesh->vert_coor_all=vert_coor_all;
    tmesh->vert_ref_all = vert_ref_all;
    tmesh->tetra_vert_all = tetra_vert_all;
    tmesh->tetra_ref_all = tetra_ref_all;
    tmesh->tria_vert_all = tria_vert_all;
    tmesh->tria_ref_all = tria_ref_all;
    tmesh->met_all = met_all;

    tmesh->E2F = E2F;
    tmesh->F2E = F2E;
    tmesh->F2N = F2N;
    tmesh->E2N = E2N;
    
//     TestMesh Structure has the following members:
//    int nv;
//    int ne;
//    int nt;
//
//    Array<int>* ien_glob;
//    ParArray<int>* ien_loc;
//    ParallelState* pstate;
//
//    double* vert_coor_all;
//    int* vert_ref_all;
//    int* tetra_vert_all;
//    int* tetra_ref_all;
//    int* tria_vert_all;
//    int* tria_ref_all;
//    double* met_all;
//
//    std::map<int,std::vector<int> > E2F;
//    std::map<int,std::vector<int> > F2E;
//    std::map<int,std::vector<int> > E2N;
//    std::map<int,std::vector<int> > F2N;
    
    return tmesh;
}

TestMesh* LoadTestMesh(MPI_Comm comm)
{
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    TestMesh* tmesh;
    
    tmesh->nv = 12;
    tmesh->ne = 12;
    tmesh->nt = 34;

/** 2) Build a global mesh in MMG5 format */
/** a) give the vertices (12 vertices with 3 coor = array of size 36) */
    double vert_coor_all[36] = { 0.0, 0.0, 0.0,
                                   0.5, 0.0, 0.0,
                                   0.5, 0.0, 1.0,
                                   0.0, 0.0, 1.0,
                                   0.0, 1.0, 0.0,
                                   0.5, 1.0, 0.0,
                                   0.5, 1.0, 1.0,
                                   0.0, 1.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   1.0, 1.0, 0.0,
                                   1.0, 0.0, 1.0,
                                   1.0, 1.0, 1.0};
    
   double* vert_coor_all_s = new double[36];
   for(int i=0;i<36;i++)
   {
       vert_coor_all_s[i] = vert_coor_all[i];
   }
    
    int  vert_ref_all[12] = {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  };
    int* vert_ref_all_s = new int[12];
    for(int i=0;i<12;i++)
    {
        vert_ref_all_s[i] = vert_ref_all[i];
    }
    /** b) give the tetrahedras (12 tetra with 4 vertices = array of size 48) */
    int tetra_vert_all[48] = { 1,  4,  2,  8,
                               8,  3,  2,  7,
                               5,  2,  6,  8,
                               5,  8,  1,  2,
                               7,  2,  8,  6,
                               2,  4,  3,  8,
                               9,  2,  3,  7,
                               7, 11,  9, 12,
                               6,  9, 10,  7,
                               6,  7,  2,  9,
                               12, 9,  7, 10,
                               9,  3, 11,  7  };
    
    int* tetra_vert_all_s = new int[48];
    for(int i=0;i<48;i++)
    {
        tetra_vert_all_s[i] = tetra_vert_all[i];
    }
    
    int tetra_ref_all[12] = {1  ,1  ,1  ,1  ,1  ,1  ,2  ,2  ,2  ,2  ,2  ,2  };
    int* tetra_ref_all_s = new int[12];
    for(int i=0;i<12;i++)
    {
        tetra_ref_all_s[i] = tetra_ref_all[i];
    }
    
    
    
    
    tmesh->ien_glob    = new Array<int>(tmesh->ne,4);
    tmesh->ien_loc  = new ParArray<int>(tmesh->ne,4,comm);
    ParallelState* pstate   = new ParallelState(tmesh->ien_loc->getNglob(),comm);
//
    int offset = pstate->getOffset(world_rank);
    std::map<int, std::set<int> > fid2t;
    std::map<std::set<int>,int> t2fid;
    std::set<std::set<int> > trian;
    std::set<int> tri;
    std::map<int,std::vector<int> > E2F;
    std::map<int,std::vector<int> > F2E;
    std::map<int,std::vector<int> > E2N;
    std::map<int,std::vector<int> > F2N;

    std::set<int> tri2;
    int fid=1;
    for(int i=0;i<tmesh->ien_glob->getNrow();i++)
    {
        for(int j=0;j<4;j++)
        {
            tmesh->ien_glob->setVal(i,j,tetra_vert_all[i*4+j]);
            E2N[i+1].push_back(tetra_vert_all[i*4+j]);
        }

        tri2.insert(tmesh->ien_glob->getVal(i,0));
        tri2.insert(tmesh->ien_glob->getVal(i,1));
        tri2.insert(tmesh->ien_glob->getVal(i,2));
        if(t2fid.find(tri2)==t2fid.end())
        {
            t2fid[tri2] = fid;
            fid2t[fid] = tri2;
            fid++;
        }
        tri2.clear();


        tri2.insert(tmesh->ien_glob->getVal(i,1));
        tri2.insert(tmesh->ien_glob->getVal(i,2));
        tri2.insert(tmesh->ien_glob->getVal(i,3));
        if(t2fid.find(tri2)==t2fid.end())
        {
            t2fid[tri2] = fid;
            fid2t[fid] = tri2;
            fid++;
        }
        tri2.clear();

        tri2.insert(tmesh->ien_glob->getVal(i,2));
        tri2.insert(tmesh->ien_glob->getVal(i,3));
        tri2.insert(tmesh->ien_glob->getVal(i,0));
        if(t2fid.find(tri2)==t2fid.end())
        {
            t2fid[tri2] = fid;
            fid2t[fid] = tri2;
            fid++;
        }
        tri2.clear();

        tri2.insert(tmesh->ien_glob->getVal(i,3));
        tri2.insert(tmesh->ien_glob->getVal(i,0));
        tri2.insert(tmesh->ien_glob->getVal(i,1));
        if(t2fid.find(tri2)==t2fid.end())
        {
            t2fid[tri2] = fid;
            fid2t[fid] = tri2;
            fid++;
        }
        tri2.clear();
    }

    std::map<int, std::set<int> >::iterator fit;
    std::map<int,int> loc2glob_face;
    int i=0;
    for(fit=fid2t.begin();fit!=fid2t.end();fit++)
    {
        std::set<int>::iterator fis;
        int j=0;
        for(fis=fit->second.begin();fis!=fit->second.end();fis++)
        {
            F2N[fit->first].push_back(*fis);
        }

        loc2glob_face[i]=fit->first;
        i++;
    }
    int nume = tmesh->ien_glob->getNrow();
    int numt = F2N.size();
    int* tria_vert_all = new int[numt*3];
    int* tria_ref_all = new int[numt];

    for(int i=0;i<numt;i++)
    {
        for(int j=0;j<3;j++)
        {
            tria_vert_all[i*3+j] = F2N[loc2glob_face[i]][j];
        }

        tria_ref_all[i] = 0;

    }

    std::map<std::set<int>,int>::iterator itts;

    for(int i=0;i<tmesh->ien_glob->getNrow();i++)
    {

        tri2.insert(tmesh->ien_glob->getVal(i,0));
        tri2.insert(tmesh->ien_glob->getVal(i,1));
        tri2.insert(tmesh->ien_glob->getVal(i,2));
        E2F[i+1].push_back(t2fid[tri2]);
        F2E[t2fid[tri2]].push_back(i+1);
        tri2.clear();

        tri2.insert(tmesh->ien_glob->getVal(i,1));
        tri2.insert(tmesh->ien_glob->getVal(i,2));
        tri2.insert(tmesh->ien_glob->getVal(i,3));
        E2F[i+1].push_back(t2fid[tri2]);
        F2E[t2fid[tri2]].push_back(i+1);
        tri2.clear();


        tri2.insert(tmesh->ien_glob->getVal(i,2));
        tri2.insert(tmesh->ien_glob->getVal(i,3));
        tri2.insert(tmesh->ien_glob->getVal(i,0));
        E2F[i+1].push_back(t2fid[tri2]);
        F2E[t2fid[tri2]].push_back(i+1);
        tri2.clear();

        tri2.insert(tmesh->ien_glob->getVal(i,3));
        tri2.insert(tmesh->ien_glob->getVal(i,0));
        tri2.insert(tmesh->ien_glob->getVal(i,1));
        E2F[i+1].push_back(t2fid[tri2]);
        F2E[t2fid[tri2]].push_back(i+1);
        tri2.clear();
    }

    /** d) give solutions values and positions */

    for(int i=0;i<tmesh->ien_loc->getNrow();i++)
    {
        for(int j=0;j<4;j++)
        {
            tmesh->ien_loc->setVal(i,j,tetra_vert_all[(i+offset)*4+j]);
        }
    }

    double* met_all_s = new double[tmesh->nv];
    for(int i=0;i<tmesh->nv;i++)
    {
        met_all_s[i] = 0.01;
    }
    
    tmesh->E2F = E2F;
    tmesh->F2E = F2E;
    tmesh->F2N = F2N;
    tmesh->E2N = E2N;
    tmesh->pstate = pstate;
    tmesh->vert_coor_all=vert_coor_all_s;
    tmesh->vert_ref_all = vert_ref_all_s;
    tmesh->tetra_vert_all = tetra_vert_all_s;
    tmesh->tetra_ref_all = tetra_ref_all_s;
    tmesh->tria_vert_all = tria_vert_all;
    tmesh->tria_ref_all = tria_ref_all;
    tmesh->met_all = met_all_s;
    
    
        // TestMesh Structure has the following members:
    //    int nv;
    //    int ne;
    //    int nt;
    //
    //    Array<int>* ien_glob;
    //    ParArray<int>* ien_loc;
    //    ParallelState* pstate;
    //
    //    double* vert_coor_all;
    //    int* vert_ref_all;
    //    int* tetra_vert_all;
    //    int* tetra_ref_all;
    //    int* tria_vert_all;
    //    int* tria_ref_all;
    //    double* met_all;
    //
    //    std::map<int,std::vector<int> > E2F;
    //    std::map<int,std::vector<int> > F2E;
    //    std::map<int,std::vector<int> > E2N;
    //    std::map<int,std::vector<int> > F2N;
        
    
    return tmesh;
    
}


int main(int argc, char** argv) {
    
    MPI_Init(NULL, NULL);
    FILE*   inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    //
    int nr = 0;
    int icomm,pos;
    int niter = 3;
    int API_mode = 0;
    int opt = 1;
    int ier;
    int ierlib;
    TestMesh* tmesh2 = LoadTestMesh(comm);
    TestMesh* tmesh = GetStructureBlockMesh(10,10,10,comm);
    PMMG_pParMesh parmesh = NULL;

    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);

    PartTmp* pTmp = getPartionIDS(tmesh->ien_loc, tmesh->pstate, comm, 4);

    std::map<int, std::set<int> > proc2face;
    std::map<int, std::set<int> > proc2nodes;
    std::map<int, int> face2proc;
    std::vector<int> partfaces;
    std::set<int> tetra_mask_set;
    std::set<int> tria_mask_set;
    std::set<int> vert_mask_set;
    int fadj;
    int offset=tmesh->pstate->getOffset(world_rank);
    for(int i=0;i<pTmp->lPart->getNrow();i++)
    {
        tetra_mask_set.insert(offset+i+1);
//
        std::set<int> owned_faces;
        for(int n=0;n<4;n++)
        {
            owned_faces.insert(tmesh->E2F[offset+i+1][n]);
            if(tmesh->F2E[tmesh->E2F[offset+i+1][n]].size()==1)
            {
                tria_mask_set.insert(tmesh->E2F[offset+i+1][n]);
            }
        }

        for(int n=0;n<4;n++)
        {
            vert_mask_set.insert(tmesh->E2N[offset+i+1][n]);
        }
//
        int start       = pTmp->xadj[i];
        int end         = pTmp->xadj[i+1];

        for(int j=start;j<end;j++)
        {
            int adjEl_id = pTmp->adjcny[j];
            int p_id = pTmp->gPart->getVal(adjEl_id,0);

            if(p_id!=world_rank) // find the adjacent element that is not on this partition.
            {
                for(int s=0;s<4;s++)
                {
                    fadj = tmesh->E2F[adjEl_id+1][s];

                    if(owned_faces.find(fadj)!=owned_faces.end())
                    {
                        partfaces.push_back(fadj);
                        proc2face[p_id].insert(fadj);
                        tria_mask_set.insert(fadj);
//
                        face2proc[fadj]=p_id;

                        for(int g=0;g<tmesh->F2N[fadj].size();g++)
                        {
                            proc2nodes[p_id].insert(tmesh->F2N[fadj][g]);
                        }
                    }
                }
            }
        }
        owned_faces.clear();
    }

    std::set<int>::iterator mask;
    int nv = vert_mask_set.size();
    int nt = tria_mask_set.size();
    int ne = tetra_mask_set.size();
int **faceNodes;
    int* vert_mask = new int[nv];
    int v=0;
    int rp = 0;
    for(mask=vert_mask_set.begin();mask!=vert_mask_set.end();mask++)
    {
        vert_mask[v] = *mask;

//        if(world_rank == rp)
//        {
//            std::cout << "vert mask " << *mask << std::endl;
//        }

        v++;
    }
//    if(world_rank == rp)
//    {
//        std::cout << std::endl;
//    }
    int* tetra_mask = new int[ne];
    v=0;
    for(mask=tetra_mask_set.begin();mask!=tetra_mask_set.end();mask++)
    {
        tetra_mask[v] = *mask;
//        if(world_rank == rp)
//        {
//            std::cout << "tetra_mask " << *mask << std::endl;
//        }
        v++;
    }
//    if(world_rank == rp)
//    {
//        std::cout << std::endl;
//    }
    int* tria_mask = new int[nt];
    v=0;
    for(mask=tria_mask_set.begin();mask!=tria_mask_set.end();mask++)
    {
        tria_mask[v] = *mask;
//        if(world_rank == rp)
//        {
//            std::cout << "tri_mask " << *mask << std::endl;
//        }
        v++;
    }
//    if(world_rank == rp)
//    {
//        std::cout << std::endl;
//    }

    std::cout << "mask values = " << nt << " " << ne << " " << nv << std::endl;
    int* inv_vert_mask = new int[tmesh->ne];

    int* inv_tria_mask = new int[tmesh->nt];

    double* vert_coor = new double[nv*3];
    int* vert_ref = new int[nv];
    double* met = new double[nv];
    int* tetra_vert = new int[ne*4];
    int* tetra_ref = new int[ne];

    int* tria_vert = new int[nt*3];
    int* tria_ref = new int[nt];


    int ncomm = proc2face.size();
    int* color_node = new int[ncomm];
    int* color_face = new int[ncomm];
    int* ntifc = new int[ncomm];
    int* npifc = new int[ncomm];
    std::map<int,std::set<int> >::iterator its;
    int t = 0;
    int u = 0;

    int** ifc_tria_glob = (int **) malloc(ncomm*sizeof(int *));
    int** ifc_tria_loc  = (int **) malloc(ncomm*sizeof(int *));

    //int *ifc_tria_glob[ncomm];
    //int *ifc_tria_loc[ncomm];
    for(its=proc2face.begin();its!=proc2face.end();its++)
    {
        color_face[t]=its->first;
        ntifc[t] = its->second.size();

        ifc_tria_glob[t] = (int *) malloc(ntifc[t]*sizeof(int));
        ifc_tria_loc[t]  = (int *) malloc(ntifc[t]*sizeof(int));

        std::set<int>::iterator itsn;
        u = 0;
        //std::cout << "number " << ntifc[t] << " for rank = " << world_rank << " -> ";
        for(itsn=its->second.begin();itsn!=its->second.end();itsn++)
        {
            ifc_tria_glob[t][u] = *itsn;

               //std::cout << ifc_tria_glob[t][u] << " ";

            //
            u++;
        }

            //std::cout << std::endl;

        t++;
    }
    int** ifc_nodes_glob = (int **) malloc(ncomm*sizeof(int *));
    int** ifc_nodes_loc  = (int **) malloc(ncomm*sizeof(int *));
//    int **ifc_nodes_glob[ncomm];
//    int **ifc_nodes_loc[ncomm];
    t = 0;
    for(its=proc2nodes.begin();its!=proc2nodes.end();its++)
    {
        color_node[t]=its->first;
        npifc[t] = its->second.size();

        ifc_nodes_glob[t] = (int *) malloc(npifc[t]*sizeof(int));
        ifc_nodes_loc[t]  = (int *) malloc(npifc[t]*sizeof(int));

        std::set<int>::iterator itsn;
        u = 0;
        std::cout << "number " << npifc[t] << " for rank = " << world_rank << " -> ";
        for(itsn=its->second.begin();itsn!=its->second.end();itsn++)
        {
            ifc_nodes_glob[t][u] = *itsn;

                std::cout << ifc_nodes_glob[t][u] << " ";

            u++;
        }

            std::cout << std::endl;

        t++;
    }

//    for(int i=0;i<ncomm;i++)
//    {
//        std::cout << "colornodes = " << color_node[i] << " " << ntifc[i] << " " << npifc[i] << std::endl;
//    }

    int k;

    get_local_mesh(nv, ne, nt, vert_mask, inv_vert_mask,
                   tetra_mask, tria_mask, inv_tria_mask,
                    vert_coor,tmesh->vert_coor_all,
                    vert_ref,tmesh->vert_ref_all,
                    tetra_vert,tmesh->tetra_vert_all,
                    tetra_ref,tmesh->tetra_ref_all,
                    tria_vert,tmesh->tria_vert_all,
                    tria_ref,tmesh->tria_ref_all,
                    met,tmesh->met_all,ncomm,
                    ntifc,ifc_tria_loc,ifc_tria_glob,
                    npifc,ifc_nodes_loc,ifc_nodes_glob);
//
   /** 1) Manually set your mesh using the PMMG_Set* functions */

     /** a) give the size of the mesh */
     int nVertices       = nv;
     int nTetrahedra     = ne;
     int nPrisms         = 0;
     int nTriangles      = nt;
     int nQuadrilaterals = 0;
     int nEdges          = 0;

    //std::cout << "stats for rank " << world_rank << " :: " << nv << " " << ne << " " << nt << std::endl;
     if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                               nQuadrilaterals,nEdges) != 1 ) {
       MPI_Finalize();
       exit(EXIT_FAILURE);
     }

     /** b) set vertices, tetrahedra, triangles */
//    if(world_rank == 0)
//    {
//        std::cout << "v = [";
//        for(int i=0;i<nVertices;i++)
//        {
//            std::cout << vert_coor[i*3+0] << ", " << vert_coor[i*3+1] << ", " << vert_coor[i*3+2] << "," << std::endl;
//        }
//        std::cout << "]"<<std::endl;
//        std::cout << "tet = [";
//        for(int i=0;i<nTetrahedra;i++)
//        {
//            std::cout << tetra_vert[i*4+0] << ", " << tetra_vert[i*4+1] << ", " << tetra_vert[i*4+2] << "," << tetra_vert[i*4+3] << "," << std::endl;
//        }
//        std::cout << "]"<<std::endl;
//        std::cout << "tria = [";
//        for(int i=0;i<nTriangles;i++)
//        {
//            std::cout << tria_vert[i*3+0] << ", " << tria_vert[i*3+1] << ", " << tria_vert[i*3+2] << "," << std::endl;
//        }
//        std::cout << "]"<<std::endl;
//    }



     if ( !opt ) {
       /* By array: give the array of the vertices coordinates such as the
        * coordinates of the k^th point are stored in vert_coor[3*(k-1)],
        * vert_coor[3*(k-1)+1] and vert_coor[3*(k-1)+2] */
       if ( PMMG_Set_vertices(parmesh,vert_coor,vert_ref) != 1 ) {
         MPI_Finalize();
         exit(EXIT_FAILURE);
       }
     }
     else {
       /* Vertex by vertex: for each vertex, give the coordinates, the reference
          and the position in mesh of the vertex */
       for ( k=0; k<nVertices; ++k ) {
         pos = 3*k;
           std::cout << "pos " << pos << std::endl;
         if ( PMMG_Set_vertex(parmesh,vert_coor[pos],vert_coor[pos+1],vert_coor[pos+2],
                              vert_ref[k], k+1) != 1 ) {
           MPI_Finalize();
           exit(EXIT_FAILURE);
         }
       }
     }


     if ( !opt ) {
       /* By array: give the array of the tetra vertices and the array of the tetra
        * references. The array of the tetra vertices is such as the four
        * vertices of the k^th tetra are respectively stored in
        * tetra_vert[4*(k-1)],tetra_vert[4*(k-1)+1],tetra_vert[4*(k-1)+2] and
        * tetra_vert[4*(k-1)+3]. */
       if ( PMMG_Set_tetrahedra(parmesh,tetra_vert,tetra_ref) != 1 ) {
         MPI_Finalize();
         exit(EXIT_FAILURE);
       }
     }
     else {
       /* Vertex by vertex: for each tetrahedra,
         give the vertices index, the reference and the position of the tetra */
       for ( k=0; k<nTetrahedra; ++k ) {
         pos = 4*k;
         if ( PMMG_Set_tetrahedron(parmesh,tetra_vert[pos],tetra_vert[pos+1],
                                   tetra_vert[pos+2],tetra_vert[pos+3],tetra_ref[k],k+1) != 1 ) {
           MPI_Finalize();
           exit(EXIT_FAILURE);
         }
       }
     }


     if ( !opt ) {
       /* By array: give the array of the tria vertices and the array of the tria
        * references. The array of the tria vertices is such as the three
        * vertices of the k^th tria are stored in
        * tria_vert[3*(k-1)], tria_vert[3*(k-1)+1] and tria_vert[4*(k-1)+2] */
       if ( PMMG_Set_triangles(parmesh,tria_vert,tria_ref) != 1 ) {
         MPI_Finalize();
         exit(EXIT_FAILURE);
       }
     }
     else {
       /* Vertex by vertex: for each triangle, give the vertices index, the
        * reference and the position of the triangle */
       for ( k=0; k<nTriangles; ++k ) {
         pos = 3*k;
         if ( PMMG_Set_triangle(parmesh,
                                tria_vert[pos],tria_vert[pos+1],tria_vert[pos+2],
                                tria_ref[k],k+1) != 1 ) {
           MPI_Finalize();
           exit(EXIT_FAILURE);
         }
       }
     }


     /** 2) Build metric in ParMmg format */
     /** Two solutions: just use the PMMG_loadMet_centralized function that will
         read a .sol(b) file formatted or manually set your sol using the
         PMMG_Set* functions */

     /** Manually set of the metric */
     /** a) give info for the metric structure: metric applied on vertex entities,
         number of vertices, the metric is scalar*/
     if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) {
       MPI_Finalize();
       exit(EXIT_FAILURE);
     }


//     if ( !opt ) {
//       /* by array */
//       if ( PMMG_Set_scalarMets(parmesh,met) != 1 ) {
//         MPI_Finalize();
//         exit(EXIT_FAILURE);
//       }
//     }
//     else {
//       /* vertex by vertex */
//       for ( k=0; k<nVertices ; k++ ) {
//         if ( PMMG_Set_scalarMet(parmesh,met[k],k+1) != 1 ) {
//           MPI_Finalize();
//           exit(EXIT_FAILURE);
//         }
//       }
//     }


     /** 3) Build solutions in PMMG format */
     /** Two solutions: just use the PMMG_loadAllSols_centralized function that
         will read a .sol(b) file formatted or manually set your solutions using
         the PMMG_Set* functions */

     /** With parmmg setters: 3 sols per vertex, the first is scalar,
         the second vectorial, the third tensorial  */
     const int nSolsAtVertices = 1; // 3 solutions per vertex
     int solType[1] = {MMG5_Tensor};

//     if ( PMMG_Set_solsAtVerticesSize(parmesh,nSolsAtVertices,nVertices,solType) != 1 ) {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//     }


     /** a) give solutions values and positions:
      - First solution (scalar) is equal to x^2 + y^2 + z^2
      - Second (vector) is (x,y,z)
      - Third (Tensor) is (100,0,0,100/(z+1),0,100/(z*z+1))
     */
     double scalar_sol[tmesh->nv],vector_sol[tmesh->nv*3],tensor_sol[tmesh->nv*6];
    double hmax = 0.1;
     for ( k=0; k<nVertices; k++ ) {
       pos = 3*k;

       /* First solution */
       scalar_sol[k] = vert_coor[pos]*vert_coor[pos]
         + vert_coor[pos+1]*vert_coor[pos+1]
         + vert_coor[pos+2]*vert_coor[pos+2];

       /* Second */
       vector_sol[3*k]   = vert_coor[pos];
       vector_sol[3*k+1] = vert_coor[pos+1];
       vector_sol[3*k+2] = vert_coor[pos+2];

        double x = vert_coor[pos];
        double y = vert_coor[pos+1];
        double z = vert_coor[pos+2];
        double a = 10*fabs(0.75-sqrt(x*x+y*y));
        double hx = min(0.002*pow(5,a),hmax);
        double hy = min(0.05*pow(2,a),hmax);
        double hz = hmax;
        double rat;
        if(x==0)
        {
            rat = 1.0;
        }
        else
        {
            rat = y/x;
        }

        double theta = atan(rat);
        //std::cout << hx << " " << hy << " " << hz << " " << theta << " " << rat << std::endl;
        double m11 = (1.0/(hx*hx))*cos(theta)*cos(theta)+(1.0/(hy*hy))*sin(theta)*sin(theta);
        double m12 = (1.0/(hx*hx)-1.0/(hy*hy))*cos(theta)*sin(theta);
        double m13 = 0.0;
        double m22 = (1.0/(hx*hx))*sin(theta)*sin(theta)+(1.0/(hy*hy))*cos(theta)*cos(theta);
        double m23 = 0.0;
        double m33 = hz;
        //std::cout << m11 << " " << m12 << "  " << m13 << " " << m22 << " " << m23 << "  " << m33 << std::endl;

        /* Third */
        tensor_sol[6*k]   = m11;
        tensor_sol[6*k+1] = m12;
        tensor_sol[6*k+2] = m13;
        tensor_sol[6*k+3] = m22;
        tensor_sol[6*k+4] = m23;
        tensor_sol[6*k+5] = m33;
        PMMG_Set_tensorMet(parmesh,m11,m12,m13,m22,m23,m33,k+1);
       /* Third */
//       tensor_sol[6*k]   = 100.;
//       tensor_sol[6*k+1] = 0.;
//       tensor_sol[6*k+2] = 0.;
//       tensor_sol[6*k+3] = 100./(vert_coor[pos+2]+1.);
//       tensor_sol[6*k+4] = 0.;
//       tensor_sol[6*k+5] = 100./(vert_coor[pos+2]*vert_coor[pos+2]+1.);
     }

     if ( !opt ) {
//       /* Give the solution by array */
//       /* First solution */
//       if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,1,scalar_sol) != 1 ) {
//         MPI_Finalize();
//         exit(EXIT_FAILURE);
//       }
//       /* Second */
//       if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,2,vector_sol) != 1 ) {
//         MPI_Finalize();
//         exit(EXIT_FAILURE);
//       }
       /* Third */
//       if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,1,tensor_sol) != 1 ) {
//         MPI_Finalize();
//         exit(EXIT_FAILURE);
//       }
     }
     else {
       /* Vertex by vertex */
       for ( k=0; k<nVertices; k++ ) {
//         /* First solution */
//         if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,1,&(scalar_sol[k]),k+1) != 1 ) {
//           MPI_Finalize();
//           exit(EXIT_FAILURE);
//         }
//         /* Second */
//         pos = 3*k;
//         if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,2,&(vector_sol[pos]),k+1) != 1 ) {
//           MPI_Finalize();
//           exit(EXIT_FAILURE);
//         }
         /* Third */
//         pos = 6*(k-1);
//         if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,1,&(tensor_sol[pos]),k+1) != 1 ) {
//           MPI_Finalize();
//           exit(EXIT_FAILURE);
//         }
       }
     }


     /** 4) Initialization of interface communicators in ParMMG.
      *     The user can choose between providing triangles (faces) interface
      *     information (through the PMMG_APIDISTRIB_faces parameter), or nodes
      *     interface information (through the PMMG_APIDISTRIB_nodes parameter).
      */

     /* Set API mode */
     if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
       MPI_Finalize();
       exit(EXIT_FAILURE);
     };

     /* Select face or node interface API */
     switch( API_mode ) {

       case PMMG_APIDISTRIB_faces :
         if( !world_rank ) printf("\n--- API mode: Setting face communicators\n");

         /* Set the number of interfaces */
         ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);

         /* Loop on each interface (proc pair) seen by the current rank) */
         for( icomm=0; icomm<ncomm; icomm++ ) {

           /* Set nb. of entities on interface and rank of the outward proc */
           ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                                  color_face[icomm],
                                                  ntifc[icomm]);

           /* Set local and global index for each entity on the interface */
           ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                                    ifc_tria_loc[icomm],
                                                    ifc_tria_glob[icomm], 1 );
         }
         break;

       case PMMG_APIDISTRIB_nodes :
         if( !world_rank ) printf("\n--- API mode: Setting node communicators\n");

         /* Set the number of interfaces */
         ier = PMMG_Set_numberOfNodeCommunicators(parmesh, ncomm);

         /* Loop on each interface (proc pair) seen by the current rank) */
         for( icomm=0; icomm<ncomm; icomm++ ) {

           /* Set nb. of entities on interface and rank of the outward proc */
           ier = PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm,
                                                  color_node[icomm],
                                                  npifc[icomm]);

           /* Set local and global index for each entity on the interface */
           ier = PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm,
                                                    ifc_nodes_loc[icomm],
                                                    ifc_nodes_glob[icomm], 1 );
         }
         break;
     }

   //  /** save mesh and interfaces **/
   //  char filemesh[48];
   //  sprintf(filemesh,"cube_in.%d.mesh",parmesh->myrank);
   //  MMG3D_saveMesh(parmesh->listgrp[0].mesh,filemesh);
   //
   //  FILE *fid;
   //  sprintf(filemesh,"cube_in.%d.mesh_parFaces",parmesh->myrank);
   //  fid = fopen(filemesh,"w");
   //  fprintf(fid,"# Number of communicators\n%d\n",ncomm);
   //  for( icomm = 0; icomm < ncomm; icomm++ ) {
   //    fprintf(fid,"\n# Color\n%d\n# Number of items\n%d\n# Local and global enumeration\n",color_face[icomm],ntifc[icomm]);
   //    for( i = 0; i < ntifc[icomm]; i++ )
   //      fprintf(fid,"%d %d\n",ifc_tria_loc[icomm][i],ifc_tria_glob[icomm][i]);
   //  }
   //  fclose(fid);
   //
   //  sprintf(filemesh,"cube_in.%d.mesh_parNodes",parmesh->myrank);
   //  fid = fopen(filemesh,"w");
   //  fprintf(fid,"# Number of communicators\n");
   //  fprintf(fid,"%d\n",ncomm);
   //  for( icomm = 0; icomm < ncomm; icomm++ ) {
   //    fprintf(fid,"\n# Color\n%d\n# Number of items\n%d\n# Local and global enumeration\n",color_node[icomm],npifc[icomm]);
   //    for( i = 0; i < npifc[icomm]; i++ )
   //      fprintf(fid,"%d %d\n",ifc_nodes_loc[icomm][i],ifc_nodes_glob[icomm][i]);
   //  }
   //  fclose(fid);


     /** ------------------------------ STEP III -------------------------- */
     /** remesh step */
    int i;
     /* Set number of iterations */
     if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) ) {
       MPI_Finalize();
       exit(EXIT_FAILURE);
     };

     /* remeshing function */
     ierlib = PMMG_parmmglib_distributed( parmesh );

        if ( ierlib == PMMG_SUCCESS )
        {
            std::cout << "succes!" << std::endl;
        }
        std::cout << "ierlib " << ierlib << std::endl;

      //if ( ierlib == PMMG_SUCCESS ) {

        /** If no remeshing is performed (zero remeshing iterations), check set
         * parallel interfaces against input data. */
        if( !niter ) {

          /* Check matching of input interface nodes with the set ones */
          if( !PMMG_Check_Set_NodeCommunicators(parmesh,ncomm,npifc,
                                                color_node,ifc_nodes_loc) ) {
            printf("### Wrong set node communicators!\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }

          /* Get input triangle nodes */
          faceNodes = (int **) malloc(ncomm*sizeof(int *));
          for( icomm = 0; icomm < ncomm; icomm++ ) {
            faceNodes[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
            for( i = 0; i < ntifc[icomm]; i++ ) {
              pos = ifc_tria_loc[icomm][i];
              faceNodes[icomm][3*i]   = tria_vert[3*(pos-1)];
              faceNodes[icomm][3*i+1] = tria_vert[3*(pos-1)+1];
              faceNodes[icomm][3*i+2] = tria_vert[3*(pos-1)+2];
            }
          }

          /* Check matching of input interface triangles with the set ones */
          if( !PMMG_Check_Set_FaceCommunicators(parmesh,ncomm,ntifc,
                                             color_face,faceNodes) ) {
            printf("### Wrong set face communicators!\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }
        }


        /** ------------------------------ STEP  IV -------------------------- */
        /** recover parallel interfaces */

        int next_node_comm,next_face_comm,*nitem_node_comm,*nitem_face_comm;
        int *color_node_out,*color_face_out;
        int **out_tria_loc, **out_node_loc;

        /* Get number of node interfaces */
        ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);

        /* Get outward proc rank and number of nodes on each interface */
        color_node_out  = (int *) malloc(next_node_comm*sizeof(int));
        nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
        for( icomm=0; icomm<next_node_comm; icomm++ )
          ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                                 &color_node_out[icomm],
                                                 &nitem_node_comm[icomm]);

        /* Get IDs of nodes on each interface */
        out_node_loc = (int **) malloc(next_node_comm*sizeof(int *));
        for( icomm=0; icomm<next_node_comm; icomm++ )
          out_node_loc[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
        ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);

        /* Get number of face interfaces */
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);

        /* Get outward proc rank and number of faces on each interface */
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);

        /* Get IDs of triangles on each interface */
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);

    /*
        for( icomm=0; icomm<next_node_comm; icomm++ )
          for( i=0; i < nitem_node_comm[icomm]; i++ )
            printf("rank %d comm %d node %d\n",parmesh->myrank,icomm,out_node_loc[icomm][i]);
        for( icomm=0; icomm<next_face_comm; icomm++ )
          for( i=0; i < nitem_face_comm[icomm]; i++ )
            printf("rank %d comm %d tria %d\n",parmesh->myrank,icomm,out_tria_loc[icomm][i]);
    */

        /** If no remeshing is performed (zero remeshing iterations), check
         *  retrieved parallel interfaces against input data. */
        if( !niter ) {

          /* Check matching of input interface nodes with the output ones */
          if( !PMMG_Check_Get_NodeCommunicators(parmesh,
                                                ncomm,npifc,
                                                color_node,ifc_nodes_loc,
                                                next_node_comm,nitem_node_comm,
                                                color_node_out,out_node_loc) ) {
            printf("### Wrong retrieved node communicators!\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }

          /* Get output triangles (as the boundary is re-generated, triangle IDs
           * have changed) */
          nVertices   = 0;
          nTetrahedra = 0;
          nTriangles  = 0;
          nEdges      = 0;
          if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                                 &nEdges) !=1 ) {
          ier = PMMG_STRONGFAILURE;
          }

          int *ref       = (int*)calloc(nTriangles,sizeof(int));
          int *required  = (int*)calloc(nTriangles,sizeof(int));
          int *triaNodes = (int*)calloc(3*nTriangles,sizeof(int));

          if ( PMMG_Get_triangles(parmesh,triaNodes,ref,required) != 1 ) {
            fprintf(stderr,"Unable to get mesh triangles\n");
            ier = PMMG_STRONGFAILURE;
          }

          int** faceNodes_out = (int **) malloc(next_face_comm*sizeof(int *));
          for( icomm = 0; icomm < next_face_comm; icomm++ ) {
            faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
            for( i = 0; i < nitem_face_comm[icomm]; i++ ) {
              pos = out_tria_loc[icomm][i];
              faceNodes_out[icomm][3*i]   = triaNodes[3*(pos-1)];
              faceNodes_out[icomm][3*i+1] = triaNodes[3*(pos-1)+1];
              faceNodes_out[icomm][3*i+2] = triaNodes[3*(pos-1)+2];
            }
          }

          /* Check matching of input interface triangles with the output ones */
          if( !PMMG_Check_Get_FaceCommunicators(parmesh,ncomm,ntifc,
                                                color_face,faceNodes,
                                                next_face_comm,nitem_face_comm,
                                                color_face_out,faceNodes_out) ) {
            printf("### Wrong retrieved face communicators!\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
          }

          free(ref);
          free(required);
          free(triaNodes);
          for( icomm = 0; icomm < ncomm; icomm++ )
            free(faceNodes[icomm]);
          free(faceNodes);
          for( icomm = 0; icomm < next_face_comm; icomm++ )
            free(faceNodes_out[icomm]);
          free(faceNodes_out);
        }

        free(color_node);
        free(color_face);
        free(npifc);
        free(ntifc);
        for( icomm = 0; icomm < ncomm; icomm++ ) {
          free(ifc_nodes_loc[icomm]);
          free(ifc_nodes_glob[icomm]);
          free(ifc_tria_loc[icomm]);
          free(ifc_tria_glob[icomm]);
        }

        free(color_node_out);
        free(color_face_out);
        free(nitem_node_comm);
        free(nitem_face_comm);
        for( icomm = 0; icomm < next_node_comm; icomm++ )
          free(out_node_loc[icomm]);
        free(out_node_loc);
        for( icomm = 0; icomm < next_face_comm; icomm++ )
          free(out_tria_loc[icomm]);
        free(out_tria_loc);

        if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                               &nEdges) !=1 ) {
          ier = PMMG_STRONGFAILURE;
        }
    
    std::cout << nVertices << " " << nTetrahedra << " " << nTriangles << " " << nEdges << std::endl;
    double* vert = new double[nVertices*3];
    int* tetra = new int[nTetrahedra*4];
    int* tria = new int[nTriangles*3];
    int* edge = new int[nEdges*2];
    int* stats = new int[4];
    stats[0]=nVertices;
    stats[1]=nTetrahedra;
    stats[2]=nTriangles;
    stats[3]=nEdges;
    int maks = largest(stats,4);
    /* Table to store the vertices/tetra/triangles/edges references */
    int *ref = new int[maks];
    int *corner = new int[nVertices];
    int *required = new int[maks];
    int *ridge = new int[nEdges];
    int nreq=0;
    int nc = 0;
    
    
    string filename = "ParMmgMesh_rank_" + std::to_string(world_rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"ParMMGbaby_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << nVertices << ", E = " << nTetrahedra << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
    for ( k=0; k<nVertices; k++ )
    {
        pos = 3*k;
        if ( PMMG_Get_vertex(parmesh,&(vert[pos]),&(vert[pos+1]),&(vert[pos+2]),
                             &(ref[k]),&(corner[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh vertex %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( corner && corner[k] )  nc++;
        if ( required && required[k] )  nreq++;
    }
    for ( k=0; k<nVertices; k++ )
    {
        pos = 3*k;
        myfile << vert[pos] << " " << vert[pos+1] << " " << vert[pos+2] << std::endl;
    }
    
    for ( k=0; k<nTetrahedra; k++ )
    {
        pos = 4*k;
        if ( PMMG_Get_tetrahedron(parmesh,
                                  &(tetra[pos  ]),&(tetra[pos+1]),
                                  &(tetra[pos+2]),&(tetra[pos+3]),
                                  &(ref[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh tetra %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( required && required[k] )  nreq++;
      
    }
    for ( k=0; k<nTetrahedra; k++ ) {
      pos = 4*k;
        myfile << tetra[pos] << " " << tetra[pos+1] << " " << tetra[pos+2] << " " << tetra[pos+3] << std::endl;
    }
        
        
        /** ------------------------------ STEP V  ---------------------------- */
        /** get results */
        /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
            that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
            using the PMMG_getMesh/PMMG_getSol functions */

        /** 1) Get the mesh with ParMmg getters and save it at the Medit file format */
//        if( !(inm = fopen(fileout,"w")) ) {
//          fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
//          exit(EXIT_FAILURE);
//        }
//        fprintf(inm,"MeshVersionFormatted 2\n");
//        fprintf(inm,"\nDimension 3\n");
//
//        /** a) get the size of the mesh: vertices, tetra, triangles, edges and
//         * allocate the arrays to receive data */
//        nVertices   = 0;
//        nTetrahedra = 0;
//        nTriangles  = 0;
//        nEdges      = 0;
//        if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
//                               &nEdges) !=1 ) {
//          ier = PMMG_STRONGFAILURE;
//        }
//
//        /* Table to store the vertices */
//        double *vert = (double*)calloc((nVertices)*3,sizeof(double));
//        if ( !vert ) {
//          perror("  ## Memory problem: point calloc");
//          nVertices = 0;
//          ier = PMMG_STRONGFAILURE;
//        }
//
//        /* Table to store the tetra */
//        int *tetra = (int*)calloc((nTetrahedra)*4,sizeof(int));
//        if ( !tetra ) {
//          perror("  ## Memory problem: tetra calloc");
//          nTetrahedra = 0;
//          ier = PMMG_STRONGFAILURE;
//        }
//
//        /* Table to store the tria */
//        int *tria = (int*)calloc((nTriangles)*3,sizeof(int));
//        if ( !tria ) {
//          perror("  ## Memory problem: tria calloc");
//          nTriangles = 0;
//          ier = PMMG_STRONGFAILURE;
//        }
//
//        /* Table to store the edges */
//        int *edge = (int*)calloc((nEdges)*2,sizeof(int));
//        if ( !edge ) {
//          perror("  ## Memory problem: edge calloc");
//          nEdges = 0;
//          ier = PMMG_STRONGFAILURE;
//        }
//
//          int* stats = new int[4];
//          stats[0]=nVertices;
//          stats[1]=nTetrahedra;
//          stats[2]=nTriangles;
//          stats[3]=nEdges;
//         int maks = largest(stats,4);
//        /* Table to store the vertices/tetra/triangles/edges references */
//        int *ref = (int*)calloc(maks,sizeof(int));
//        if ( !ref ) {
//          perror("  ## Memory problem: ref calloc");
//          MPI_Finalize();
//          exit(EXIT_FAILURE);
//        }
//
//        /* Table to know if a vertex is corner */
//        int *corner = (int*)calloc(nVertices,sizeof(int));
//        if ( !corner ) {
//          perror("  ## Memory problem: corner calloc");
//          MPI_Finalize();
//          exit(EXIT_FAILURE);
//        }
//
//        /* Table to know if a vertex/tetra/tria/edge is required */
//        int *required = (int*)calloc(maks,sizeof(int));
//        if ( !required ) {
//          perror("  ## Memory problem: required calloc");
//          MPI_Finalize();
//          exit(EXIT_FAILURE);
//        }
//
//        /* Table to know if an edge delimits a sharp angle */
//        int *ridge = (int*)calloc(nEdges ,sizeof(int));
//        if ( !ridge ) {
//          perror("  ## Memory problem: ridge calloc");
//          MPI_Finalize();
//          exit(EXIT_FAILURE);
//        }
//
//        /** b) Vertex recovering */
//        nreq = nc = 0;
//        fprintf(inm,"\nVertices\n%d\n",nVertices);
//
//        if ( !opt ) {
//          /* By array */
//          if ( PMMG_Get_vertices(parmesh,vert,ref,corner,required) != 1 ) {
//            fprintf(inm,"Unable to get mesh vertices \n");
//            ier = PMMG_STRONGFAILURE;
//          }
//          for ( k=0; k<nVertices; k++ ) {
//            if ( corner && corner[k] )  nc++;
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        else {
//          /* Vertex by vertex */
//          for ( k=0; k<nVertices; k++ ) {
//            pos = 3*k;
//            if ( PMMG_Get_vertex(parmesh,&(vert[pos]),&(vert[pos+1]),&(vert[pos+2]),
//                                 &(ref[k]),&(corner[k]),&(required[k])) != 1 ) {
//              fprintf(inm,"Unable to get mesh vertex %d \n",k);
//              ier = PMMG_STRONGFAILURE;
//            }
//            if ( corner && corner[k] )  nc++;
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        for ( k=0; k<nVertices; k++ ) {
//          pos = 3*k;
//          fprintf(inm,"%.15lg %.15lg %.15lg %d \n",vert[pos],vert[pos+1],vert[pos+2],ref[k]);
//        }
//
//        fprintf(inm,"\nCorners\n%d\n",nc);
//        for ( k=0; k<nVertices; k++ ) {
//          if ( corner && corner[k] )  fprintf(inm,"%d \n",k);
//        }
//        fprintf(inm,"\nRequiredVertices\n%d\n",nreq);
//        for ( k=0; k<nVertices; k++ ) {
//          if ( required && required[k] )  fprintf(inm,"%d \n",k);
//        }
//        free(corner);
//        corner = NULL;
//
//        /** d) Triangles recovering */
//        nreq = 0;
//        fprintf(inm,"\nTriangles\n%d\n",nTriangles);
//
//        if ( !opt ) {
//          /* By array */
//          if ( PMMG_Get_triangles(parmesh,tria,ref,required) != 1 ) {
//            fprintf(inm,"Unable to get mesh triangles\n");
//            ier = PMMG_STRONGFAILURE;
//          }
//          for ( k=0; k<nTriangles; k++ ) {
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        else {
//          /* Triangle by triangle */
//          for ( k=0; k<nTriangles; k++ ) {
//            pos = 3*k;
//            if ( PMMG_Get_triangle(parmesh,&(tria[pos]),&(tria[pos+1]),&(tria[pos+2]),
//                                   &(ref[k]),&(required[k])) != 1 ) {
//              fprintf(inm,"Unable to get mesh triangle %d \n",k);
//              ier = PMMG_STRONGFAILURE;
//            }
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        for ( k=0; k<nTriangles; k++ ) {
//          pos = 3*k;
//          fprintf(inm,"%d %d %d %d \n",tria[pos],tria[pos+1],tria[pos+2],ref[k]);
//        }
////
////
//        fprintf(inm,"\nRequiredTriangles\n%d\n",nreq);
//        for ( k=0; k<nTriangles; k++ ) {
//          if ( required && required[k] )  fprintf(inm,"%d \n",k);
//        }
//
//        /** e) Edges recovering */
//        nreq = 0;nr = 0;
//        fprintf(inm,"\nEdges\n%d\n",nEdges);
//
//        if ( !opt ) {
//          /* By array */
//          if ( PMMG_Get_edges(parmesh,edge,ref,ridge,required) != 1 ) {
//            fprintf(inm,"Unable to get mesh edges\n");
//            ier = PMMG_STRONGFAILURE;
//          }
//          for ( k=0; k<nEdges; k++ ) {
//            if ( ridge && ridge[k] )  nr++;
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        else {
//          /* Edge by edge */
//          for ( k=0; k<nEdges; k++ ) {
//            pos = 2*k;
//            if ( PMMG_Get_edge(parmesh,&(edge[pos]),&(edge[pos+1]),
//                               &(ref[k]),&(ridge[k]),&(required[k])) != 1 ) {
//              fprintf(inm,"Unable to get mesh edge %d \n",k);
//              ier = PMMG_STRONGFAILURE;
//            }
//            if ( ridge && ridge[k] )  nr++;
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        for ( k=0; k<nEdges; k++ ) {
//          pos = 2*k;
//          fprintf(inm,"%d %d %d \n",edge[pos],edge[pos+1],ref[k]);
//        }
//
//        fprintf(inm,"\nRequiredEdges\n%d\n",nreq);
//        for ( k=0; k<nEdges; k++ ) {
//          if ( required && required[k] )  fprintf(inm,"%d \n",k);
//        }
//        fprintf(inm,"\nRidges\n%d\n",nr);
//        for ( k=0; k<nEdges; k++ ) {
//          if ( ridge && ridge[k] )  fprintf(inm,"%d \n",k);
//        }
//
//        /** c) Tetra recovering */
//        nreq = 0;
//        fprintf(inm,"\nTetrahedra\n%d\n",nTetrahedra);
//
//        if ( !opt ) {
//          /* By array */
//          if ( PMMG_Get_tetrahedra(parmesh,tetra,ref,required) != 1 ) {
//            fprintf(inm,"Unable to get mesh tetra\n");
//            ier = PMMG_STRONGFAILURE;
//          }
//          for ( k=0; k<nTetrahedra; k++ ) {
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        else {
//          /* Tetra by tetra */
//          for ( k=0; k<nTetrahedra; k++ ) {
//            pos = 4*k;
//            if ( PMMG_Get_tetrahedron(parmesh,
//                                      &(tetra[pos  ]),&(tetra[pos+1]),
//                                      &(tetra[pos+2]),&(tetra[pos+3]),
//                                      &(ref[k]),&(required[k])) != 1 ) {
//              fprintf(inm,"Unable to get mesh tetra %d \n",k);
//              ier = PMMG_STRONGFAILURE;
//            }
//            if ( required && required[k] )  nreq++;
//          }
//        }
//        for ( k=0; k<nTetrahedra; k++ ) {
//          pos = 4*k;
//          fprintf(inm,"%d %d %d %d %d \n",
//                  tetra[pos],tetra[pos+1],tetra[pos+2],tetra[pos+3],ref[k]);
//        }
//
//        fprintf(inm,"\nRequiredTetrahedra\n%d\n",nreq);
//        for ( k=0; k<nTetrahedra; k++ ) {
//          if ( required && required[k] )  fprintf(inm,"%d \n",k);
//        }
//
//        fprintf(inm,"\nEnd\n");
//        fclose(inm);
//
//        free(vert)    ; vert     = NULL;
//        free(tetra)   ; tetra    = NULL;
//        free(tria)    ; tria     = NULL;
//        free(edge)    ; edge     = NULL;
//        free(ref)     ; ref      = NULL;
//        free(required); required = NULL;
//        free(ridge)   ; ridge    = NULL;
//
//        /** 3) Get the metric with ParMmg getters */
//        if ( !(inm = fopen(metout,"w")) ) {
//          fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT METRIC FILE.\n");
//          exit(EXIT_FAILURE);
//        }
//        fprintf(inm,"MeshVersionFormatted 2\n");
//        fprintf(inm,"\nDimension 3\n");
//
//        /** a) get the size of the metric: type of entity to which apply the
//            metric(SolAtVertices,...), number of entities to which apply the metric,
//            type of solution (scalar, tensor...) */
//        nVertices = 0;
//        int typEntity,typSol;
//        if ( PMMG_Get_metSize(parmesh,&typEntity,&nVertices,&typSol) != 1 ) {
//          printf("Unagle to get metric size\n");
//          nVertices = 0;
//          ier = PMMG_LOWFAILURE;
//        }
//
//        /* We set a scalar metric so the output metric must be scalar */
//        if ( ( typEntity != MMG5_Vertex )  || ( typSol != MMG5_Scalar ) ) {
//          MPI_Finalize();
//          exit(EXIT_FAILURE);
//        }
//
//        /** b) Vertex recovering */
//        double *sol = (double*)calloc(nVertices+1 ,sizeof(double));
//
//        fprintf(inm,"\nSolAtVertices\n%d\n",nVertices);
//        fprintf(inm,"1 1 \n\n");
//        if ( !opt ) {
//          /* by array */
//          if ( PMMG_Get_scalarMets(parmesh,sol) != 1 ) {
//            fprintf(inm,"Unable to get metrics\n");
//            ier = PMMG_LOWFAILURE;
//          }
//        }
//        else {
//          for ( k=0; k<nVertices; k++ ) {
//            /* Vertex by vertex */
//            if ( PMMG_Get_scalarMet(parmesh,&sol[k]) != 1 ) {
//              fprintf(inm,"Unable to get metrics %d \n",k);
//              ier = PMMG_LOWFAILURE;
//            }
//          }
//        }
//        for ( k=0; k<nVertices; k++ ) {
//          fprintf(inm,"%.15lg \n",sol[k]);
//        }
//
//        fprintf(inm,"\nEnd\n");
//        fclose(inm);
//
//        free(sol);
//
//        /** 4) Get the solutions with ParMmg getters */
//        // To implement when ParMmg will be ready
////      }
////      else if ( ierlib == PMMG_STRONGFAILURE ) {
////        fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
////      }

      /** 5) Free the PMMG5 structures */
//      PMMG_Free_all(PMMG_ARG_start,
//                    PMMG_ARG_ppParMesh,&parmesh,
//                    PMMG_ARG_end);

//      free(fileout); fileout = NULL;
//      free(solout) ; solout  = NULL;
//      free(metout) ; metout  = NULL;

    MPI_Finalize();
    
    return 0;
     
}
