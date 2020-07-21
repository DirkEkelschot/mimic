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
    int offsety_0 = 0;
    int offsety_1 = N;
    int offsetz = N*N;
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


    if ( MMG3D_Set_meshSize(mmgMesh,nvertices,tets.size(),0,triangles.size(),0,0) != 1 )  exit(EXIT_FAILURE);
    for(int i=0;i<nvertices;i++)
    {
        mmgMesh->point[i+1].c[0]  = verts[i].x;  mmgMesh->point[i+1].c[1]  = verts[i].y; mmgMesh->point[i+1].c[2]  = verts[i].z; mmgMesh->point[i+1].ref  = 0;

    }
    
    std::cout << "Test mesh stats:" << std::endl;
    std::cout << " Number of tetrahedra = " << tets.size() << std::endl;
    std::cout << " Number of triangles = " << triangles.size() << std::endl;
    std::cout << " Number of triangles map = " << Tr.size() << std::endl;
    std::cout << " size F2E = " << F2E.size() << std::endl;
    std::cout << " Number of vertices = " << nvertices << std::endl;
    
    
    for(int i=1;i<=tets.size();i++)
    {
        mmgMesh->tetra[i].v[0] = tets[i-1][0]+1;
        mmgMesh->tetra[i].v[1] = tets[i-1][1]+1;
        mmgMesh->tetra[i].v[2] = tets[i-1][2]+1;
        mmgMesh->tetra[i].v[3] = tets[i-1][3]+1;
       // mmgMesh->tetra[i].ref  = 1;

    }

    MMG3D_saveMesh(mmgMesh,"meshname.mesh");
    
    MmgTdata->mmgMesh = mmgMesh;
    MmgTdata->E2N = E2N;
    MmgTdata->F2N = F2N;
    MmgTdata->E2F = E2F;
    MmgTdata->F2E = F2E;
    MmgTdata->E2N = E2N;
    return MmgTdata;
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
    std::clock_t startr;
    startr = std::clock();
    //  GetXadjandAdjcyArrays(iee,ien,comm);
    //  Example3DPartitioningWithParVarParMetis();
    //  ExampleUS3DPartitioningWithParVarParMetis();
    //Example3DPartitioningWithParVarParMetis();
//============================================================

    //const char* fn_conn="grids/piston/conn.h5";
    const char* fn_conn="grids/piston/conn.h5";
    const char* fn_grid="grids/piston/grid.h5";
    const char* fn_data="grids/adept/data.h5";
    const char* fn_adept="grids/adept/conn.h5";

//
//    Array<int>*    zdefs  = ReadDataSetFromGroupFromFile<int>(fn_adept,"zones","zdefs");
//    Array<char>*  znames  = ReadDataSetFromGroupFromFile<char>(fn_adept,"zones","znames");
//    PlotBoundaryData(znames,zdefs,comm);
//    ParArray<int>* ien    = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
//    ParArray<int>* ief    = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
//    //ParArray<int>* ife  = ReadDataSetFromFileInParallel<int>(fn_conn,"ife",comm,info);
//    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
//
//    UnitTestEigenDecomp();




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
    
//    PMMG_pParMesh parmesh;
//    PMMG_Init_parMesh(PMMG_ARG_start,
//    PMMG_ARG_ppParMesh, &parmesh,
//    PMMG_ARG_pMesh, PMMG_ARG_pMet,
//    PMMG_ARG_dim, 3,
//    PMMG_ARG_MPIComm, MPI_COMM_WORLD,
//    PMMG_ARG_end);
    
    //int nvertices = 12;
    //if ( MMG3D_Set_meshSize(mmgMesh,12,12,0,20,0,0) != 1 )  exit(EXIT_FAILURE);

    double N = 2.0;
    
    MmgTestData* mmgdata = GetStructuredBlockMmgMesh(N, mmgMesh, mmgSol);

    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,mmgMesh->np,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
    double hmax = 0.2;
    /** b) give solutions values and positions */
    for(k=1 ; k<=mmgMesh->np ; k++)
    {
      //mmgSol->m[k] = 0.5;
        double x = mmgMesh->point[k].c[0];
        double y = mmgMesh->point[k].c[1];
        double z = mmgMesh->point[k].c[2];
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
        if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,k) != 1 ) exit(EXIT_FAILURE);
        //if ( MMG3D_Set_tensorSol(mmgSol, m11,m12,m13,m22,m23,m33,k) != 1 ) exit(EXIT_FAILURE);
        //if ( MMG3D_Set_tensorSol(mmgSol, 1.0,0.0,0.0,0.1,0.0,1.0,k) != 1 ) exit(EXIT_FAILURE);
      /* or with the api function :
         if ( MMG3D_Set_scalarSol(mmgSol,0.5,k) != 1 ) exit(EXIT_FAILURE); */
    }
    /** 4) If you don't use the API functions, you MUST call
        the MMG3D_Set_handGivenMesh() function. Don't call it if you use
        the API functions */
    //MMG3D_Set_handGivenMesh(mmgMesh);

    /** 5) (not mandatory): check if the number of given entities match with mesh size */
    if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

//    * ------------------------------ STEP  II --------------------------
//    * remesh function
//     WARNING: the MMG3D_mmg3dlib function returns 1 if success, 0 if fail.
//     The MMG3D4 library was working opposite.
    //ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

    ofstream myfile2;
    myfile2.open("mmgMesh_v3.dat");
    myfile2 << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile2 <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    //std::cout << " verts check " << LVerts.size() << " " << hx.size() << std::endl;
    myfile2 <<"ZONE N = " << mmgMesh->np << ", E = " << mmgMesh->ne << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
    for(int i=0;i<mmgMesh->np;i++)
    {
        myfile2 << mmgMesh->point[i+1].c[0] << " " <<mmgMesh->point[i+1].c[1] << " " << mmgMesh->point[i+1].c[2] << std::endl;
    }

    for(int i=1;i<=mmgMesh->ne;i++)
    {
        myfile2 << mmgMesh->tetra[i].v[0] << " " << mmgMesh->tetra[i].v[1] << " " << mmgMesh->tetra[i].v[2] << " " << mmgMesh->tetra[i].v[3] << std::endl;
    }

    myfile2.close();

//    for(int i=1;i<=mmgMesh->ne;i++)
//    {
//       std::cout << i << " :: " << mmgMesh->tetra[i].v[0] << " " << mmgMesh->tetra[i].v[1] << " " << mmgMesh->tetra[i].v[2] << " " << mmgMesh->tetra[i].v[3] << std::endl;
//    }

    ParArray<int>* ien_tet = new ParArray<int>(mmgMesh->ne,4,comm);
    int Nel = ien_tet->getNglob();
    ParallelState* pstate  = new ParallelState(Nel,comm);
    int nloc_tet = pstate->getNloc(world_rank);
    int offset = pstate->getOffset(world_rank);

    for(int i=0;i<nloc_tet;i++)
    {
        for(int j=0;j<4;j++)
        {
            ien_tet->setVal(i,j,(mmgMesh->tetra[offset+i+1].v[j]-1));
        }
    }

    PartTmp* pTmp = getPartionIDS(ien_tet, pstate, comm, 4);
//    std::map<int,std::vector<int> > adj_elements;
////    for(int i=0;i<mmgdata->F2E.size();i++)
////    {
////        for(int j=0;j<mmgdata->F2E[i].size();j++)
////        {
////            std::cout << mmgdata->F2E[i][j] << " ";
////        }
////
////        std::cout << std::endl;
////    }
    
    std::map<int,std::vector<int> > E2F1 = mmgdata->E2F;
    std::map<int,std::vector<int> >::iterator maps1;
    if(world_rank == 0)
    {
        std::cout << "E2F1 " << E2F1.size() << std::endl;

        std::cout << "before = " << std::endl;
        for(maps1=E2F1.begin();maps1!=E2F1.end();maps1++)
        {
            int l = maps1->second.size();
            
            for(int i=0;i<l;i++)
            {
                std::cout << maps1->second[i] << " ";
            }
            
            std::cout << std::endl;
            
        }
    }
    UpdateConnectivityMmgMesh(mmgdata,comm);
    
    std::map<int,std::vector<int> > E2F2 = mmgdata->E2F;
    std::map<int,std::vector<int> >::iterator maps2;
    if(world_rank == 0)
    {
        std::cout << "E2F2 " << E2F2.size() << std::endl;

        std::cout << "after = " << std::endl;
        for(maps2=E2F2.begin();maps2!=E2F2.end();maps2++)
        {
            int l = maps2->second.size();

            for(int i=0;i<l;i++)
            {
                std::cout << maps2->second[i] << " ";
            }

            std::cout << std::endl;

        }
    }
    
    std::map<int,std::vector<int> > E2F = mmgdata->E2F;
    std::map<int,std::vector<int> > F2N = mmgdata->F2N;
    int fadj = 0;
    std::map<int, std::set<int> > proc2face;
    std::map<int, std::set<int> > proc2nodes;
    std::map<int, int> face2proc;
    std::vector<int> partfaces;
    for(int i=0;i<pTmp->lPart->getNrow();i++)
    {
        std::set<int> owned_faces;
        for(int n=0;n<4;n++)
        {
           owned_faces.insert(E2F[offset+i][n]);
        }
        //
        int start = pTmp->xadj[i];
        int end   = pTmp->xadj[i+1];
        //
        for(int j=start;j<end;j++)
        {
            int adjEl_id = pTmp->adjcny[j];
            int p_id = pTmp->gPart->getVal(adjEl_id,0);
            if(p_id!=world_rank) // find the adjacent element that is not on this partition.
            {
                //adj_elements[p_id].push_back(adjEl_id);

                for(int s=0;s<4;s++)
                {
                    fadj = E2F[adjEl_id][s];

                    if(owned_faces.find(fadj)!=owned_faces.end())
                    {
                        partfaces.push_back(fadj);
                        proc2face[p_id].insert(fadj);
                        face2proc[fadj]=p_id;
                        for(int g=0;g<3;g++)
                        {
                            proc2nodes[p_id].insert(F2N[fadj][g]);
                        }
                        
                    }
                }
            }
        }
        owned_faces.clear();
    }
    
    int ncomm = proc2face.size();
    int* color_node = new int[ncomm];
    int* color_face = new int[ncomm];
    int* ntifc = new int[ncomm];
    int* npifc = new int[ncomm];
    std::map<int,std::set<int> >::iterator its;
    int t = 0;
    int u = 0;
    
    int *ifc_tria_glob[ncomm];
    int *ifc_tria_loc[ncomm];
    for(its=proc2face.begin();its!=proc2face.end();its++)
    {
        color_face[t]=its->first;
        ntifc[t] = its->second.size();
        
        ifc_tria_glob[t] = new int[ntifc[t]];
        ifc_tria_loc[t]  = new int[ntifc[t]];
        
        std::set<int>::iterator itsn;
        u = 0;
        for(itsn=its->second.begin();itsn!=its->second.end();itsn++)
        {
            ifc_tria_glob[t][u] = *itsn;
            u++;
        }
        t++;
    }
    
    int *ifc_nodes_glob[ncomm];
    int *ifc_nodes_loc[ncomm];
    t = 0;
    for(its=proc2nodes.begin();its!=proc2nodes.end();its++)
    {
        color_node[t]=its->first;
        npifc[t] = its->second.size();
        
        ifc_nodes_glob[t] = new int[npifc[t]];
        ifc_nodes_loc[t] = new int[npifc[t]];
        
        std::set<int>::iterator itsn;
        u = 0;
        for(itsn=its->second.begin();itsn!=its->second.end();itsn++)
        {
            ifc_nodes_glob[t][u] = *itsn;
            u++;
        }
        t++;
    }
    
    
//    if(world_rank == 0)
//    {
//        std::cout << "rank = " << world_rank << std::endl;
//        for(int i=0;i<ncomm;i++)
//        {
//            std::cout << "faces SHARED with proc " << color_face[i] << " -> ";
//            for(int j=0;j<ntifc[i];j++)
//            {
//                std::cout << ifc_tria_glob[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
//
//        for(int i=0;i<ncomm;i++)
//        {
//            std::cout << "nodes SHARED with proc " << color_node[i] << " -> ";
//
//            for(int j=0;j<npifc[i];j++)
//            {
//                std::cout << ifc_nodes_glob[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    
    
    MPI_Finalize();
    
    return 0;
     
}
