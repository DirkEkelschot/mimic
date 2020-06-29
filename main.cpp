#include "adapt_io.h"
#include "adapt_compute.h"
#include "adapt_operations.h"
#include "adapt_math.h"
#include <math.h>
#include "adapt_output.h"
#include "adapt_partition.h"
//#include "adapt.h"
#include <iomanip>
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
    int nloc = P->getParallelState()->getNloc(world_rank);
    int of = P->getParallelState()->getOffset(world_rank);
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

Array<double>* ComputeHessianNew(Partition* P, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    Array<double>* hessian = new Array<double>(Nel,3);
    
    std::map<int,std::vector<int> > gE2lV       = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gV       = P->getGlobElem2GlobVerts();
    std::vector<Vert> locVerts                  = P->getLocalVerts();
    int loc_vid = 0;int glob_vid=0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    std::map<int,std::map<int,int> >::iterator itadj;
    std::vector<double> Hi;
    double u_ip1jk = 0.0;
    double u_im1jk = 0.0;
    double u_ijp1k = 0.0;
    double u_ijm1k = 0.0;
    double u_ijkp1 = 0.0;
    double u_ijkm1 = 0.0;
    int offset = P->getParallelState()->getOffset(rank);
    int* xadj = P->getXadj();
    int* adjcny = P->getAdjcny();
    int nloc = P->getParallelState()->getNloc(rank);
    Vert Vip1jk;
    Vert Vim1jk;
    Vert Vijp1k;
    Vert Vijm1k;
    Vert Vijkp1;
    Vert Vijkm1;
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
        
        Vip1jk.x = 0.0;Vip1jk.y=0.0;Vip1jk.z=0.0;
        Vim1jk.x = 0.0;Vim1jk.y=0.0;Vim1jk.z=0.0;
        Vijp1k.x = 0.0;Vijp1k.y=0.0;Vijp1k.z=0.0;
        Vijm1k.x = 0.0;Vijm1k.y=0.0;Vijm1k.z=0.0;
        Vijkp1.x = 0.0;Vijkp1.y=0.0;Vijkp1.z=0.0;
        Vijkm1.x = 0.0;Vijkm1.y=0.0;Vijkm1.z=0.0;
        Vert Vijk = ComputeCenterCoord(Pijk,8);

        
        /*
        if(Vijk.y<1.0e-03)
	{
		for(int i=0;i<8;i++)
		{
			std::cout << Pijk[i*3+0] << " " << Pijk[i*3+1] << " " << Pijk[i*3+2] << std::endl;

		}
		
	}
	*/
	delete[] Pijk; 
        int tel   = 0;
        int f = 0;
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            
            double u_po = P->getU0atGlobalElem(adjEl_id);
            double* Po = new double[np*3];
            
            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                glob_vid  = gE2gV[adjEl_id][k];
                Po[k*3+0] = locVerts[loc_vid].x;
                Po[k*3+1] = locVerts[loc_vid].y;
                Po[k*3+2] = locVerts[loc_vid].z;
		if(adjEl_id == 5233975 && rank == 0)
		{
			f = 1;
			std::cout <<"check dan " <<  glob_vid << " :: " << Po[k*3+0] << " " << Po[k*3+1] << " " << Po[k*3+2] << std::endl;
		}
            }
            
            Vert Vpo = ComputeCenterCoord(Po,8);
            if(f == 1)
	    {

		std::cout << "check dan 2 " << Vpo.x << " " << Vpo.y << " " << Vpo.z << std::endl;
	    	f = 0;
	    }
            Vrt->setVal(tel,0,Vpo.x-Vijk.x);
            Vrt->setVal(tel,1,Vpo.y-Vijk.y);
            Vrt->setVal(tel,2,Vpo.z-Vijk.z);
            b->setVal(tel,0,u_po-u_ijk);
            
            delete[] Po;
            
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
	
	if(i==2359155 && rank == 0)
	{
		std::cout << "element = " << i << std::endl;
		std::cout << "Vijk = " << Vijk.x << " " << Vijk.y << " " << Vijk.z << std::endl;
		std::cout << "nadj = " << nadj << std::endl;	
		for(int s=0;s<nadj;s++)
		{
			for(int j=0;j<3;j++)
			{
				std::cout << "dXnew["<<s<<","<<j<<"]=" << Vrt->getVal(s,j) << "; ";
			}
			std::cout << std::endl;

		}
		std::cout << std::endl;
		for(int s=0;s<nadj;s++)
		{

			std::cout << "b["<<s<<"]="<< b->getVal(s,0) << std::endl;
		}

		for(int s=0;s<3;s++)
		{
			std::cout << "x=["<<s<<"]="<< x->getVal(s,0) << std::endl;
		}
	}

        
        hessian->setVal(i,0,x->getVal(0,0));
        hessian->setVal(i,1,x->getVal(1,0));
        hessian->setVal(i,2,x->getVal(2,0));
        
        vijkIDs.clear();
        delete Vrt_T;
        delete Vrt;
        delete b;
        delete[] A_cm;
    }
    
    delete[] xadj;
    delete[] adjcny;
        
    return hessian;
}


Array<double>* ComputeHessian(Partition* P, int Nel, std::map<int,std::map<int,int> > E2Etopo, MPI_Comm comm)
{
    
    // d^2u/dx^2 = (u_(i+1,j,k)-2u_(i,j,k)+u_(i-1,j,k))/(dx)^2 // faces 0-1
    // d^2u/dy^2 = (u_(i,j+1,k)-2u_(i,j,k)+u_(i,j-1,k))/(dy)^2 // faces 2-3
    // d^2u/dz^2 = (u_(i,j,k+1)-2u_(i,j,k)+u_(i,j,k-1))/(dz)^2 // faces 4-5
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    
    Array<double>* hessian = new Array<double>(Nel,3);
    
    std::map<int,std::vector<int> > gE2lV       = P->getGlobElem2LocVerts();

    std::vector<Vert> locVerts                  = P->getLocalVerts();
    //
    double d2udx2 = 0.0;
    double d2udy2 = 0.0;
    double d2udz2 = 0.0;
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    Array<double>* Vrt_T = new Array<double>(3,6);
    Array<double>* Vrt = new Array<double>(6,3);
    Array<double>* b = new Array<double>(6,1);
    
    std::map<int,std::map<int,int> >::iterator itadj;
    std::vector<double> Hi;
	        double u_ip1jk = 0.0; 
        double u_im1jk = 0.0; 
        double u_ijp1k = 0.0; 
        double u_ijm1k = 0.0; 
        double u_ijkp1 = 0.0; 
        double u_ijkm1 = 0.0; 
    for(itadj=E2Etopo.begin();itadj!=E2Etopo.end();itadj++)
    {
        for(int i=0;i<6;i++)
        {
            
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,i,0.0);
                Vrt->setVal(i,j,0.0);
            }
        }
        
        b->setVal(0,0,0.0);
        b->setVal(1,0,0.0);
        b->setVal(2,0,0.0);
        b->setVal(3,0,0.0);
        b->setVal(4,0,0.0);
        b->setVal(5,0,0.0);
        
        std::map<int,int>::iterator itmap;
        double u_ijk = P->getU0atGlobalElem(itadj->first);
        std::vector<int> vijkIDs = gE2lV[itadj->first];
        
        double* Pijk = new double[np*3];
        //std::cout << "========================" << std::endl;
        for(int i=0;i<vijkIDs.size();i++)
        {
            loc_vid     = vijkIDs[i];
            Pijk[i*3+0] = locVerts[loc_vid].x;
            Pijk[i*3+1] = locVerts[loc_vid].y;
            Pijk[i*3+2] = locVerts[loc_vid].z;
        }
        Vert Vip1jk;
	Vip1jk.x = 0.0;Vip1jk.y=0.0;Vip1jk.z=0.0;
	Vert Vim1jk;
 	Vim1jk.x = 0.0;Vim1jk.y=0.0;Vim1jk.z=0.0;
	Vert Vijp1k;
	Vijp1k.x = 0.0;Vijp1k.y=0.0;Vijp1k.z=0.0;
	Vert Vijm1k;
	Vijm1k.x = 0.0;Vijm1k.y=0.0;Vijm1k.z=0.0;
	Vert Vijkp1;
	Vijkp1.x = 0.0;Vijkp1.y=0.0;Vijkp1.z=0.0;
	Vert Vijkm1;
	Vijkm1.x = 0.0;Vijkm1.y=0.0;Vijkm1.z=0.0;
        Vert Vijk = ComputeCenterCoord(Pijk,8);


	//std::cout << "Verts " << std::endl;
	//std::cout << " ijk " << Vijk.x << " " << Vijk.y << " " << Vijk.z << std::endl;
        delete[] Pijk;
        if(itadj->second.find(0) != itadj->second.end() && itadj->second.find(1) != itadj->second.end())
        {
            //std::cout << " it " << std::endl; 
            std::vector<int> vip1jkIDs = gE2lV[itadj->second[0]];
            double* Pip1jk = new double[np*3];
            for(int i=0;i<vip1jkIDs.size();i++)
            {
                loc_vid  = vip1jkIDs[i];
                Pip1jk[i*3+0] = locVerts[loc_vid].x;
                Pip1jk[i*3+1] = locVerts[loc_vid].y;
                Pip1jk[i*3+2] = locVerts[loc_vid].z;
            }
            
            Vip1jk = ComputeCenterCoord(Pip1jk,8);
            
            Vrt->setVal(0,0,Vip1jk.x-Vijk.x);
            Vrt->setVal(0,1,Vip1jk.y-Vijk.y);
            Vrt->setVal(0,2,Vip1jk.z-Vijk.z);
            
            //-------------------------------------------------------------------
            
            std::vector<int> vim1jkIDs = gE2lV[itadj->second[1]];
            double* Pim1jk = new double[np*3];
            for(int i=0;i<vim1jkIDs.size();i++)
            {
                loc_vid  = vim1jkIDs[i];
                Pim1jk[i*3+0] = locVerts[loc_vid].x;
                Pim1jk[i*3+1] = locVerts[loc_vid].y;
                Pim1jk[i*3+2] = locVerts[loc_vid].z;
            }
            //std::cout << std::endl;
            Vim1jk = ComputeCenterCoord(Pim1jk,8);
            
            Vrt->setVal(1,0,Vim1jk.x-Vijk.x);
            Vrt->setVal(1,1,Vim1jk.y-Vijk.y);
            Vrt->setVal(1,2,Vim1jk.z-Vijk.z);

            
            /*
            double dim1jk = sqrt((Vim1jk.x-Vijk.x)*(Vim1jk.x-Vijk.x)+
                                 (Vim1jk.y-Vijk.y)*(Vim1jk.y-Vijk.y)+
                                 (Vim1jk.z-Vijk.z)*(Vim1jk.z-Vijk.z));
            
            double dip1jk = sqrt((Vip1jk.x-Vijk.x)*(Vip1jk.x-Vijk.x)+
                                 (Vip1jk.y-Vijk.y)*(Vip1jk.y-Vijk.y)+
                                 (Vip1jk.z-Vijk.z)*(Vip1jk.z-Vijk.z));
            */
            
            
            u_ip1jk = P->getU0atGlobalElem(itadj->second[0]);
            u_im1jk = P->getU0atGlobalElem(itadj->second[1]);

	    //std::cout << " 0 " << Vip1jk.x << " " << Vip1jk.y << " " << Vip1jk.z << std::endl;
	    //std::cout << " 1 " << Vim1jk.x << " " << Vim1jk.y << " " << Vim1jk.z << std::endl;
            b->setVal(0,0,u_ip1jk-u_ijk);
            b->setVal(1,0,u_im1jk-u_ijk);
            
            delete[] Pip1jk;
            delete[] Pim1jk;
            //std::cout << "b0 = " << u_ip1jk-u_ijk << " " << u_ip1jk << " " << u_ijk << std::endl;
            //std::cout << "b1 = " << u_im1jk-u_ijk << " " << u_im1jk << " " << u_ijk << std::endl;
            //d2udx2 = (u_ip1jk-2*u_ijk+u_im1jk)/(dim1jk+dip1jk);
            
            
        }
        
        
        if(itadj->second.find(2) != itadj->second.end() && itadj->second.find(3) != itadj->second.end())
        {
            std::vector<int> vijp1kIDs = gE2lV[itadj->second[2]];
            double* Pijp1k = new double[np*3];
            for(int i=0;i<vijp1kIDs.size();i++)
            {
                loc_vid  = vijp1kIDs[i];
                Pijp1k[i*3+0] = locVerts[loc_vid].x;
                Pijp1k[i*3+1] = locVerts[loc_vid].y;
                Pijp1k[i*3+2] = locVerts[loc_vid].z;
            }
            Vijp1k = ComputeCenterCoord(Pijp1k,8);
            
            Vrt->setVal(2,0,Vijp1k.x-Vijk.x);
            Vrt->setVal(2,1,Vijp1k.y-Vijk.y);
            Vrt->setVal(2,2,Vijp1k.z-Vijk.z);
            
            //-------------------------------------------------------------------
            
            std::vector<int> vijm1kIDs = gE2lV[itadj->second[3]];
            double* Pijm1k = new double[np*3];
            for(int i=0;i<vijm1kIDs.size();i++)
            {
                loc_vid  = vijm1kIDs[i];
                Pijm1k[i*3+0] = locVerts[loc_vid].x;
                Pijm1k[i*3+1] = locVerts[loc_vid].y;
                Pijm1k[i*3+2] = locVerts[loc_vid].z;
            }
            //std::cout << std::endl;
            Vijm1k = ComputeCenterCoord(Pijm1k,8);
            
            Vrt->setVal(3,0,Vijm1k.x-Vijk.x);
            Vrt->setVal(3,1,Vijm1k.y-Vijk.y);
            Vrt->setVal(3,2,Vijm1k.z-Vijk.z);
        
            /*
            
            double dijm1k = sqrt((Vijm1k.x-Vijk.x)*(Vijm1k.x-Vijk.x)+
                                 (Vijm1k.y-Vijk.y)*(Vijm1k.y-Vijk.y)+
                                 (Vijm1k.z-Vijk.z)*(Vijm1k.z-Vijk.z));
            
            double dijp1k = sqrt((Vijp1k.x-Vijk.x)*(Vijp1k.x-Vijk.x)+
                                 (Vijp1k.y-Vijk.y)*(Vijp1k.y-Vijk.y)+
                                 (Vijp1k.z-Vijk.z)*(Vijp1k.z-Vijk.z));
            */
            
            u_ijp1k = P->getU0atGlobalElem(itadj->second[2]);
            u_ijm1k = P->getU0atGlobalElem(itadj->second[3]);
            //std::cout << u_ijp1k << " " << u_ijm1k << "=========> " << " " << Vijm1k.x << " " << Vijk.x << " " << Vijp1k.x << " :: " << " " << Vijm1k.y << " " << Vijk.y << " " << Vijp1k.y << " ::  " << Vijm1k.z << " " << Vijk.z << " " << Vijp1k.z << std::endl;        
            b->setVal(2,0,u_ijp1k-u_ijk);
            b->setVal(3,0,u_ijm1k-u_ijk);
	    //std::cout << " 2 " << Vijp1k.x << " " << Vijp1k.y << " " << Vijp1k.z << std::endl;
            //std::cout << " 3 " << Vijm1k.x << " " << Vijm1k.y << " " << Vijm1k.z<< std::endl;
            delete[] Pijp1k;
            delete[] Pijm1k;
            //std::cout << "b2 = " << u_ijp1k-u_ijk << " " << u_ijp1k<< " " << u_ijk << std::endl;
            //std::cout << "b3 = " << u_ijm1k-u_ijk << " " << u_ijm1k<< " " << u_ijk << std::endl;
//
//            d2udy2 = (u_ijp1k-2*u_ijk+u_ijm1k)/(dijm1k+dijp1k);
//
            
        }
        
        
        if(itadj->second.find(4) != itadj->second.end() && itadj->second.find(5) != itadj->second.end())
        {
            
            //d2udz2 = (u_ijkp1-2*u_ijk+u_ijkm1);
            
            std::vector<int> vijkp1IDs = gE2lV[itadj->second[4]];
            double* Pijkp1 = new double[np*3];
            for(int i=0;i<vijkp1IDs.size();i++)
            {
                loc_vid  = vijkp1IDs[i];
                Pijkp1[i*3+0] = locVerts[loc_vid].x;
                Pijkp1[i*3+1] = locVerts[loc_vid].y;
                Pijkp1[i*3+2] = locVerts[loc_vid].z;
            }
            Vijkp1 = ComputeCenterCoord(Pijkp1,8);
            
            Vrt->setVal(4,0,Vijkp1.x-Vijk.x);
            Vrt->setVal(4,1,Vijkp1.y-Vijk.y);
            Vrt->setVal(4,2,Vijkp1.z-Vijk.z);
            
            
            //-------------------------------------------------------------------
            
            
            std::vector<int> vijkm1IDs = gE2lV[itadj->second[5]];
            double* Pijkm1 = new double[np*3];
            for(int i=0;i<vijkm1IDs.size();i++)
            {
                loc_vid  = vijkm1IDs[i];
                Pijkm1[i*3+0] = locVerts[loc_vid].x;
                Pijkm1[i*3+1] = locVerts[loc_vid].y;
                Pijkm1[i*3+2] = locVerts[loc_vid].z;
            }
            Vijkm1 = ComputeCenterCoord(Pijkm1,8);
            
            Vrt->setVal(5,0,Vijkm1.x-Vijk.x);
            Vrt->setVal(5,1,Vijkm1.y-Vijk.y);
            Vrt->setVal(5,2,Vijkm1.z-Vijk.z);
            
            /*
            double dijkm1 = sqrt((Vijkm1.x-Vijk.x)*(Vijkm1.x-Vijk.x)+
                                 (Vijkm1.y-Vijk.y)*(Vijkm1.y-Vijk.y)+
                                 (Vijkm1.z-Vijk.z)*(Vijkm1.z-Vijk.z));
            
            double dijkp1 = sqrt((Vijkp1.x-Vijk.x)*(Vijkp1.x-Vijk.x)+
                                 (Vijkp1.y-Vijk.y)*(Vijkp1.y-Vijk.y)+
                                 (Vijkp1.z-Vijk.z)*(Vijkp1.z-Vijk.z));
            */
            
            u_ijkp1 = P->getU0atGlobalElem(itadj->second[4]);
            u_ijkm1 = P->getU0atGlobalElem(itadj->second[5]);
            
            b->setVal(4,0,u_ijkp1-u_ijk);
            b->setVal(5,0,u_ijkm1-u_ijk);
            //std::cout << " 4 " << Vijkp1.x << " " << Vijkp1.y << " " << Vijkp1.z  << std::endl;
            //std::cout << " 5 " << Vijkm1.x << " " << Vijkm1.y << " " << Vijkm1.z  << std::endl;
            //std::cout << "b4 = " << u_ijkp1-u_ijk << std::endl;
            //std::cout << "b5 = " << u_ijkm1-u_ijk << std::endl;

            //d2udy2 = (u_ijkp1-2*u_ijk+u_ijkm1)/(dijkm1+dijkp1);
            delete[] Pijkp1;
            delete[] Pijkm1;
        }

        //std::cout << std::endl;
        

//        for(int i=0;i<6;i++)
//        {
//            std::cout << b->getVal(i,0) << " ";
//        }
//        std::cout << std::endl;
	/*
	std::cout << "=====" << std::endl;
	 for(int i=0;i<6;i++)
        {    
            for(int j=0;j<3;j++)
            {    
                std::cout << Vrt->getVal(i,j) << " ";
            }    
	    std::cout << std::endl;
        }    
        std::cout << "=====" << std::endl;
	for(int i=0;i<3;i++)
        {
            for(int j=0;j<6;j++)
            {
                Vrt_T->setVal(i,j,Vrt->getVal(j,i));
            }
        }

        Array<double>* R    = MatMul(Vrt_T,Vrt);
        bool isdiag = isDiagonalMatrix(R);
        if(isdiag==true)
        {
            Array<double>*Rinv = new Array<double>(3,3);
            for(int i=0;i<Rinv->getNrow();i++)
            {
                for(int j=0;j<Rinv->getNcol();j++)
                {
                    Rinv->setVal(i,j,0.0);
                }
            }
            
            if(R->getVal(0,0)!=0.0)
            {
                Rinv->setVal(0,0,1.0/R->getVal(0,0));
            }
            else
            {
                Rinv->setVal(0,0,0.0);
            }
            if(R->getVal(1,1)!=0.0)
            {
                Rinv->setVal(1,1,1.0/R->getVal(1,1));
            }
            else
            {
                Rinv->setVal(1,1,0.0);
            }
            if(R->getVal(2,2)!=0.0)
            {
                Rinv->setVal(2,2,1.0/R->getVal(2,2));
            }
            else
            {
                Rinv->setVal(2,2,0.0);
            }
            Array<double>* Rn   = MatMul(Rinv,Vrt_T);
            Array<double>* x    = MatMul(Rn,b);

            hessian->setVal(e,0,x->getVal(0,0));
            hessian->setVal(e,1,x->getVal(1,0));
            hessian->setVal(e,2,x->getVal(2,0));

            if(std::isnan(x->getVal(0,0)) || fabs(x->getVal(0,0))>1.0e7)
	    {
		x->setVal(0,0,0.0);
		x->setVal(1,0,0.0);
		x->setVal(2,0,0.0);
		std::cout << rank << " diag -> " << " NaN" << std::endl;
		for(int u=0;u<Rinv->getNrow();u++)
		{
			for(int w=0;w<Rinv->getNcol();w++)
			{
				std::cout << R->getVal(u,w) << " ";
			}
			std::cout << std::endl;
		}
            }

	    delete Rinv;
            delete[] R;
            delete[] Rn;
        }
        else
        {
	    
            Array<double>*Rinv  = MatInv(R);
            Array<double>* Rn   = MatMul(Rinv,Vrt_T);
            Array<double>* x    = MatMul(Rn,b);
            if(std::isnan(x->getVal(0,0)) || x->getVal(0,0)>1.0e16)
	    {
		std::cout << rank << " non diag -> " << " NaN" << std::endl;
		for(int u=0;u<Vrt->getNrow();u++)
		{
			for(int w=0;w<Vrt->getNcol();w++)
			{
				std::cout << Vrt->getVal(u,w) << " ";
			}
			std::cout << std::endl;
		}
	    }
	    
	    double* A_cm = new double[6*3];
	    for(int i=0;i<6;i++)
    	    {
        	for(int j=0;j<3;j++)
        	{
           		A_cm[j*6+i] = Vrt->getVal(i,j);
        	}
    	    }
	    Array<double>* x = SolveQR(A_cm,6,3,b);
            if(std::isnan(x->getVal(0,0)) || x->getVal(0,0)>1.0e02)
	    {
		std::cout << rank << " non diag -> " << " NaN" << std::endl;
                for(int u=0;u<Vrt->getNrow();u++)
                {
                        for(int w=0;w<Vrt->getNcol();w++)
                        {
                                std::cout << Vrt->getVal(u,w) << " ";
                        }
                        std::cout << std::endl;
                }
		for(int u=0;u<x->getNrow();u++)
                {    
                        for(int w=0;w<x->getNcol();w++)
                        {    
                                std::cout << x->getVal(u,w) << " "; 
                        }    
                        std::cout << std::endl;
                } 
		x->setVal(0,0,0.0);
		x->setVal(1,0,0.0);
		x->setVal(2,0,0.0); 

	    }
	*/

	if(itadj->second.find(0) != itadj->second.end() && itadj->second.find(1) != itadj->second.end()
		&& itadj->second.find(2) != itadj->second.end() && itadj->second.find(3) != itadj->second.end() && itadj->second.find(4) != itadj->second.end()  && itadj->second.find(5) != itadj->second.end())
        {  
            double* A_cm = new double[6*3];
            for(int i=0;i<6;i++)
            {
                for(int j=0;j<3;j++)
                {
                    A_cm[j*6+i] = Vrt->getVal(i,j);
                }
            }
            b->setVal(0,0,u_ip1jk-u_ijk);
            b->setVal(1,0,u_im1jk-u_ijk);
            b->setVal(2,0,u_ijp1k-u_ijk);
            b->setVal(3,0,u_ijm1k-u_ijk);
            b->setVal(4,0,u_ijkp1-u_ijk);
            b->setVal(5,0,u_ijkm1-u_ijk);

            Array<double>* x = SolveQR(A_cm,6,3,b);
	 
            hessian->setVal(e,0,x->getVal(0,0));
            hessian->setVal(e,1,x->getVal(1,0));
            hessian->setVal(e,2,x->getVal(2,0));
            
            delete[] A_cm;
	    }
	    else
	    {
            hessian->setVal(e,0,1.0);
            hessian->setVal(e,1,1.0);
            hessian->setVal(e,2,1.0);
            
        }

        e++;
    }
    
    return hessian;
}


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
    const char* fn_conn="grids/adept/conn.h5";
    const char* fn_grid="grids/adept/grid.h5";
    const char* fn_data="grids/adept/data.h5";
    const char* fn_adept="grids/adept/conn.h5";
    
    
    Array<int>*    zdefs = ReadDataSetFromGroupFromFile<int>(fn_adept,"zones","zdefs");
    Array<char>*  znames = ReadDataSetFromGroupFromFile<char>(fn_adept,"zones","znames");
    PlotBoundaryData(znames,zdefs,comm);
    ParArray<int>* iee    = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);
    Array<double>* xcn2    = ReadDataSetFromFile<double>(fn_grid,"xcn");
    if(world_rank == 0)
    {
    for(int i=0;i<iee->getNcol();i++)
    {
	std::cout << iee->getVal(2359155,i) << " ";
    }
    }    
    ParArray<int>* ien    = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    if(world_rank == 1)
    {
	std::vector<int> vertvec;
	for(int i=0;i<ien->getNcol()-1;i++)
	{
		//std::cout << ien->getVal(5233975-2839081,i) << " ";
		vertvec.push_back(ien->getVal(5233975-2839081,i+1)-1);
	}

	std::cout << std::endl;
	int np = 8;int loc_vid; 
        double* Pijk = new double[np*3];
        for(int k=0;k<vertvec.size();k++)
        {    
            loc_vid     = vertvec[k];
	    std::cout << loc_vid << std::endl;
            Pijk[k*3+0] = xcn2->getVal(loc_vid,0);
            Pijk[k*3+1] = xcn2->getVal(loc_vid,1);
            Pijk[k*3+2] = xcn2->getVal(loc_vid,2);
        }    
     
        Vert Vijk = ComputeCenterCoord(Pijk,8);

        std::cout <<"result=" << Vijk.x << " " << Vijk.y << " " << Vijk.z << std::endl;
    }    
    std::cout << std::endl;
    
    ParArray<int>* ief    = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    
    
    ParArray<int>* ife    = ReadDataSetFromFileInParallel<int>(fn_conn,"ife",comm,info);

    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
//    ParArray<double>* ifn = ReadDataSetFromFileInParallel<double>(fn_grid,"ifn",comm,info);
//    ParArray<double>* boundaries = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_6","boundaries",comm,info);
    
    int Nel = ien->getNglob();
    int Nel_part = ien->getNrow();
    ParArray<double>* interior   = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",Nel,comm,info);
    ParArray<double>* xcn_new   = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","xcn",Nel,comm,info);
//===================================================================================

    ParallelState* pstate = new ParallelState(ien->getNglob(),comm);
    //
    int nglob = ien->getNglob();
    int nrow  = ien->getNrow();
    int ncol  = ien->getNcol()-1;
    //
    ParArray<int>* ien_copy = new ParArray<int>(nglob,ncol,comm);
    //
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            ien_copy->setVal(i,j,ien->getVal(i,j+1)-1);
        }
    }
    //
    int ncol_ief = ief->getNcol()-1;
    ParArray<int>* ief_copy = new ParArray<int>(nglob,ncol_ief,comm);

    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol_ief;j++)
        {
            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
        }
    }

    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(ien_copy,comm,8);

    ParallelState* xcn_parstate = new ParallelState(xcn_new->getNglob(),comm);
    ParArray<double>* var = new ParArray<double>(Nel,1,comm);
    
    for(int i=0;i<Nel_part;i++)
    {
        double v = (i+ien->getOffset(world_rank))*0.2420;
        //var->setVal(i,0,v);
        var->setVal(i,0,interior->getVal(i,0));
    }
    
    double t0 = MPI_Wtime();
    Partition* P = new Partition(ien_copy, ief_copy, parmetis_pstate,pstate, xcn_new,xcn_parstate,var,comm);
    double t00 = MPI_Wtime();
    double timing = t00-t0;
    double max_time = 0.0;
    MPI_Allreduce(&timing, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    
    if (world_rank == 0)
    {
        std::cout << "t_max := " << max_time  << std::endl;
    }
    
    
//    double  t1 = MPI_Wtime();
//    std::map<int,map<int,int> > E2Etopo = getElement2ElementTopology(P, comm);
//    double t2 = MPI_Wtime();
//    double timing2 = t2-t1;
//      //std::cout << "time spent scheming " << t2-t1 << std::endl;
//    double max_time2 = 0.0;
//    MPI_Allreduce(&timing2, &max_time2, 1, MPI_DOUBLE, MPI_MAX, comm);
//    if (world_rank == 0)
//    {
//        std::cout << "t_max2 := " << max_time2  << std::endl;
//    }

    std::map<int,std::map<int,int> >::iterator itadj;
//    for(itadj=E2Etopo.begin();itadj!=E2Etopo.end();itadj++)
//    {
//        std::map<int,int>::iterator itm;
//
//        for(itm=itadj->second.begin();itm!=itadj->second.end();itm++)
//        {
//            if(itm->first == 1 || itm->first == 0)
//            {
//                std::cout << itadj->first << " => " << itm->first << " " << itm->second;
//                std::cout << std::endl;
//
//            }
//        }
//    }
    
       
     double  t3 = MPI_Wtime();
     Array<double>* hessian = ComputeHessianNew(P,ien_copy->getNrow(),comm);
     double  t4 = MPI_Wtime();
      double timing3 = t4-t3;
    double max_time3 = 0.0;
    MPI_Allreduce(&timing3, &max_time3, 1, MPI_DOUBLE, MPI_MAX, comm);
    if (world_rank == 0)
    {
        std::cout << "t_max3 := " << max_time3  << std::endl;
    }
    
    
    std::vector<Vert> LocalVs = P->getLocalVerts();
    int e = 0;
    double dudx = 0.0;
    double dudy = 0.0;
    double dudz = 0.0;
    int loc_v_id = 0;
    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
    std::vector<Vert> LVerts =  P->getLocalVerts();
    //Array<double>* hx=new Array<double>(ien_copy->getNrow(),ien_copy->getNcol());
    //Array<double>* hy=new Array<double>(ien_copy->getNrow(),ien_copy->getNcol());
    //Array<double>* hz=new Array<double>(ien_copy->getNrow(),ien_copy->getNcol());
    std::map<int,std::vector<double> > collect_drhodx;
    std::map<int,std::vector<double> > collect_drhody;
    std::map<int,std::vector<double> > collect_drhodz;
    int loc_v = 0;
    std::set<int> unique_verts;
    std::vector<int> unique_verts_vec;
    for(int i=0;i<ien_copy->getNrow();i++)
    {
        
        double drhodx_e = hessian->getVal(i,0);
        double drhody_e = hessian->getVal(i,1);
        double drhodz_e = hessian->getVal(i,2);
        //if(world_rank == 0)
        //{
	//    if(drhodx_e==0.0 || drhody_e == 0.0 || drhodz_e == 0.0)
	//    {
        //    std::cout << i << " " << drhodx_e << " " << drhody_e  << " " << drhodz_e << std::endl;
        //	}
        //}

        for(int j=0;j<8;j++)
        {
            loc_v = loc_elem2verts_loc[i][j];
	    collect_drhodx[loc_v].push_back(drhodx_e);
	    collect_drhody[loc_v].push_back(drhody_e);
	    collect_drhodz[loc_v].push_back(drhodz_e);
        }
    }
     
    std::map<int,std::vector<double> >::iterator it_rhos;
    
    double sumdx = 0;
    double sumdy = 0;
    double sumdz = 0;
    
    int c = 0;
    
    std::vector<double> hx;
    std::vector<double> hy;
    std::vector<double> hz;
    
    for(it_rhos=collect_drhodx.begin();it_rhos!=collect_drhodx.end();it_rhos++)
    {
        sumdx = 0;
        sumdy = 0;
        sumdz = 0;

        unique_verts_vec.push_back(it_rhos->first);

        for(int q = 0;q<it_rhos->second.size();q++)
        {
            sumdx = sumdx + it_rhos->second[q];
            sumdy = sumdy + collect_drhody[it_rhos->first][q];
            sumdz = sumdz + collect_drhodz[it_rhos->first][q];
        }
 
        hx.push_back(sumdx/it_rhos->second.size());
        hy.push_back(sumdy/it_rhos->second.size());
        hz.push_back(sumdz/it_rhos->second.size());
        c++;
    }
    
    //std::cout << hx.size() << " " << LVerts.size() << " " << unique_verts_vec.size() << std::endl;
    
    string filename = "quantity_rank_" + std::to_string(world_rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"drhox\", \"drhoy\", \"drhoz\"" << std::endl;
    int nvert = unique_verts_vec.size();
    myfile <<"ZONE N = " << nvert << ", E = " << ien_copy->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    for(int i=0;i<nvert;i++)
    {
       myfile << LVerts[unique_verts_vec[i]].x << "   " << LVerts[unique_verts_vec[i]].y << "   " << LVerts[unique_verts_vec[i]].z << "   " << hx[i] << " " << hy[i] << " " << hz[i] << std::endl;
    }

    for(int i=0;i<ien_copy->getNrow();i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    myfile.close();
     
    MPI_Finalize();
    
    return 0;
     
}
