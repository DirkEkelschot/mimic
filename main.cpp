#include "adapt_io.h"
#include "adapt_compute.h"
#include "adapt_operations.h"
#include "adapt_math.h"
#include <math.h>
#include "adapt_output.h"
#include "adapt_partition.h"
//#include "adapt_datatype.h"
#include "adapt_geometry.h"
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
    M[0] = 0.25;M[1]=-0.3;M[2]=0.4;
    M[3] = -0.3;M[4]=1.25;M[5]=0.1;
    M[6] = 0.4;M[7]=0.1;M[8]=0.25;
    
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
    int nloc = P->getLocalPartition()->getNloc(world_rank);
    int of = P->getLocalPartition()->getOffset(world_rank);
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
*/

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

Array<double>* ComputedUdXi(Partition* P, std::vector<double> U, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::vector<Vert> locVerts            = P->getLocalVerts();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    int offset  = P->getLocalPartition()->getOffset(rank);
    int* xadj   = P->getXadj();
    int* adjcny = P->getAdjcny();
    int nloc    = P->getLocalPartition()->getNrow();
    int lid;
    Array<double>* grad = new Array<double>(nloc,3);
    std::vector<std::vector<double> > store_coords;
    std::vector<double> upo_stored;
    std::vector<double> upijk_stored;
    std::vector<double> fac_stored;
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
        lid = gE2lE[i+offset];
        double u_ijk = U[lid];
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
        //double Volijk = ComputeVolumeHexCell(Pijk);           
        int tel   = 0;
        std::vector<double> tmp;
        if(i==1000)
        {
            std::cout << "adjEl_id ";
        }
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            if(i==1000)
            {
                std::cout << adjEl_id << " ";
            }
            
            lid = gE2lE[adjEl_id];
            double u_po = U[lid];
            double* Po = new double[np*3];

            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                Po[k*3+0] = locVerts[loc_vid].x;
                Po[k*3+1] = locVerts[loc_vid].y;
                Po[k*3+2] = locVerts[loc_vid].z;
            }


            Vert* Vpo = ComputeCenterCoord(Po,8);
            //double Volpo = ComputeVolumeHexCell(Po);

            double wi = sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                             (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                             (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));

            double fac = (wi);

            Vrt->setVal(tel,0,(Vpo->x-Vijk->x));
            Vrt->setVal(tel,1,(Vpo->y-Vijk->y));
            Vrt->setVal(tel,2,(Vpo->z-Vijk->z));
//
//          //b->setVal(tel,0,fac*(u_po-u_ijk));
            upo_stored.push_back(u_po);
            upijk_stored.push_back(u_ijk);
            fac_stored.push_back(fac);
            delete[] Po;
            delete Vpo;
            tel++;
        }
        if(i==1000)
        {
            std::cout << std::endl;
        }

        double min = *min_element(fac_stored.begin(), fac_stored.end());
        double max = *max_element(fac_stored.begin(), fac_stored.end());

        std::vector<double> fac_stored_update(nadj);

        for(int s=0;s<nadj;s++)
        {
            fac_stored_update[s] = (1.0/fac_stored[s]);
        }

        double* A_cm = new double[nadj*3];
        double sum =0.0;
        for(int s=0;s<nadj;s++)
        {
            for(int j=0;j<3;j++)
            {
                A_cm[j*nadj+s] = fac_stored_update[s]*Vrt->getVal(s,j);
            }
        }


        for(int s=0;s<nadj;s++)
        {
            b->setVal(s,0,fac_stored_update[s]*(upo_stored[s]-upijk_stored[s]));
        }
        Array<double>* x = SolveQR(A_cm,nadj,3,b);


        grad->setVal(i,0,x->getVal(0,0));
        grad->setVal(i,1,x->getVal(1,0));
        grad->setVal(i,2,x->getVal(2,0));


        upo_stored.clear();
        upijk_stored.clear();
        fac_stored.clear();
        fac_stored_update.clear();
        store_coords.clear();
        delete[] Pijk;
        vijkIDs.clear();
        delete Vrt_T;
        delete Vrt;
        delete b;
        delete[] A_cm;
        delete x;
        delete Vijk;
        
    }
    std::cout << std::endl;
     
      
    
    //delete[] xadj;
    //delete[] adjcny;
    //gE2lV.clear();
    //locVerts.clear();
    
    return grad;
}



Array<double>* ComputedUdXi_v2(Partition* P, std::vector<double> U, std::vector<double> Uvert, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::vector<Vert> locVerts            = P->getLocalVerts();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    int offset  = P->getLocalPartition()->getOffset(rank);
    int* xadj   = P->getXadj();
    int* adjcny = P->getAdjcny();
    int nloc    = P->getLocalPartition()->getNrow();
    int lid;
    Array<double>* grad = new Array<double>(nloc,3);
    std::vector<std::vector<double> > store_coords;
    std::vector<double> upo_stored;
    std::vector<double> upijk_stored;
    std::vector<double> fac_stored;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<double> uz;
    std::vector<double> uu;
    std::set<int> unique_verts;

    for(int i = 0;i<nloc;i++)
    {
        int start = xadj[i];
        int end   = xadj[i+1];
        
        
//        Array<double>* Vrt_T = new Array<double>(3,nadj);
//        Array<double>* Vrt   = new Array<double>(nadj,3);
//        Array<double>* b     = new Array<double>(nadj,1);
//
//        for(int q=0;q<nadj;q++)
//        {
//            for(int j=0;j<3;j++)
//            {
//                Vrt_T->setVal(j,q,0.0);
//                Vrt->setVal(q,j,0.0);
//            }
//        }
        lid = gE2lE[i+offset];
        double u_ijk = U[lid];
        std::vector<int> vijkIDs = gE2lV[i+offset];
        double* Pijk = new double[np*3];
        for(int k=0;k<vijkIDs.size();k++)
        {
            loc_vid     = vijkIDs[k];
            unique_verts.insert(loc_vid);
            ux.push_back(locVerts[loc_vid].x);
            uy.push_back(locVerts[loc_vid].y);
            uz.push_back(locVerts[loc_vid].z);
            uu.push_back(Uvert[loc_vid]);
            
            
        }
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            lid = gE2lE[adjEl_id];
           
            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                
                if(unique_verts.find(loc_vid)==unique_verts.end())
                {
                    unique_verts.insert(loc_vid);
                    ux.push_back(locVerts[loc_vid].x);
                    uy.push_back(locVerts[loc_vid].y);
                    uz.push_back(locVerts[loc_vid].z);
                    uu.push_back(Uvert[loc_vid]);
                }
            }
        }
        
        int nadj = uu.size();
        unique_verts.clear();
        uu.clear();
        ux.clear();
        uy.clear();
        uz.clear();
        
    }
     
      
    
    //delete[] xadj;
    //delete[] adjcny;
    //gE2lV.clear();
    //locVerts.clear();
    
    return grad;
}


/*
Array<double>* ComputeHessian(Partition* P, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int nloc = P->getLocalPartition()->getNrow();
    Array<double>* hessian = new Array<double>(Nel,3);
    
    std::map<int,std::vector<int> > gE2lV       = P->getGlobalElement2LocalVert();

    std::vector<Vert> locVerts                  = P->getLocalVerts();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    std::map<int,std::map<int,int> >::iterator itadj;

    int offset = P->getLocalPartition()->getOffset(rank);
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
        
        delete[] Pijk;
        
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
        delete x;
    }
    
    //delete[] xadj;
    //delete[] adjcny;
        
    return hessian;
}

*/


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


std::vector<std::vector<double> > ComputeDistances(Partition* P, ParallelState* pstate, ParArray<int>* iee, Array<int>* ifn, int Nel, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert> LocalVs = P->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    
    std::vector<int> ElemPart = P->getLocElem();
    
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;
    double* Pijk = new double[8*3];
    double* Padj = new double[8*3];
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    int offset = pstate->getOffset(world_rank);
    Vert* Vc = new Vert;
    for(int i=0;i<iee->getNrow();i++)
    {
        int elID = i+offset;//ElemPart[i];
        //std::vector<int> vijkIDs = gE2lV[elID];
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid].x;
            Pijk[k*3+1] = LocalVs[loc_vid].y;
            Pijk[k*3+2] = LocalVs[loc_vid].z;
        }
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        for(int j=0;j<6;j++)
        {
            adjID = iee->getVal(i,j);

            if(adjID<Nel)
            {
                //std::cout << gE2lV[adjID].size() << std::endl;
                for(int k=0;k<gE2lV[adjID].size();k++)
                {
                    loc_vid     = gE2lV[adjID][k];
                    Padj[k*3+0] = LocalVs[loc_vid].x;
                    Padj[k*3+1] = LocalVs[loc_vid].y;
                    Padj[k*3+2] = LocalVs[loc_vid].z;
                }
                Vert* Vadj = ComputeCenterCoord(Padj,8);

                d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                         (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                         (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));

                dist.push_back(d);
            }
            else
            {
                int fid = gE2gF[elID][j];

                std::vector<int> faceverts;
                //double* face_adj = new double[4*3];
                Vc->x = 0.0;Vc->y = 0.0;Vc->z = 0.0;
                for(int s=0;s<4;s++)
                {
                    int gvid = ifn->getVal(fid,s);
                    int lvid = gV2lV[gvid];

//                  face_adj[s*3+0] = LocalVs[lvid].x;
//                  face_adj[s*3+1] = LocalVs[lvid].y;
//                  face_adj[s*3+2] = LocalVs[lvid].z;
                    
                    Vc->x = Vc->x+LocalVs[lvid].x;
                    Vc->y = Vc->y+LocalVs[lvid].y;
                    Vc->z = Vc->z+LocalVs[lvid].z;
                }
//
                Vc->x = Vc->x/4.0;
                Vc->y = Vc->y/4.0;
                Vc->z = Vc->z/4.0;

                cou++;

                faceverts.clear();

                d = 2.0*sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                             (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                             (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

                dist.push_back(d);
            }
        }
        iee_dist.push_back(dist);
        dist.clear();
    }
    return iee_dist;
}


Array<double>* ComputedUdx(Partition* P, ParallelState* pstate, ParArray<int>* iee, std::map<int,std::vector<int> > iee_vec,std::map<int,std::vector<int> > ief_vec, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, Array<double>* bound, MPI_Comm comm, Array<int>* ife)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert> LocalVs = P->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::vector<int> ElemPart             = P->getLocElem();
    std::vector<int> Loc_Elem          = P->getLocElem();
    int* xadj                             = P->getXadj();
    int* adjcny                           = P->getAdjcny();
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;
    double* Pijk = new double[8*3];
    double* Padj = new double[8*3];
    
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    int offset = pstate->getOffset(world_rank);
    Vert* Vc = new Vert;
    int lid = 0;
    double u_ijk, u_po;

    int nLocElem = Loc_Elem.size();
    Array<double>* grad = new Array<double>(nLocElem,3);

    for(int i=0;i<nLocElem;i++)
    {
        int nadj = 6;
    
       // std::vector<int> adj_el_real;
//        for(int q=0;q<6;q++)
//        {
//            if(iee->getVal(i,q)<Nel)
//            {
//                adj_el_real.push_back(iee->getVal(i,q));
//                nadj++;
//            }
//        }
                
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        std::vector<std::vector<int> > MatVrt;
        
        for(int q=0;q<nadj;q++)
        {
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        int elID = Loc_Elem[i];
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid].x;
            Pijk[k*3+1] = LocalVs[loc_vid].y;
            Pijk[k*3+2] = LocalVs[loc_vid].z;
        }
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        
        lid = gE2lE[elID];
        u_ijk = U[lid];
        int t = 0;
        int start = xadj[i];
        int end   = xadj[i+1];
        for(int j=0;j<6;j++)
        {
            
            int fid = gE2gF[elID][j];
            adjID   = iee_vec[elID][j];
            //adjID = adjcny[start+j];
            int fid2 = ief_vec[elID][j];
//            if(world_rank == 0)
//            {
//            std::cout << "---" << world_rank  << " " << offset << " fid = " << fid << " -> (" << ife->getVal(fid,0) << " " << ife->getVal(fid,1) << ") --- (" << adjID << ", " << elID <<") ??? fid2 = " << fid2 << " -> (" << ife->getVal(fid2,0) << " " << ife->getVal(fid2,1) << ")" << std::endl;
//            }
            
            if(adjID<Nel)
            {
                lid = gE2lE[adjID];
                u_po = U[lid];

                for(int k=0;k<gE2lV[adjID].size();k++)
                {
                    loc_vid     = gE2lV[adjID][k];
                    Padj[k*3+0] = LocalVs[loc_vid].x;
                    Padj[k*3+1] = LocalVs[loc_vid].y;
                    Padj[k*3+2] = LocalVs[loc_vid].z;
                }
                
                Vert* Vadj = ComputeCenterCoord(Padj,8);

                d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                         (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                         (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));
                
                Vrt->setVal(t,0,(Vadj->x-Vijk->x));
                Vrt->setVal(t,1,(Vadj->y-Vijk->y));
                Vrt->setVal(t,2,(Vadj->z-Vijk->z));
                                
                b->setVal(t,0,u_po-u_ijk);
                delete Vadj;
                dist.push_back(d);
                t++;

                
            }
            
            else
            {
                int fid = gE2gF[elID][j];

                Vc->x = 0.0;Vc->y = 0.0;Vc->z = 0.0;
                
                for(int s=0;s<4;s++)
                {
                    int gvid = ifn->getVal(fid,s);
                    int lvid = gV2lV[gvid];

                    Vc->x = Vc->x+LocalVs[lvid].x;
                    Vc->y = Vc->y+LocalVs[lvid].y;
                    Vc->z = Vc->z+LocalVs[lvid].z;
                }
                
                Vc->x = Vc->x/4.0;
                Vc->y = Vc->y/4.0;
                Vc->z = Vc->z/4.0;

                d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                             (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                             (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

                u_po = ghost->getVal(adjID-Nel,0);
                
                double u_fpo = bound->getVal(adjID-Nel,0);
 
                Vrt->setVal(t,0,(Vc->x-Vijk->x));
                Vrt->setVal(t,1,(Vc->y-Vijk->y));
                Vrt->setVal(t,2,(Vc->z-Vijk->z));
//
                b->setVal(t,0,u_po-u_ijk);
                t++;
                dist.push_back(d);
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
            //std::cout << x->getVal(0,0) << " " << x->getVal(1,0) << " " << x->getVal(2,0) << std::endl;
            grad->setVal(i,0,x->getVal(0,0));
            grad->setVal(i,1,x->getVal(1,0));
            grad->setVal(i,2,x->getVal(2,0));
            
            delete[] A_cm;
            delete x;
             
             
        }
//        if(world_rank == 1)
//        {
//            std::cout << std::endl;
//        }
        //adj_el_real.clear();
        
        delete Vrt_T;
        delete Vrt;
        delete b;
        
        iee_dist.push_back(dist);
        dist.clear();
    }
    return grad;
}


Array<double>* ComputedUdx2(Partition* P, ParallelState* pstate, ParArray<int>* iee, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, Array<double>* bound, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert> LocalVs = P->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::vector<int> ElemPart             = P->getLocElem();
    
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;
    double* Pijk = new double[8*3];
    double* Padj = new double[8*3];
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    int offset = pstate->getOffset(world_rank);
    Vert* Vc = new Vert;
    int lid = 0;
    double u_ijk, u_po;

    
    Array<double>* grad = new Array<double>(iee->getNrow(),3);
    
    for(int i=0;i<iee->getNrow();i++)
    {
        int nadj = 6;
        
       // std::vector<int> adj_el_real;
//        for(int q=0;q<6;q++)
//        {
//            if(iee->getVal(i,q)<Nel)
//            {
//                adj_el_real.push_back(iee->getVal(i,q));
//                nadj++;
//            }
//        }
                
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        std::vector<std::vector<int> > MatVrt;
        
        for(int q=0;q<nadj;q++)
        {
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        int elID = i+offset;//ElemPart[i];
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid].x;
            Pijk[k*3+1] = LocalVs[loc_vid].y;
            Pijk[k*3+2] = LocalVs[loc_vid].z;
        }
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        lid = gE2lE[i+offset];
        u_ijk = U[lid];
        int t = 0;
        for(int j=0;j<6;j++)
        {
            adjID = iee->getVal(i,j);

            if(adjID<Nel)
            {
                lid = gE2lE[adjID];
                u_po = U[lid];

                for(int k=0;k<gE2lV[adjID].size();k++)
                {
                    loc_vid     = gE2lV[adjID][k];
                    Padj[k*3+0] = LocalVs[loc_vid].x;
                    Padj[k*3+1] = LocalVs[loc_vid].y;
                    Padj[k*3+2] = LocalVs[loc_vid].z;
                }
                
                Vert* Vadj = ComputeCenterCoord(Padj,8);

                d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                         (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                         (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));
                
                Vrt->setVal(t,0,(Vadj->x-Vijk->x));
                Vrt->setVal(t,1,(Vadj->y-Vijk->y));
                Vrt->setVal(t,2,(Vadj->z-Vijk->z));
                
                b->setVal(t,0,u_po-u_ijk);
                delete Vadj;
                dist.push_back(d);
                t++;
//                Vecb.push_back(u_po-u_ijk);
//                MatVrt.push_back(tmp);
                
            }
            else
            {
                //int fid = gE2gF[elID][j];
                int fid = ief->getVal(elID,j);
                //double* face_adj = new double[4*3];
                Vc->x = 0.0;Vc->y = 0.0;Vc->z = 0.0;
                
                for(int s=0;s<4;s++)
                {
                    int gvid = ifn->getVal(fid,s);
                    int lvid = gV2lV[gvid];
                    
                    Vc->x = Vc->x+LocalVs[lvid].x;
                    Vc->y = Vc->y+LocalVs[lvid].y;
                    Vc->z = Vc->z+LocalVs[lvid].z;
                }
                
                Vc->x = Vc->x/4.0;
                Vc->y = Vc->y/4.0;
                Vc->z = Vc->z/4.0;
                
                d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                             (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                             (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

                u_po = ghost->getVal(adjID-Nel,0);
                
                double u_fpo = bound->getVal(adjID-Nel,0);
 
                Vrt->setVal(t,0,(Vc->x-Vijk->x));
                Vrt->setVal(t,1,(Vc->y-Vijk->y));
                Vrt->setVal(t,2,(Vc->z-Vijk->z));
                b->setVal(t,0,u_po-u_ijk);
                t++;
                dist.push_back(d);
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
            
//            if(nadj<7)
//            {
//                std::cout << "NADJ = " << nadj << std::endl;
//                for(int s=0;s<nadj;s++)
//                {
//                    for(int j=0;j<3;j++)
//                    {
//                        std::cout << Vrt->getVal(s,j) << " ";
//                    }
//
//                    std::cout << std::endl;
//                }
//                std::cout << std::endl;
//            }
            
            
            Array<double>* x = SolveQR(A_cm,nadj,3,b);
        
            grad->setVal(i,0,x->getVal(0,0));
            grad->setVal(i,1,x->getVal(1,0));
            grad->setVal(i,2,x->getVal(2,0));
            
            delete[] A_cm;
            delete x;
        }
        
        //adj_el_real.clear();
        
        delete Vrt_T;
        delete Vrt;
        delete b;
        
        iee_dist.push_back(dist);
        dist.clear();
    }
    return grad;
}

std::vector<double> ReduceToVertices(Partition* P, Array<double>* Uelem, MPI_Comm comm)
{
    std::vector<double> Uelem_all = P->PartitionAuxilaryData(Uelem, comm);
    std::map<int,std::vector<double> > collect_Ui;
    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
    for(int i=0;i<Uelem_all.size();i++)
    {
        double uinew = Uelem_all[i];
        int loc_v;
        for(int j=0;j<8;j++)
        {
            loc_v = loc_elem2verts_loc[i][j];
            collect_Ui[loc_v].push_back(uinew);
        }
    }
    std::map<int,std::vector<double> >::iterator it_rhos;
    std::vector<double> uivert;
    
    for(it_rhos=collect_Ui.begin();it_rhos!=collect_Ui.end();it_rhos++)
    {
        double sum_u = 0;
        
        for(int q = 0;q<it_rhos->second.size();q++)
        {
            sum_u    = sum_u + it_rhos->second[q];
        }
        uivert.push_back(sum_u/it_rhos->second.size());
        
    }
        
    return uivert;
    
}


struct ModifiedGreenGauss
{
    Array<double>* elem_center;
    std::map<int, int> ief_map;
    std::map<int,vector<Vec3D*> > face_n;
    std::map<int,vector<Vec3D*> > face_r;
    std::map<int,vector<Vec3D*> > xfxc;
    std::map<int,vector<double> > dr;
    std::map<int,double> vol;
    std::map<int,vector<double> > dS;
    
};


//void CheckSurfaceOrientation(double* Plane)
//{
//    double*E0 = new double[2];
//}


ModifiedGreenGauss* ComputeMofiedGreenGaussData(Partition* Pa, ParArray<int>* iee, std::map<int,std::vector<int> > iee_vec, Array<int>* ifn, ParArray<int>* ief, std::map<int,std::vector<int> > ief_vec, Array<int>* ife, Array<double>* ghost, std::map<int,double> U, MPI_Comm comm)
{
    int nlocElem, start, end, offset, nloc, np, loc_vid, size, rank, lid;
    int vf0, vf1, vf2, vf3, vf4, vf5, vf6, vf7, fid;
    double wi, ds0, ds1 ,ds2, ds3, ds4, ds5, u_po,orient0,orient1,orient2,orient3,orient4,orient5,L0,L1,L2,L3,L4,L5;
    
    double* v0=new double[3];
    double* v1=new double[3];
    double* v2=new double[3];
    int Nel = iee->getNglob();
    
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> > gE2lV = Pa->getGlobElem2LocVerts();
    std::vector<Vert> locVerts            = Pa->getLocalVerts();
    offset                                = Pa->getLocalPartition()->getOffset(rank);
    nloc                                  = Pa->getLocalPartition()->getNloc(rank);
    std::map<int,int> gE2lE               = Pa->getGlobalElement2LocalElement();
    std::map<int,int> lE2gE               = Pa->getLocalElement2GlobalElement();
    int* xadj                             = Pa->getXadj();
    int* adjcny                           = Pa->getAdjcny();
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    np                                    = 8;
    double* Pijk                          = new double[np*3];
    std::vector<int> Loc_Elem             = Pa->getLocElem();
    int nLocElem                          = Loc_Elem.size();
    Array<double>* cc                     = new Array<double>(nLocElem,3);
    std::map<int,int> gV2lV               = Pa->getGlobalVert2LocalVert();
    
    std::vector<int> ElemPart               = Pa->getLocElem();
    int nLocElemInc = ElemPart.size();

    std::vector<int> vijkIDs;
    std::map<int,vector<Vec3D*> > normals;
    std::map<int,vector<Vec3D*> > rvector;
    std::map<int,vector<Vec3D*> > dxfxc;
    std::map<int,double> vol;
    std::map<int,vector<double> > dr;
    std::map<int,vector<double> > dS;
    double* Po  = new double[np*3];
    int tel     = 0;
    Vert* fc0 = new Vert;
    Vert* fc1 = new Vert;
    Vert* fc2 = new Vert;
    Vert* fc3 = new Vert;
    Vert* fc4 = new Vert;
    Vert* fc5 = new Vert;
        
    ModifiedGreenGauss* mggData = new ModifiedGreenGauss;
    std::vector<Vert*> face;
    for(int i=0;i<nLocElem;i++)
    {
        int gEl = Loc_Elem[i];

        vijkIDs = gE2lV[gEl];

        for(int k=0;k<vijkIDs.size();k++)
        {
           loc_vid     = vijkIDs[k];
           Pijk[k*3+0] = locVerts[loc_vid].x;
           Pijk[k*3+1] = locVerts[loc_vid].y;
           Pijk[k*3+2] = locVerts[loc_vid].z;
        }

        Vert* Vijk = ComputeCenterCoord(Pijk, np);
        double Vol = ComputeVolumeHexCell(Pijk);
        vol[gEl]   = Vol;
        
        for(int s=0;s<6;s++)
        {
            int faceid = ief_vec[gEl][s];

            Vert* Vface = new Vert;
            
            for(int r=0;r<4;r++)
            {
                int gvid = ifn->getVal(faceid,r);
                int lvid = gV2lV[gvid];
                
                Vert* V = new Vert;
                V->x    = locVerts[lvid].x;
                V->y    = locVerts[lvid].y;
                V->z    = locVerts[lvid].z;
                
                Vface->x = Vface->x+locVerts[lvid].x;
                Vface->y = Vface->y+locVerts[lvid].y;
                Vface->z = Vface->z+locVerts[lvid].z;
                
                face.push_back(V);
                
            }

            Vface->x = Vface->x/4.0;
            Vface->y = Vface->y/4.0;
            Vface->z = Vface->z/4.0;
            
            L0 = sqrt((Vface->x-Vijk->x)*(Vface->x-Vijk->x)
                     +(Vface->y-Vijk->y)*(Vface->y-Vijk->y)
                     +(Vface->z-Vijk->z)*(Vface->z-Vijk->z));
            
            Vec3D* r0 = new Vec3D;
            
            r0->c0 = (Vface->x-Vijk->x);
            r0->c1 = (Vface->y-Vijk->y);
            r0->c2 = (Vface->z-Vijk->z);
            
            v0[0] = face[1]->x-face[0]->x;
            v0[1] = face[1]->y-face[0]->y;
            v0[2] = face[1]->z-face[0]->z;

            v1[0] = face[3]->x-face[0]->x;
            v1[1] = face[3]->y-face[0]->y;
            v1[2] = face[3]->z-face[0]->z;
            
            Vec3D* n0 = ComputeSurfaceNormal(v0,v1);
            orient0   = DotVec3D(r0,n0);
            
            if(orient0<0.0)
            {
                NegateVec3D(n0);
            }
            
            double orientaft = DotVec3D(r0,n0);
            
            ds0 = ComputeSurfaceArea(v0,v1);
            dS[gEl].push_back(ds0);
            normals[gEl].push_back(n0);
            dxfxc[gEl].push_back(r0);
            face.clear();
        }
        
        tel = 0;
        start = xadj[i];
        end = xadj[i+1];
       
        for(int j=0;j<6;j++)
        {
           int fid   = ief_vec[gEl][j];
           int adjID = iee_vec[gEl][j];

           if(adjID<Nel)// If internal element;
           {
               lid  = gE2lE[adjID];
               u_po = U[adjID];
               
               for(int k=0;k<gE2lV[adjID].size();k++)
               {
                   loc_vid     = gE2lV[adjID][k];
                   Po[k*3+0] = locVerts[loc_vid].x;
                   Po[k*3+1] = locVerts[loc_vid].y;
                   Po[k*3+2] = locVerts[loc_vid].z;
               }
               
               Vert* Vpo = ComputeCenterCoord(Po,8);

               double d = sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                               (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                               (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));
            
               Vec3D* rf = new Vec3D;
               rf->c0    = (Vpo->x-Vijk->x)/d;
               rf->c1    = (Vpo->y-Vijk->y)/d;
               rf->c2    = (Vpo->z-Vijk->z)/d;
               
               rvector[gEl].push_back(rf);
               dr[gEl].push_back(d);
           }
           else // If boundary face then search data in the correct ghost cell;
           {
               fid = ief_vec[gEl][j];

               //double* face_adj = new double[4*3];
               Vert* Vpo = new Vert;
               Vpo->x = 0.0;Vpo->y = 0.0;Vpo->z = 0.0;
               for(int s=0;s<4;s++)
               {
                   int gvid = ifn->getVal(fid,s);
                   int lvid = gV2lV[gvid];

                   Vpo->x = Vpo->x+locVerts[lvid].x;
                   Vpo->y = Vpo->y+locVerts[lvid].y;
                   Vpo->z = Vpo->z+locVerts[lvid].z;
               }

               Vpo->x = Vpo->x/4.0;
               Vpo->y = Vpo->y/4.0;
               Vpo->z = Vpo->z/4.0;

               double d = 2.0*sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                            (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                            (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));
               
               u_po = ghost->getVal(adjID-Nel,0);
               
               Vec3D* rf = new Vec3D;
               rf->c0    = (Vpo->x-Vijk->x)/d;
               rf->c1    = (Vpo->y-Vijk->y)/d;
               rf->c2    = (Vpo->z-Vijk->z)/d;

               rvector[gEl].push_back(rf);
               dr[gEl].push_back(d);
               
               
               delete Vpo;
           }
           tel++;
        }
    }
    

    mggData->elem_center  = cc;
    mggData->face_n = normals;
    mggData->face_r = rvector;
    mggData->xfxc = dxfxc;
    mggData->dr = dr;
    mggData->dS = dS;
    mggData->vol = vol;
    return mggData;
}


Array<double>* ReconstructGradient(Partition* Pa, std::map<int,std::vector<int> > iee_vec, std::map<int,std::vector<int> > ief_vec, std::map<int,double> U, ModifiedGreenGauss* mgg, Array<double>* ghost, Array<int>* ife,MPI_Comm comm)
{
    int lid, gEl, adjID, l_adjid, size, rank;
    double u_c, u_nb, gu_c_vx, gu_c_vy, gu_c_vz, gu_nb_vx, gu_nb_vy, gu_nb_vz,sum_phix,sum_phiy,sum_phiz,dphi_dn,Vol;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    int Nel = Pa->getLocalPartition()->getNglob();
    
    std::map<int,int> gE2lE                 = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem               = Pa->getLocElem();
    int nLoc_Elem                           = Loc_Elem.size();
    
    Array<double>* gu_c_x      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_y      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_z      = new Array<double>(nLoc_Elem,1);
    
    Array<double>* gu_c_old    = new Array<double>(nLoc_Elem,3);
    
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    
    for(int i=0;i<nLoc_Elem;i++)
    {
        gu_c_x->setVal(i,0,0.0);
        gu_c_y->setVal(i,0,0.0);
        gu_c_z->setVal(i,0,0.0);
        
        gu_c_old->setVal(i,0,0.0);
        gu_c_old->setVal(i,1,0.0);
        gu_c_old->setVal(i,2,0.0);
    }
    
    std::map<int,vector<Vec3D*> > normals   = mgg->face_n;
    std::map<int,vector<Vec3D*> > rvector   = mgg->face_r;
    std::map<int,vector<Vec3D*> > dxfxc     = mgg->xfxc;
    std::map<int,vector<double> > dS        = mgg->dS;
    std::map<int,vector<double> > dr        = mgg->dr;
    std::map<int,double > vol               = mgg->vol;
    int it = 0;
    double alpha   = 0.0;
    double L2normx = 0.0;
    double L2normy = 0.0;
    double L2normz = 0.0;
    std::vector<Vec3D*> n_grads;


    for(int it=0;it<10;it++)
    {
        for(int i=0;i<nLoc_Elem;i++)
        {
            gu_c_old->setVal(i,0,gu_c_x->getVal(i,0));
            gu_c_old->setVal(i,1,gu_c_y->getVal(i,0));
            gu_c_old->setVal(i,2,gu_c_z->getVal(i,0));
        }
        
         //communicate grad phi!!!
        
        std::map<int,double> dUdx_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_x, comm);
        std::map<int,double> dUdy_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_y, comm);
        std::map<int,double> dUdz_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_z, comm);

        for(int i=0;i<nLoc_Elem;i++)
         {
             gEl = Loc_Elem[i];
             lid = gE2lE[gEl];
             u_c = U[gEl];

             gu_c_vx = gu_c_x->getVal(lid,0);
             gu_c_vy = gu_c_y->getVal(lid,0);
             gu_c_vz = gu_c_z->getVal(lid,0);

             sum_phix = 0.0;
             sum_phiy = 0.0;
             sum_phiz = 0.0;
             
             for(int j=0;j<6;j++)
             {
                 adjID   = iee_vec[gEl][j];
                 int fid = ief_vec[gEl][j];
                 
                 if(adjID<Nel)
                 {
                     //l_adjid = gE2lE[adjID];
                     gu_nb_vx = dUdx_p_bnd[adjID];
                     gu_nb_vy = dUdy_p_bnd[adjID];
                     gu_nb_vz = dUdz_p_bnd[adjID];
                     u_nb = U[adjID];
                 }
                 else
                 {
                     u_nb     = ghost->getVal(adjID-Nel,0);
                     
                     gu_nb_vx = gu_c_x->getVal(lid,0);
                     gu_nb_vy = gu_c_y->getVal(lid,0);
                     gu_nb_vz = gu_c_z->getVal(lid,0);
                     
                 }
                 
                 Vec3D* nj          = normals[gEl][j];
                 Vec3D* rj          = rvector[gEl][j];
                 
                 double alpha     = DotVec3D(nj,rj);
                 
                 Vec3D* nf_m_arf    = new Vec3D;

                 nf_m_arf->c0=nj->c0-alpha*rj->c0;
                 nf_m_arf->c1=nj->c1-alpha*rj->c1;
                 nf_m_arf->c2=nj->c2-alpha*rj->c2;
                 //std::cout << alpha << std::endl;
                 dphi_dn = alpha * (u_nb - u_c)/dr[gEl][j]; +  0.5 * ((gu_nb_vx + gu_c_vx) * nf_m_arf->c0 +  (gu_nb_vy + gu_c_vy) * nf_m_arf->c1 +  (gu_nb_vz + gu_c_vz) * nf_m_arf->c2);
                 
                 sum_phix = sum_phix+dphi_dn*dxfxc[gEl][j]->c0*dS[gEl][j];
                 sum_phiy = sum_phiy+dphi_dn*dxfxc[gEl][j]->c1*dS[gEl][j];
                 sum_phiz = sum_phiz+dphi_dn*dxfxc[gEl][j]->c2*dS[gEl][j];
             }
             
             Vol = vol[gEl];
             gu_c_x->setVal(i,0,1.0/Vol*sum_phix);
             gu_c_y->setVal(i,0,1.0/Vol*sum_phiy);
             gu_c_z->setVal(i,0,1.0/Vol*sum_phiz);
         }
        dUdx_p_bnd.clear();
        dUdy_p_bnd.clear();
        dUdz_p_bnd.clear();
        L2normx = 0.0;
        L2normy = 0.0;
        L2normz = 0.0;
        for(int n=0;n<nLoc_Elem;n++)
        {
            L2normx = L2normx+ sqrt((gu_c_x->getVal(n,0)-gu_c_old->getVal(n,0))*(gu_c_x->getVal(n,0)-gu_c_old->getVal(n,0)));

            L2normy = L2normy+ sqrt((gu_c_y->getVal(n,0)-gu_c_old->getVal(n,1))*(gu_c_y->getVal(n,0)-gu_c_old->getVal(n,1)));

            L2normz = L2normz+ sqrt((gu_c_z->getVal(n,0)-gu_c_old->getVal(n,2))*(gu_c_z->getVal(n,0)-gu_c_old->getVal(n,2)));

        }
        std::cout << it << " -> convergence " <<  L2normx <<","<<L2normy<<","<<L2normz<<","<< std::endl;
        if(L2normx<1.0e-07 && L2normy<1.0e-07 && L2normz<1.0e-07)
        {
            break;
        }
        
    }
    
    return gu_c_old;
    
}


void ComputeFaceValues(Partition* P, Array<double>* U, MPI_Comm comm)
{
    int nface, start, end, rank, size;
    int nloc, offset, adjEl_id, leid, geid, i, t;
    double u_o,u_c;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    offset  = P->getLocalPartition()->getOffset(rank);
    int* xadj   = P->getXadj();
    int* adjcny = P->getAdjcny();
    nloc    = P->getLocalPartition()->getNrow();
    
    nface = 6; // # hardcoded for hexes for now
    
    std::vector<double> Uelem_all         = P->PartitionAuxilaryData(U, comm);
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::map<int,int> lE2gE               = P->getLocalElement2GlobalElement();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    Array<double>* Uf                     = new Array<double>(nloc,nface);
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();

    
    
    for(i=0;i<nloc;i++)
    {
        start  = xadj[i];
        end    = xadj[i+1];
        geid   = lE2gE[i];
        u_c    = U->getVal(i,0);
        t = 0;
        for(int j=start;j<end;j++)
        {
            adjEl_id = adjcny[j];
            leid     = gE2lE[adjEl_id];
            u_o      = Uelem_all[leid];
            
            Uf->setVal(i,t,u_c-u_o);
            
            t++;
        }
    }
}


Array<double>* ComputeVolumes(Partition* Pa)
{
    int loc_vid;
    std::map<int,int> gE2lE                 = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem               = Pa->getLocElem();
    int nLocElem                            = Loc_Elem.size();
    std::map<int,std::vector<int> > gE2lV   = Pa->getGlobElem2LocVerts();
    std::vector<Vert> locVerts              = Pa->getLocalVerts();
    Array<double>* Volumes = new Array<double>(nLocElem,1);
    std::vector<int> vijkIDs;
    double* Pijk = new double[8*3];
    for(int i=0;i<nLocElem;i++)
    {
       int gEl = Loc_Elem[i];

       vijkIDs = gE2lV[gEl];

       for(int k=0;k<vijkIDs.size();k++)
       {
          loc_vid     = vijkIDs[k];
          Pijk[k*3+0] = locVerts[loc_vid].x;
          Pijk[k*3+1] = locVerts[loc_vid].y;
          Pijk[k*3+2] = locVerts[loc_vid].z;
       }

       double Vol = ComputeVolumeHexCell(Pijk);
       Volumes->setVal(i,0,Vol);
    }
    return Volumes;
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
    const char* fn_conn="cases/cylinder/anisotropic_16k/conn_aniso_16k.h5";
    const char* fn_grid="cases/cylinder/anisotropic_16k/grid_aniso_16k.h5";
    const char* fn_data="cases/cylinder/anisotropic_16k/data_aniso_16k.h5";
    const char* fn_adept="cases/cylinder/anisotropic_16k/conn_aniso_16k.h5";
    
    Array<int>*    zdefs = ReadDataSetFromGroupFromFile<int>(fn_adept,"zones","zdefs");
    Array<char>*  znames = ReadDataSetFromGroupFromFile<char>(fn_adept,"zones","znames");
    ParArray<int>* ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    ParArray<int>* iee = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);
    ParArray<double>* xcn_def = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);
    Array<int>* ife = ReadDataSetFromFile<int>(fn_conn,"ife");
    //ParArray<double>* xcn_def = ReadDataSetFromFileInParallel<double>(fn_data,"xcn",comm,info);
    //ParArray<double>* ifn = ReadDataSetFromFileInParallel<double>(fn_grid,"ifn",comm,info);
    Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
    int nFglob = ifn->getNrow();
    //ParArray<double>* boundaries = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","boundaries",comm,info);
    PlotBoundaryData(znames,zdefs,comm);
    
    int Nel = ien->getNglob();
    int Nel_part = ien->getNrow();
    //ParArray<double>* xcn_def        = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","xcn",Nel,0,comm,info);
    ParArray<double>* interior = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
    Array<double>* ghost = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);
    Array<double>* bound = ReadDataSetFromRunInFile<double>(fn_data,"run_1","boundaries");
    
    std::cout << "boundary faces = " << bound->getNrow() << " " << ghost->getNrow() << std::endl;
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
    delete ien;
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
    delete ief;
    
    int ncol_iee = iee->getNcol()-1;
    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
    
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol_iee;j++)
        {
            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
        }
    }
    delete iee;
    
    int nrow_ifn = ifn->getNrow();
    int ncol_ifn = 4;
    Array<int>* ifn_copy = new Array<int>(nFglob,ncol_ifn);
    for(int i=0;i<nrow_ifn;i++)
    {
        for(int j=0;j<ncol_ifn;j++)
        {
            ifn_copy->setVal(i,j,ifn->getVal(i,j+1)-1);
        }
    }
    delete ifn;
    
    int nrow_ife = ife->getNrow();
    int ncol_ife = 2;
    Array<int>* ife_copy = new Array<int>(nrow_ife,ncol_ife);
    for(int i=0;i<nrow_ife;i++)
    {
        for(int j=0;j<ncol_ife;j++)
        {
            ife_copy->setVal(i,j,ife->getVal(i,j)-1);
        }
    }
    delete ife;

    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(ien_copy,comm,8);

    ParallelState* xcn_parstate = new ParallelState(xcn_def->getNglob(),comm);
    Array<double>* Uivar = new Array<double>(Nel_part,1);
    for(int i=0;i<Nel_part;i++)
    {
        Uivar->setVal(i,0,interior->getVal(i,0));
    }
    delete interior;
   
    Partition* P = new Partition(ien_copy, iee_copy, ief_copy, parmetis_pstate,pstate, xcn_def,xcn_parstate,Uivar,comm);
    
    //std::map<int,std::vector<int> > iee_loc = P->getElement2EntityPerPartition(iee_copy,comm);
    
//    std::map<int,std::vector<int> >::iterator it;
//    for(it=iee_loc.begin();it!=iee_loc.end();it++)
//    {
//        std::cout << it->first << " ";
//        for(int q=0;q<it->second.size();q++)
//        {
//            std::cout  << it->second[q] << " ";
//        }
//        std::cout << std::endl;
//    }
    
    
    std::vector<double> Uaux = P->PartitionAuxilaryData(Uivar, comm);
    
    std::map<int,std::vector<int> > iee_loc = P->getElement2EntityPerPartition(iee_copy,comm);
    std::map<int,std::vector<int> > ief_loc = P->getElement2EntityPerPartition(ief_copy,comm);
    
    Array<double>* dUdXi_v2 = ComputedUdx(P, pstate, iee_copy, iee_loc, ief_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm, ife_copy);

   
//    if(world_rank == 1)
//    {
//        std::map<int,std::vector<int> >::iterator itmp;
//        for(itmp=iee_loc.begin();itmp!=iee_loc.end();itmp++)
//        {
//            std::cout << "rank = " << world_rank << " "  << itmp->first << std::endl;
//        }
//    }
    
    std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uivar,comm);

    ModifiedGreenGauss* mggData = ComputeMofiedGreenGaussData(P,iee_copy,iee_loc,ifn_copy,ief_copy,ief_loc,ife_copy,ghost,UauxNew,comm);
////
   Array<double>* dUdXi = ReconstructGradient(P,iee_loc,ief_loc,UauxNew,mggData,ghost,ife_copy,comm);
//
    Array<double>* dUidxi = new Array<double>(dUdXi->getNrow(),1);
    Array<double>* dUidyi = new Array<double>(dUdXi->getNrow(),1);
    Array<double>* dUidzi = new Array<double>(dUdXi->getNrow(),1);

    for(int i=0;i<dUdXi->getNrow();i++)
    {
        dUidxi->setVal(i,0,dUdXi->getVal(i,0));
        dUidyi->setVal(i,0,dUdXi->getVal(i,1));
        dUidzi->setVal(i,0,dUdXi->getVal(i,2));
    }
//    std::vector<double> dUdxaux = P->PartitionAuxilaryData(dUidxi, comm);
//
//    Array<double>* dU2dXi2 = ReconstructGradient(P,iee_loc,dUdxaux,mggData,ghost,ife_copy);
//
//    Array<double>* dudx = new Array<double>(dUdXi->getNrow(),1);
//    Array<double>* dudy = new Array<double>(dUdXi->getNrow(),1);
//    Array<double>* dudz = new Array<double>(dUdXi->getNrow(),1);
//
//    Array<double>* du2dx2 = new Array<double>(dUdXi->getNrow(),1);
//    Array<double>* du2dxy = new Array<double>(dUdXi->getNrow(),1);
//    Array<double>* du2dxz = new Array<double>(dUdXi->getNrow(),1);
//
//    for(int i=0;i<dU2dXi2->getNrow();i++)
//    {
//        dudx->setVal(i,0,dUdXi->getVal(i,0));
//        dudy->setVal(i,0,dUdXi->getVal(i,1));
//        dudz->setVal(i,0,dUdXi->getVal(i,2));
//
//        du2dx2->setVal(i,0,dU2dXi2->getVal(i,0));
//        du2dxy->setVal(i,0,dU2dXi2->getVal(i,1));
//        du2dxz->setVal(i,0,dU2dXi2->getVal(i,2));
//    }
////
////
////
//    int nloc_e = P->getNlocElem();
//
//    int nElemLoc = P->getNlocElem();
//    //std::vector<double> Uaux = P->PartitionAuxilaryData(Uivar, comm);
//
//    // ===================================================================
//    // ======================TEST AUX PARTITIONING========================
//    // ===================================================================
//    //std::vector<double> Uaux = P->PartitionAuxilaryData(Ui, comm);
////    std::vector<double> Uelem = P->getUelem();
////    //std::cout << "Are sizes equal? -> " << Uaux.size() << " " << Uelem.size() << std::endl;
////    for(int i=0;i<Uaux.size();i++)
////    {
////        if(Uaux[i] - Uelem[i] > 1.0e-09)
////        {
////            std::cout << Uaux[i] - Uelem[i] << std::endl;
////        }
////    }
////    std::cout << "Are sizes equal? -> " << Uaux.size() << " " << Uelem.size() << std::endl;
//    // ===================================================================
//    // ===================================================================
//    // ===================================================================
//
//    //std::vector<std::vector<double> > dist = ComputeDistances(P, pstate, iee_copy, ifn_copy, Nel, comm);
//
//
////    int offset = pstate->getOffsets()[world_rank];
////
////    //Array<double>* dUdXi = ComputedUdx(P, pstate, iee_copy, iee_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm);
////    //Array<double>* dUdXi = ComputedUdx(P, pstate, iee_loc, iee_copy, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm);
////    //Array<double>* dUdXi = ComputedUdx2(P, pstate, iee_copy, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm);
////
////
////

//
    
    Array<double>* dudx_v2 = new Array<double>(dUdXi_v2->getNrow(),1);
    Array<double>* dudy_v2 = new Array<double>(dUdXi_v2->getNrow(),1);
    Array<double>* dudz_v2 = new Array<double>(dUdXi_v2->getNrow(),1);

    for(int i=0;i<dUdXi_v2->getNrow();i++)
    {
        dudx_v2->setVal(i,0,dUdXi_v2->getVal(i,0));
        dudy_v2->setVal(i,0,dUdXi_v2->getVal(i,1));
        dudz_v2->setVal(i,0,dUdXi_v2->getVal(i,2));

    }
    std::vector<double> dUdxaux_old = P->PartitionAuxilaryData(dudx_v2, comm);
//
//    Array<double>* dU2dXi2_old = ComputedUdx(P, pstate, iee_copy, iee_loc, ifn_copy, ief_copy, Nel, dUdxaux_old, ghost, bound, comm);
////
//    Array<double>* du2dx2_old = new Array<double>(dU2dXi2_old->getNrow(),1);
//    Array<double>* du2dxy_old = new Array<double>(dU2dXi2_old->getNrow(),1);
//    Array<double>* du2dxz_old = new Array<double>(dU2dXi2_old->getNrow(),1);
//
//    for(int i=0;i<dU2dXi2_old->getNrow();i++)
//    {
//        du2dx2_old->setVal(i,0,dU2dXi2_old->getVal(i,0));
//        du2dxy_old->setVal(i,0,dU2dXi2_old->getVal(i,1));
//        du2dxz_old->setVal(i,0,dU2dXi2_old->getVal(i,2));
//    }
////
    std::vector<double> u_v = ReduceToVertices(P,Uivar,comm);
////
    std::vector<double> dudx_v = ReduceToVertices(P,dUidxi,comm);
    std::vector<double> dudy_v = ReduceToVertices(P,dUidyi,comm);
    std::vector<double> dudz_v = ReduceToVertices(P,dUidzi,comm);
//
//    std::vector<double> du2dx2_v = ReduceToVertices(P,du2dx2,comm);
//    std::vector<double> du2dxy_v = ReduceToVertices(P,du2dxy,comm);
//    std::vector<double> du2dxz_v = ReduceToVertices(P,du2dxz,comm);
//
//    std::vector<double> du2dx2_v_old = ReduceToVertices(P,du2dx2_old,comm);
////
//
//////
    std::vector<double> dudx_v_v2 = ReduceToVertices(P,dudx_v2,comm);
    std::vector<double> dudy_v_v2 = ReduceToVertices(P,dudy_v2,comm);
    std::vector<double> dudz_v_v2 = ReduceToVertices(P,dudz_v2,comm);
////
    std::vector<Vert> Verts =  P->getLocalVerts();
    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
//
//    //std::cout << "gradsizes " << dUdXi->getNrow() << " " << dUdXi_v2->getNrow() << " " << Verts.size() << std::endl;
//    //=============================================================================
//    //=====================Output the data in Tecplot format=======================
//    //=============================================================================
        string filename = "gradients_" + std::to_string(world_rank) + ".dat";
        ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
        //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"dx\", \"dy\", \"dz\", \"dx_old\", \"dy_old\", \"dz_old\", \"dx_v2\", \"dy_v2\", \"dz_v2\", \"dx2_old\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"dx\", \"dy\", \"dz\", \"dx_old\", \"dy_old\", \"dz_old\"" << std::endl;
        int nvert = Verts.size();
        myfile <<"ZONE N = " << nvert << ", E = " << ien_copy->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;

        for(int i=0;i<nvert;i++)
        {
            myfile<< Verts[i].x << " " << Verts[i].y << " " << Verts[i].z << " " << u_v[i] << " "<< dudx_v[i] << " " << dudy_v[i] << " " << dudz_v[i] << " " << dudx_v_v2[i] << " " << dudy_v_v2[i] << " " << dudz_v_v2[i] << std::endl;
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
//
     
    //=============================================================================
    //=====================Output the data in Tecplot format=======================
    //=============================================================================
    
    
    
//    Array<double>* dUdXi = ComputedUdXi(P,Uaux,ien_copy->getNglob(),comm);
//    //Array<double>* dUdXi = ComputedUdXi_v2(P,Uaux,Uvert,ien_copy->getNglob(),comm);
//
//
//
//
//
//    double  t6 = MPI_Wtime();
//    //P->PartitionAuxilaryData(grad, comm);
//    double  t7 = MPI_Wtime();
//    double timing6 = t7-t6;
//    double max_time6 = 0.0;
//    MPI_Allreduce(&timing6, &max_time6, 1, MPI_DOUBLE, MPI_MAX, comm);
//    if (world_rank == 0)
//    {
//      std::cout << "t_max6 := " << max_time6  << std::endl;
//
//    }
//
//
//
////    if (world_rank == 0)
////    {
////        std::cout << max_time_part << " " << max_time_grad << " " << max_time_aux << std::endl;
////    }
//
//    //==============================================================================================
//    //================================Compute Eigenvalues/vectors on verts for output===============
//    //==============================================================================================
//
//    std::map<int,std::vector<double> > collect_d2Udxix;
//    std::map<int,std::vector<double> > collect_d2Udxiy;
//    std::map<int,std::vector<double> > collect_d2Udxiz;
//    std::map<int,std::vector<double> > collect_dUidX;
//    std::map<int,std::vector<double> > collect_dUidY;
//    std::map<int,std::vector<double> > collect_dUidZ;
//    std::map<int,std::vector<double> > collect_dU2idX2;
//    std::map<int,std::vector<double> > collect_dU2idXY;
//    std::map<int,std::vector<double> > collect_dU2idXZ;
//
//    std::map<int,std::vector<double> > collect_dU2idYX;
//    std::map<int,std::vector<double> > collect_dU2idY2;
//    std::map<int,std::vector<double> > collect_dU2idYZ;
//
//    std::map<int,std::vector<double> > collect_dU2idZX;
//    std::map<int,std::vector<double> > collect_dU2idZY;
//    std::map<int,std::vector<double> > collect_dU2idZ2;
//
//    std::map<int,std::vector<double> > collect_Ui;
//    std::map<int,std::vector<double> > collect_Vi;
//    std::map<int,std::vector<double> > collect_Di;
//    std::vector<double> dUdxiaux;
//    std::vector<double> dU2dxixaux;
//    std::vector<double> dU2dxiyaux;
//    std::vector<double> dU2dxizaux;
////    Array<double>* d2Udxix = new Array<double>(Nel_part,1);
////    Array<double>* d2Udxiy = new Array<double>(Nel_part,1);
////    Array<double>* d2Udxiz = new Array<double>(Nel_part,1);
//    std::map<int,std::vector<double> >::iterator it_rhos;
//    std::map<int,std::vector<double> >::iterator it_vi;
//    std::map<int,std::vector<double> >::iterator it_di;
//    std::vector<std::vector<double> > Hess;
//    std::vector<double> d2Udxix;
//    std::vector<double> d2Udxiy;
//    std::vector<double> d2Udxiz;
//    std::vector<std::vector<double> > Hess_El;
//    std::vector<Vert> LocalVs = P->getLocalVerts();
//    int e = 0;
//    double dudx = 0.0;
//    double dudy = 0.0;
//    double dudz = 0.0;
//    int loc_v_id = 0;
//    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
//    std::vector<Vert> LVerts =  P->getLocalVerts();
//
//    int loc_v = 0;
//
//    Array<double>* dUidxi = new Array<double>(Nel_part,1);
//    Array<double>* dUidyi = new Array<double>(Nel_part,1);
//    Array<double>* dUidzi = new Array<double>(Nel_part,1);
//
//    Array<double>* dU2idx2i_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idxyi_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idxzi_new = new Array<double>(Nel_part,1);
//
//    Array<double>* dU2idyxi_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idy2i_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idyzi_new = new Array<double>(Nel_part,1);
//
//    Array<double>* dU2idzxi_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idzyi_new = new Array<double>(Nel_part,1);
//    Array<double>* dU2idz2i_new = new Array<double>(Nel_part,1);
//
//
//    for(int i=0;i<Nel_part;i++)
//    {
//        dUidxi->setVal(i,0,dUdXi->getVal(i,0));
//        dUidyi->setVal(i,1,dUdXi->getVal(i,1));
//        dUidzi->setVal(i,2,dUdXi->getVal(i,2));
//    }
//
//    std::vector<double> dUdXaux = P->PartitionAuxilaryData(dUidxi, comm);
//    std::vector<double> dUdYaux = P->PartitionAuxilaryData(dUidyi, comm);
//    std::vector<double> dUdZaux = P->PartitionAuxilaryData(dUidzi, comm);
//
//    Array<double>* dU2dXi2 = ComputedUdXi(P,dUdXaux,ien_copy->getNglob(),comm);
//    Array<double>* dU2dYi2 = ComputedUdXi(P,dUdYaux,ien_copy->getNglob(),comm);
//    Array<double>* dU2dZi2 = ComputedUdXi(P,dUdZaux,ien_copy->getNglob(),comm);
//
//    for(int i=0;i<Nel_part;i++)
//    {
//        dU2idx2i_new->setVal(i,0,dU2dXi2->getVal(i,0));
//        dU2idxyi_new->setVal(i,0,dU2dXi2->getVal(i,1));
//        dU2idxzi_new->setVal(i,0,dU2dXi2->getVal(i,2));
//
//        dU2idyxi_new->setVal(i,0,dU2dYi2->getVal(i,0));
//        dU2idy2i_new->setVal(i,0,dU2dYi2->getVal(i,1));
//        dU2idyzi_new->setVal(i,0,dU2dYi2->getVal(i,2));
//
//        dU2idzxi_new->setVal(i,0,dU2dZi2->getVal(i,0));
//        dU2idzyi_new->setVal(i,0,dU2dZi2->getVal(i,1));
//        dU2idz2i_new->setVal(i,0,dU2dZi2->getVal(i,2));
//
//        std::cout << "Hessian symmetric ??? " << i << std::endl;
//        std::cout << dU2idx2i_new->getVal(i,0) << " " << dU2idxyi_new->getVal(i,0) << " " << dU2idxzi_new->getVal(i,0) << std::endl;
//        std::cout << dU2idyxi_new->getVal(i,0) << " " << dU2idy2i_new->getVal(i,0) << " " << dU2idyzi_new->getVal(i,0) << std::endl;
//        std::cout << dU2idzxi_new->getVal(i,0) << " " << dU2idzyi_new->getVal(i,0) << " " << dU2idz2i_new->getVal(i,0) << std::endl;
//        std::cout << "======================" << std::endl;
//
//    }
//
//    std::vector<double> dU2dX2aux = P->PartitionAuxilaryData(dU2idx2i_new, comm);
//    std::vector<double> dU2dXYaux = P->PartitionAuxilaryData(dU2idxyi_new, comm);
//    std::vector<double> dU2dXZaux = P->PartitionAuxilaryData(dU2idxzi_new, comm);
//
//    std::vector<double> dU2dYXaux = P->PartitionAuxilaryData(dU2idyxi_new, comm);
//    std::vector<double> dU2dY2aux = P->PartitionAuxilaryData(dU2idy2i_new, comm);
//    std::vector<double> dU2dYZaux = P->PartitionAuxilaryData(dU2idyzi_new, comm);
//
//    std::vector<double> dU2dZXaux = P->PartitionAuxilaryData(dU2idzxi_new, comm);
//    std::vector<double> dU2dZYaux = P->PartitionAuxilaryData(dU2idzyi_new, comm);
//    std::vector<double> dU2dZ2aux = P->PartitionAuxilaryData(dU2idz2i_new, comm);
//
//    for(int i=0;i<nElemLoc;i++)
//    {
//        double uinew = Uaux[i];
//        double duidxnew = dUdXaux[i];
//        double duidynew = dUdYaux[i];
//        double duidznew = dUdZaux[i];
//
//        double du2idx2new = dU2dX2aux[i];
//        double du2idxynew = dU2dXYaux[i];
//        double du2idxznew = dU2dXZaux[i];
//
//        double du2idyxnew = dU2dYXaux[i];
//        double du2idy2new = dU2dY2aux[i];
//        double du2idyznew = dU2dYZaux[i];
//
//        double du2idzxnew = dU2dZXaux[i];
//        double du2idzynew = dU2dZYaux[i];
//        double du2idz2new = dU2dZ2aux[i];
//
//
//
////        if(world_rank == 0)
////        if(world_rank == 0)
////        {
////            std::cout << "du2idx2new  " << du2idx2new << std::endl;
////        }
//        for(int j=0;j<8;j++)
//        {
//            loc_v = loc_elem2verts_loc[i][j];
//            collect_Ui[loc_v].push_back(uinew);
//
//            collect_dUidX[loc_v].push_back(duidxnew);
//            collect_dUidY[loc_v].push_back(duidynew);
//            collect_dUidZ[loc_v].push_back(duidznew);
//
//            collect_dU2idX2[loc_v].push_back(du2idx2new);
//            collect_dU2idXY[loc_v].push_back(du2idxynew);
//            collect_dU2idXZ[loc_v].push_back(du2idxznew);
//
//            collect_dU2idYX[loc_v].push_back(du2idyxnew);
//            collect_dU2idY2[loc_v].push_back(du2idy2new);
//            collect_dU2idYZ[loc_v].push_back(du2idyznew);
//
//            collect_dU2idZX[loc_v].push_back(du2idzxnew);
//            collect_dU2idZY[loc_v].push_back(du2idzynew);
//            collect_dU2idZ2[loc_v].push_back(du2idz2new);
//
//        }
//    }
//    int c= 0;
//    double sum_u;
//
//    double sum_dudx;
//    double sum_dudy;
//    double sum_dudz;
//
//    double sum_du2dx2;
//    double sum_du2dxy;
//    double sum_du2dxz;
//
//    double sum_du2dyx;
//    double sum_du2dy2;
//    double sum_du2dyz;
//
//    double sum_du2dzx;
//    double sum_du2dzy;
//    double sum_du2dz2;
//
//    std::vector<double> uivert;
//    std::vector<double> dudxivert;
//    std::vector<double> dudyivert;
//    std::vector<double> dudzivert;
//
//    std::vector<double> du2dx2vert;
//    std::vector<double> du2dxyvert;
//    std::vector<double> du2dxzvert;
//
//    std::vector<double> du2dyxvert;
//    std::vector<double> du2dy2vert;
//    std::vector<double> du2dyzvert;
//
//    std::vector<double> du2dzxvert;
//    std::vector<double> du2dzyvert;
//    std::vector<double> du2dz2vert;
//
//    std::vector<double> WRarr_00;
//    std::vector<double> WRarr_01;
//    std::vector<double> WRarr_02;
//
//    std::vector<double> WRarr_10;
//    std::vector<double> WRarr_11;
//    std::vector<double> WRarr_12;
//
//    std::vector<double> WRarr_20;
//    std::vector<double> WRarr_21;
//    std::vector<double> WRarr_22;
//
//    double*Hmet = new double[9];
//
//    for(it_rhos=collect_Ui.begin();it_rhos!=collect_Ui.end();it_rhos++)
//    {
//        sum_u = 0;
//        sum_dudx = 0;
//        sum_dudy = 0;
//        sum_dudz = 0;
//        sum_du2dx2 = 0;
//        sum_du2dxy = 0;
//        sum_du2dxz = 0;
//        sum_du2dyx = 0;
//        sum_du2dy2 = 0;
//        sum_du2dyz = 0;
//        sum_du2dzx = 0;
//        sum_du2dzy = 0;
//        sum_du2dz2 = 0;
//        for(int q = 0;q<it_rhos->second.size();q++)
//        {
//            sum_u    = sum_u + it_rhos->second[q];
//            sum_dudx = sum_dudx + collect_dUidX[it_rhos->first][q];
//            sum_dudy = sum_dudy + collect_dUidY[it_rhos->first][q];
//            sum_dudz = sum_dudz + collect_dUidZ[it_rhos->first][q];
//
//            sum_du2dx2 = sum_du2dx2 + collect_dU2idX2[it_rhos->first][q];
//            sum_du2dxy = sum_du2dxy + collect_dU2idXY[it_rhos->first][q];
//            sum_du2dxz = sum_du2dxz + collect_dU2idXZ[it_rhos->first][q];
//
//            sum_du2dyx = sum_du2dyx + collect_dU2idYX[it_rhos->first][q];
//            sum_du2dy2 = sum_du2dy2 + collect_dU2idY2[it_rhos->first][q];
//            sum_du2dyz = sum_du2dyz + collect_dU2idYZ[it_rhos->first][q];
//
//            sum_du2dzx = sum_du2dzx + collect_dU2idZX[it_rhos->first][q];
//            sum_du2dzy = sum_du2dzy + collect_dU2idZY[it_rhos->first][q];
//            sum_du2dz2 = sum_du2dz2 + collect_dU2idZ2[it_rhos->first][q];
//
//        }
////        Hmet[0] = 1.0+sum_du2dx2/it_rhos->second.size();
////        Hmet[1] = sum_du2dxy/it_rhos->second.size();
////        Hmet[2] = sum_du2dxz/it_rhos->second.size();
////
////        Hmet[3] = sum_du2dyx/it_rhos->second.size();
////        Hmet[4] = 1.0+sum_du2dy2/it_rhos->second.size();
////        Hmet[5] = sum_du2dyz/it_rhos->second.size();
////
////        Hmet[6] = sum_du2dzx/it_rhos->second.size();
////        Hmet[7] = sum_du2dzy/it_rhos->second.size();
////        Hmet[8] = 1.0+sum_du2dz2/it_rhos->second.size();
////
////        double * WR = new double[3];
////        double * WI = new double[3];
////        double * V = new double[3*3];
////        double * iV = new double[3*3];
////
////        EigenDecomp(3, Hmet,  WR,  WI, V, iV );
////
////        WRarr_0.push_back(fabs(WR[0]) );
////        WRarr_1.push_back(fabs(WR[1]) );
////        WRarr_2.push_back(fabs(WR[2]));
//
//        uivert.push_back(sum_u/it_rhos->second.size());
//
//        dudxivert.push_back(sum_dudx/it_rhos->second.size());
//        dudyivert.push_back(sum_dudy/it_rhos->second.size());
//        dudzivert.push_back(sum_dudz/it_rhos->second.size());
//
//        du2dx2vert.push_back(sum_du2dx2/it_rhos->second.size());
//        du2dxyvert.push_back(sum_du2dxy/it_rhos->second.size());
//        du2dxzvert.push_back(sum_du2dxz/it_rhos->second.size());
//
//        du2dyxvert.push_back(sum_du2dyx/it_rhos->second.size());
//        du2dy2vert.push_back(sum_du2dy2/it_rhos->second.size());
//        du2dyzvert.push_back(sum_du2dyz/it_rhos->second.size());
//
//        du2dzxvert.push_back(sum_du2dzx/it_rhos->second.size());
//        du2dzyvert.push_back(sum_du2dzy/it_rhos->second.size());
//        du2dz2vert.push_back(sum_du2dz2/it_rhos->second.size());
//
//        c++;
//    }
//    double du2dx2vert_max = *std::max_element(du2dx2vert.begin(), du2dx2vert.end());
//    double du2dxyvert_max = *std::max_element(du2dxyvert.begin(), du2dxyvert.end());
//    double du2dxzvert_max = *std::max_element(du2dxzvert.begin(), du2dxzvert.end());
//    double du2dyxvert_max = *std::max_element(du2dyxvert.begin(), du2dyxvert.end());
//    double du2dy2vert_max = *std::max_element(du2dy2vert.begin(), du2dy2vert.end());
//    double du2dyzvert_max = *std::max_element(du2dyzvert.begin(), du2dyzvert.end());
//    double du2dzxvert_max = *std::max_element(du2dzxvert.begin(), du2dzxvert.end());
//    double du2dzyvert_max = *std::max_element(du2dzyvert.begin(), du2dzyvert.end());
//    double du2dz2vert_max = *std::max_element(du2dz2vert.begin(), du2dz2vert.end());
//
//    double R = 0.056;
//    double hmin = 0.000001;
//    double hmax = 0.01;
//
//
//
//    for(int i=0;i<LVerts.size();i++)
//    {
//        double r = sqrt((LVerts[i].x-0.0)*(LVerts[i].x-0.0)
//                       +(LVerts[i].y-0.0)*(LVerts[i].y-0.0));
//
//        //std::cout << "r = " << r << " " << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << std::endl;
//        if (r < R)
//        {
//            dudxivert[i] = 0.0;
//            dudyivert[i] = 0.0;
//            dudzivert[i] = 0.0;
//
//            du2dx2vert[i] = 0.0;
//            du2dxyvert[i] = 0.0;
//            du2dxzvert[i] = 0.0;
//
//            du2dyxvert[i] = 0.0;
//            du2dy2vert[i] = 0.0;
//            du2dyzvert[i] = 0.0;
//
//            du2dzxvert[i] = 0.0;
//            du2dzyvert[i] = 0.0;
//            du2dz2vert[i] = 0.0;
//        }
//
//        if(du2dx2vert[i]<0.1*du2dx2vert_max)
//        {
//            dudxivert[i] = 0.0;
//            dudyivert[i] = 0.0;
//            dudzivert[i] = 0.0;
//
//            du2dx2vert[i] = 0.0;
//            du2dxyvert[i] = 0.0;
//            du2dxzvert[i] = 0.0;
//
//            du2dyxvert[i] = 0.0;
//            du2dy2vert[i] = 0.0;
//            du2dyzvert[i] = 0.0;
//
//            du2dzxvert[i] = 0.0;
//            du2dzyvert[i] = 0.0;
//            du2dz2vert[i] = 0.0;
//        }
//
//        Hmet[0] = du2dx2vert[i];
//        Hmet[1] = du2dxyvert[i];
//        Hmet[2] = du2dxzvert[i];
//
//        Hmet[3] = du2dyxvert[i];
//        Hmet[4] = du2dy2vert[i];
//        Hmet[5] = du2dyzvert[i];
//
//        Hmet[6] = du2dzxvert[i];
//        Hmet[7] = du2dzyvert[i];
//        Hmet[8] = du2dz2vert[i];
//
//        double * WR = new double[3];
//        double * WI = new double[3];
//        double * V = new double[3*3];
//        double * iV = new double[3*3];
//        double* WRn = new double[3];
//        Array<double>* DR  = new Array<double>(3,3);
//        Array<double>* VR  = new Array<double>(3,3);
//        Array<double>* iVR = new Array<double>(3,3);
//        for(int i=0;i<3;i++)
//        {
//            WR[i]  = 0.0;
//            WRn[i] = 0.0;
//            WI[i]  = 0.0;
//            for(int j=0;j<3;j++)
//            {
//                DR->setVal(i,j,0.0);
//                VR->setVal(i,j,0.0);
//                iVR->setVal(i,j,0.0);
//                V[i*3+j] = 0.0;
//                iV[i*3+j] = 0.0;
//            }
//
//        }
//
//        EigenDecomp(3, Hmet,  WR,  WI, V, iV );
//
//        int f = 12;
//        WRn[0] = std::min(std::max(f*fabs(WR[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
//        WRn[1] = std::min(std::max(f*fabs(WR[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
//        WRn[2] = 1.0/(0.005*0.005);
//
////
//        DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
//        DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
//        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);
//
//        VR->setVal(0,0,V[0]);VR->setVal(0,1,V[1]);VR->setVal(0,2,V[2]);
//        VR->setVal(1,0,V[3]);VR->setVal(1,1,V[4]);VR->setVal(1,2,V[5]);
//        VR->setVal(2,0,V[6]);VR->setVal(2,1,V[7]);VR->setVal(2,2,V[8]);
//
//        iVR->setVal(0,0,iV[0]);iVR->setVal(0,1,iV[1]);iVR->setVal(0,2,iV[2]);
//        iVR->setVal(1,0,iV[3]);iVR->setVal(1,1,iV[4]);iVR->setVal(1,2,iV[5]);
//        iVR->setVal(2,0,iV[6]);iVR->setVal(2,1,iV[7]);iVR->setVal(2,2,iV[8]);
//
////        VRi->setVal(2,0,iV[6]);VRi->setVal(2,1,iV[7]);VRi->setVal(2,2,iV[8]);
//
//        Array<double>* Rs = MatMul(VR,DR);
//        Array<double>* Rf = MatMul(Rs,iVR);
//        if(i==2)
//        {
//        std::cout << "========================================================" << std::endl;
//        std::cout << "H " << std::endl;
//        std::cout << Hmet[0] << " " <<  Hmet[1] << " " <<  Hmet[2] << std::endl;
//        std::cout << Hmet[3] << " " <<  Hmet[4] << " " <<  Hmet[5] << std::endl;
//        std::cout << Hmet[6] << " " <<  Hmet[7] << " " <<  Hmet[8] << std::endl;
//        std::cout << "DR " << std::endl;
//        std::cout << DR->getVal(0,0) << " " <<  DR->getVal(0,1) << " " <<  DR->getVal(0,2) << std::endl;
//        std::cout << DR->getVal(1,0) << " " <<  DR->getVal(1,1) << " " <<  DR->getVal(1,2) << std::endl;
//        std::cout << DR->getVal(2,0) << " " <<  DR->getVal(2,1) << " " <<  DR->getVal(2,2) << std::endl;
//
//        std::cout << "VR " << std::endl;
//        std::cout << VR->getVal(0,0) << " " <<  VR->getVal(0,1) << " " <<  VR->getVal(0,2) << std::endl;
//        std::cout << VR->getVal(1,0) << " " <<  VR->getVal(1,1) << " " <<  VR->getVal(1,2) << std::endl;
//        std::cout << VR->getVal(2,0) << " " <<  VR->getVal(2,1) << " " <<  VR->getVal(2,2) << std::endl;
//
//        std::cout << "iVR " << std::endl;
//        std::cout << iVR->getVal(0,0) << " " <<  iVR->getVal(0,1) << " " <<  iVR->getVal(0,2) << std::endl;
//        std::cout << iVR->getVal(1,0) << " " <<  iVR->getVal(1,1) << " " <<  iVR->getVal(1,2) << std::endl;
//        std::cout << iVR->getVal(2,0) << " " <<  iVR->getVal(2,1) << " " <<  iVR->getVal(2,2) << std::endl;
//
//        std::cout << "Rs " << std::endl;
//        std::cout << Rs->getVal(0,0) << " " <<  Rs->getVal(0,1) << " " <<  Rs->getVal(0,2) << std::endl;
//        std::cout << Rs->getVal(1,0) << " " <<  Rs->getVal(1,1) << " " <<  Rs->getVal(1,2) << std::endl;
//        std::cout << Rs->getVal(2,0) << " " <<  Rs->getVal(2,1) << " " <<  Rs->getVal(2,2) << std::endl;
//        std::cout << "Rf " << std::endl;
//        std::cout << Rf->getVal(0,0) << " " <<  Rf->getVal(0,1) << " " <<  Rf->getVal(0,2) << std::endl;
//        std::cout << Rf->getVal(1,0) << " " <<  Rf->getVal(1,1) << " " <<  Rf->getVal(1,2) << std::endl;
//        std::cout << Rf->getVal(2,0) << " " <<  Rf->getVal(2,1) << " " <<  Rf->getVal(2,2) << std::endl;
//        std::cout << "========================================================" << std::endl;
//        }
//
////        WRarr_00.push_back( Rf->getVal(0,0) );
////        WRarr_01.push_back( Rf->getVal(0,1) );
////        WRarr_02.push_back( Rf->getVal(0,2) );
////
////        WRarr_10.push_back( Rf->getVal(1,0) );
////        WRarr_11.push_back( Rf->getVal(1,1) );
////        WRarr_12.push_back( Rf->getVal(1,2) );
////
////        WRarr_20.push_back( Rf->getVal(2,0) );
////        WRarr_21.push_back( Rf->getVal(2,1) );
////        WRarr_22.push_back( Rf->getVal(2,2) );
//
//        WRarr_00.push_back( Hmet[0] );
//        WRarr_01.push_back( Hmet[1] );
//        WRarr_02.push_back( Hmet[2] );
//
//        WRarr_10.push_back( Hmet[1] );
//        WRarr_11.push_back( Hmet[4] );
//        WRarr_12.push_back( Hmet[5] );
//
//        WRarr_20.push_back( Hmet[2] );
//        WRarr_21.push_back( Hmet[5] );
//        WRarr_22.push_back( Hmet[8] );
//
//
//        delete DR;
//        delete VR;
//        delete iVR;
//        delete Rs;
//        delete Rf;
//
//        delete[] iV;
//        delete[] V;
//        delete[] WR;
//        delete[] WI;
//        delete[] WRn;
//
//    }
//
//
//    int nvert = uivert.size();
//    UnitTestEigenDecomp();
//
//    //==============================================================================================
//    //================================Output the data in Tecplot format=============================
//    //==============================================================================================
//    string filename10 = "pre_elem_rank_" + std::to_string(world_rank) + ".dat";
//    ofstream myfile10;
//    myfile10.open(filename10);
//    //string filename2 = "eig_slice_rank_" + std::to_string(world_rank) + ".dat";
//    //ofstream myfile2;
//    //myfile2.open(filename2);
//    string filename3 = "pre_metric_rank_" + std::to_string(world_rank) + "_v1.dat";
//    ofstream myfile3;
//    myfile3.open(filename3);
//    string filename = "eig_rank_" + std::to_string(world_rank) + ".dat";
//    ofstream myfile;
//    myfile.open(filename);
//    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"hxx\", \"hxy\", \"hxz\", \"hyx\", \"hyy\", \"hyz\", \"hzx\", \"hzy\", \"hzz\"" << std::endl;
//    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"hxx\", \"hxy\", \"hxz\", \"hyx\", \"hyy\", \"hyz\", \"hzx\", \"hzy\", \"hzz\"" << std::endl;
//    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"var\", \"gradx_var\", \"grady_var\", \"gradz_var\", \"v00\", \"v01\", \"v02\", \"v10\", \"v11\", \"v12\", \"v20\", \"v21\", \"v22\"" << std::endl;
//    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"var\", \"gradx_var\", \"grady_var\", \"gradz_var\",\"gxx\",\"gxy\",\"gxz\",\"gyx\",\"gyy\",\"gyz\",\"gzx\",\"gzy\",\"gzz\",\"H00\",\"H01\",\"H02\",\"H10\",\"H11\",\"H12\",\"H20\",\"H21\",\"H22\"" << std::endl;
//    //myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"var\", \"gradx_var\", \"grady_var\", \"gradz_var\"" << std::endl;
//    std::cout << " verts check " << LVerts.size() << " " << WRarr_00.size() << std::endl;
//    int nvert2 = LVerts.size();
//    myfile <<"ZONE N = " << nvert2 << ", E = " << ien_copy->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
//    for(int i=0;i<nvert;i++)
//    {
//       //myfile << LVerts[i].x << "   " << LVerts[i].y << "   " << LVerts[i].z << "   " << hxx[i] << " " << hxy[i] << " " << hxz[i] << " " << hyx[i] << " " << hyy[i] << " " << hyz[i] << " " << hzx[i] << " " << hzy[i] << " " << hzz[i] << std::endl;
//        myfile << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << uivert[i] << " " << dudxivert[i] << " " << dudyivert[i]<< " " << dudzivert[i] << " " << du2dx2vert[i] << " " << du2dxyvert[i]<< " " << du2dxzvert[i]<< " " << du2dyxvert[i]<< " " << du2dy2vert[i]<< " " << du2dyzvert[i]<< " " << du2dzxvert[i]<< " " << du2dzyvert[i]<< " " << du2dz2vert[i] << " " << WRarr_00[i]<< " " << WRarr_01[i]<< " " << WRarr_02[i] << " " << WRarr_10[i]<< " " << WRarr_11[i]<< " " << WRarr_12[i] << " " << WRarr_20[i]<< " " << WRarr_21[i]<< " " << WRarr_22[i] <<std::endl;
//        //myfile << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << uivert[i] << " " << dudxivert[i] << " " << dudyivert[i]<< " " << dudzivert[i] << std::endl;
//
//        myfile3 << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << WRarr_00[i]<< " " << WRarr_01[i]<< " " << WRarr_02[i] << " " << WRarr_10[i]<< " " << WRarr_11[i]<< " " << WRarr_12[i] << " " << WRarr_20[i]<< " " << WRarr_21[i]<< " " << WRarr_22[i]<< std::endl;
//
//
//
//        //myfile << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << uivert[i] << " " << dudxivert[i] << " " << dudyivert[i]<< " " << dudzivert[i]<< " " << Hess[0][i] << " " << Hess[1][i] << " " << Hess[2][i] << " " << Hess[3][i] << " " << Hess[4][i] << " " << Hess[5][i] << " " << Hess[6][i] << " " << Hess[7][i] << " " << Hess[8][i] << std::endl;
//        //myfile << LVerts[i].x << "   " << LVerts[i].y << "   " << LVerts[i].z << "   " << DiPlot[0][i] << " " << DiPlot[1][i] << " " << DiPlot[2][i] << " " << ViPlot[0][i] << " " << ViPlot[1][i] << " " << ViPlot[2][i] << "   " << ViPlot[3][i] << " " << ViPlot[4][i] << " " << ViPlot[5][i] << "   " << ViPlot[6][i] << " " << ViPlot[7][i] << " " << ViPlot[8][i] << std::endl;
//        //myfile << LVerts[i].x << "   " << LVerts[i].y << "   " << LVerts[i].z << "   " << ViPlot[0][i] << " " << ViPlot[1][i] << " " << ViPlot[2][i] << std::endl;
//       //myfile3 << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << Hess[0][i] << " " << Hess[1][i] << " " << Hess[2][i] << " " << Hess[3][i] << " " << Hess[4][i] << " " << Hess[5][i] << " " << Hess[6][i] << " " << Hess[7][i] << " " << Hess[8][i] << std::endl;
////       if(fabs(LVerts[i].z)<1.0e-02)
////       {
////           myfile2 << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << DiPlot->getVal(0,i) << " " << DiPlot->getVal(1,i) << " " << DiPlot->getVal(2,i) << " " << ViPlot->getVal(0,i) << " " << ViPlot->getVal(1,i) << " " << ViPlot->getVal(2,i) << " " << ViPlot->getVal(3,i) << " " << ViPlot->getVal(4,i) << " " << ViPlot->getVal(5,i) << " " << ViPlot->getVal(6,i) << " " << ViPlot->getVal(7,i) << " " << ViPlot->getVal(8,i) << std::endl;
////     //  myfile3 << LVerts[i].x << " " << LVerts[i].y << " " << LVerts[i].z << " " << Hess[0][i] << " " << Hess[1][i] << " " << Hess[2][i] << " " << Hess[3][i] << " " << Hess[4][i] << " " << Hess[5][i] << " " << Hess[6][i] << " " << Hess[7][i] << " " << Hess[8][i] << " " << DiPlot->getVal(i,0) << " " << DiPlot->getVal(i,1) << " " << DiPlot->getVal(i,2) << std::endl;
////       }
//    }
//
//    for(int i=0;i<ien_copy->getNrow();i++)
//    {
//       myfile10 << loc_elem2verts_loc[i][0]+1 << " " <<
//                 loc_elem2verts_loc[i][1]+1 << " " <<
//                 loc_elem2verts_loc[i][2]+1 << " " <<
//                 loc_elem2verts_loc[i][3]+1 << " " <<
//                 loc_elem2verts_loc[i][4]+1 << " " <<
//                 loc_elem2verts_loc[i][5]+1 << " " <<
//                 loc_elem2verts_loc[i][6]+1 << " " <<
//                 loc_elem2verts_loc[i][7]+1 << std::endl;
//
//
//       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
//                 loc_elem2verts_loc[i][1]+1 << "  " <<
//                 loc_elem2verts_loc[i][2]+1 << "  " <<
//                 loc_elem2verts_loc[i][3]+1 << "  " <<
//                 loc_elem2verts_loc[i][4]+1 << "  " <<
//                 loc_elem2verts_loc[i][5]+1 << "  " <<
//                 loc_elem2verts_loc[i][6]+1 << "  " <<
//                 loc_elem2verts_loc[i][7]+1 << std::endl;
//    }
//    myfile.close();
//    //myfile2.close();
//    myfile3.close();
//    myfile10.close();
    
    
    MPI_Finalize();
    
    return 0;
     
}
