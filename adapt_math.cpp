#include "adapt_math.h"

int geqrf(int m, int n, double* A, int lda, double *tau)
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

int ormqr(char side, char trans, int m, int n, int k, double *A, int lda, double *tau, double* C, int ldc)
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

int trtrs(char uplo, char trans, char diag, int n, int nrhs, double* A, int lda, double* B, int ldb)
{
    int info = 0;
    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs,
            A, &lda, B, &ldb, &info);
    return info;
}

Array<double>* SolveQR(double* A, int m, int n, Array<double>* b)
{
    int nrow = m;
    int ncol = n;
    Array<double>* out;
    int info = 0;
    
    double tau[ncol];
    int lwork=1;
    double iwork;
    Array<double>* b_copy = new Array<double>(m,1);
    for(int i=0;i<m;i++)
    {
    b_copy->setVal(i,0,b->getVal(i,0));
    }
    // DGEQRF for Q*R=A, i.e., A and tau hold R and Householder reflectors

    geqrf(nrow, ncol, A, nrow, tau);
    
    ormqr('L', 'T', nrow, 1, ncol, A, nrow, tau, b_copy->data, nrow);
    
    trtrs('U', 'N', 'N', ncol, 1, A, nrow, b_copy->data, nrow);
    
    out = new Array<double>(ncol,1);
    for(int i = 0;i<ncol;i++)
    {
        out->setVal(i,0,b_copy->getVal(i,0));
    }
   
    delete b_copy;
    return out;
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
    
    dgetrf_(&n, &n, R->data, &n, Pivot, &info);
    dgetri_(&n, R->data, &n, Pivot, WORK, &size, &info);
    
    return R;
}
 



/*
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
*/
//#endif
