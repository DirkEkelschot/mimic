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



Eig* ComputeEigenDecomp(int n, double * A)
{
  Eig* eig = new Eig;
  char JOBVL = 'V';
  char JOBVR = 'N';
  int size = 10*n;
  double WORK [size];
  int info;
  int i,j;
  int Pivot[n];
  double * WR = new double[n];
  double * WI = new double[n];
  double * V = new double[n*n];
  double * iV = new double[n*n];
    
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
    
    eig->Dre = WR;
    eig->Dim = WI;
    
    eig->V   = V;
    eig->iV  = iV;
    return eig;
}

SVD* ComputeSVD(int M, int N, double * A)
{
    SVD* svd = new SVD;
    double wkopt;
    int size = 100*N;
    double WORK [size];
    int info, lwork;
    int i,j;

    int ldu = M;
    int ldvt = N;
    int lda = M;
    char all = 'A';
    double* s = new double[N];
    double* u = new double[ldu*M];
    double* vt = new double[ldvt*N];
    //double u[ldu*M], vt[ldvt*N];
    dgesvd_( &all, &all, &M, &N, A, &lda, s, u, &ldu, vt, &ldvt, WORK, &size, &info );
    svd->s  = s;
    svd->u  = u;
    svd->vt = vt;
    return svd;
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
 




 
 void UnitTestSVD()
 {
     double* A = new double[9];
     A[0] = 2758;A[1]=-2477.29;A[2]=1.00706e-12;
     A[3] = -2477.29;A[4]=324.239;A[5]=2.29616e-11;
     A[6] = 1.00706e-12;A[7]=2.29616e-11;A[8]=1.55016e-11;
     SVD* svd = ComputeSVD(3, 3, A);
     std::cout << " s " << std::endl;
     std::cout << svd->s[0] << " " << svd->s[1] << " " << svd->s[2] << std::endl;
     std::cout << " u " << std::endl;
     std::cout << svd->u[0] << " " << svd->u[1] << " " << svd->u[2] << std::endl;
     std::cout << svd->u[3] << " " << svd->u[4] << " " << svd->u[5] << std::endl;
     std::cout << svd->u[6] << " " << svd->u[7] << " " << svd->u[8] << std::endl;
     std::cout << " vt " << std::endl;
     std::cout << svd->vt[0] << " " << svd->vt[1] << " " << svd->vt[2] << std::endl;
     std::cout << svd->vt[3] << " " << svd->vt[4] << " " << svd->vt[5] << std::endl;
     std::cout << svd->vt[6] << " " << svd->vt[7] << " " << svd->vt[8] << std::endl;
     
     
     delete svd;
 //    U= [[-8.48791491e-01 -5.28727723e-01  1.09121085e-14]
 //    [ 5.28727723e-01 -8.48791491e-01  1.25551127e-14]
 //        [ 2.62386865e-15  1.64262071e-14  1.00000000e+00]];
 //
 //    Vh= [[-8.48791491e-01  5.28727723e-01  2.62386865e-15]
 //    [ 5.28727723e-01  8.48791491e-01 -1.64262071e-14]
 //    [ 1.09121085e-14  1.25551127e-14  1.00000000e+00]];

 }
 
 
void UnitTestEigenDecomp()
{
    double *M = new double[3*3];
    M[0] = 0.25;M[1]=-0.3;M[2]=0.4;
    M[3] = -0.3;M[4]=1.25;M[5]=0.1;
    M[6] = 0.4;M[7]=0.1;M[8]=0.25;
    
//    double * WR = new double[3];
//    double * WI = new double[3];
//    double * V = new double[3*3];
//    double * iV = new double[3*3];
    //EigenDecomp(3, M,  WR,  WI, V, iV );
    Eig* eig = ComputeEigenDecomp(3,M);
    for(int i=0;i<3;i++)
    {
        std::cout << "eigenvalues = " << eig->Dre[i] << " + " << eig->Dim[i] << "i" << std::endl;
    }
    
    double* Vref = new double[3*3];
    double* iVref = new double[3*3];
    
    Vref[0] = 0.71598;Vref[1] = 0.643504;Vref[2] = -0.270694;
    Vref[3] = 0.193611;Vref[4] = 0.189507;Vref[5] = 0.962602;
    Vref[6] = -0.670736;Vref[7] = 0.741613;Vref[8] = -0.0110941;
    
    iVref[0] = 0.71598;iVref[1] = 0.193611;iVref[2] = -0.670736;
    iVref[3] = 0.643504;iVref[4] = 0.189507;iVref[5] = 0.741613;
    iVref[6] = -0.270694;iVref[7] = 0.962602;iVref[8] = -0.0110941;
    
    std::cout << "V ==> "  << std::endl;;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << eig->V[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
    std::cout << "Vinv ==> " << std::endl;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << eig->iV[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//#endif
