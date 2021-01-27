#include "adapt_array.h"


#ifndef ADAPT_MATH_H
#define ADAPT_MATH_H

/*
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
*/

struct SVD{
    double* s;
    double* u;
    double* vt;
};

struct Eig{
    double* Dre;
    double* Dim;
    double* V;
    double* iV;
};

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

    void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a,
                int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                double* work, int* lwork, int* info );
}

int geqrf(int m, int n, double* A, int lda, double *tau);
int ormqr(char side, char trans, int m, int n, int k, double *A, int lda, double *tau, double* C, int ldc);
int trtrs(char uplo, char trans, char diag, int n, int nrhs, double* A, int lda, double* B, int ldb);


Array<double>* SolveQR(double* A, int m, int n, Array<double>* b);

void EigenDecomp(int n, double * A,  double * WR, double * WI, double * V, double * iV );

bool isDiagonalMatrix(Array<double>* Msq);

Array<double>* MatInv(Array<double>* A);

Eig* ComputeEigenDecomp(int n, double * A);

SVD* ComputeSVD(int M, int N, double * A);

void UnitTestSVD();

void UnitTestEigenDecomp();


#endif
