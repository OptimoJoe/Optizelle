#include "linalg.h"
#include "FortranCInterface.h"
extern "C" {
// Double BLAS and LAPACK routines
#define dcopy_fortran FortranCInterface_GLOBAL (dcopy,DCOPY)
void dcopy_fortran(int* m,double* x,int* incx,double* y,int* incy);
#define daxpy_fortran FortranCInterface_GLOBAL (daxpy,DAXPY)
void daxpy_fortran(int* m,double* alpha,double* x,int* incx,double* y,
    int* incy);
#define dscal_fortran FortranCInterface_GLOBAL (dscal,DSCAL)
void dscal_fortran(int* m,double* alpha,double* x,int* incx);
#define ddot_fortran FortranCInterface_GLOBAL (ddot,DDOT)
double ddot_fortran(int* m,double* x,int* incx,double* y,int* incy);
#define dsyr2k_fortran FortranCInterface_GLOBAL (dsyr2k,DSYR2K)
void dsyr2k_fortran(char* uplo,char* trans,int* n,int* k,double* alpha,
    double* A,int* lda,double* B,int* ldb,double* beta,
    double* C,int* ldc);
#define dsyevr_fortran FortranCInterface_GLOBAL (dsyevr,DSYEVR)
void dsyevr_fortran(char* jobz,char* range,char* uplo,int* n,double *A,
    int* lda,double* vl,double* vu,int* il,int* iu,double* abstol,
    int* m,double* w,double* z,int* ldz,int* isuppz,double* work,
    int* lwork,int* iwork,int* liwork,int* info);
#define dstemr_fortran FortranCInterface_GLOBAL (dstemr,DSTEMR)
void dstemr_fortran(char* jobz,char* range,int* n,double *D,double *E,
    double* vl,double* vu,int* il,int* iu,int* m,double* w,double* z,
    int* ldz,int* nzc,int* isuppz,int* tryrac,double* work,int* lwork,
    int* iwork,int* liwork,int* info);
#define dstevr_fortran FortranCInterface_GLOBAL (dstevr,DSTEVR)
void dstevr_fortran(char* jobz,char* range,int* n,double *D,double *E,
    double* vl,double* vu,int* il,int* iu,double* abstol,int* m,
    double* w,double* z,int* ldz,int* isuppz,double* work,int* lwork,
    int* iwork,int* liwork,int* info);
#define dlamch_fortran FortranCInterface_GLOBAL (dlamch,DLAMCH)
double dlamch_fortran(char* cmach);
#define dgemm_fortran FortranCInterface_GLOBAL (dgemm,DGEMM)
void dgemm_fortran(char* transa,char* transb,int* m,int* n,int* k,
    double* alpha, double* A, int* lda,double* B,int* ldb,
    double* beta,double* C,int *ldc);
#define dsymm_fortran FortranCInterface_GLOBAL (dsymm,DSYMM)
void dsymm_fortran(char* side,char* uplo,int* m,int* n,double* alpha,
    double* A, int* lda,double* B,int* ldb,double* beta,double* C,
    int *ldc);
#define dsymv_fortran FortranCInterface_GLOBAL (dsymv,DSYMV)
void dsymv_fortran(char* uplo,int* n,double* alpha,double* A, int* lda,
    double* x,int* incx,double* beta,double* y,int *incy);
#define dpotrf_fortran FortranCInterface_GLOBAL (dpotrf,DPOTRF)
void dpotrf_fortran(char* uplo,int *n,double* A,int* lda,int* info);
#define dtrtri_fortran FortranCInterface_GLOBAL (dtrtri,DTRTRI)
void dtrtri_fortran(char* uplo,char* diag,int* n,double* A,int* lda,int* info);

// Float BLAS and LAPACK routines
#define scopy_fortran FortranCInterface_GLOBAL (scopy,SCOPY)
void scopy_fortran(int* m,float* x,int* incx,float* y,int* incy);
#define saxpy_fortran FortranCInterface_GLOBAL (saxpy,SAXPY)
void saxpy_fortran(int* m,float* alpha,float* x,int* incx,float* y,
    int* incy);
#define sscal_fortran FortranCInterface_GLOBAL (sscal,SSCAL)
void sscal_fortran(int* m,float* alpha,float* x,int* incx);
#define sdot_fortran FortranCInterface_GLOBAL (sdot,SDOT)
float sdot_fortran(int* m,float* x,int* incx,float* y,int* incy);
#define ssyr2k_fortran FortranCInterface_GLOBAL (ssyr2k,SSYR2K)
void ssyr2k_fortran(char* uplo,char* trans,int* n,int* k,float* alpha,
    float* A,int* lda,float* B,int* ldb,float* beta,
    float* C,int* ldc);
#define ssyevr_fortran FortranCInterface_GLOBAL (ssyevr,SSYEVR)
void ssyevr_fortran(char* jobz,char* range,char* uplo,int* n,float *A,
    int* lda,float* vl,float* vu,int* il,int* iu,float* abstol,
    int* m,float* w,float* z,int* ldz,int* isuppz,float* work,
    int* lwork,int* iwork,int* liwork,int* info);
#define sstemr_fortran FortranCInterface_GLOBAL (sstemr,SSTEMR)
void sstemr_fortran(char* jobz,char* range,int* n,float *D,float *E,
    float* vl,float* vu,int* il,int* iu,int* m,float* w,float* z,
    int* ldz,int* nzc,int* isuppz,int* tryrac,float* work,int* lwork,
    int* iwork,int* liwork,int* info);
#define sstevr_fortran FortranCInterface_GLOBAL (sstevr,SSTEVR)
void sstevr_fortran(char* jobz,char* range,int* n,float *D,float *E,
    float* vl,float* vu,int* il,int* iu,float* abstol,int* m,
    float* w,float* z,int* ldz,int* isuppz,float* work,int* lwork,
    int* iwork,int* liwork,int* info);
#define slamch_fortran FortranCInterface_GLOBAL (slamch,SLAMCH)
float slamch_fortran(char* cmach);
#define sgemm_fortran FortranCInterface_GLOBAL (sgemm,SGEMM)
void sgemm_fortran(char* transa,char* transb,int* m,int* n,int* k,
    float* alpha, float* A, int* lda,float* B,int* ldb,
    float* beta,float* C,int *ldc);
#define ssymm_fortran FortranCInterface_GLOBAL (ssymm,SSYMM)
void ssymm_fortran(char* side,char* uplo,int* m,int* n,float* alpha,
    float* A, int* lda,float* B,int* ldb,float* beta,float* C,
    int *ldc);
#define ssymv_fortran FortranCInterface_GLOBAL (ssymv,SSYMV)
void ssymv_fortran(char* uplo,int* n,float* alpha,float* A, int* lda,
    float* x,int* incx,float* beta,float* y,int *incy);
#define spotrf_fortran FortranCInterface_GLOBAL (spotrf,SPOTRF)
void spotrf_fortran(char* uplo,int *n,float* A,int* lda,int* info);
#define strtri_fortran FortranCInterface_GLOBAL (strtri,STRTRI)
void strtri_fortran(char* uplo,char* diag,int* n,float* A,int* lda,int* info);
}

namespace peopt {

    template <>
    void copy <double> (int n,const double* x,int incx,double* y,int incy) {
        dcopy_fortran(&n,const_cast <double*> (x),&incx,y,&incy);
    }

    template <>
    void copy <float> (int n,const float* x,int incx,float* y,int incy) {
        scopy_fortran(&n,const_cast <float*> (x),&incx,y,&incy);
    }

    template <>
    void axpy <double> (
        int n,double alpha,const double* x,int incx,double* y,int incy
    ) {
        daxpy_fortran(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
    }
    
    template <>
    void axpy <float> (
        int n,float alpha,const float* x,int incx,float* y,int incy
    ) {
        saxpy_fortran(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
    }
    
    template <>
    void scal <double> (int n,double alpha,double* x,int incx) {
        dscal_fortran(&n,&alpha,x,&incx);
    }

    template <>
    void scal <float> (int n,float alpha,float* x,int incx) {
        sscal_fortran(&n,&alpha,x,&incx);
    }

    template <>
    double dot <double> (
        int n,const double* x,int incx,const double* y,int incy
    ) {
        return ddot_fortran(&n,const_cast <double*> (x),&incx,
            const_cast <double*> (y),&incy);
    }

    template <>
    float dot <float> (int n,const float* x,int incx,const float* y,int incy) {
        return sdot_fortran(&n,const_cast <float*> (x),&incx,
            const_cast <float*> (y),&incy);
    }

    template <>
    void syr2k(char uplo,char trans,int n,int k,double alpha,
        const double* A,int lda,const double* B,int ldb,double beta,
        double* C,int ldc
    ) {
        dsyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
            const_cast <double*> (B),&ldb,&beta,C,&ldc);
    }

    template <>
    void syr2k(char uplo,char trans,int n,int k,float alpha,
        const float* A,int lda,const float* B,int ldb,float beta,
        float* C,int ldc
    ) {
        ssyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
            const_cast <float*> (B),&ldb,&beta,C,&ldc);
    }

    template <>
    void syevr(char jobz,char range,char uplo,int n,double *A,int lda,
        double vl,double vu,int il,int iu,double abstol,int& m,
        double* w,double* z,int ldz,int* isuppz,double* work,int lwork,
        int* iwork,int liwork,int& info
    ) {
        dsyevr_fortran(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
            &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
    }

    template <>
    void syevr(char jobz,char range,char uplo,int n,float *A,int lda,
        float vl,float vu,int il,int iu,float abstol,int& m,
        float* w,float* z,int ldz,int* isuppz,float* work,int lwork,
        int* iwork,int liwork,int& info
    ) {
        ssyevr_fortran(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
            &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stemr(char jobz,char range,int n,double *D,double *E,double vl,
        double vu,int il,int iu,int& m,double* w,double* z,int ldz,int nzc,
        int* isuppz,int& tryrac,double* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        dstemr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
            isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stemr(char jobz,char range,int n,float *D,float *E,float vl,
        float vu,int il,int iu,int& m,float* w,float* z,int ldz,int nzc,
        int* isuppz,int& tryrac,float* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        sstemr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
            isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stevr(char jobz,char range,int n,double *D,double *E,double vl,
        double vu,int il,int iu,double abstol,int& m,double* w,double* z,
        int ldz,int* isuppz,double* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        dstevr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,z,&ldz,
            isuppz,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stevr(char jobz,char range,int n,float *D,float *E,float vl,
        float vu,int il,int iu,float abstol,int& m,float* w,float* z,
        int ldz,int* isuppz,float* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        sstevr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,z,&ldz,
            isuppz,work,&lwork,iwork,&liwork,&info);
    }

    template <>
    double lamch(char cmach) {
        return dlamch_fortran(&cmach);
    }

    template <>
    float lamch(char cmach) {
        return slamch_fortran(&cmach);
    }

    template <>
    void gemm(char transa,char transb,int m,int n,int k,double alpha,
        const double* A,int lda,const double* B,int ldb,double beta,
        double* C,int ldc
    ) {
        dgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }
    
    template <>
    void gemm(char transa,char transb,int m,int n,int k,float alpha,
        const float* A,int lda,const float* B,int ldb,float beta,
        float* C,int ldc
    ) {
        sgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void symm(char side,char uplo,int m,int n,double alpha,const double* A,
        int lda,const double* B,int ldb,double beta,double* C,int ldc
    ) {
        dsymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void symm(char side,char uplo,int m,int n,float alpha,const float* A,
        int lda,const float* B,int ldb,float beta,float* C,int ldc
    ) {
        ssymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }
    
    template <>
    void symv(char uplo,int n,double alpha,const double* A,int lda,
        const double* x,int incx,double beta,double* y,int incy
    ) {
        dsymv_fortran(&uplo,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
    }
    
    template <>
    void symv(char uplo,int n,float alpha,const float* A,int lda,
        const float* x,int incx,float beta,float* y,int incy
    ) {
        ssymv_fortran(&uplo,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
    }

    template <>
    void potrf(char uplo,int n,double* A,int lda,int& info) {
        dpotrf_fortran(&uplo,&n,A,&lda,&info);
    }
    
    template <>
    void potrf(char uplo,int n,float* A,int lda,int& info) {
        spotrf_fortran(&uplo,&n,A,&lda,&info);
    }

    template <>
    void trtri(char uplo,char diag,int n,double* A,int lda,int& info) {
        dtrtri_fortran(&uplo,&diag,&n,A,&lda,&info);
    }
    
    template <>
    void trtri(char uplo,char diag,int n,float* A,int lda,int& info) {
        strtri_fortran(&uplo,&diag,&n,A,&lda,&info);
    }
    
    // Indexing function for matrices.  This assumes that the upper right
    // corner has index (1,1).  Basically, we're using the indexing scheme
    // of math and Fortran and not of C.
    unsigned int ijtok(unsigned int i,unsigned int j,unsigned int m){
        return (i-1)+(j-1)*m;
    }


}
