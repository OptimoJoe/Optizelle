#include "linalg.h"
#include "FortranCInterface.h"
using peopt::Integer;
extern "C" {
// Double BLAS and LAPACK routines
#define dcopy_fortran FortranCInterface_GLOBAL (dcopy,DCOPY)
void dcopy_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy);
#define daxpy_fortran FortranCInterface_GLOBAL (daxpy,DAXPY)
void daxpy_fortran(Integer* m,double* alpha,double* x,Integer* incx,double* y,
    Integer* incy);
#define dscal_fortran FortranCInterface_GLOBAL (dscal,DSCAL)
void dscal_fortran(Integer* m,double* alpha,double* x,Integer* incx);
#define ddot_fortran FortranCInterface_GLOBAL (ddot,DDOT)
double ddot_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy);
#define dsyr2k_fortran FortranCInterface_GLOBAL (dsyr2k,DSYR2K)
void dsyr2k_fortran(char* uplo,char* trans,Integer* n,Integer* k,double* alpha,
    double* A,Integer* lda,double* B,Integer* ldb,double* beta,
    double* C,Integer* ldc);
#define dsyevr_fortran FortranCInterface_GLOBAL (dsyevr,DSYEVR)
void dsyevr_fortran(char* jobz,char* range,char* uplo,Integer* n,double *A,
    Integer* lda,double* vl,double* vu,Integer* il,Integer* iu,double* abstol,
    Integer* m,double* w,double* z,Integer* ldz,Integer* isuppz,double* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define dstemr_fortran FortranCInterface_GLOBAL (dstemr,DSTEMR)
void dstemr_fortran(char* jobz,char* range,Integer* n,double *D,double *E,
    double* vl,double* vu,Integer* il,Integer* iu,Integer* m,double* w,
    double* z,Integer* ldz,Integer* nzc,Integer* isuppz,Integer* tryrac,
    double* work,Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define dstevr_fortran FortranCInterface_GLOBAL (dstevr,DSTEVR)
void dstevr_fortran(char* jobz,char* range,Integer* n,double *D,double *E,
    double* vl,double* vu,Integer* il,Integer* iu,double* abstol,Integer* m,
    double* w,double* z,Integer* ldz,Integer* isuppz,double* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define dlamch_fortran FortranCInterface_GLOBAL (dlamch,DLAMCH)
double dlamch_fortran(char* cmach);
#define dgemm_fortran FortranCInterface_GLOBAL (dgemm,DGEMM)
void dgemm_fortran(char* transa,char* transb,Integer* m,Integer* n,Integer* k,
    double* alpha, double* A, Integer* lda,double* B,Integer* ldb,
    double* beta,double* C,Integer *ldc);
#define dsymm_fortran FortranCInterface_GLOBAL (dsymm,DSYMM)
void dsymm_fortran(char* side,char* uplo,Integer* m,Integer* n,double* alpha,
    double* A, Integer* lda,double* B,Integer* ldb,double* beta,double* C,
    Integer *ldc);
#define dsymv_fortran FortranCInterface_GLOBAL (dsymv,DSYMV)
void dsymv_fortran(char* uplo,Integer* n,double* alpha,double* A, Integer* lda,
    double* x,Integer* incx,double* beta,double* y,Integer *incy);
#define dpotrf_fortran FortranCInterface_GLOBAL (dpotrf,DPOTRF)
void dpotrf_fortran(char* uplo,Integer *n,double* A,Integer* lda,Integer* info);
#define dtrtri_fortran FortranCInterface_GLOBAL (dtrtri,DTRTRI)
void dtrtri_fortran(char* uplo,char* diag,Integer* n,double* A,Integer* lda,
    Integer* info);
#define drotg_fortran FortranCInterface_GLOBAL (drotg,DROTG)
void drotg_fortran(double* a,double* b,double* c,double *s);
#define drot_fortran FortranCInterface_GLOBAL (drot,DROT)
void drot_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy,
    double* c,double *s);
#define dtpsv_fortran FortranCInterface_GLOBAL (dtpsv,DTPSV)
void dtpsv_fortran(char* uplo,char* trans,char* diag,Integer* n,double* Ap,
    double* x,Integer* incx);

// Float BLAS and LAPACK routines
#define scopy_fortran FortranCInterface_GLOBAL (scopy,SCOPY)
void scopy_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy);
#define saxpy_fortran FortranCInterface_GLOBAL (saxpy,SAXPY)
void saxpy_fortran(Integer* m,float* alpha,float* x,Integer* incx,float* y,
    Integer* incy);
#define sscal_fortran FortranCInterface_GLOBAL (sscal,SSCAL)
void sscal_fortran(Integer* m,float* alpha,float* x,Integer* incx);
#define sdot_fortran FortranCInterface_GLOBAL (sdot,SDOT)
float sdot_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy);
#define ssyr2k_fortran FortranCInterface_GLOBAL (ssyr2k,SSYR2K)
void ssyr2k_fortran(char* uplo,char* trans,Integer* n,Integer* k,float* alpha,
    float* A,Integer* lda,float* B,Integer* ldb,float* beta,
    float* C,Integer* ldc);
#define ssyevr_fortran FortranCInterface_GLOBAL (ssyevr,SSYEVR)
void ssyevr_fortran(char* jobz,char* range,char* uplo,Integer* n,float *A,
    Integer* lda,float* vl,float* vu,Integer* il,Integer* iu,float* abstol,
    Integer* m,float* w,float* z,Integer* ldz,Integer* isuppz,float* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define sstemr_fortran FortranCInterface_GLOBAL (sstemr,SSTEMR)
void sstemr_fortran(char* jobz,char* range,Integer* n,float *D,float *E,
    float* vl,float* vu,Integer* il,Integer* iu,Integer* m,float* w,float* z,
    Integer* ldz,Integer* nzc,Integer* isuppz,Integer* tryrac,float* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define sstevr_fortran FortranCInterface_GLOBAL (sstevr,SSTEVR)
void sstevr_fortran(char* jobz,char* range,Integer* n,float *D,float *E,
    float* vl,float* vu,Integer* il,Integer* iu,float* abstol,Integer* m,
    float* w,float* z,Integer* ldz,Integer* isuppz,float* work,Integer* lwork,
    Integer* iwork,Integer* liwork,Integer* info);
#define slamch_fortran FortranCInterface_GLOBAL (slamch,SLAMCH)
float slamch_fortran(char* cmach);
#define sgemm_fortran FortranCInterface_GLOBAL (sgemm,SGEMM)
void sgemm_fortran(char* transa,char* transb,Integer* m,Integer* n,Integer* k,
    float* alpha, float* A, Integer* lda,float* B,Integer* ldb,
    float* beta,float* C,Integer *ldc);
#define ssymm_fortran FortranCInterface_GLOBAL (ssymm,SSYMM)
void ssymm_fortran(char* side,char* uplo,Integer* m,Integer* n,float* alpha,
    float* A, Integer* lda,float* B,Integer* ldb,float* beta,float* C,
    Integer *ldc);
#define ssymv_fortran FortranCInterface_GLOBAL (ssymv,SSYMV)
void ssymv_fortran(char* uplo,Integer* n,float* alpha,float* A, Integer* lda,
    float* x,Integer* incx,float* beta,float* y,Integer *incy);
#define spotrf_fortran FortranCInterface_GLOBAL (spotrf,SPOTRF)
void spotrf_fortran(char* uplo,Integer *n,float* A,Integer* lda,Integer* info);
#define strtri_fortran FortranCInterface_GLOBAL (strtri,STRTRI)
void strtri_fortran(char* uplo,char* diag,Integer* n,float* A,Integer* lda,
    Integer* info);
#define srotg_fortran FortranCInterface_GLOBAL (srotg,SROTG)
void srotg_fortran(float* a,float* b,float* c,float *s);
#define srot_fortran FortranCInterface_GLOBAL (srot,SROT)
void srot_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy,
    float* c,float *s);
#define stpsv_fortran FortranCInterface_GLOBAL (stpsv,STPSV)
void stpsv_fortran(char* uplo,char* trans,char* diag,Integer* n,float* Ap,
    float* x,Integer* incx);
}

namespace peopt {

    template <>
    void copy <double> (Integer n,const double* x,Integer incx,double* y,
        Integer incy
    ) {
        dcopy_fortran(&n,const_cast <double*> (x),&incx,y,&incy);
    }

    template <>
    void copy <float> (Integer n,const float* x,Integer incx,float* y,
        Integer incy
    ) {
        scopy_fortran(&n,const_cast <float*> (x),&incx,y,&incy);
    }

    template <>
    void axpy <double> (
        Integer n,double alpha,const double* x,Integer incx,double* y,
        Integer incy
    ) {
        daxpy_fortran(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
    }
    
    template <>
    void axpy <float> (
        Integer n,float alpha,const float* x,Integer incx,float* y,
        Integer incy
    ) {
        saxpy_fortran(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
    }
    
    template <>
    void scal <double> (Integer n,double alpha,double* x,Integer incx) {
        dscal_fortran(&n,&alpha,x,&incx);
    }

    template <>
    void scal <float> (Integer n,float alpha,float* x,Integer incx) {
        sscal_fortran(&n,&alpha,x,&incx);
    }

    template <>
    double dot <double> (
        Integer n,const double* x,Integer incx,const double* y,Integer incy
    ) {
        return ddot_fortran(&n,const_cast <double*> (x),&incx,
            const_cast <double*> (y),&incy);
    }

    template <>
    float dot <float> (
        Integer n,const float* x,Integer incx,const float* y,Integer incy
    ) {
        return sdot_fortran(&n,const_cast <float*> (x),&incx,
            const_cast <float*> (y),&incy);
    }

    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,double alpha,
        const double* A,Integer lda,const double* B,Integer ldb,double beta,
        double* C,Integer ldc
    ) {
        dsyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
            const_cast <double*> (B),&ldb,&beta,C,&ldc);
    }

    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,float alpha,
        const float* A,Integer lda,const float* B,Integer ldb,float beta,
        float* C,Integer ldc
    ) {
        ssyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
            const_cast <float*> (B),&ldb,&beta,C,&ldc);
    }

    template <>
    void syevr(char jobz,char range,char uplo,Integer n,double *A,Integer lda,
        double vl,double vu,Integer il,Integer iu,double abstol,Integer& m,
        double* w,double* z,Integer ldz,Integer* isuppz,double* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info
    ) {
        dsyevr_fortran(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
            &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
    }

    template <>
    void syevr(char jobz,char range,char uplo,Integer n,float *A,Integer lda,
        float vl,float vu,Integer il,Integer iu,float abstol,Integer& m,
        float* w,float* z,Integer ldz,Integer* isuppz,float* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info
    ) {
        ssyevr_fortran(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
            &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stemr(char jobz,char range,Integer n,double *D,double *E,double vl,
        double vu,Integer il,Integer iu,Integer& m,double* w,double* z,
        Integer ldz,Integer nzc,Integer* isuppz,Integer& tryrac,double* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info
    ) {
        dstemr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
            isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stemr(char jobz,char range,Integer n,float *D,float *E,float vl,
        float vu,Integer il,Integer iu,Integer& m,float* w,float* z,
        Integer ldz,Integer nzc,Integer* isuppz,Integer& tryrac,float* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info
    ) {
        sstemr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
            isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stevr(char jobz,char range,Integer n,double *D,double *E,double vl,
        double vu,Integer il,Integer iu,double abstol,Integer& m,double* w,
        double* z,Integer ldz,Integer* isuppz,double* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info
    ) {
        dstevr_fortran(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,z,&ldz,
            isuppz,work,&lwork,iwork,&liwork,&info);
    }
    
    template <>
    void stevr(char jobz,char range,Integer n,float *D,float *E,float vl,
        float vu,Integer il,Integer iu,float abstol,Integer& m,float* w,
        float* z,Integer ldz,Integer* isuppz,float* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info
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
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,double alpha,
        const double* A,Integer lda,const double* B,Integer ldb,double beta,
        double* C,Integer ldc
    ) {
        dgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }
    
    template <>
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,float alpha,
        const float* A,Integer lda,const float* B,Integer ldb,float beta,
        float* C,Integer ldc
    ) {
        sgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,double alpha,const double* A,
        Integer lda,const double* B,Integer ldb,double beta,double* C,
        Integer ldc
    ) {
        dsymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,float alpha,const float* A,
        Integer lda,const float* B,Integer ldb,float beta,float* C,
        Integer ldc
    ) {
        ssymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }
    
    template <>
    void symv(char uplo,Integer n,double alpha,const double* A,Integer lda,
        const double* x,Integer incx,double beta,double* y,Integer incy
    ) {
        dsymv_fortran(&uplo,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
    }
    
    template <>
    void symv(char uplo,Integer n,float alpha,const float* A,Integer lda,
        const float* x,Integer incx,float beta,float* y,Integer incy
    ) {
        ssymv_fortran(&uplo,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
    }

    template <>
    void potrf(char uplo,Integer n,double* A,Integer lda,Integer& info) {
        dpotrf_fortran(&uplo,&n,A,&lda,&info);
    }
    
    template <>
    void potrf(char uplo,Integer n,float* A,Integer lda,Integer& info) {
        spotrf_fortran(&uplo,&n,A,&lda,&info);
    }

    template <>
    void trtri(
        char uplo,char diag,Integer n,double* A,Integer lda,Integer& info
    ) {
        dtrtri_fortran(&uplo,&diag,&n,A,&lda,&info);
    }
    
    template <>
    void trtri(
        char uplo,char diag,Integer n,float* A,Integer lda,Integer& info
    ) {
        strtri_fortran(&uplo,&diag,&n,A,&lda,&info);
    }

    template <>
    void rotg <double> (double a,double b,double& c,double& s) {
        drotg_fortran(&a,&b,&c,&s);
    }
    
    template <>
    void rotg <float> (float a,float b,float& c,float& s) {
        srotg_fortran(&a,&b,&c,&s);
    }
    
    template <>
    void rot <double> (
        Integer n,const double* x,Integer incx,double* y,Integer incy,
        double c,double s
    ){
        drot_fortran(&n,const_cast <double*> (x),&incx,y,&incy,&c,&s);
    }
    
    template <>
    void rot <float> (
        Integer n,const float* x,Integer incx,float* y,Integer incy,
        float c,float s
    ){
        srot_fortran(&n,const_cast <float*> (x),&incx,y,&incy,&c,&s);
    }
    
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,const double* Ap,double* x,
        Integer incx
    ) {
        dtpsv_fortran(&uplo,&trans,&diag,&n,const_cast <double*> (Ap),x,&incx);
    }
    
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,const float* Ap,float* x,
        Integer incx
    ) {
        stpsv_fortran(&uplo,&trans,&diag,&n,const_cast <float*> (Ap),x,&incx);
    }
    
    // Indexing function for matrices.  This assumes that the upper left 
    // corner has index (1,1).  Basically, we're using the indexing scheme
    // of math and Fortran and not of C.
    Natural ijtok(Natural i,Natural j,Natural m){
        return (i-1)+(j-1)*m;
    }

    // Indexing for packed storage.  This assumes the upper left corner has
    // index (1,1).
    Natural ijtokp(Natural i,Natural j) {
        return i+j*(j-1)/Natural(2)-Natural(1);
    }

    // Indexing for vectors.  Assumes the first index is 1.
    Natural itok(Natural i) {
        return i-Natural(1);
    }


}
