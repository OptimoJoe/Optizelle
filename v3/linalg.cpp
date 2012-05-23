#include "linalg.h"
extern "C" {
#ifdef CAPS_FORTRAN
    #ifdef NOUNDER_FORTRAN
        void DCOPY(int* m,double* x,int* incx,double* y,int* incy);
        void DAXPY(
            int* m,double* alpha,double* x,int* incx,double* y,int* incy);
        void DSCAL(int* m,double* alpha,double* x,int* incx);
        double DDOT(int* m,double* x,int* incx,double* y,int* incy);
        void DSYR2K(char* uplo,char* trans,int* n,int* k,double* alpha,
            double* A,int* lda,double* B,int* ldb,double* beta,
            double* C,int* ldc);
        void DSYEVR(char* jobz,char* range,char* uplo,int* n,double *A,
            int* lda,double* vl,double* vu,int* il,int* iu,double* abstol,
            int* m,double* w,double* z,int* ldz,int* isuppz,double* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void DSTEMR(char* jobz,char* range,int* n,double *D,double *E,
            double* vl,double* vu,int* il,int* iu,int* m,double* w,double* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,double* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void DGEMM(char* transa,char* transb,int* m,int* n,int* k,
            double* alpha, double* A, int* lda,double* B,int* ldb,
            double* beta,double* C,int *ldc);
        void DSYMM(char* side,char* uplo,int* m,int* n,double* alpha,
            double* A, int* lda,double* B,int* ldb,double* beta,double* C,
            int *ldc);
        void DSYMV(char* uplo,int* n,double* alpha,double* A, int* lda,
            double* x,int* incx,double* beta,double* y,int *incy);
        void DPOTRF(char* uplo,int *n,double* A,int* lda,int* info);
        void DTRTRI(char* uplo,char* diag,int* n,double* A,int* lda,int* info);
        
        void SCOPY(int* m,float* x,int* incx,float* y,int* incy);
        void SAXPY(
            int* m,float* alpha,float* x,int* incx,float* y,int* incy);
        void SSCAL(int* m,float* alpha,float* x,int* incx);
        float SDOT(int* m,float* x,int* incx,float* y,int* incy);
        void SSYR2K(char* uplo,char* trans,int* n,int* k,float* alpha,
            float* A,int* lda,float* B,int* ldb,float* beta,
            float* C,int* ldc);
        void SSYEVR(char* jobz,char* range,char* uplo,int* n,float *A,
            int* lda,float* vl,float* vu,int* il,int* iu,float* abstol,
            int* m,float* w,float* z,int* ldz,int* isuppz,float* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void SSTEMR(char* jobz,char* range,int* n,float *D,float *E,
            float* vl,float* vu,int* il,int* iu,int* m,float* w,float* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,float* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void SGEMM(char* transa,char* transb,int* m,int* n,int* k,
            float* alpha, float* A, int* lda,float* B,int* ldb,
            float* beta,float* C,int *ldc);
        void SSYMM(char* side,char* uplo,int* m,int* n,float* alpha,
            float* A, int* lda,float* B,int* ldb,float* beta,float* C,
            int *ldc);
        void SSYMV(char* uplo,int* n,float* alpha,float* A, int* lda,
            float* x,int* incx,float* beta,float* y,int *incy);
        void SPOTRF(char* uplo,int *n,float* A,int* lda,int* info);
        void STRTRI(char* uplo,char* diag,int* n,float* A,int* lda,int* info);
    #else
        void DCOPY_(int* m,double* x,int* incx,double* y,int* incy);
        void DAXPY_(
            int* m,double* alpha,double* x,int* incx,double* y,int* incy);
        void DSCAL_(int* m,double* alpha,double* x,int* incx);
        double DDOT_(int* m,double* x,int* incx,double* y,int* incy);
        void DSYR2K_(char* uplo,char* trans,int* n,int* k,double* alpha,
            double* A,int* lda,double* B,int* ldb,double* beta,
            double* C,int* ldc);
        void DSYEVR_(char* jobz,char* range,char* uplo,int* n,double *A,
            int* lda,double* vl,double* vu,int* il,int* iu,double* abstol,
            int* m,double* w,double* z,int* ldz,int* isuppz,double* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void DSTEMR_(char* jobz,char* range,int* n,double *D,double *E,
            double* vl,double* vu,int* il,int* iu,int* m,double* w,double* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,double* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void DGEMM_(char* transa,char* transb,int* m,int* n,int* k,
            double* alpha, double* A, int* lda,double* B,int* ldb,
            double* beta,double* C,int *ldc);
        void DSYMM_(char* side,char* uplo,int* m,int* n,double* alpha,
            double* A, int* lda,double* B,int* ldb,double* beta,double* C,
            int *ldc);
        void DSYMV_(char* uplo,int* n,double* alpha,double* A, int* lda,
            double* x,int* incx,double* beta,double* y,int *incy);
        void DPOTRF_(char* uplo,int *n,double* A,int* lda,int* info);
        void DTRTRI_(char* uplo,char* diag,int* n,double* A,int* lda,int* info);

        void SCOPY_(int* m,float* x,int* incx,float* y,int* incy);
        void SAXPY_(
            int* m,float* alpha,float* x,int* incx,float* y,int* incy);
        void SSCAL_(int* m,float* alpha,float* x,int* incx);
        float SDOT_(int* m,float* x,int* incx,float* y,int* incy);
        void SSYR2K_(char* uplo,char* trans,int* n,int* k,float* alpha,
            float* A,int* lda,float* B,int* ldb,float* beta,
            float* C,int* ldc);
        void SSYEVR_(char* jobz,char* range,char* uplo,int* n,float *A,
            int* lda,float* vl,float* vu,int* il,int* iu,float* abstol,
            int* m,float* w,float* z,int* ldz,int* isuppz,float* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void SSTEMR_(char* jobz,char* range,int* n,float *D,float *E,
            float* vl,float* vu,int* il,int* iu,int* m,float* w,float* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,float* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void SGEMM_(char* transa,char* transb,int* m,int* n,int* k,
            float* alpha, float* A, int* lda,float* B,int* ldb,
            float* beta,float* C,int *ldc);
        void SSYMM_(char* side,char* uplo,int* m,int* n,float* alpha,
            float* A, int* lda,float* B,int* ldb,float* beta,float* C,
            int *ldc);
        void SSYMV_(char* uplo,int* n,float* alpha,float* A, int* lda,
            float* x,int* incx,float* beta,float* y,int *incy);
        void SPOTRF_(char* uplo,int *n,float* A,int* lda,int* info);
        void STRTRI_(char* uplo,char* diag,int* n,float* A,int* lda,int* info);
    #endif
#else
    #ifdef NOUNDER_FORTRAN
        void dcopy(int* m,double* x,int* incx,double* y,int* incy);
        void daxpy(
            int* m,double* alpha,double* x,int* incx,double* y,int* incy);
        void dscal(int* m,double* alpha,double* x,int* incx);
        double ddot(int* m,double* x,int* incx,double* y,int* incy);
        void dsyr2k(char* uplo,char* trans,int* n,int* k,double* alpha,
            double* A,int* lda,double* B,int* ldb,double* beta,
            double* C,int* ldc);
        void dsyevr(char* jobz,char* range,char* uplo,int* n,double *A,
            int* lda,double* vl,double* vu,int* il,int* iu,double* abstol,
            int* m,double* w,double* z,int* ldz,int* isuppz,double* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void dstemr(char* jobz,char* range,int* n,double *D,double *E,
            double* vl,double* vu,int* il,int* iu,int* m,double* w,double* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,double* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void dgemm(char* transa,char* transb,int* m,int* n,int* k,
            double* alpha, double* A, int* lda,double* B,int* ldb,
            double* beta,double* C,int *ldc);
        void dsymm(char* side,char* uplo,int* m,int* n,double* alpha,
            double* A, int* lda,double* B,int* ldb,double* beta,double* C,
            int *ldc);
        void dsymv(char* uplo,int* n,double* alpha,double* A, int* lda,
            double* x,int* incx,double* beta,double* y,int *incy);
        void dpotrf(char* uplo,int *n,double* A,int* lda,int* info);
        void dtrtri(char* uplo,char* diag,int* n,double* A,int* lda,int* info);

        void scopy(int* m,float* x,int* incx,float* y,int* incy);
        void saxpy(
            int* m,float* alpha,float* x,int* incx,float* y,int* incy);
        void sscal(int* m,float* alpha,float* x,int* incx);
        float sdot(int* m,float* x,int* incx,float* y,int* incy);
        void ssyr2k(char* uplo,char* trans,int* n,int* k,float* alpha,
            float* A,int* lda,float* B,int* ldb,float* beta,
            float* C,int* ldc);
        void ssyevr(char* jobz,char* range,char* uplo,int* n,float *A,
            int* lda,float* vl,float* vu,int* il,int* iu,float* abstol,
            int* m,float* w,float* z,int* ldz,int* isuppz,float* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void sstemr(char* jobz,char* range,int* n,float *D,float *E,
            float* vl,float* vu,int* il,int* iu,int* m,float* w,float* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,float* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void sgemm(char* transa,char* transb,int* m,int* n,int* k,
            float* alpha, float* A, int* lda,float* B,int* ldb,
            float* beta,float* C,int *ldc);
        void ssymm(char* side,char* uplo,int* m,int* n,float* alpha,
            float* A, int* lda,float* B,int* ldb,float* beta,float* C,
            int *ldc);
        void ssymv(char* uplo,int* n,float* alpha,float* A, int* lda,
            float* x,int* incx,float* beta,float* y,int *incy);
        void spotrf(char* uplo,int *n,float* A,int* lda,int* info);
        void strtri(char* uplo,char* diag,int* n,float* A,int* lda,int* info);
    #else
        void dcopy_(int* m,double* x,int* incx,double* y,int* incy);
        void daxpy_(
            int* m,double* alpha,double* x,int* incx,double* y,int* incy);
        void dscal_(int* m,double* alpha,double* x,int* incx);
        double ddot_(int* m,double* x,int* incx,double* y,int* incy);
        void dsyr2k_(char* uplo,char* trans,int* n,int* k,double* alpha,
            double* A,int* lda,double* B,int* ldb,double* beta,
            double* C,int* ldc);
        void dsyevr_(char* jobz,char* range,char* uplo,int* n,double *A,
            int* lda,double* vl,double* vu,int* il,int* iu,double* abstol,
            int* m,double* w,double* z,int* ldz,int* isuppz,double* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void dstemr_(char* jobz,char* range,int* n,double *D,double *E,
            double* vl,double* vu,int* il,int* iu,int* m,double* w,double* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,double* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void dgemm_(char* transa,char* transb,int* m,int* n,int* k,
            double* alpha, double* A, int* lda,double* B,int* ldb,
            double* beta,double* C,int *ldc);
        void dsymm_(char* side,char* uplo,int* m,int* n,double* alpha,
            double* A, int* lda,double* B,int* ldb,double* beta,double* C,
            int *ldc);
        void dsymv_(char* uplo,int* n,double* alpha,double* A, int* lda,
            double* x,int* incx,double* beta,double* y,int *incy);
        void dpotrf_(char* uplo,int *n,double* A,int* lda,int* info);
        void dtrtri_(char* uplo,char* diag,int* n,double* A,int* lda,int* info);
        
        void scopy_(int* m,float* x,int* incx,float* y,int* incy);
        void saxpy_(
            int* m,float* alpha,float* x,int* incx,float* y,int* incy);
        void sscal_(int* m,float* alpha,float* x,int* incx);
        float sdot_(int* m,float* x,int* incx,float* y,int* incy);
        void ssyr2k_(char* uplo,char* trans,int* n,int* k,float* alpha,
            float* A,int* lda,float* B,int* ldb,float* beta,
            float* C,int* ldc);
        void ssyevr_(char* jobz,char* range,char* uplo,int* n,float *A,
            int* lda,float* vl,float* vu,int* il,int* iu,float* abstol,
            int* m,float* w,float* z,int* ldz,int* isuppz,float* work,
            int* lwork,int* iwork,int* liwork,int* info);
        void sstemr_(char* jobz,char* range,int* n,float *D,float *E,
            float* vl,float* vu,int* il,int* iu,int* m,float* w,float* z,
            int* ldz,int* nzc,int* isuppz,int* tryrac,float* work,int* lwork,
            int* iwork,int* liwork,int* info);
        void sgemm_(char* transa,char* transb,int* m,int* n,int* k,
            float* alpha, float* A, int* lda,float* B,int* ldb,
            float* beta,float* C,int *ldc);
        void ssymm_(char* side,char* uplo,int* m,int* n,float* alpha,
            float* A, int* lda,float* B,int* ldb,float* beta,float* C,
            int *ldc);
        void ssymv_(char* uplo,int* n,float* alpha,float* A, int* lda,
            float* x,int* incx,float* beta,float* y,int *incy);
        void spotrf_(char* uplo,int *n,float* A,int* lda,int* info);
        void strtri_(char* uplo,char* diag,int* n,float* A,int* lda,int* info);
    #endif
#endif
}

namespace peopt {

    template <>
    void copy <double> (int n,const double* x,int incx,double* y,int incy) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DCOPY(&n,const_cast <double*> (x),&incx,y,&incy);
            #else
                DCOPY_(&n,const_cast <double*> (x),&incx,y,&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dcopy(&n,const_cast <double*> (x),&incx,y,&incy);
            #else
                dcopy_(&n,const_cast <double*> (x),&incx,y,&incy);
            #endif
        #endif
    }

    template <>
    void copy <float> (int n,const float* x,int incx,float* y,int incy) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SCOPY(&n,const_cast <float*> (x),&incx,y,&incy);
            #else
                SCOPY_(&n,const_cast <float*> (x),&incx,y,&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                scopy(&n,const_cast <float*> (x),&incx,y,&incy);
            #else
                scopy_(&n,const_cast <float*> (x),&incx,y,&incy);
            #endif
        #endif
    }

    template <>
    void axpy <double> (
        int n,double alpha,const double* x,int incx,double* y,int incy
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DAXPY(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
            #else
                DAXPY_(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                daxpy(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
            #else
                daxpy_(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
            #endif
        #endif
    }
    
    template <>
    void axpy <float> (
        int n,float alpha,const float* x,int incx,float* y,int incy
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SAXPY(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
            #else
                SAXPY_(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                saxpy(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
            #else
                saxpy_(&n,&alpha,const_cast <float*> (x),&incx,y,&incy);
            #endif
        #endif
    }
    
    template <>
    void scal <double> (int n,double alpha,double* x,int incx) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSCAL(&n,&alpha,x,&incx);
            #else
                DSCAL_(&n,&alpha,x,&incx);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dscal(&n,&alpha,x,&incx);
            #else
                dscal_(&n,&alpha,x,&incx);
            #endif
        #endif
    }

    template <>
    void scal <float> (int n,float alpha,float* x,int incx) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSCAL(&n,&alpha,x,&incx);
            #else
                SSCAL_(&n,&alpha,x,&incx);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                sscal(&n,&alpha,x,&incx);
            #else
                sscal_(&n,&alpha,x,&incx);
            #endif
        #endif
    }

    template <>
    double dot <double> (
        int n,const double* x,int incx,const double* y,int incy
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                return DDOT(&n,const_cast <double*> (x),&incx,
                    const_cast <double*> (y),&incy);
            #else
                return DDOT_(&n,const_cast <double*> (x),&incx,
                    const_cast <double*> (y),&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                return ddot(&n,const_cast <double*> (x),&incx,
                    const_cast <double*> (y),&incy);
            #else
                return ddot_(&n,const_cast <double*> (x),&incx,
                    const_cast <double*> (y),&incy);
            #endif
        #endif
    }

    template <>
    float dot <float> (int n,const float* x,int incx,const float* y,int incy) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                return SDOT(&n,const_cast <float*> (x),&incx,
                    const_cast <float*> (y),&incy);
            #else
                return SDOT_(&n,const_cast <float*> (x),&incx,
                    const_cast <float*> (y),&incy);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                return sdot(&n,const_cast <float*> (x),&incx,
                    const_cast <float*> (y),&incy);
            #else
                return sdot_(&n,const_cast <float*> (x),&incx,
                    const_cast <float*> (y),&incy);
            #endif
        #endif
    }

    template <>
    void syr2k(char uplo,char trans,int n,int k,double alpha,
        const double* A,int lda,const double* B,int ldb,double beta,
        double* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSYR2K(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
                    const_cast <double*> (B),&ldb,&beta,C,&ldc);
            #else
                DSYR2K_(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
                    const_cast <double*> (B),&ldb,&beta,C,&ldc);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dsyr2k(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
                    const_cast <double*> (B),&ldb,&beta,C,&ldc);
            #else
                dsyr2k_(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
                    const_cast <double*> (B),&ldb,&beta,C,&ldc);
            #endif
        #endif
    }

    template <>
    void syr2k(char uplo,char trans,int n,int k,float alpha,
        const float* A,int lda,const float* B,int ldb,float beta,
        float* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSYR2K(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
                    const_cast <float*> (B),&ldb,&beta,C,&ldc);
            #else
                SSYR2K_(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
                    const_cast <float*> (B),&ldb,&beta,C,&ldc);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                ssyr2k(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
                    const_cast <float*> (B),&ldb,&beta,C,&ldc);
            #else
                ssyr2k_(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
                    const_cast <float*> (B),&ldb,&beta,C,&ldc);
            #endif
        #endif
    }

    template <>
    void syevr(char jobz,char range,char uplo,int n,double *A,int lda,
        double vl,double vu,int il,int iu,double abstol,int& m,
        double* w,double* z,int ldz,int* isuppz,double* work,int lwork,
        int* iwork,int liwork,int& info
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSYEVR(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #else
                DSYEVR_(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dsyevr(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #else
                dsyevr_(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #endif
        #endif
    }

    template <>
    void syevr(char jobz,char range,char uplo,int n,float *A,int lda,
        float vl,float vu,int il,int iu,float abstol,int& m,
        float* w,float* z,int ldz,int* isuppz,float* work,int lwork,
        int* iwork,int liwork,int& info
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSYEVR(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #else
                SSYEVR_(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                ssyevr(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #else
                ssyevr_(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
                    &m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            #endif
        #endif
    }
    
    template <>
    void stemr(char jobz,char range,int n,double *D,double *E,double vl,
        double vu,int il,int iu,int& m,double* w,double* z,int ldz,int nzc,
        int* isuppz,int& tryrac,double* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSTEMR(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #else
                DSTEMR_(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dstemr(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #else
                dstemr_(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #endif
        #endif
    }
    
    template <>
    void stemr(char jobz,char range,int n,float *D,float *E,float vl,
        float vu,int il,int iu,int& m,float* w,float* z,int ldz,int nzc,
        int* isuppz,int& tryrac,float* work,int lwork,int* iwork,int liwork,
        int& info
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSTEMR(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #else
                SSTEMR_(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                sstemr(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #else
                sstemr_(&jobz,&range,&n,D,E,&vl,&vu,&il,&iu,&m,w,z,&ldz,&nzc,
                    isuppz,&tryrac,work,&lwork,iwork,&liwork,&info);
            #endif
        #endif
    }

    template <>
    void gemm(char transa,char transb,int m,int n,int k,double alpha,
        const double* A,int lda,const double* B,int ldb,double beta,
        double* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DGEMM(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #else
                DGEMM_(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dgemm(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #else
                dgemm_(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #endif
    }
    
    template <>
    void gemm(char transa,char transb,int m,int n,int k,float alpha,
        const float* A,int lda,const float* B,int ldb,float beta,
        float* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SGEMM(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #else
                SGEMM_(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                sgemm(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #else
                sgemm_(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #endif
    }

    template <>
    void symm(char side,char uplo,int m,int n,double alpha,const double* A,
        int lda,const double* B,int ldb,double beta,double* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSYMM(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #else
                DSYMM_(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dsymm(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #else
                dsymm_(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #endif
    }

    template <>
    void symm(char side,char uplo,int m,int n,float alpha,const float* A,
        int lda,const float* B,int ldb,float beta,float* C,int ldc
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSYMM(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #else
                SSYMM_(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                ssymm(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #else
                ssymm_(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
            #endif
        #endif
    }
    
    template <>
    void symv(char uplo,int n,double alpha,const double* A,int lda,
        const double* x,int incx,double beta,double* y,int incy
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DSYMV(&uplo,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
            #else
                DSYMV_(&uplo,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dsymv(&uplo,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
            #else
                dsymv_(&uplo,&n,&alpha,const_cast <double*> (A),
                    &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
            #endif
        #endif
    }
    
    template <>
    void symv(char uplo,int n,float alpha,const float* A,int lda,
        const float* x,int incx,float beta,float* y,int incy
    ) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SSYMV(&uplo,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
            #else
                SSYMV_(&uplo,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                ssymv(&uplo,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
            #else
                ssymv_(&uplo,&n,&alpha,const_cast <float*> (A),
                    &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
            #endif
        #endif
    }

    template <>
    void potrf(char uplo,int n,double* A,int lda,int& info) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DPOTRF(&uplo,&n,A,&lda,&info);
            #else
                DPOTRF_(&uplo,&n,A,&lda,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dpotrf(&uplo,&n,A,&lda,&info);
            #else
                dpotrf_(&uplo,&n,A,&lda,&info);
            #endif
        #endif
    }
    
    template <>
    void potrf(char uplo,int n,float* A,int lda,int& info) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                SPOTRF(&uplo,&n,A,&lda,&info);
            #else
                SPOTRF_(&uplo,&n,A,&lda,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                spotrf(&uplo,&n,A,&lda,&info);
            #else
                spotrf_(&uplo,&n,A,&lda,&info);
            #endif
        #endif
    }

    template <>
    void trtri(char uplo,char diag,int n,double* A,int lda,int& info) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                DTRTRI(&uplo,&diag,&n,A,&lda,&info);
            #else
                DTRTRI_(&uplo,&diag,&n,A,&lda,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                dtrtri(&uplo,&diag,&n,A,&lda,&info);
            #else
                dtrtri_(&uplo,&diag,&n,A,&lda,&info);
            #endif
        #endif
    }
    
    template <>
    void trtri(char uplo,char diag,int n,float* A,int lda,int& info) {
        #ifdef CAPS_FORTRAN
            #ifdef NOUNDER_FORTRAN
                STRTRI(&uplo,&diag,&n,A,&lda,&info);
            #else
                STRTRI_(&uplo,&diag,&n,A,&lda,&info);
            #endif
        #else
            #ifdef NOUNDER_FORTRAN
                strtri(&uplo,&diag,&n,A,&lda,&info);
            #else
                strtri_(&uplo,&diag,&n,A,&lda,&info);
            #endif
        #endif
    }
    
    // Absolute value
    template <>
    double abs(double alpha) {
        return std::fabs(alpha);
    }
    template <>
    float abs(float alpha) {
        return std::fabs(alpha);
    }
    
    // Square root
    template <>
    double sqrt(double alpha) {
        return std::sqrt(alpha);
    }
    template <>
    float sqrt(float alpha) {
        return std::sqrt(alpha);
    }
    
    // Logarithm 
    template <>
    double log(double alpha) {
        return std::log(alpha);
    }
    template <>
    float log(float alpha) {
        return std::log(alpha);
    }
    
    // Indexing function for matrices.  This assumes that the upper right
    // corner has index (1,1).  Basically, we're using the indexing scheme
    // of math and Fortran and not of C.
    unsigned int ijtok(unsigned int i,unsigned int j,unsigned int m){
        return (i-1)+(j-1)*m;
    }
}
