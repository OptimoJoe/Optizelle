/*
Copyright 2013-2014 OptimoJoe.

For the full copyright notice, see LICENSE.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Joseph Young (joe@optimojoe.com)
*/

#include "optizelle/linalg.h"
#include "FortranCInterface.h"
using Optizelle::Integer;
extern "C" {
// Double and float BLAS and LAPACK routines
#define drotg_fortran FortranCInterface_GLOBAL (drotg,DROTG)
void drotg_fortran(double* a,double* b,double* c,double *s);
#define srotg_fortran FortranCInterface_GLOBAL (srotg,SROTG)
void srotg_fortran(float* a,float* b,float* c,float *s);

#define drot_fortran FortranCInterface_GLOBAL (drot,DROT)
void drot_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy,
    double* c,double *s);
#define srot_fortran FortranCInterface_GLOBAL (srot,SROT)
void srot_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy,
    float* c,float *s);

#define dcopy_fortran FortranCInterface_GLOBAL (dcopy,DCOPY)
void dcopy_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy);
#define scopy_fortran FortranCInterface_GLOBAL (scopy,SCOPY)
void scopy_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy);

#define daxpy_fortran FortranCInterface_GLOBAL (daxpy,DAXPY)
void daxpy_fortran(Integer* m,double* alpha,double* x,Integer* incx,double* y,
    Integer* incy);
#define saxpy_fortran FortranCInterface_GLOBAL (saxpy,SAXPY)
void saxpy_fortran(Integer* m,float* alpha,float* x,Integer* incx,float* y,
    Integer* incy);

#define dscal_fortran FortranCInterface_GLOBAL (dscal,DSCAL)
void dscal_fortran(Integer* m,double* alpha,double* x,Integer* incx);
#define sscal_fortran FortranCInterface_GLOBAL (sscal,SSCAL)
void sscal_fortran(Integer* m,float* alpha,float* x,Integer* incx);

#define ddot_fortran FortranCInterface_GLOBAL (ddot,DDOT)
double ddot_fortran(Integer* m,double* x,Integer* incx,double* y,Integer* incy);
#define sdot_fortran FortranCInterface_GLOBAL (sdot,SDOT)
float sdot_fortran(Integer* m,float* x,Integer* incx,float* y,Integer* incy);

#define dgemv_fortran FortranCInterface_GLOBAL (dgemv,DGEMV)
void dgemv_fortran(char* trans,Integer* m,Integer* n,double* alpha,double* A,
    Integer* lda,double* x,Integer* incx,double* beta,double* y,Integer *incy);
#define sgemv_fortran FortranCInterface_GLOBAL (sgemv,SGEMV)
void sgemv_fortran(char* trans,Integer* m,Integer* n,float* alpha,float* A,
    Integer* lda,float* x,Integer* incx,float* beta,float* y,Integer *incy);

#define dsymv_fortran FortranCInterface_GLOBAL (dsymv,DSYMV)
void dsymv_fortran(char* uplo,Integer* n,double* alpha,double* A, Integer* lda,
    double* x,Integer* incx,double* beta,double* y,Integer *incy);
#define ssymv_fortran FortranCInterface_GLOBAL (ssymv,SSYMV)
void ssymv_fortran(char* uplo,Integer* n,float* alpha,float* A, Integer* lda,
    float* x,Integer* incx,float* beta,float* y,Integer *incy);

#define dspmv_fortran FortranCInterface_GLOBAL (dspmv,DSPMV)
void dspmv_fortran(char* uplo,Integer* n,double* alpha,double* Ap,
    double* x,Integer* incx,double* beta,double* y,Integer *incy);
#define sspmv_fortran FortranCInterface_GLOBAL (sspmv,SSPMV)
void sspmv_fortran(char* uplo,Integer* n,float* alpha,float* Ap,
    float* x,Integer* incx,float* beta,float* y,Integer *incy);

#define dtrsv_fortran FortranCInterface_GLOBAL (dtrsv,DTPSV)
void dtrsv_fortran(char* uplo,char* trans,char* diag,Integer* n,double* A,
    Integer* lda,double* x,Integer* incx);
#define strsv_fortran FortranCInterface_GLOBAL (strsv,STPSV)
void strsv_fortran(char* uplo,char* trans,char* diag,Integer* n,float* A,
    Integer* lda,float* x,Integer* incx);

#define dtpsv_fortran FortranCInterface_GLOBAL (dtpsv,DTPSV)
void dtpsv_fortran(char* uplo,char* trans,char* diag,Integer* n,double* Ap,
    double* x,Integer* incx);
#define stpsv_fortran FortranCInterface_GLOBAL (stpsv,STPSV)
void stpsv_fortran(char* uplo,char* trans,char* diag,Integer* n,float* Ap,
    float* x,Integer* incx);

#define dgemm_fortran FortranCInterface_GLOBAL (dgemm,DGEMM)
void dgemm_fortran(char* transa,char* transb,Integer* m,Integer* n,Integer* k,
    double* alpha, double* A, Integer* lda,double* B,Integer* ldb,
    double* beta,double* C,Integer *ldc);
#define sgemm_fortran FortranCInterface_GLOBAL (sgemm,SGEMM)
void sgemm_fortran(char* transa,char* transb,Integer* m,Integer* n,Integer* k,
    float* alpha, float* A, Integer* lda,float* B,Integer* ldb,
    float* beta,float* C,Integer *ldc);

#define dsymm_fortran FortranCInterface_GLOBAL (dsymm,DSYMM)
void dsymm_fortran(char* side,char* uplo,Integer* m,Integer* n,double* alpha,
    double* A, Integer* lda,double* B,Integer* ldb,double* beta,double* C,
    Integer *ldc);
#define ssymm_fortran FortranCInterface_GLOBAL (ssymm,SSYMM)
void ssymm_fortran(char* side,char* uplo,Integer* m,Integer* n,float* alpha,
    float* A, Integer* lda,float* B,Integer* ldb,float* beta,float* C,
    Integer *ldc);

#define dsyr2k_fortran FortranCInterface_GLOBAL (dsyr2k,DSYR2K)
void dsyr2k_fortran(char* uplo,char* trans,Integer* n,Integer* k,double* alpha,
    double* A,Integer* lda,double* B,Integer* ldb,double* beta,
    double* C,Integer* ldc);
#define ssyr2k_fortran FortranCInterface_GLOBAL (ssyr2k,SSYR2K)
void ssyr2k_fortran(char* uplo,char* trans,Integer* n,Integer* k,float* alpha,
    float* A,Integer* lda,float* B,Integer* ldb,float* beta,
    float* C,Integer* ldc);

#define dlamch_fortran FortranCInterface_GLOBAL (dlamch,DLAMCH)
double dlamch_fortran(char* cmach);
#define slamch_fortran FortranCInterface_GLOBAL (slamch,SLAMCH)
float slamch_fortran(char* cmach);

#define dsyevr_fortran FortranCInterface_GLOBAL (dsyevr,DSYEVR)
void dsyevr_fortran(char* jobz,char* range,char* uplo,Integer* n,double *A,
    Integer* lda,double* vl,double* vu,Integer* il,Integer* iu,double* abstol,
    Integer* m,double* w,double* z,Integer* ldz,Integer* isuppz,double* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define ssyevr_fortran FortranCInterface_GLOBAL (ssyevr,SSYEVR)
void ssyevr_fortran(char* jobz,char* range,char* uplo,Integer* n,float *A,
    Integer* lda,float* vl,float* vu,Integer* il,Integer* iu,float* abstol,
    Integer* m,float* w,float* z,Integer* ldz,Integer* isuppz,float* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);

#define dstemr_fortran FortranCInterface_GLOBAL (dstemr,DSTEMR)
void dstemr_fortran(char* jobz,char* range,Integer* n,double *D,double *E,
    double* vl,double* vu,Integer* il,Integer* iu,Integer* m,double* w,
    double* z,Integer* ldz,Integer* nzc,Integer* isuppz,Integer* tryrac,
    double* work,Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define sstemr_fortran FortranCInterface_GLOBAL (sstemr,SSTEMR)
void sstemr_fortran(char* jobz,char* range,Integer* n,float *D,float *E,
    float* vl,float* vu,Integer* il,Integer* iu,Integer* m,float* w,float* z,
    Integer* ldz,Integer* nzc,Integer* isuppz,Integer* tryrac,float* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);

#define dstevr_fortran FortranCInterface_GLOBAL (dstevr,DSTEVR)
void dstevr_fortran(char* jobz,char* range,Integer* n,double *D,double *E,
    double* vl,double* vu,Integer* il,Integer* iu,double* abstol,Integer* m,
    double* w,double* z,Integer* ldz,Integer* isuppz,double* work,
    Integer* lwork,Integer* iwork,Integer* liwork,Integer* info);
#define sstevr_fortran FortranCInterface_GLOBAL (sstevr,SSTEVR)
void sstevr_fortran(char* jobz,char* range,Integer* n,float *D,float *E,
    float* vl,float* vu,Integer* il,Integer* iu,float* abstol,Integer* m,
    float* w,float* z,Integer* ldz,Integer* isuppz,float* work,Integer* lwork,
    Integer* iwork,Integer* liwork,Integer* info);

#define dspevx_fortran FortranCInterface_GLOBAL (dspevx,DSPEVX)
void dspevx_fortran(char* jobz,char* range,char* uplo,Integer* n,double* Ap,
    double* vl,double* vu,Integer* il,Integer* iu,double* abstol,
    Integer* m,double* w,double* z,Integer* ldz,double* work,
    Integer* iwork,Integer* ifail,Integer* info);
#define sspevx_fortran FortranCInterface_GLOBAL (sspevx,SSPEVX)
void sspevx_fortran(char* jobz,char* range,char* uplo,Integer* n,float* Ap,
    float* vl,float* vu,Integer* il,Integer* iu,float* abstol,
    Integer* m,float* w,float* z,Integer* ldz,float* work,
    Integer* iwork,Integer* ifail,Integer* info);

#define dtptrs_fortran FortranCInterface_GLOBAL (dtptrs,DTPTRS)
void dtptrs_fortran(char* uplo,char* trans,char* diag,Integer *n,Integer* nrhs,
    double* Ap,double* B,Integer* ldb,Integer* info);
#define stptrs_fortran FortranCInterface_GLOBAL (stptrs,STPTRS)
void stptrs_fortran(char* uplo,char* trans,char* diag,Integer *n,Integer* nrhs,
    float* Ap,float* B,Integer* ldb,Integer* info);

#define dtprfs_fortran FortranCInterface_GLOBAL (dtprfs,DTPRFS)
void dtprfs_fortran(char* uplo,char* trans,char* diag,Integer *n,Integer* nrhs,
    double* Ap,double* B,Integer* ldb,double* X,Integer* ldx,double* ferr,
    double* berr,double* work,Integer* iwork,Integer* info);
#define stprfs_fortran FortranCInterface_GLOBAL (stprfs,STPRFS)
void stprfs_fortran(char* uplo,char* trans,char* diag,Integer *n,Integer* nrhs,
    float* Ap,float* B,Integer* ldb,float* X,Integer* ldx,float* ferr,
    float* berr,float* work,Integer* iwork,Integer* info);

#define dtrcon_fortran FortranCInterface_GLOBAL (dtrcon,DTRCON)
void dtrcon_fortran(char* norm,char* uplo,char* diag,Integer *n,
    double* A,Integer* lda,double* rcond,double* work,Integer* iwork,
    Integer* info);
#define strcon_fortran FortranCInterface_GLOBAL (strcon,STRCON)
void strcon_fortran(char* norm,char* uplo,char* diag,Integer *n,
    float* A,Integer* lda,float* rcond,float* work,Integer* iwork,
    Integer* info);

#define dtpcon_fortran FortranCInterface_GLOBAL (dtpcon,DTPCON)
void dtpcon_fortran(char* norm,char* uplo,char* diag,Integer *n,
    double* Ap,double* rcond,double* work,Integer* iwork,Integer* info);
#define stpcon_fortran FortranCInterface_GLOBAL (stpcon,STPCON)
void stpcon_fortran(char* norm,char* uplo,char* diag,Integer *n,
    float* Ap,float* rcond,float* work,Integer* iwork,Integer* info);

#define dpotrf_fortran FortranCInterface_GLOBAL (dpotrf,DPOTRF)
void dpotrf_fortran(char* uplo,Integer *n,double* A,Integer* lda,Integer* info);
#define spotrf_fortran FortranCInterface_GLOBAL (spotrf,SPOTRF)
void spotrf_fortran(char* uplo,Integer *n,float* A,Integer* lda,Integer* info);

#define dpotri_fortran FortranCInterface_GLOBAL (dpotri,DPOTRF)
void dpotri_fortran(char* uplo,Integer *n,double* A,Integer* lda,Integer* info);
#define spotri_fortran FortranCInterface_GLOBAL (spotri,SPOTRF)
void spotri_fortran(char* uplo,Integer *n,float* A,Integer* lda,Integer* info);

#define dpftrf_fortran FortranCInterface_GLOBAL (dpftrf,DPFTRF)
void dpftrf_fortran(char* transr,char* uplo,Integer* n,double* Arf,
    Integer* info);
#define spftrf_fortran FortranCInterface_GLOBAL (spftrf,SPFTRF)
void spftrf_fortran(char* transr,char* uplo,Integer* n,float* Arf,
    Integer* info);

#define dtrtri_fortran FortranCInterface_GLOBAL (dtrtri,DTRTRI)
void dtrtri_fortran(char* uplo,char* diag,Integer* n,double* A,Integer* lda,
    Integer* info);
#define strtri_fortran FortranCInterface_GLOBAL (strtri,STRTRI)
void strtri_fortran(char* uplo,char* diag,Integer* n,float* A,Integer* lda,
    Integer* info);

#define dspgst_fortran FortranCInterface_GLOBAL (dspgst,DSPGST)
void dspgst_fortran(Integer* itype,char* uplo,Integer* n,double* Ap,double* Bp,
    Integer* info);
#define sspgst_fortran FortranCInterface_GLOBAL (sspgst,SSPGST)
void sspgst_fortran(Integer* itype,char* uplo,Integer* n,float* Ap,float* Bp,
    Integer* info);

#define dgeqrf_fortran FortranCInterface_GLOBAL (dgeqrf,DGEQRF)
void dgeqrf_fortran(Integer* m,Integer* n,double* A,Integer* lda,double* tau,
    double* work,Integer* lwork,Integer* info);
#define sgeqrf_fortran FortranCInterface_GLOBAL (sgeqrf,SGEQRF)
void sgeqrf_fortran(Integer* m,Integer* n,float* A,Integer* lda,float* tau,
    float* work,Integer* lwork,Integer* info);

#define dorgqr_fortran FortranCInterface_GLOBAL (dorgqr,DORGQR)
void dorgqr_fortran(Integer* m,Integer* n,Integer* k,double* A,Integer* lda,
    double* tau,double* work,Integer* lwork,Integer* info);
#define sorgqr_fortran FortranCInterface_GLOBAL (sorgqr,SORGQR)
void sorgqr_fortran(Integer* m,Integer* n,Integer* k,float* A,Integer* lda,
    float* tau,float* work,Integer* lwork,Integer* info);

#define dtrttf_fortran FortranCInterface_GLOBAL (dtrttf,DTRTTF)
void dtrttf_fortran(char* transr,char* uplo,Integer* n,double* A,Integer* lda,
    double* Arf,Integer* info);
#define strttf_fortran FortranCInterface_GLOBAL (strttf,STRTTF)
void strttf_fortran(char* transr,char* uplo,Integer* n,float* A,Integer* lda,
    float* Arf,Integer* info);

#define dtrttp_fortran FortranCInterface_GLOBAL (dtrttp,DTRTTP)
void dtrttp_fortran(char* uplo,Integer* n,double* A,Integer* lda,double* Ap,
    Integer* info);
#define strttp_fortran FortranCInterface_GLOBAL (strttp,STRTTP)
void strttp_fortran(char* uplo,Integer* n,float* A,Integer* lda,float* Ap,
    Integer* info);

#define dtfttr_fortran FortranCInterface_GLOBAL (dtfttr,DTFTTR)
void dtfttr_fortran(char* transr,char* uplo,Integer* n,double* Arf,
    double* A,Integer* lda,Integer* info);
#define stfttr_fortran FortranCInterface_GLOBAL (stfttr,STFTTR)
void stfttr_fortran(char* transr,char* uplo,Integer* n,float* Arf,
    float* A,Integer* lda,Integer* info);

#define dtfttp_fortran FortranCInterface_GLOBAL (dtfttp,DTFTTP)
void dtfttp_fortran(char* transr,char* uplo,Integer* n,double* Arf,
    double* Ap,Integer* info);
#define stfttp_fortran FortranCInterface_GLOBAL (stfttp,STFTTP)
void stfttp_fortran(char* transr,char* uplo,Integer* n,float* Arf,
    float* Ap,Integer* info);

#define dtpttr_fortran FortranCInterface_GLOBAL (dtpttr,DTPTTR)
void dtpttr_fortran(char* uplo,Integer* n,double* Ap,double* A,Integer* lda,
    Integer* info);
#define stpttr_fortran FortranCInterface_GLOBAL (stpttr,STPTTR)
void stpttr_fortran(char* uplo,Integer* n,float* Ap,float* A,Integer* lda,
    Integer* info);
}

namespace Optizelle {

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
        Integer n,double* x,Integer incx,double* y,Integer incy,
        double c,double s
    ){
        drot_fortran(&n,x,&incx,y,&incy,&c,&s);
    }
    template <>
    void rot <float> (
        Integer n,float* x,Integer incx,float* y,Integer incy,
        float c,float s
    ){
        srot_fortran(&n,x,&incx,y,&incy,&c,&s);
    }

    template <>
    void copy <double> (Integer n,double const * const x,Integer incx,double* y,
        Integer incy
    ) {
        dcopy_fortran(&n,const_cast <double*> (x),&incx,y,&incy);
    }
    template <>
    void copy <float> (Integer n,float const * const x,Integer incx,float* y,
        Integer incy
    ) {
        scopy_fortran(&n,const_cast <float*> (x),&incx,y,&incy);
    }

    template <>
    void axpy <double> (
        Integer n,double alpha,double const * const x,Integer incx,double* y,
        Integer incy
    ) {
        daxpy_fortran(&n,&alpha,const_cast <double*> (x),&incx,y,&incy);
    }
    template <>
    void axpy <float> (
        Integer n,float alpha,float const * const x,Integer incx,float* y,
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
        Integer n,double const * const x,Integer incx,double const * const y,
        Integer incy
    ) {
        return ddot_fortran(&n,const_cast <double*> (x),&incx,
            const_cast <double*> (y),&incy);
    }
    template <>
    float dot <float> (
        Integer n,float const * const x,Integer incx,float const * const y,
        Integer incy
    ) {
        return sdot_fortran(&n,const_cast <float*> (x),&incx,
            const_cast <float*> (y),&incy);
    }

    template <>
    void gemv(char trans,Integer m,Integer n,double alpha,double const *const A,
        Integer lda,double const * const x,Integer incx,double beta,double* y,
        Integer incy
    ) {
        dgemv_fortran(&trans,&m,&n,&alpha,const_cast <double*> (A),&lda,
            const_cast <double*> (x),&incx,&beta,y,&incy);
    }
    template <>
    void gemv(char trans,Integer m,Integer n,float alpha,float const * const A,
        Integer lda,float const * const x,Integer incx,float beta,float* y,
        Integer incy
    ) {
        sgemv_fortran(&trans,&m,&n,&alpha,const_cast <float*> (A),&lda,
            const_cast <float*> (x),&incx,&beta,y,&incy);
    }
    
    template <>
    void symv(char uplo,Integer n,double alpha,double const*const A,Integer lda,
        double const * const x,Integer incx,double beta,double* y,Integer incy
    ) {
        dsymv_fortran(&uplo,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (x),&incx,&beta,y,&incy); 
    }
    template <>
    void symv(char uplo,Integer n,float alpha,float const * const A,Integer lda,
        float const * const x,Integer incx,float beta,float* y,Integer incy
    ) {
        ssymv_fortran(&uplo,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (x),&incx,&beta,y,&incy); 
    }

    template <>
    void spmv(char uplo,Integer n,double alpha,double const * const Ap,
        double const * const x,Integer incx,double beta,double* y,Integer incy
    ) {
        dspmv_fortran(&uplo,&n,&alpha,const_cast <double*> (Ap),
            const_cast <double*> (x),&incx,&beta,y,&incy);
    }
    template <>
    void spmv(char uplo,Integer n,float alpha,float const * const Ap,
        float const * const x,Integer incx,float beta,float* y,Integer incy
    ) {
        sspmv_fortran(&uplo,&n,&alpha,const_cast <float*> (Ap),
            const_cast <float*> (x),&incx,&beta,y,&incy);
    }
    
    template <>
    void trsv(
        char uplo,char trans,char diag,Integer n,double const * const A,
        Integer lda,double* x,Integer incx
    ) {
        dtrsv_fortran(&uplo,&trans,&diag,&n,const_cast <double*> (A),&lda,x,
            &incx);
    }
    template <>
    void trsv(
        char uplo,char trans,char diag,Integer n,float const * const A,
        Integer lda,float* x,Integer incx
    ) {
        strsv_fortran(&uplo,&trans,&diag,&n,const_cast <float*> (A),&lda,x,
            &incx);
    }
    
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,double const * const Ap,
        double* x,Integer incx
    ) {
        dtpsv_fortran(&uplo,&trans,&diag,&n,const_cast <double*> (Ap),x,&incx);
    }
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,float const * const Ap,
        float* x,Integer incx
    ) {
        stpsv_fortran(&uplo,&trans,&diag,&n,const_cast <float*> (Ap),x,&incx);
    }

    template <>
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc
    ) {
        dgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }
    template <>
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc
    ) {
        sgemm_fortran(&transa,&transb,&m,&n,&k,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc
    ) {
        dsymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <double*> (A),
            &lda,const_cast <double*> (B),&ldb,&beta,C,&ldc); 
    }
    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc
    ) {
        ssymm_fortran(&side,&uplo,&m,&n,&alpha,const_cast <float*> (A),
            &lda,const_cast <float*> (B),&ldb,&beta,C,&ldc); 
    }

    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc
    ) {
        dsyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <double*> (A),&lda,
            const_cast <double*> (B),&ldb,&beta,C,&ldc);
    }
    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc
    ) {
        ssyr2k_fortran(&uplo,&trans,&n,&k,&alpha,const_cast <float*> (A),&lda,
            const_cast <float*> (B),&ldb,&beta,C,&ldc);
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
    void spevx(char jobz,char range,char uplo,Integer n,double* Ap,
        double vl,double vu,Integer il,Integer iu,double abstol,Integer& m,
        double* w,double* z,Integer ldz,double* work,Integer* iwork,
        Integer* ifail,Integer& info
    ) {
        dspevx_fortran(&jobz,&range,&uplo,&n,Ap,&vl,&vu,&il,&iu,&abstol,&m,w,z,
            &ldz,work,iwork,ifail,&info);
    }
    template <>
    void spevx(char jobz,char range,char uplo,Integer n,float* Ap,
        float vl,float vu,Integer il,Integer iu,float abstol,Integer& m,
        float* w,float* z,Integer ldz,float* work,Integer* iwork,
        Integer* ifail,Integer& info
    ) {
        sspevx_fortran(&jobz,&range,&uplo,&n,Ap,&vl,&vu,&il,&iu,&abstol,&m,w,z,
            &ldz,work,iwork,ifail,&info);
    }
    
    template <>
    void tptrs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        double const * const Ap,double* B,Integer ldb,Integer& info
    ) {
        dtptrs_fortran(&uplo,&trans,&diag,&n,&nrhs,const_cast <double*>(Ap),
            B,&ldb,&info);
    }
    template <>
    void tptrs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        float const * const Ap,float* B,Integer ldb,Integer& info
    ) {
        stptrs_fortran(&uplo,&trans,&diag,&n,&nrhs,const_cast <float*>(Ap),
            B,&ldb,&info);
    }
    
    template <>
    void tprfs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        double const * const Ap,double const * const B,Integer ldb,
        double const * const X,Integer ldx,double* ferr,double* berr,
        double* work,Integer* iwork,Integer& info
    ) {
        dtprfs_fortran(&uplo,&trans,&diag,&n,&nrhs,const_cast <double*>(Ap),
            const_cast <double*> (B),&ldb,const_cast <double*> (X),&ldx,
            ferr,berr,work,iwork,&info);
    }
    template <>
    void tprfs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        float const * const Ap,float const * const B,Integer ldb,
        float const * const X,Integer ldx,float* ferr,float* berr,
        float* work,Integer* iwork,Integer& info
    ) {
        stprfs_fortran(&uplo,&trans,&diag,&n,&nrhs,const_cast <float*>(Ap),
            const_cast <float*> (B),&ldb,const_cast <float*> (X),&ldx,
            ferr,berr,work,iwork,&info);
    }
    
    template <>
    void trcon(char norm,char uplo,char diag,Integer n,
        double const * const A,Integer lda,double& rcond,double* work,
        Integer* iwork,Integer& info
    ) {
        dtrcon_fortran(&norm,&uplo,&diag,&n,const_cast <double*>(A),&lda,&rcond,
            work,iwork,&info); 
    }
    template <>
    void trcon(char norm,char uplo,char diag,Integer n,
        float const * const A,Integer lda,float& rcond,float* work,
        Integer* iwork,Integer& info
    ) {
        strcon_fortran(&norm,&uplo,&diag,&n,const_cast <float*>(A),&lda,&rcond,
            work,iwork,&info); 
    }
    
    template <>
    void tpcon(char norm,char uplo,char diag,Integer n,
        double const * const Ap,double& rcond,double* work,Integer* iwork,
        Integer& info
    ) {
        dtpcon_fortran(&norm,&uplo,&diag,&n,const_cast <double*>(Ap),&rcond,
            work,iwork,&info); 
    }
    template <>
    void tpcon(char norm,char uplo,char diag,Integer n,
        float const * const Ap,float& rcond,float* work,Integer* iwork,
        Integer& info
    ) {
        stpcon_fortran(&norm,&uplo,&diag,&n,const_cast <float*>(Ap),&rcond,
            work,iwork,&info); 
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
    void potri(char uplo,Integer n,double* A,Integer lda,Integer& info) {
        dpotri_fortran(&uplo,&n,A,&lda,&info);
    }
    template <>
    void potri(char uplo,Integer n,float* A,Integer lda,Integer& info) {
        spotri_fortran(&uplo,&n,A,&lda,&info);
    }

    template <>
    void pftrf(char transr,char uplo,Integer n,double* Arf,Integer& info) {
        dpftrf_fortran(&transr,&uplo,&n,Arf,&info);
    }
    template <>
    void pftrf(char transr,char uplo,Integer n,float* Arf,Integer& info) {
        spftrf_fortran(&transr,&uplo,&n,Arf,&info);
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
    void spgst(Integer itype,char uplo,Integer n,double* Ap,
        double const * const Bp,Integer& info
    ) {
        dspgst_fortran(&itype,&uplo,&n,Ap,const_cast <double*> (Bp),&info);
    }
    template <>
    void spgst(Integer itype,char uplo,Integer n,float* Ap,
        float const * const Bp,Integer& info
    ) {
        sspgst_fortran(&itype,&uplo,&n,Ap,const_cast <float*> (Bp),&info);
    }

    template <>
    void geqrf(Integer m,Integer n,double* A,Integer lda,double* tau,
        double* work,Integer lwork,Integer& info
    ) {
        dgeqrf_fortran(&m,&n,A,&lda,tau,work,&lwork,&info);
    }
    template <>
    void geqrf(Integer m,Integer n,float* A,Integer lda,float* tau,
        float* work,Integer lwork,Integer& info
    ) {
        sgeqrf_fortran(&m,&n,A,&lda,tau,work,&lwork,&info);
    }

    template <>
    void orgqr(Integer m,Integer n,Integer k,double* A,Integer lda,double* tau,
        double* work,Integer lwork,Integer& info
    ) {
        dorgqr_fortran(&m,&n,&k,A,&lda,tau,work,&lwork,&info);
    }
    template <>
    void orgqr(Integer m,Integer n,Integer k,float* A,Integer lda,float* tau,
        float* work,Integer lwork,Integer& info
    ) {
        sorgqr_fortran(&m,&n,&k,A,&lda,tau,work,&lwork,&info);
    }

    template <>
    void trttf(char transr,char uplo,Integer n,double const * const A,
        Integer lda,double* Arf,Integer& info
    ) {
        dtrttf_fortran(&transr,&uplo,&n,const_cast <double*> (A),&lda,Arf,
            &info);
    }
    template <>
    void trttf(char transr,char uplo,Integer n,float const * const A,
        Integer lda,float* Arf,Integer& info
    ) {
        strttf_fortran(&transr,&uplo,&n,const_cast <float*> (A),&lda,Arf,
            &info);
    }

    template <>
    void trttp(char uplo,Integer n,double const * const A,Integer lda,
        double* Ap,Integer& info
    ) {
        dtrttp_fortran(&uplo,&n,const_cast <double*> (A),&lda,Ap,&info);
    }
    template <>
    void trttp(char uplo,Integer n,float const * const A,Integer lda,
        float* Ap,Integer& info
    ) {
        strttp_fortran(&uplo,&n,const_cast <float*> (A),&lda,Ap,&info);
    }

    template <>
    void tfttr(char transr,char uplo,Integer n,double const * const Arf,
        double* A,Integer lda,Integer& info
    ) {
        dtfttr_(&transr,&uplo,&n,const_cast <double*> (Arf),A,&lda,&info);
    }
    template <>
    void tfttr(char transr,char uplo,Integer n,float const * const Arf,
        float* A,Integer lda,Integer& info
    ) {
        stfttr_fortran(&transr,&uplo,&n,const_cast <float*> (Arf),A,&lda,&info);
    }

    template <>
    void tfttp(char transr,char uplo,Integer n,double const * const Arf,
        double* Ap,Integer& info
    ) {
        dtfttp_fortran(&transr,&uplo,&n,const_cast <double*> (Arf),Ap,&info);
    }
    template <>
    void tfttp(char transr,char uplo,Integer n,float const * const Arf,
        float* Ap,Integer& info
    ) {
        stfttp_fortran(&transr,&uplo,&n,const_cast <float*> (Arf),Ap,&info);
    }

    template <>
    void tpttr(char uplo,Integer n,double const * const Ap,double* A,
        Integer lda,Integer& info
    ) {
        dtpttr_fortran(&uplo,&n,const_cast <double*> (Ap),A,&lda,&info);
    }
    template <>
    void tpttr(char uplo,Integer n,float const * const Ap,float* A,
        Integer lda,Integer& info
    ) {
        stpttr_fortran(&uplo,&n,const_cast <float*> (Ap),A,&lda,&info);
    }
    
    // Indexing function for matrices.  This assumes that the upper left 
    // corner has index (1,1).  Basically, we're using the indexing scheme
    // of math and Fortran and not of C.
    Natural ijtok(Natural const & i,Natural const & j,Natural const & m){
        return i-1+(j-1)*m;
    }

    // Indexing for packed storage where i<=j.  This assumes the upper left
    // corner has index (1,1).
    Natural ijtokp(Natural const & i,Natural const & j) {
        return i-1 + j*(j-1)/2;
    }
    
    // Indexing formula for symmetric matrices in RPF format where i<=j.
    Natural ijtokrf(Natural const & i,Natural const & j,Natural const & m) {
        // Determine k where m = 2k when m is even or m = 2k+1 when m is odd
        const Natural k = m/2; 

        // Return the index
        return i<=k && j<=k ? j-1 + k+1 + (2*k+1)*(i-1) : i-1 + (j-1-k)*m;
    }

    // Indexing for vectors.  Assumes the first index is 1.
    Natural itok(Natural const & i) {
        return i-Natural(1);
    }
    
    namespace TruncatedStop{
        // Converts the truncated CG stopping condition to a string 
        std::string to_string(t const & trunc_stop){
            switch(trunc_stop){
            case NotConverged:
                return "NotConverged";
            case NegativeCurvature:
                return "NegativeCurvature";
            case RelativeErrorSmall:
                return "RelativeErrorSmall";
            case MaxItersExceeded:
                return "MaxItersExceeded";
            case TrustRegionViolated:
                return "TrustRegionViolated";
            case NanDetected:
                return "NanDetected";
            case NonProjector:
                return "NonProjector";
            case NonSymmetric:
                return "NonSymmetric";
            case LossOfOrthogonality:
                return "LossOfOrthogonality";
            case InvalidTrustRegionOffset:
                return "InvalidTrustRegionOffset";
            case TooManyFailedSafeguard:
                return "TooManyFailedSafeguard";
            case ObjectiveIncrease:
                return "ObjectiveIncrease";
            default:
                throw;
            }
        }
        
        // Converts a string to a truncated CG stopping condition
        t from_string(std::string const & trunc_stop){
            if(trunc_stop=="NotConverged")
                return NotConverged;
            else if(trunc_stop=="NegativeCurvature")
                return NegativeCurvature;
            else if(trunc_stop=="RelativeErrorSmall")
                return RelativeErrorSmall;
            else if(trunc_stop=="MaxItersExceeded")
                return MaxItersExceeded;
            else if(trunc_stop=="TrustRegionViolated")
                return TrustRegionViolated;
            else if(trunc_stop=="NanDetected")
                return NanDetected;
            else if(trunc_stop=="NonProjector")
                return NonProjector;
            else if(trunc_stop=="NonSymmetric")
                return NonSymmetric;
            else if(trunc_stop=="LossOfOrthogonality")
                return LossOfOrthogonality;
            else if(trunc_stop=="InvalidTrustRegionOffset")
                return InvalidTrustRegionOffset;
            else if(trunc_stop=="TooManyFailedSafeguard")
                return TooManyFailedSafeguard;
            else if(trunc_stop=="ObjectiveIncrease")
                return ObjectiveIncrease;
            else
                throw;
        }

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name) {
            if( name=="NotConverged" ||
                name=="NegativeCurvature" ||
                name=="RelativeErrorSmall" ||
                name=="MaxItersExceeded" ||
                name=="TrustRegionViolated" ||
                name=="NanDetected" ||
                name=="NonProjector" ||
                name=="NonSymmetric" ||
                name=="LossOfOrthogonality" ||
                name=="InvalidTrustRegionOffset" ||
                name=="TooManyFailedSafeguard" ||
                name=="ObjectiveIncrease"
            )
                return true;
            else
                return false;
        }
    }
}
