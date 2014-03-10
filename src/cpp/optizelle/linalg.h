/*
Copyright 2013 OptimoJoe.

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

#ifndef LINALG_H 
#define LINALG_H 

#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <utility>
#include <string>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <random>

// Putting this into a class prevents its construction.  Essentially, we use
// this trick in order to create modules like in ML.  It also allows us to
// created templated namespaces.
#define NO_CONSTRUCTORS(Name) \
    Name() = delete; \
    Name(Name const &) = delete; \
    Name & operator = (Name const &) = delete; \
    ~Name() = delete;

// Disallows copying or assigning.  This is useful for things like the
// StateManipulator or FunctionModification classes.
#define NO_COPY_ASSIGNMENT(Name) \
    Name(Name const &) = delete; \
    Name & operator = (Name const &) = delete;

// Disallows copying, assigning, or default construction.  This is useful for
// for things like the State::t types.
#define NO_DEFAULT_COPY_ASSIGNMENT(Name) \
    Name() = delete; \
    Name(Name const &) = delete; \
    Name& operator = (Name const &) = delete;

//---Optizelle0---
namespace Optizelle {
//---Optizelle1---

    typedef size_t Natural;
    typedef ptrdiff_t Integer;

    template <typename Real>
    Real sq(Real const & x) {
        return x*x;
    }

    template <typename Real>
    int sgn(Real const & val) {
        return (Real(0) < val) - (val < Real(0));
    }
    
    template <typename Real>
    void xerbla(std::string srname,Integer info) {
        std::cout << " ** On entry to " << srname << " parameter number "
            << info << " had an illegal value" << std::endl;
    }

    template <typename Real>
    void rotg(Real a,Real b,Real& c,Real& s) {
        Real r,z;
        Real roe = b;
        if(fabs(a) > fabs(b)) roe = a;
        Real scale = fabs(a) + fabs(b);
        if(scale == Real(0.)) {
            c = Real(1.0);
            s = Real(0.0);
            r = Real(0.0);
            z = Real(0.0);
        } else {
            r = scale*sqrt(sq <Real> (a/scale) + sq <Real> (b/scale));
            r = sgn <Real> (roe)*r;
            c = a/r;
            s = b/r;
            z = Real(1.0);
            if(fabs(a) > fabs(b)) z=s;
            if(fabs(b) >= fabs(a) && c!=Real(0.)) z=Real(1.0)/c;
        }
        a=r;
        b=z;
    }
    template <>
    void rotg(double a,double b,double& c,double& s);
    template <>
    void rotg(float a,float b,float& c,float& s);

    template <typename Real>
    void rot(Integer n,Real* x,Integer incx,Real* y,Integer incy,Real c,Real s){
        if(n<=0) return;
        for(Integer i=0,ix=0,iy=0;
            i<n;
            i++,ix+=incx,iy+=incy
        ) {
            Real temp = c*x[ix] + s*y[iy];
            y[iy] = c*y[iy] - s*x[ix];
            x[ix] = temp;
        }
    }
    template <>
    void rot <double> (
        Integer n,double* x,Integer incx,double* y,Integer incy,
        double c,double s);
    template <>
    void rot <float> (
        Integer n,float* x,Integer incx,float* y,Integer incy,
        float c,float s);

    template <typename Real>
    void copy(Integer n,Real const * const x,Integer incx,Real* y,Integer incy){
        if(n<=0) return;
        for(Integer i=0,ix=0,iy=0;
            i<n;
            i++,ix+=incx,iy+=incy
        )
            y[iy]=x[ix];
    }
    template <>
    void copy <double> (Integer n,double const * const x,Integer incx,double* y,
        Integer incy);
    template <>
    void copy <float> (Integer n,float const * const x,Integer incx,float* y,
        Integer incy);
    
    template <typename Real>
    void axpy(Integer n,Real alpha,Real const * const x,Integer incx,
        Real* y,Integer incy);
    template <>
    void axpy <double> (
        Integer n,double alpha,double const * const x,Integer incx,double* y,
        Integer incy);
    template <>
    void axpy <float> (
        Integer n,float alpha,float const * const x,Integer incx,float* y,
        Integer incy);
    
    template <typename Real>
    void scal(Integer n,const Real alpha,Real* x,Integer incx) {
        if(n<=0) return;
        for(int i=0,ix=0;
            i<n;
            i++,ix+=incx
        )
            x[ix]*=alpha;
    }
    template <>
    void scal <double> (Integer n,double alpha,double* x,Integer incx);
    template <>
    void scal <float> (Integer n,float alpha,float* x,Integer incx);
    
    template <typename Real>
    Real dot(Integer n,Real const * const x,Integer incx,Real const * const y,
        Integer incy);
    template <>
    double dot <double> (
        Integer n,double const * const x,Integer incx,double const * const y,
        Integer incy);
    template <>
    float dot <float> (
        Integer n,float const * const x,Integer incx,float const * const y,
        Integer incy);
    

    template <typename Real>
    void gemv(char trans,Integer m,Integer n,Real alpha,Real const * const A,
        Integer lda,Real const * const x,Integer incx,Real beta,Real* y,
        Integer incy);
    template <>
    void gemv(char trans,Integer m,Integer n,double alpha,double const *const A,
        Integer lda,double const * const x,Integer incx,double beta,double* y,
        Integer incy);
    template <>
    void gemv(char trans,Integer m,Integer n,float alpha,float const * const A,
        Integer lda,float const * const x,Integer incx,float beta,float* y,
        Integer incy);
    
    template <typename Real>
    void symv(char uplo,Integer n,Real alpha,Real const * const A,Integer lda,
        Real const * const x,Integer incx,Real beta,Real* y,Integer incy);
    template <>
    void symv(char uplo,Integer n,double alpha,double const*const A,Integer lda,
        double const * const x,Integer incx,double beta,double* y,Integer incy);
    template <>
    void symv(char uplo,Integer n,float alpha,float const * const A,Integer lda,
        float const * const x,Integer incx,float beta,float* y,Integer incy);

    template <typename Real>
    void spmv(char uplo,Integer n,Real alpha,Real const * const Ap,
        Real const * const x,Integer incx,Real beta,Real* y,Integer incy);
    template <>
    void spmv(char uplo,Integer n,double alpha,double const * const Ap,
        double const * const x,Integer incx,double beta,double* y,Integer incy);
    template <>
    void spmv(char uplo,Integer n,float alpha,float const * const Ap,
        float const * const x,Integer incx,float beta,float* y,Integer incy);
    
    template <typename Real>
    void trsv(char uplo,char trans,char diag,Integer n,Real const * const A,
        Integer lda,Real* x,Integer incx); 
    template <>
    void trsv(
        char uplo,char trans,char diag,Integer n,double const * const A,
        Integer lda,double* x,Integer incx);
    template <>
    void trsv(
        char uplo,char trans,char diag,Integer n,float const * const A,
        Integer lda,float* x,Integer incx);
    
    // NOTE: this routine is not fully general.  It only implements what we
    // need.
    template <typename Real>
    void tpsv(char uplo,char trans,char diag,Integer n,Real const * const Ap,
        Real* x,Integer incx
    ) {
        // Test the input parameters.
        Integer info=0;
        if(uplo != 'U')
            info = 1;
        else if(trans != 'N')
            info = 2;
        else if(diag != 'U' && diag!= 'N')
            info = 3;
        else if(n<0)
            info = 4;
        else if(incx==0)
            info = 7;
        if(info !=0) {
            xerbla <Real> ("TPSV",info);
            return;
        }

        // Quick return if possible.
        if(n==0) return;

        bool nounit = diag == 'N';

        // Set up the start point in X if the increment is not unity. This
        // will be  ( N - 1 )*INCX  too small for descending loops.
        Integer kx = 0;
        if(incx < 0) kx = 0 - (n-1)*incx;
        else if(incx !=1) kx=0;

        // Start the operations. In this version the elements of AP are
        // accessed sequentially with one pass through AP.
        if(trans == 'N') {

            // Form  x := inv( A )*x.
            if(uplo=='U') {
                Integer kk = (n*(n+1))/2;
                if(incx==1) {
                    for(Integer j=n;j>=1;j--) {
                        if(x[j-1]!=Real(0.)) {
                            if(nounit) x[j-1] = x[j-1]/Ap[kk-1];
                            Real temp = x[j-1];
                            Integer k = kk-1;
                            for(Integer i=j-1;i>=1;i--) {
                                x[i-1] = x[i-1] - temp * Ap[k-1];
                                k--;
                            }
                        }
                        kk -= j;
                    }
                } else {
                    Integer jx = kx + (n-1)*incx;
                    for(Integer j=1;j>=1;j--) {
                        if(x[jx-1]!=Real(0.)) {
                            if(nounit) x[jx-1] = x[jx-1]/Ap[kk-1];
                            Real temp = x[jx-1];
                            Integer ix=jx;
                            for(Integer k=kk-1;k>=kk-j+1;k--) {
                                ix = ix - incx;
                                x[ix-1] = x[ix-1] - temp*Ap[k-1];
                            }
                        }
                        jx -= incx;
                        kk -= j;
                    }
                }
            }
        }
    }
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,double const * const Ap,
        double* x,Integer incx);
    template <>
    void tpsv(
        char uplo,char trans,char diag,Integer n,float const * const Ap,
        float* x,Integer incx);

    template <typename Real>
    void gemm(char transa,char transb,Integer m,Integer n,Integer k,Real alpha,
        Real const * const A,Integer lda,Real const * const B,Integer ldb,
        Real beta,Real* C,Integer ldc);
    template <>
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc);
    template <>
    void gemm(
        char transa,char transb,Integer m,Integer n,Integer k,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc);
    
    template <typename Real>
    void symm(char side,char uplo,Integer m,Integer n,Real alpha,
        Real const * const A,Integer lda,Real const * const B,Integer ldb,
        Real beta,Real* C,Integer ldc);
    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc);
    template <>
    void symm(
        char side,char uplo,Integer m,Integer n,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc);

    
    template <typename Real>
    void syr2k(char uplo,char trans,Integer n,Integer k,Real alpha,
        Real const * const A,Integer lda,Real const * const B,Integer ldb,
        Real beta,Real* C,Integer ldc);
    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,double alpha,
        double const * const A,Integer lda,double const * const B,Integer ldb,
        double beta,double* C,Integer ldc);
    template <>
    void syr2k(char uplo,char trans,Integer n,Integer k,float alpha,
        float const * const A,Integer lda,float const * const B,Integer ldb,
        float beta,float* C,Integer ldc);

    template <typename Real>
    Real lamch(char cmach);
    template <>
    double lamch(char cmach); 
    template <>
    float lamch(char cmach);

    template <typename Real>
    void syevr(char jobz,char range,char uplo,Integer n,Real *A,Integer lda,
        Real vl,Real vu,Integer il,Integer iu,Real abstol,Integer& m,
        Real* w,Real* z,Integer ldz,Integer* isuppz,Real* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info);
    template <>
    void syevr(char jobz,char range,char uplo,Integer n,double *A,Integer lda,
        double vl,double vu,Integer il,Integer iu,double abstol,Integer& m,
        double* w,double* z,Integer ldz,Integer* isuppz,double* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    template <>
    void syevr(char jobz,char range,char uplo,Integer n,float *A,Integer lda,
        float vl,float vu,Integer il,Integer iu,float abstol,Integer& m,
        float* w,float* z,Integer ldz,Integer* isuppz,float* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    
    template <typename Real>
    void stemr(char jobz,char range,Integer n,Real *D,Real *E,Real vl,Real vu,
        Integer il,Integer iu,Integer& m,Real* w,Real* z,Integer ldz,
        Integer nzc,Integer* isuppz,Integer& tryrac,Real* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    template <>
    void stemr(char jobz,char range,Integer n,double *D,double *E,double vl,
        double vu,Integer il,Integer iu,Integer& m,double* w,double* z,
        Integer ldz,Integer nzc,Integer* isuppz,Integer& tryrac,double* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    template <>
    void stemr(char jobz,char range,Integer n,float *D,float *E,float vl,
        float vu,Integer il,Integer iu,Integer& m,float* w,float* z,
        Integer ldz,Integer nzc,Integer* isuppz,Integer& tryrac,float* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    
    template <typename Real>
    void stevr(char jobz,char range,Integer n,Real *D,Real *E,Real vl,Real vu,
        Integer il,Integer iu,Real abstol, Integer& m,Real* w,Real* z,
        Integer ldz,Integer* isuppz,Real* work,Integer lwork,Integer* iwork,
        Integer liwork,Integer& info);
    template <>
    void stevr(char jobz,char range,Integer n,double *D,double *E,double vl,
        double vu,Integer il,Integer iu,double abstol,Integer& m,double* w,
        double* z,Integer ldz,Integer* isuppz,double* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info);
    template <>
    void stevr(char jobz,char range,Integer n,float *D,float *E,float vl,
        float vu,Integer il,Integer iu,float abstol,Integer& m,float* w,
        float* z,Integer ldz,Integer* isuppz,float* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info);

    template <typename Real>
    void spevx(char jobz,char range,char uplo,Integer n,Real* Ap,
        Real vl,Real vu,Integer il,Integer iu,Real abstol,Integer& m,
        Real* w,Real* z,Integer ldz,Real* work,Integer* iwork,
        Integer* ifail,Integer& info);
    template <>
    void spevx(char jobz,char range,char uplo,Integer n,double* Ap,
        double vl,double vu,Integer il,Integer iu,double abstol,Integer& m,
        double* w,double* z,Integer ldz,double* work,Integer* iwork,
        Integer* ifail,Integer& info);
    template <>
    void spevx(char jobz,char range,char uplo,Integer n,float* Ap,
        float vl,float vu,Integer il,Integer iu,float abstol,Integer& m,
        float* w,float* z,Integer ldz,float* work,Integer* iwork,
        Integer* ifail,Integer& info);

    template <typename Real>
    void tptrs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        Real const * const Ap,Real* B,Integer ldb,Integer& info);
    template <>
    void tptrs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        double const * const Ap,double* B,Integer ldb,Integer& info);
    template <>
    void tptrs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        float const * const Ap,float* B,Integer ldb,Integer& info);

    template <typename Real>
    void tprfs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        Real const * const Ap,Real const * const B,Integer ldb,
        Real const * const X,Integer ldx,Real* ferr,Real* berr,Real* work,
        Integer* iwork,Integer& info);
    template <>
    void tprfs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        double const * const Ap,double const * const B,Integer ldb,
        double const * const X,Integer ldx,double* ferr,double* berr,
        double* work,Integer* iwork,Integer& info);
    template <>
    void tprfs(char uplo,char trans,char diag,Integer n,Integer nrhs,
        float const * const Ap,float const * const B,Integer ldb,
        float const * const X,Integer ldx,float* ferr,float* berr,
        float* work,Integer* iwork,Integer& info);

    template <typename Real>
    void trcon(char norm,char uplo,char diag,Integer n,
        Real const * const A,Integer lda,Real& rcond,Real* work,Integer* iwork,
        Integer& info);
    template <>
    void trcon(char norm,char uplo,char diag,Integer n,
        double const * const A,Integer lda,double& rcond,double* work,
        Integer* iwork,Integer& info);
    template <>
    void trcon(char norm,char uplo,char diag,Integer n,
        float const * const A,Integer lda,float& rcond,float* work,
        Integer* iwork,Integer& info);

    template <typename Real>
    void tpcon(char norm,char uplo,char diag,Integer n,
        Real const * const Ap,Real& rcond,Real* work,Integer* iwork,
        Integer& info);
    template <>
    void tpcon(char norm,char uplo,char diag,Integer n,
        double const * const Ap,double& rcond,double* work,Integer* iwork,
        Integer& info);
    template <>
    void tpcon(char norm,char uplo,char diag,Integer n,
        float const * const Ap,float& rcond,float* work,Integer* iwork,
        Integer& info);

    template <typename Real>
    void potrf(char uplo,Integer n,Real* A,Integer lda,Integer& info);
    template <>
    void potrf(char uplo,Integer n,double* A,Integer lda,Integer& info);
    template <>
    void potrf(char uplo,Integer n,float* A,Integer lda,Integer& info);
    
    template <typename Real>
    void potri(char uplo,Integer n,Real* A,Integer lda,Integer& info);
    template <>
    void potri(char uplo,Integer n,double* A,Integer lda,Integer& info);
    template <>
    void potri(char uplo,Integer n,float* A,Integer lda,Integer& info);
    
    template <typename Real>
    void pftrf(char transr,char uplo,Integer n,Real* Arf,Integer& info);
    template <>
    void pftrf(char transr,char uplo,Integer n,double* Arf,Integer& info);
    template <>
    void pftrf(char transr,char uplo,Integer n,float* Arf,Integer& info);

    template <typename Real>
    void trtri(char uplo,char diag,Integer n,Real* A,Integer lda,Integer& info);
    template <>
    void trtri(
        char uplo,char diag,Integer n,double* A,Integer lda,Integer& info);
    template <>
    void trtri(
        char uplo,char diag,Integer n,float* A,Integer lda,Integer& info);

    template <typename Real>
    void spgst(Integer itype,char uplo,Integer n,Real* Ap,Real const * const Bp,
        Integer& info);
    template <>
    void spgst(Integer itype,char uplo,Integer n,double* Ap,
        double const * const Bp,Integer& info);
    template <>
    void spgst(Integer itype,char uplo,Integer n,float* Ap,
        float const * const Bp,Integer& info);

    template <typename Real>
    void geqrf(Integer m,Integer n,Real* A,Integer lda,Real* tau,
        Real* work,Integer lwork,Integer& info);
    template <>
    void geqrf(Integer m,Integer n,double* A,Integer lda,double* tau,
        double* work,Integer lwork,Integer& info);
    template <>
    void geqrf(Integer m,Integer n,float* A,Integer lda,float* tau,
        float* work,Integer lwork,Integer& info);

    template <typename Real>
    void orgqr(Integer m,Integer n,Integer k,Real* A,Integer lda,Real* tau,
        Real* work,Integer lwork,Integer& info);
    template <>
    void orgqr(Integer m,Integer n,Integer k,double* A,Integer lda,double* tau,
        double* work,Integer lwork,Integer& info);
    template <>
    void orgqr(Integer m,Integer n,Integer k,float* A,Integer lda,float* tau,
        float* work,Integer lwork,Integer& info);

    template <typename Real>
    void trttf(char transr,char uplo,Integer n,Real const * const A,Integer lda,
        Real* Arf,Integer& info);
    template <>
    void trttf(char transr,char uplo,Integer n,double const * const A,
        Integer lda,double* Arf,Integer& info);
    template <>
    void trttf(char transr,char uplo,Integer n,float const * const A,
        Integer lda,float* Arf,Integer& info);

    template <typename Real>
    void trttp(char uplo,Integer n,Real const * const A,Integer lda,
        Real* Ap,Integer& info);
    template <>
    void trttp(char uplo,Integer n,double const * const A,Integer lda,
        double* Ap,Integer& info);
    template <>
    void trttp(char uplo,Integer n,float const * const A,Integer lda,
        float* Ap,Integer& info);

    template <typename Real>
    void tfttr(char transr,char uplo,Integer n,Real const * const Arf,Real* A,
        Integer lda,Integer& info);
    template <>
    void tfttr(char transr,char uplo,Integer n,double const * const Arf,
        double* A,Integer lda,Integer& info);
    template <>
    void tfttr(char transr,char uplo,Integer n,float const * const Arf,
        float* A,Integer lda,Integer& info);

    template <typename Real>
    void tfttp(char transr,char uplo,Integer n,Real const * const Arf,Real* Ap,
        Integer& info);
    template <>
    void tfttp(char transr,char uplo,Integer n,double const * const Arf,
        double* Ap,Integer& info);
    void tfttp(char transr,char uplo,Integer n,float const * const Arf,
        float* Ap,Integer& info);

    template <typename Real>
    void tpttr(char uplo,Integer n,Real const * const A,Real* Ap,Integer lda,
        Integer& info);
    template <>
    void tpttr(char uplo,Integer n,double const * const Ap,double* A,
        Integer lda,Integer& info);
    template <>
    void tpttr(char uplo,Integer n,float const * const Ap,float* A,
        Integer lda,Integer& info);

    // Indexing function for matrices
    Natural ijtok(Natural const & i,Natural const & j,Natural const & m);

    // Indexing for packed storage where i<=j
    Natural ijtokp(Natural const & i,Natural const & j); 

    // Indexing formula for symmetric matrices in RPF format where i<=j.
    Natural ijtokrf(Natural const & i,Natural const & j,Natural const & m);

    // Indexing for vectors 
    Natural itok(Natural const & i);

    //---Operator0---
    // A linear operator specification, A : X->Y
    template <
        typename Real,
        template <typename> class X,
        template <typename> class Y
    >
    struct Operator {
        // Create some type shortcuts
        typedef typename X <Real>::Vector X_Vector;
        typedef typename Y <Real>::Vector Y_Vector;

        // y = A(x)
        virtual void eval(X_Vector const & x,Y_Vector &y) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~Operator() {}
    };
    //---Operator1---

    /* Given a Schur decomposition of A, A=V D V', solve the Sylvester equation
    
       A X + X A = B

    */
    template <typename Real>
    void sylvester(
        Natural const & m,
        Real const * const V,
        Real const * const D,
        Real const * const B,
        Real * X
    ) {

        // Find V' B V
        std::vector <Real> tmp(m*m);
        std::vector <Real> VtBV(m*m);
        // tmp <- B V
        symm <Real> ('L','U',m,m,Real(1.),&(B[0]),m,&(V[0]),m,Real(0.),
            &(tmp[0]),m); 
        // VtBV <- V' B V
        gemm <Real> ('T','N',m,m,m,Real(1.),&(V[0]),m,&(tmp[0]),m,Real(0.),
            &(VtBV[0]),m);

        // Solve for each column of X.  In theory, we only need half of these
        // elements since X is symmetric.
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for(Natural j=1;j<=m;j++) {
            for(Natural i=1;i<=j;i++) 
                X[ijtok(i,j,m)]=VtBV[ijtok(i,j,m)]/(D[i-1]+D[j-1]);
        }

        // Realransform the solution back, X = V X V'
        // tmp <- V X
        symm <Real> ('R','U',m,m,Real(1.),&(X[0]),m,&(V[0]),m,Real(0.),
            &(tmp[0]),m);
        // X <- V X V'
        gemm <Real> ('N','T',m,m,m,Real(1.),&(tmp[0]),m,&(V[0]),m,Real(0.),
            &(X[0]),m);
    }

    // Find a bound on the smallest eigenvalue of the given matrix A such
    // that lambda_min(A) < alpha where alpha is returned from this function.
    template <typename Real>
    Real lanczos(
        Natural const & m,
        Real const * const A,
        Natural const & max_iter,
        Real const & tol
    ) {
        // Create the initial Krylov vector
        std::vector <Real> v(m,Real(1./std::sqrt(Real(m))));

        // Get the next Krylov vector and orthgonalize it
        std::vector <Real> w(m);
        // w <- A v
        symv <Real> ('U',m,Real(1.),&(A[0]),m,&(v[0]),1,Real(0.),&(w[0]),1);
        // alpha[0] <- <Av,v>
        std::vector <Real> alpha;
        alpha.emplace_back(dot <Real> (m,&(w[0]),1,&(v[0]),1));
        // w <- Av - <Av,v> v
        axpy <Real> (m,-alpha[0],&(v[0]),1,&(w[0]),1);

        // Store the norm of the Arnoldi vector w in the off diagonal part of T.
        // By T, we mean the triagonal matrix such that A = Q T Q'.
        std::vector <Real> beta;
        beta.emplace_back(std::sqrt(dot <Real> (m,&(w[0]),1,&(w[0]),1)));

        // Allocate memory for solving an eigenvalue problem for the Ritz
        // values and vectors later.
        std::vector <Integer> isuppz;
        std::vector <Real> work(1);
        std::vector <Integer> iwork(1);
        Integer lwork=-1;
        Integer liwork=-1;
        Integer info;
        Integer nevals;
        std::vector <Real> W;
        std::vector <Real> Z;
        std::vector <Real> D;
        std::vector <Real> E;

        // Start Lanczos
        std::vector <Real> v_old(m);
        for(Natural i=0;i<max_iter;i++) {
            // Save the current Arnoldi vector
            copy <Real> (m,&(v[0]),1,&(v_old[0]),1);

            // Copy the candidate Arnoldi vector to the current Arnoldi vector
            copy <Real> (m,&(w[0]),1,&(v[0]),1);

            // Get the normalized version of this vector.  This is now a real
            // Arnoldi vector.
            scal <Real> (m,Real(1.)/beta[i],&(v[0]),1);

            // Get the new Arnoldi vector, w <- A v
            symv <Real> ('U',m,Real(1.),&(A[0]),m,&(v[0]),1,Real(0.),&(w[0]),1);

            // Orthogonalize against v_old and v using modified Gram-Schdmit.

            // First, we orthogonalize against v_old
            // w <- Av - <Av,v_old> v_old.  Due to symmetry, <Av,v_old>=beta.
            axpy <Real> (m,-beta[i],&(v_old[0]),1,&(w[0]),1);

            // Now, we orthogonalize against v
            // Find the Gram-Schmidt coefficient
            alpha.emplace_back(dot <Real> (m,&(w[0]),1,&(v[0]),1));
            // Orthogonlize w to v
            axpy <Real> (m,-alpha[i+1],&(v[0]),1,&(w[0]),1);

            // Store the norm of the Arnoldi vector w in the off diagonal part
            // of T.
            beta.emplace_back(std::sqrt(dot <Real> (m,&(w[0]),1,&(w[0]),1)));
   
            // Figure out the workspaces for the eigenvalues and eigenvectors
            Natural k=alpha.size();  // Size of the eigenvalue subproblem
            D.resize(alpha.size());
            copy <Real> (k,&(alpha[0]),1,&(D[0]),1);
            E.resize(beta.size());
            copy <Real> (k,&(beta[0]),1,&(E[0]),1);
            isuppz.resize(2*k);
            lwork=20*k;
            work.resize(lwork);
            liwork=10*k;
            iwork.resize(liwork);
            W.resize(k);
            Z.resize(k*k);
            Optizelle::stevr <Real> ('V','A',k,&(D[0]),&(E[0]),Real(0.),
                Real(0.),0,0,Optizelle::lamch <Real> ('S'),
                nevals,&(W[0]),&(Z[0]),k,&(isuppz[0]),&(work[0]),
                lwork,&(iwork[0]),liwork,info);

            // Find beta_i |s_{i1}| where s_{i1} is the last element
            // of the 1st Ritz vector, which corresponds to the smallest
            // Ritz value.
            Real err_est = fabs(Z[ijtok(k,1,k)])*beta[i+1];

            // Stop of the error estimates are small
            if(err_est < tol)
                break;
        }

        // Return the smallest Ritz value
        return W[0];
    }

    // Solve the symmetric eigenvalue problem A x = lambda x for the leftmost
    // eigenvalue using the implicitely restarted Arnoldi method.  Here, A is in
    // rectangular packed format (RPF).
    template <typename Real>
    std::pair <Real,Real> syiram(
        Natural const & m,
        Real const * const Ap,
        Natural const & iter_innr_max,
        Natural const & iter_outr_max,
        Real const & tol
    ) {
        // Allocate memory for solving a dense eigenvalue problem in LAPACK.
        Integer eig_size = m<=iter_innr_max ? m : iter_innr_max;
        std::vector <Real> Aeig(eig_size*(eig_size+1)/2);
        Integer nevals(0);
        std::vector <Real> W(eig_size);
        std::vector <Real> Z(0);
        std::vector <Real> work(8*eig_size);
        std::vector <Integer> ifail(0);
        std::vector <Integer> iwork(5*eig_size);
        Integer info=0;

        // If the size of the matrix is equal to or smaller than the maximum 
        // number of inner iterations, just use LAPACK to find the smallest
        // eigenvalue
        if(m<=iter_innr_max) {
            // Copy in the operator for the eigenvalue problem 
            copy <Real> (m*(m+1)/2,&(Ap[0]),1,&(Aeig[0]),1);

            // Solve the eigenvalue problem for the smallest eigenvalue
            spevx <Real> ('N','I','U',m,&(Aeig[0]),Real(0.),
                Real(0.),1,1,Real(2.)*lamch <Real> ('S'),nevals,
                &(W[0]),&(Z[0]),m,&(work[0]),&(iwork[0]),&(ifail[0]),info);
            return std::pair <Real,Real> (W[0],Real(0.));
        }

        // Initialize memory for all of the Krylov vectors
        std::vector <Real> V(m*(iter_innr_max+1));

        // Allocate memory for the norm of the current Krylov vector.  This is
        // also the element below the diagonal in H (if we stored it
        // explicitely).
        Real norm_v;
        
        // Initialize the first Krylov vector 
        std::mt19937 gen(1);
        std::uniform_real_distribution<> dis(0, 1);
        for(Natural i=1;i<=m;i++)
            V[ijtok(i,1,m)]=dis(gen);
        norm_v=sqrt(dot <Real> (m,&(V[ijtok(1,1,m)]),1,&(V[ijtok(1,1,m)]),1));
        scal <Real> (m,Real(1.)/norm_v,&(V[ijtok(1,1,m)]),1);

        // Allocate memory for the upper Hessenberg matrix that arises from
        // Gram-Schmidt.  Given that we're working with a symmetric matrix, this
        // is really tridiagonal.  At this point, we're just going to leave the
        // matrix since we have the operator Q' H Q later and I don't currently
        // know of a great way to compute this if we have two vectors
        // representing the tridiagonal form of H.  Also, we have Hp to
        // represent the packed version of this matrix and H to represent the
        // symmetric version of this.  We prefer the use of Hp, but convert to
        // H when required.
        std::vector <Real> Hp(iter_innr_max*(iter_innr_max+1)/2,Real(0.));
        std::vector <Real> H(iter_innr_max*iter_innr_max,Real(0.));

        // Allocate memory for the QR factorizations used in the restart
        std::vector <Real> Q(iter_innr_max*iter_innr_max);
        std::vector <Real> Q_all(iter_innr_max*iter_innr_max);
        std::vector <Real> tau(iter_innr_max);
        Integer lwork_qr=iter_innr_max;
        std::vector <Real> work_qr(lwork_qr);

        // Allocate memory for use in the QR-shifts 
        std::vector <Real> QR_tmp(iter_innr_max*iter_innr_max);

        // Allocate memory for the vector of all ones.  This makes it easier
        // to define the identiy matrix later
        std::vector <Real> e(iter_innr_max,Real(1.));

        // Allocate memory for the starting Krylov vectors after the restart
        std::vector <Real> v1(m);
        std::vector <Real> v2(m);

        // Continue to compute until we converge
        for(Natural iter_outr=1;iter_outr<=iter_outr_max;iter_outr++) {
            // If we restart, we already have a second Krylov vector.
            // Otherwise, we compute the second Krylov vector.
            Natural gs_start= iter_outr==1 ? 1 : 2; 

            // Do the Arnoldi iteration
            for(Natural k=gs_start;k<=iter_innr_max;k++) {
                // Gram-Schmidt with DGKS correction
                spmv <Real> ('U',m,Real(1.),&(Ap[0]),&(V[ijtok(1,k,m)]),1,
                    Real(0.),&(V[ijtok(1,k+1,m)]),1);
                for(Natural dgks_iter=1;dgks_iter<=2;dgks_iter++) {
                    for(Natural i=1;i<=k;i++) {
                        Real alpha = dot <Real> (
                            m,&(V[ijtok(1,k+1,m)]),1,&(V[ijtok(1,i,m)]),1);
                        Hp[ijtokp(i,k)]+=alpha;
                        axpy <Real> (m,-alpha,&(V[ijtok(1,i,m)]),1,
                            &(V[ijtok(1,k+1,m)]),1);
                    }
                }
                norm_v=sqrt(
                    dot<Real>(m,&(V[ijtok(1,k+1,m)]),1,&(V[ijtok(1,k+1,m)]),1));
                scal <Real> (m,Real(1.)/norm_v,&(V[ijtok(1,k+1,m)]),1);
            }

            // Find the Ritz values of H 
            copy <Real> (iter_innr_max*(iter_innr_max+1)/2,&(Hp[0]),1,
                &(Aeig[0]),1);
            spevx <Real> ('N','A','U',iter_innr_max,&(Aeig[0]),Real(0.),
                Real(0.),0,0,Real(2.)*lamch <Real> ('S'),nevals,
                &(W[0]),&(Z[0]),iter_innr_max,&(work[0]),&(iwork[0]),
                &(ifail[0]),info);

            // Q_all <- I
            scal <Real> (iter_innr_max*iter_innr_max,Real(0.),&(Q_all[0]),1);
            axpy <Real> (iter_innr_max,Real(1.),&(e[0]),1,&(Q_all[0]),
                iter_innr_max+1);

            // Next, do a number of QR-iterations using the Ritz values to the
            // right of the leftmost Ritz value as shifts
            for(Natural i=iter_innr_max;i>=2;i--) {

                // Convert Hp to H.  We need H for the QR factorization routine,
                // geqrf, below.
                tpttr <Real> ('U',iter_innr_max,&(Hp[0]),&(H[0]),iter_innr_max,
                    info);
                for(Natural j=1;j<=iter_innr_max;j++)
                    copy <Real> (iter_innr_max-j,
                        &(H[ijtok(j,j+1,iter_innr_max)]),iter_innr_max,
                        &(H[ijtok(j+1,j,iter_innr_max)]),1);

                // Q <- H - theta_i I
                copy <Real> (iter_innr_max*iter_innr_max,&(H[0]),1,&(Q[0]),1);
                axpy <Real> (iter_innr_max,-W[ijtok(i,1,iter_innr_max)],
                    &(e[0]),1,&(Q[0]),iter_innr_max+1);

                // Q <- QR(H - theta_i I)
                geqrf <Real> (iter_innr_max,iter_innr_max,&(Q[0]),iter_innr_max,
                    &(tau[0]),&(work_qr[0]),lwork_qr,info);
                orgqr <Real> (iter_innr_max,iter_innr_max,iter_innr_max,
                    &(Q[0]),iter_innr_max,&(tau[0]),&(work_qr[0]),lwork_qr,
                    info);

                // QR_tmp <- H Q
                symm <Real> ('L','U',iter_innr_max,iter_innr_max,Real(1.),
                    &(H[0]),iter_innr_max,&(Q[0]),iter_innr_max,Real(0.),
                    &(QR_tmp[0]),iter_innr_max);

                // H <- Q' H Q
                gemm <Real> ('T','N',iter_innr_max,iter_innr_max,iter_innr_max,
                    Real(1.),&(Q[0]),iter_innr_max,&(QR_tmp[0]),iter_innr_max,
                    Real(0.),&(H[0]),iter_innr_max);

                // Convert back to Hp
                trttp <Real> ('U',iter_innr_max,&(H[0]),iter_innr_max,&(Hp[0]),
                    info); 

                // QR_tmp <- Q Q_all
                gemm <Real> ('N','N',iter_innr_max,iter_innr_max,iter_innr_max,
                    Real(1.),&(Q_all[0]),iter_innr_max,&(Q[0]),iter_innr_max,
                    Real(0.),&(QR_tmp[0]),iter_innr_max);

                // Q_all <- Q Q_all
                copy <Real> (iter_innr_max*iter_innr_max,&(QR_tmp[0]),1,
                    &(Q_all[0]),1);
            }
            
            // Find the starting Krylov vectors

            // v1 <- V Q_all(*,1)
            gemv <Real> ('N',m,iter_innr_max,Real(1.),&(V[0]),m,
                &(Q_all[ijtok(1,1,iter_innr_max)]),1,Real(0.),&(v1[0]),1);

            // v2 <- V Q_all(*,2)
            gemv <Real> ('N',m,iter_innr_max,Real(1.),&(V[0]),m,
                &(Q_all[ijtok(1,2,iter_innr_max)]),1,Real(0.),&(v2[0]),1);

            // v2 <- V Q_all(*,2) H(2,1)
            scal <Real> (m,Hp[ijtokp(2,1)],&(v2[0]),1);

            // v2 <- V Q_all(*,2) H(2,1) 
            //           + norm_v Q_all(iter_innr_max,1) V(*,iter_innr_max+1)
            axpy <Real> (m,norm_v*Q_all[ijtok(iter_innr_max,1,iter_innr_max)],
                &(V[ijtok(1,iter_innr_max+1,m)]),1,&(v2[0]),1);

            // Setup the upper Hessenberg matrix that stores the Gram-Schmidt
            // coefficients and the norms

            // H(1,1) stays the same
            Real H11 = Hp[ijtokp(1,1)];

            // norm_v (H(2,1)) <- || v2 ||, but we're not storing the lower
            // diagonal pieces
            norm_v=sqrt (dot <Real> (m,&(v2[0]),1,&(v2[0]),1));

            // Zero out the upper Hessenberg matrix since we just add elements
            // into this matrix since we're using the DGKS correction
            scal <Real> (iter_innr_max*(iter_innr_max+1)/2,Real(0.),&(Hp[0]),1);

            // Insert H11
            Hp[ijtokp(1,1)]=H11;

            // Insert the new Krylov vectors into V
            scal <Real> (m,Real(1.)/norm_v,&(v2[0]),1);
            copy <Real> (m,&(v1[0]),1,&(V[ijtok(1,1,m)]),1);
            copy <Real> (m,&(v2[0]),1,&(V[ijtok(1,2,m)]),1);

            // Check the stopping condition.  Basically, we can look at the Ritz
            // values of the new H and multiply it by the norm of the last
            // Krylov vector.  However, since this is 1x1, we can just look at
            // H(2,1), which is || v2 ||.  If we're converged, the current Ritz
            // value is just H(1,1).
            if(norm_v < tol)
                break;
        }

        // Return the smallest current Ritz value and the estimated error
        return std::pair <Real,Real> (Hp[ijtokp(1,1)],norm_v);
    }

    // Solve the generalized, symmetric eigenvalue problem A x = lambda x for
    // the leftmost eigenvalue using the implicitely restarted Arnoldi method.
    // Here, A and B is in rectangular packed format (RPF) and we assume that
    // B is positive definite.
    template <typename Real>
    std::pair <Real,Real> gsyiram(
        Natural const & m,
        Real const * const Arf,
        Real const * const Brf,
        Natural const & iter_innr_max,
        Natural const & iter_outr_max,
        Real const & tol
    ) {
        // First, make a copy of Brf since the Choleski factorization is
        // destructive
        std::vector <Real> Brf0 (m*(m+1)/2);
        copy <Real> (m*(m+1)/2,Brf,1,&(Brf0[0]),1);

        // Then, find the Choleski factorization of Brf0
        Integer info;
        pftrf('N','U',m,&(Brf0[0]),info);

        // Next, find the packed version of Arf and Brf0
        std::vector <Real> Ap(m*(m+1)/2);
        tfttp <Real> ('N','U',m,Arf,&(Ap[0]),info);
        std::vector <Real> Bp(m*(m+1)/2);
        tfttp <Real> ('N','U',m,&(Brf0[0]),&(Bp[0]),info);

        // Ap <- inv(U') A inv(U) where U comes from the above Choleski
        // fractorization of Brf.
        spgst(1,'U',m,&(Ap[0]),&(Bp[0]),info);

        // Now, find the smallest eigenvalue of Ap = inv(U') A inv(U) 
        return syiram <Real> (m,&(Ap[0]),iter_innr_max,iter_outr_max,tol);
    }

    // Solves a quadratic equation
    //
    // a x^2 + b x + c = 0
    //
    // Here, we assume that a, b, and c are not all zero.
    //
    // (input) a, b, c : Coefficients of the quadratic
    // (output) nroots : Number of roots
    // (output) r1 : First root, if it exists
    // (output) r2 : Second root, if it exists
    template <typename Real>
    void quad_equation(
        Real const & a,
        Real const & b,
        Real const & c,
        Natural & nroots,
        Real & r1,
        Real & r2
    ) {

        // It's sort of hard to tell if we have a quadratic or a linear since
        // we don't have a good way with the information provided to tell if
        // the quadratic coefficient is small.  As such, we do a bad, hard
        // check if the leading coefficient is zero.  If it is not the case,
        // we assume that we have a quadratic and we use the most stable
        // equation that we can for the root.
        if( a != Real(0.) ) { 
            if(b < Real(0.)) {
                r1 = (-b + sqrt(b*b-Real(4.)*a*c)) / (Real(2.)*a);
                r2 = (Real(2.)*c) / (-b + sqrt(b*b-Real(4.)*a*c));
            } else {
                r1 = (Real(2.)*c) / (-b - sqrt(b*b-Real(4.)*a*c));
                r2 = (-b - sqrt(b*b-Real(4.)*a*c)) / (Real(2.)*a);
            }
            nroots = 2;

        // Now, in the case that a is zero, but b is not, we have a linear
        // function and we can solve for the root.
        } else if( b != Real(0.)) {
            r1 = -c/b;
            nroots = 1;

        // Here, we have a constant function.  Now, we could have no roots
        // if c is zero.  Alternatively, we could have an infinity number of
        // roots if c is zero.  There's not really a good way to denote all
        // of these cases, so we just assume that c is not zero and return
        // zero roots.
        } else {
            nroots = 0;
        }
    }

    // Reasons we stop the Krylov method
    namespace KrylovStop{
        enum t{
            NegativeCurvature,        // Negative curvature detected
            RelativeErrorSmall,       // Relative error is small
            MaxItersExceeded,         // Maximum number of iterations exceeded
            TrustRegionViolated,      // Trust-region radius violated
            Instability,              // Some sort of instability detected.
                                      // Normally, a NaN.
            InvalidTrustRegionCenter  // Trust-region center is chosen such that
                                      // || C(x_cntr) || > delta where delta is
                                      // the trust-region radius.  This means
                                      // that our starting solution of 0
                                      // violates the trust-region.
        };

        // Converts the Krylov stopping condition to a string 
        std::string to_string(t const & krylov_stop);
        
        // Converts a string to a Krylov stopping condition
        t from_string(std::string const & krylov_stop);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // A orthogonalizes a vector Bx to a list of other Bxs.  
    template <
        typename Real,
        template <typename> class XX
    >
    void Aorthogonalize(
        std::list <typename XX <Real>::Vector> const & Bvs,
        std::list <typename XX <Real>::Vector> const & ABvs,
        typename XX <Real>::Vector & Bx,
        typename XX <Real>::Vector & ABx
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Orthogonalize the vectors
        for(typename std::list <X_Vector>::const_iterator
                Bv=Bvs.begin(),
                ABv=ABvs.begin();
            Bv!=Bvs.end();
            Bv++,ABv++
        ) {
            Real beta=X::innr(*ABv,Bx);
            X::axpy(Real(-1.)*beta,*Bv,Bx);
            X::axpy(Real(-1.)*beta,*ABv,ABx);
        }
    }

    // Computes the truncated projected conjugate direction algorithm in order
    // to solve Ax=b where we restrict x to be in the range of B and that
    // || C (x - x_cntr) || <= delta.  The parameters are as follows.
    // 
    // (input) A : Operator in the system A B x = b.
    // (input) b : Right hand side in the system A B x = b.
    // (input) B : Projection in the system A B x = b.
    // (input) C : Operator that modifies the shape of the trust-region.
    // (input) eps : Stopping tolerance.
    // (input) iter_max :  Maximum number of iterations.
    // (input) orthog_max : Maximum number of orthgonalizations.  If this
    //     number is 1, then we do the conjugate gradient algorithm.
    // (input) delta : Trust region radius.  If this number is infinity, we
    //     do not scale the final step if we detect negative curvature.
    // (input) x_cntr : Center of the trust-region. 
    // (input) do_orthog_check : Orthogonality check for projected algorithms 
    // (output) x : Final solution x.
    // (output) x_cp : The Cauchy-Point, which is defined as the solution x
    //     after a single iteration.
    // (output) norm_Br : The norm ||B r|| of the final residual.
    // (output) iter : The number of iterations required to converge. 
    // (output) krylov_stop : The reason why the Krylov method was terminated.
    template <
        typename Real,
        template <typename> class XX
    >
    void truncated_cd(
        Operator <Real,XX,XX> const & A,
        typename XX <Real>::Vector const & b,
        Operator <Real,XX,XX> const & B,
        Operator <Real,XX,XX> const & C,
        Real const & eps,
        Natural const & iter_max,
        Natural const & orthog_max,
        Real const & delta,
        typename XX <Real>::Vector const & x_cntr,
	bool const & do_orthog_check,
        typename XX <Real>::Vector & x,
        typename XX <Real>::Vector & x_cp,
        Real & norm_Br0,
        Real & norm_Br,
        Natural & iter,
        KrylovStop::t & krylov_stop
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Set the tolerance for our orthogonality check
        const Real eps_orthog(0.5);

        // Allocate memory for the projected search direction and the
        // the operator applied to the projection
        X_Vector Bp(X::init(x));
        X_Vector ABp(X::init(x));

        // Allocate memory for the previous search directions
        std::list <X_Vector> Bps;
        std::list <X_Vector> ABps;
        
        // Allocate memory for the residuals and projected residuals 
        std::list <X_Vector> rs;
        std::list <X_Vector> Brs;
        std::list <Real> norm_Brs;

        // Allocate memory for the orthogonality check matrix.  This is
        // inv(D) ( Brs' rs - D^2 ) inv(D) = inv(D) Brs' rs inv(D) - I,
        // where D = diag(||Br_i||).  In theory, the projected/preconditioned
        // residuals should be orthogonal to the residuals.  In practice,
        // this may not be the case,  Hence, we look at how antisymmetric
        // this matrix is and bail if it becomes too bad.  Note, we really
        // need to know ( i, j, <orthog_check>_ij ).  We need the indices in
        // order to help eliminate elements when we don't do full 
        // orthogonalization.
        typedef typename
            std::list < std::pair< std::pair <Natural,Natural> , Real > >
            indexed_matrix;
        indexed_matrix orthog_check;

        // Allocate memory for the norm of the trust-region operator
        // applied at the trial step, || C( (x-x_cntr) + alpha Bp) ||.  
        Real norm_C_x_m_xcntr_p_a_Bp(std::numeric_limits <Real>::quiet_NaN());

        // Allocate memory for the line-search used in the algorithm.
        Real alpha(std::numeric_limits<Real>::quiet_NaN());

        // Initialize x to zero. 
        X::zero(x);

        // Allocate memory for and find the initial residual, A*x-b = -b.
        // Also, allocate memory for the projected residual.
        X_Vector r(X::init(x));
        X_Vector Br(X::init(x));
        X::copy(b,r);
        X::scal(Real(-1.),r);

        // Allocate memory for two additional work vectors.  We use these
        // when calcuating the norm of the operator C applied to the
        // trial step.  In theory, we could get rid of at least one of these
        // work vectors.  However, by doing so, we have to calculate this
        // norm implicitely, which involves additional inner products and adding
        // the results together.  This can have some numerical difficulties,
        // so we just sacrifice the additional memory.
        X_Vector x_tmp1(X::init(x));
        X_Vector x_tmp2(X::init(x));

        // Allocate memory for the quantity x-x_cntr
        X_Vector x_m_xcntr(X::init(x));
        X::copy(x_cntr,x_m_xcntr);
        X::scal(Real(-1.),x_m_xcntr);

        // Verify that ||C(x-x_cntr)|| <= delta.  This insures that our initial
        // iterate lies inside the trust-region.  If it does not, we exit.
        C.eval(x_m_xcntr,x_tmp1);
        Real norm_C_x_m_xcntr = sqrt(X::innr(x_tmp1,x_tmp1));
        if(norm_C_x_m_xcntr > delta) {
            X::zero(x_cp);
            iter=0;
            krylov_stop = KrylovStop::InvalidTrustRegionCenter;
            return;
        }

        // Find the Br application.  We use this to find the B-norm of the
        // residual.  Then, we save the original B-norm of the residual.
        B.eval(r,Br);
        norm_Br0 = sqrt(X::innr(Br,Br));
        norm_Br = norm_Br0; 

        // Find the projected search direction 
        X::copy(Br,Bp);
        X::scal(Real(-1.),Bp);

        // Loop until the maximum iteration
        for(iter=1;iter<=iter_max;iter++){
        
            // If the norm of the residual is small relative to the starting
            // residual, exit
            if(norm_Br <= eps*norm_Br0) {
                iter--;
                krylov_stop = KrylovStop::RelativeErrorSmall;
                break;
            }

            // Find the ABp application
            A.eval(Bp,ABp);

            // Orthogonalize this direction to the previous directions
            Aorthogonalize <Real,XX> (Bps,ABps,Bp,ABp); 

            // Check if this direction is a descent direction.  If it is not,
            // flip it so that it is.
            if(X::innr(Bp,r) > Real(0.)) {
                X::scal(Real(-1.),Bp);
                X::scal(Real(-1.),ABp);
            }

            // Find || Bp ||_A^2.  Note, this is how we detect negative
            // curvature as well as instability in the algorithm due to things
            // like NaNs in the operator calculations.
            Real Anorm_Bp_2 = X::innr(Bp,ABp);

            // If we detect negative curvature or instability, in this case
            // that || Bp ||_A=NaN, then we don't need to compute the
            // following steps.
            if( Anorm_Bp_2 > Real(0.) && Anorm_Bp_2 == Anorm_Bp_2 ) {
                // Find || Bp ||_A.
                Real Anorm_Bp = sqrt(Anorm_Bp_2);

                // Check if we need to eliminate any vectors for
                // orthogonalization.
                if(Bps.size()==orthog_max) {
                    Bps.pop_front();
                    ABps.pop_front();
                    rs.pop_front();
                    Brs.pop_front();
                    norm_Brs.pop_front();

                    // Don't remove elements from the orthogonality check 
                    // matrix until we have enough elements to remove.
                    if(iter>orthog_max) {
                        // Figure out the index to remove from orthog_check
                        Natural ii=iter-orthog_max;

                        // Remove all elements associated with this index
                        typename indexed_matrix::iterator
                            ele=orthog_check.begin();
                        while(ele!=orthog_check.end()) {
                            if(ele->first.first==ii || ele->first.second==ii)
                                ele=orthog_check.erase(ele);
                            else
                                ele++;
                        }
                    }
                }

                // Store the previous directions
                Bps.emplace_back(std::move(X::init(x)));
                X::copy(Bp,Bps.back());
                X::scal(Real(1.)/Anorm_Bp,Bps.back());
                
                ABps.emplace_back(std::move(X::init(x)));
                X::copy(ABp,ABps.back());
                X::scal(Real(1.)/Anorm_Bp,ABps.back());

                // Store the previous residuals
                rs.emplace_back(std::move(X::init(x)));
                X::copy(r,rs.back());

                Brs.emplace_back(std::move(X::init(x)));
                X::copy(Br,Brs.back());

                norm_Brs.emplace_back(norm_Br);

                // Build new pieces of the orthogonality check matrix
                typename std::list <X_Vector>::reverse_iterator
                    Br_star=Brs.rbegin();
                typename std::list <Real>::reverse_iterator
                    norm_Bri=norm_Brs.rbegin(); 

                // Start out at the iteration number and count backwards
                // the number of elements equal to the number of residual
                // vectors.
                for(Natural ii=iter;ii>iter-rs.size();ii--) {
                    typename std::list <X_Vector>::reverse_iterator
                        r_star=rs.rbegin();
                    typename std::list <Real>::reverse_iterator
                        norm_Brj=norm_Brs.rbegin(); 

                    // Use the same counting scheme as ii
                    for(Natural jj=iter;jj>iter-rs.size();jj--) {

                        // Don't recompute elements that we've already stored
                        if(ii >= iter || jj >= iter) {
                            // < Br_i, r_j > / || Br_i || || Br_j ||
                            orthog_check.emplace_back(
                                std::pair <std::pair <Natural,Natural> , Real> (
                                    std::pair <Natural,Natural> (ii,jj),
                                    X::innr(*Br_star,*r_star)
                                        / ((*norm_Brj)*(*norm_Bri))
                                )
                            );

                            // If we're on a diagonal element, make sure to
                            // remove the piece of the identity, 1.
                            if(ii==jj)
                                orthog_check.back().second -= Real(1.);
                        }

                        r_star++;
                        norm_Brj++;
                    }
                    Br_star++;
                    norm_Bri++;
                }

                // Do an exact linesearch in the computed direction
                alpha = -X::innr(r,Bp) / Anorm_Bp_2;

                // Find the trial step and the norm || C(x+alpha Bp) || 

                // x_tmp1 <- (x-x_cntr) + alpha Bp
                X::copy(x_m_xcntr,x_tmp1);
                X::axpy(alpha,Bp,x_tmp1);

                // x_tmp2 <- C( (x-x_cntr) + alpha Bp)
                C.eval(x_tmp1,x_tmp2);

                // norm_C_x_m_xcntr_p_a_Bp <- || C( (x-x_cntr) + alpha Bp) ||
                norm_C_x_m_xcntr_p_a_Bp = sqrt(X::innr(x_tmp2,x_tmp2));
            }

            // Calculate the Frobenius norm of the orgonality check.
            Real norm_orthogcheck(0.);
            for(typename indexed_matrix::iterator ele=orthog_check.begin();
                ele!=orthog_check.end();
                ele++
            ) 
                norm_orthogcheck += (ele->second)*(ele->second);
            norm_orthogcheck = sqrt(norm_orthogcheck);

            // Check if we have instability.  This is either a NaN in our
            // normalization factor for Gram-Schmidt or if the norm
            // of our orthogonality check matrix is greater than our
            // tolerance.
            bool instability = (Anorm_Bp_2!=Anorm_Bp_2)
                || (do_orthog_check && norm_orthogcheck > eps_orthog);

            // If we have negative curvature or our trial point is outside the
            // trust-region radius, terminate truncated-PCD and find our final
            // step by scaling the search direction until we reach the
            // trust-region radius.  In the case of a NaN, we can't trust
            // our current search direction, Bp, so we take our last search
            // direction, which is stored in the Bps.
            if( Anorm_Bp_2 <= Real(0.) ||
                norm_C_x_m_xcntr_p_a_Bp >= delta ||
                instability 
            ) {

                // If we have instability, we terminate without taking a step
                if(instability)
                    ;

                // If we're paying attention to the trust-region, scale the
                // step appropriately.
                else if(delta < std::numeric_limits <Real>::infinity()) {
                    // Find sigma so that
                    // || C((x-x_cntr) + sigma Bp) || = delta.
                    // This can be found by finding the positive root of the
                    // quadratic
                    // || C(x-x_cntr) + sigma C(Bp) ||^2 = delta^2.
                    // Specifically, we want the positive root of
                    // sigma^2 <C(Bp),C(Bp)>
                    //     + sigma (2 <C(Bp),C(x-x_cntr)>)
                    //     + (<C(x-x_cntr),C(x-x_cntr)>-delta^2).

                    // x_tmp1 <- C(Bp)
                    C.eval(Bp,x_tmp1);

                    // x_tmp2 <- C(x-x_cntr)
                    C.eval(x_m_xcntr,x_tmp2);

                    // Solve the quadratic equation for the positive root 
                    Real aa = X::innr(x_tmp1,x_tmp1);
                    Real bb = Real(2.)*X::innr(x_tmp1,x_tmp2);
                    Real cc = norm_C_x_m_xcntr*norm_C_x_m_xcntr-delta*delta;
                    Natural nroots;
                    Real r1;
                    Real r2;
                    quad_equation(aa,bb,cc,nroots,r1,r2);
                    Real sigma = r1 > r2 ? r1 : r2;

                    // Take the step and find its residual. Then, compute
                    // B-norm of the residual. 
                    X::axpy(sigma,Bp,x);
                    X::axpy(sigma,Bp,x_m_xcntr);
                    X::axpy(sigma,ABp,r);
                    B.eval(r,Br);
                    norm_Br=sqrt(X::innr(Br,Br));

                // In the case that we're ignoring the trust-region, we
                // take a step of unit scale. 
                } else {
                    X::axpy(Real(1.),Bp,x);
                    X::axpy(Real(1.),Bp,x_m_xcntr);
                    X::axpy(Real(1.),ABp,r);
                    B.eval(r,Br);
                    norm_Br=sqrt(X::innr(Br,Br));
                }

                // Determine why we stopped
                if(instability)
                    krylov_stop = KrylovStop::Instability;
                else if(Anorm_Bp_2 <= Real(0.)) 
                    krylov_stop = KrylovStop::NegativeCurvature;
                else
                    krylov_stop = KrylovStop::TrustRegionViolated;
 
                // If this is the first iteration, save the Cauchy-Point
                if(iter==1) X::copy(x,x_cp);
                break;
            }

            // Take a step in this direction
            X::axpy(alpha,Bp,x);
            X::axpy(alpha,Bp,x_m_xcntr);

            // If this is the first iteration, save the Cauchy-Point
            if(iter==1) X::copy(x,x_cp);

            // Update the norm of x
            norm_C_x_m_xcntr = norm_C_x_m_xcntr_p_a_Bp;

            // Find the new residual and projected residual
            X::axpy(alpha,ABp,r);
            B.eval(r,Br);
            norm_Br = sqrt(X::innr(Br,Br));

            // Find the projected steepest descent direction
            X::copy(Br,Bp);
            X::scal(Real(-1.),Bp);	
        }

        // If we've exceeded the maximum iteration, make sure to denote this
        if(iter > iter_max)
            krylov_stop=KrylovStop::MaxItersExceeded;

        // Adjust the iteration number if we ran out of iterations
        iter = iter > iter_max ? iter_max : iter;
    }

    // Solve a 2x2 linear system in packed storage.  This is done through
    // Gaussian elimination with complete pivoting.  In addition, this assumes
    // that the system is nonsingular.
    //
    // (input) A : 2x2 matrix in packed storage.  That's a length 3 vector.
    //    Note, we pass this in by value, which initiates a copy.  This is
    //    because we modify the matrix, so we were going to need a copy
    //    anyway.
    // (input) b : Vector of length 2.
    // (output) x : Solution to the linear system.
    template <typename Real>
    void solve2x2(
        std::vector <Real> A,
        std::vector <Real> b,
        std::vector <Real> & x
    ) {
        // Find the largest element of A in absolute value.  Store in i.
        Natural i=0;  Real val=fabs(A[0]);
        for(Natural j=1;j<=2;j++) {
            if(fabs(A[j]) < val) {
                i=j;
                val=fabs(A[j]);
            }
        }

        // Determine the row and column pivots
        std::vector <Natural> p(3);
        std::vector <Natural> q(3);
        if(i==0) {
            p[1]=1; p[2]=2;
            q[1]=1; q[2]=2;
        } else if(i==1) {
            p[1]=2; p[2]=1;
            q[1]=1; q[2]=2;
        } else {
            p[1]=2; p[2]=1;
            q[1]=2; q[2]=1;
        }

        // Do a step of Gaussian elimination
        Real alpha = -A[ijtokp(p[2],q[1])] / A[ijtokp(p[1],q[1])];
        A[ijtokp(p[2],q[2])] = A[ijtokp(p[2],q[2])]
            + alpha*A[ijtokp(p[1],q[2])];
        b[itok(p[2])] = b[itok(p[2])] + alpha * b[itok(p[1])];

        // Do back subsitutition
        x.resize(2);
        x[itok(p[2])] = b[itok(p[2])] / A[ijtokp(p[2],q[2])];
        x[itok(p[1])] =
            (b[itok(p[1])] - A[ijtokp(p[1],q[2])] * x[itok(p[2])])
            / A[ijtokp(p[1],q[1])];
    }

    // Evaluate a two variable objective function of the form
    //
    // f(x) = x'*A*x + a'*x
    //
    // where A is held in packed storage.
    //
    // (input) A : 2x2 matrix in packed storage.  That's a length 3 vector.
    // (input) a : Vector of length 2.
    // (input) x : Vector of length 2.
    // (return) objective value
    template <typename Real>
    Real obj2x2(
        std::vector <Real> const & A,
        std::vector <Real> const & a,
        std::vector <Real> const & x
    ) {
        return (A[0]*x[0]+a[0])*x[0] + (A[2]*x[1]+a[1])*x[1]
            + Real(2.)*A[1]*x[0]*x[1];
    }

    // Optimize a two variable, box constrained, quadratic program of the form
    //
    // min <Ax,x> + <a,x> st lb <= x <= ub
    //
    // This is accomplished through brute force.  Namely, an active set method
    // where we just check all possible combinations of active sets.
    //
    // (input) A : 2x2 matrix in packed storage.  That's a length 3 vector.
    // (input) a : Vector of length 2.
    // (input) x : Vector of length 2.
    // (input) lb : Vector of length 2.
    // (input) ub : Vector of length 2.
    // (output) x : The optimal solution.
    template <typename Real>
    void quad2x2(
        std::vector <Real> const & A,
        std::vector <Real> const & a,
        std::vector <Real> const & lb,
        std::vector <Real> const & ub,
        std::vector <Real> & x
    ) {

        // Store the best objective function up until this point
        Real f_x = std::numeric_limits <Real>::infinity();
        
        // List all of the points to check for feasibility and optimality
        std::list <std::vector <Real> > zs;
        std::vector <Real> z(2);

        // Unconstrained minimum
        std::vector <Real> minus_a(2); minus_a[0]=-a[0]; minus_a[1]=-a[1];
        solve2x2 <Natural,Real> (A,minus_a,z);
        zs.emplace_back(z);  

        // z1 to the lower bound
        z[0] = lb[0];
        z[1] = -(a[1]+Real(2.)*A[0]*A[1]*z[0])/(Real(2.)*A[2]);
        zs.emplace_back(z);
       
        // z2 to the lower bound
        z[1] = lb[1];
        z[0] = -(a[0]+Real(2.)*A[0]*A[1]*z[1])/(2*A[0]); 
        zs.emplace_back(z);
        
        // z1 to the upper bound 
        z[0] = ub[0];
        z[1] = -(a[1]+Real(2.)*A[0]*A[1]*z[0])/(2*A[2]);
        zs.emplace_back(z);
        
        // z2 to the upper bound
        z[1] = ub[1];
        z[0] = -(a[0]+Real(2.)*A[0]*A[1]*z[1])/(2*A[0]); 
        zs.emplace_back(z);
       
        // Lower left corner
        z[0] = lb[0];
        z[1] = lb[1]; 
        zs.emplace_back(z);
        
        // Lower right corner
        z[0] = ub[0];
        z[1] = lb[1]; 
        zs.emplace_back(z);

        // Upper right corner
        z[0] = ub[0];
        z[1] = ub[1]; 
        zs.emplace_back(z);

        // Upper left corner
        z[0] = lb[0];
        z[1] = ub[1]; 
        zs.emplace_back(z);

        // Find the feasible point with the lowest objective value
        for(typename std::list <std::vector <Real> >::iterator zp=zs.begin();
            zp!=zs.end();
            zp++
        ) {
            Real f_z = obj2x2 <Natural,Real> (A,a,*zp);
            if((*zp)[0]>=lb[0] && (*zp)[1]>=lb[1] &&
               (*zp)[0]<=ub[0] && (*zp)[1]<=ub[1] &&
               f_z<f_x
           ){
                x = *zp;
                f_x = f_z;
            }
        }
    }

    // Orthogonalizes a vector x to a list of other xs.  
    template <
        typename Real,
        template <typename> class XX
    >
    void orthogonalize(
        std::list <typename XX <Real>::Vector> const & vs,
        typename XX <Real>::Vector & x,
        Real * R
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Orthogonalize the vectors
        Natural i=0;
        for(typename std::list <X_Vector>::const_iterator v=vs.begin();
            v!=vs.end();
            v++
        ) {
            Real beta=X::innr(*v,x);
            X::axpy(Real(-1.)*beta,*v,x);
            R[i] = beta;
            i++;
        }
    }

    // Solves for the linear solve iterate update dx in the current Krylov space
    template <
        typename Real,
        template <typename> class XX
    >
    void solveInKrylov(
        Natural const & m,
        Real const * const R,
        Real const * const Qt_e1,
        std::list <typename XX <Real>::Vector> const & vs,
        Operator <Real,XX,XX> const & Mr_inv,
        typename XX <Real>::Vector const & x,
        typename XX <Real>::Vector & dx
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        
        // Allocate memory for the solution of the triangular solve 
        std::vector <Real> y(m);

        // Create one temporary element required to solve for the iterate
        X_Vector V_y(X::init(x));

        // Solve the system for y
        copy <Real> (m,&(Qt_e1[0]),1,&(y[0]),1);
        tpsv <Real> ('U','N','N',m,&(R[0]),&(y[0]),1);

#if 0
        // Check the condition number on R.  If it's large, something
        // wrong has occured and we need to bail.
        Integer info(0);
        Real rcond(0.);
        std::vector <Real> work(3*m);
        std::vector <Integer> iwork(m);
        tpcon('I','U','N',m,&(R[0]),rcond,&(work[0]),&(iwork[0]),info);
        if(rcond < std::numeric_limits <Real>::epsilon()*1e3)
            return false;
#endif

        // Compute tmp = V y
        X::zero(V_y);
        typename std::list <X_Vector>::const_iterator vv=vs.begin();
        for(Natural j=0;j<m;j++) {
            X::axpy(Real(y[j]),*vv,V_y);
            vv++;
        }

        // Right recondition the above linear combination
        Mr_inv.eval(V_y,dx);
    }

    // Resets the GMRES method.  This does a number of things
    // 1.  Calculates the preconditioned residual.
    // 2.  Finds the norm of the preconditioned residual.
    // 3.  Finds the initial Krylov vector.
    // 4.  Initializes the list of Krylov vectors.
    // 5.  Finds the initial RHS for the least squares system, Q' norm(w1) e1.
    // 6.  Clears out all of the old Givens rotations
    // These steps are required during initialization as well as during a
    // restart of GMRES
    template <
        typename Real,
        template <typename> class XX
    >
    void resetGMRES(
        typename XX <Real>::Vector const & rtrue,
        Operator <Real,XX,XX> const & Ml_inv,
        Natural const & rst_freq,
        typename XX <Real>::Vector & v,
        std::list <typename XX <Real>::Vector> & vs,
        typename XX <Real>::Vector & r,
        Real & norm_r,
        std::vector <Real> & Qt_e1,
        std::list <std::pair<Real,Real> > & Qts
    ){
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Apply the left preconditioner to the true residual.  This
        // completes #1
        Ml_inv.eval(rtrue,r);

        // Store the norm of the preconditioned residual.  This completes #2.
        norm_r = sqrt(X::innr(r,r));

        // Find the initial Krylov vector.  This completes #3.
        X::copy(r,v);
        X::scal(Real(1.)/norm_r,v);

        // Clear memory for the list of Krylov vectors and insert the first
        // vector.  This completes #4.
        vs.clear();
        vs.emplace_back(std::move(X::init(rtrue)));
        X::copy(v,vs.back());

        // Find the initial right hand side for the vector Q' norm(w1) e1.  This
        // completes #5.
        scal <Real> (rst_freq+1,Real(0.),&(Qt_e1[0]),1);
        Qt_e1[0] = norm_r;

        // Clear out the Givens rotations.  This completes #6.
        Qts.clear();
    }

    // A function that has free reign to manipulate and change the stopping
    // tolerance for GMRES.  This should be used cautiously.
    template <
        typename Real,
        template <typename> class XX
    >
    struct GMRESManipulator {
        // Disallow constructors
        NO_COPY_ASSIGNMENT(GMRESManipulator);
       
        // Give an empty default constructor
        GMRESManipulator() {}

        // Application
        virtual void eval(
            Natural const & iter,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & b,
            Real & eps
        ) const = 0; 

        // Allow the derived class to deallocate memory
        virtual ~GMRESManipulator() {}
    };
    
    // An empty manipulator that does nothing 
    template <
        typename Real,
        template <typename> class XX
    >
    struct EmptyGMRESManipulator : public GMRESManipulator <Real,XX> {
        // Disallow constructors
        NO_COPY_ASSIGNMENT(EmptyGMRESManipulator);
       
        // Give an empty default constructor
        EmptyGMRESManipulator() {}

        // Application
        virtual void eval(
            Natural const & iter,
            typename XX <Real>::Vector const & x,
            typename XX <Real>::Vector const & b,
            Real & eps
        ) const { } 
    };

    // Computes the GMRES algorithm in order to solve A(x)=b.
    // (input) A : Operator that computes A(x)
    // (input) b : Right hand side
    // (input) eps : Relative stopping tolerance.  We check the relative 
    //    difference between the current and original preconditioned
    //    norm of the residual.
    // (input) iter_max : Maximum number of iterations
    // (input) rst_freq : Restarts GMRES every rst_freq iterations.  If we don't
    //    want restarting, set this to zero. 
    // (input) Ml_inv : Operator that computes the left preconditioner
    // (input) Mr_inv : Operator that computes the right preconditioner
    // (input/output) x : Initial guess of the solution.  Returns the final
    //    solution.
    // (return) (norm_rtrue,iter) : Final norm of the true residual and
    //    the number of iterations computed.  They are returned in a STL pair.
    template <
        typename Real,
        template <typename> class XX
    >
    std::pair <Real,Natural> gmres(
        Operator <Real,XX,XX> const & A,
        typename XX <Real>::Vector const & b,
        Real eps,
        Natural iter_max,
        Natural rst_freq,
        Operator <Real,XX,XX> const & Ml_inv,
        Operator <Real,XX,XX> const & Mr_inv,
        GMRESManipulator <Real,XX> const & gmanip,
        typename XX <Real>::Vector & x
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Adjust the restart frequency if it is too big
        rst_freq = rst_freq > iter_max ? iter_max : rst_freq;

        // Adjust the restart frequency if none is desired.
        rst_freq = rst_freq == 0 ? iter_max : rst_freq;

        // Allocate memory for the residual
        X_Vector r(X::init(x));
        
        // Allocate memory for the iterate update 
        X_Vector dx(X::init(x));
        
        // Allocate memory for x + dx 
        X_Vector x_p_dx(X::init(x));
        
        // Allocate memory for the true residual
        X_Vector rtrue(X::init(x));
        
        // Allocate memory for the norm of the true, preconditioned, and
        // original true norm of the residual
        Real norm_rtrue;
        Real norm_r;

        // Allocate memory for the R matrix in the QR factorization of H where
        // A V = V H + e_m' w_m
        // Note, this size is restricted to be no larger than the restart
        // frequency
        std::vector <Real> R(rst_freq*(rst_freq+1)/2);

        // Allocate memory for the normalized Krylov vector
        X_Vector v(X::init(x));

        // Allocate memory for w, the orthogonalized, but not normalized vector
        X_Vector w(X::init(x));

        // Allocate memory for the list of Krylov vectors
        std::list <X_Vector> vs;

        // Allocate memory for right hand side of the linear system, the vector
        // Q' norm(w1) e1.  Since we have a problem overdetermined by a single
        // index at each step, the size of this vector is the restart frequency
        // plus 1.
        std::vector <Real> Qt_e1(rst_freq+1);

        // Allocoate memory for the Givens rotations
        std::list <std::pair<Real,Real> > Qts;

        // Allocate a temporary work element
        X_Vector A_Mrinv_v(X::init(x));

        // Allocate memory for the subiteration number of GMRES taking into
        // account restarting
        Natural i(0);

        // Find the true residual and its norm
        A.eval(x,rtrue);
        X::scal(Real(-1.),rtrue);
        X::axpy(Real(1.),b,rtrue);
        norm_rtrue = sqrt(X::innr(rtrue,rtrue));

        // Initialize the GMRES algorithm
        resetGMRES<Real,XX> (rtrue,Ml_inv,rst_freq,v,vs,r,norm_r,
            Qt_e1,Qts);
            
        // If for some bizarre reason, we're already optimal, don't do any work 
        gmanip.eval(0,x,b,eps);
        if(norm_rtrue <= eps) iter_max=0;	

        // Iterate until the maximum iteration
        Natural iter;
        for(iter = 1; iter <= iter_max;iter++) {

            // Find the current iterate taking into account restarting
            i = iter % rst_freq;

            // We the above remainder is zero, we're on our final iteration
            // before restarting.  However, the iterate in this case is equal to
            // the restart frequency and not zero since our factorization has
            // size rst_freq x rst_freq.
            if(i == 0) i = rst_freq;

            // Find the next Krylov vector
            Mr_inv.eval(v,w);
            A.eval(w,A_Mrinv_v);
            Ml_inv.eval(A_Mrinv_v,w);

            // Orthogonalize this Krylov vector with respect to the rest
            orthogonalize <Real,XX> (vs,w,&(R[(i-1)*i/2]));

            // Find the norm of the remaining, orthogonalized vector
            Real norm_w = sqrt(X::innr(w,w));

            // Normalize the orthogonalized Krylov vector and insert it into the
            // list of Krylov vectros
            X::copy(w,v);
            X::scal(Real(1.)/norm_w,v);
            vs.emplace_back(std::move(X::init(x)));
            X::copy(v,vs.back());

            // Apply the existing Givens rotations to the new column of R
            Natural j=1;
            for(typename std::list <std::pair<Real,Real> >::iterator
                    Qt=Qts.begin();
                Qt!=Qts.end();
                Qt++
            ) { 
                rot <Real> (1,&(R[(j-1)+(i-1)*i/2]),1,&(R[j+(i-1)*i/2]),1,
                    Qt->first,Qt->second);
                j++;
            }

            // Form the new Givens rotation
            Qts.emplace_back(std::pair <Real,Real> ());
            rotg <Real> (R[(i-1)+i*(i-1)/2],norm_w,
                Qts.back().first,Qts.back().second);

            // Apply this new Givens rotation to the last element of R and 
            // norm(w).  This fixes our system R.
            rot <Real> (1,&(R[(i-1)+i*(i-1)/2]),1,
                &(norm_w),1,Qts.back().first,Qts.back().second);

            // Apply the new givens rotation to the RHS.  This also determines
            // the new norm of the preconditioned residual.
            rot <Real> (1,&(Qt_e1[i-1]),1,&(Qt_e1[i]),
                1,Qts.back().first,Qts.back().second);
            norm_r = fabs(Qt_e1[i]);
                
            // Solve for the new iterate update
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,Mr_inv,x,dx);

            // Find the current iterate, its residual, the residual's norm
            X::copy(x,x_p_dx);
            X::axpy(Real(1.),dx,x_p_dx);
            A.eval(x_p_dx,rtrue);
            X::scal(Real(-1.),rtrue);
            X::axpy(Real(1.),b,rtrue);
            norm_rtrue = sqrt(X::innr(rtrue,rtrue));

            // Adjust the stopping tolerance
            gmanip.eval(i,x_p_dx,b,eps);

            // Determine if we should exit since the norm of the true residual
            // is small
            if(norm_rtrue <= eps) break;	

            // If we've hit the restart frequency, reset the Krylov spaces and
            // factorizations
            if(i%rst_freq==0) {

                // Move to the new iterate
                X::copy(x_p_dx,x);

                // Reset the GMRES algorithm
                resetGMRES<Real,XX> (rtrue,Ml_inv,rst_freq,v,vs,r,norm_r,
                    Qt_e1,Qts);

                // Make sure to correctly indicate that we're now working on
                // iteration 0 of the next round of GMRES.  If we exit
                // immediately thereafter, we use this check to make sure we
                // don't do any additional solves for x.
                i = 0;
            }
        }

        // Adjust the iteration number if we ran out of iterations
        iter = iter > iter_max ? iter_max : iter;

        // As long as we didn't just solve for our new ierate, go ahead and
        // solve for it now.
        if(i > 0){ 
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,Mr_inv,x,dx);
            X::axpy(Real(1.),dx,x);
        }

        // Return the norm and the residual
        return std::pair <Real,Natural> (norm_rtrue,iter);
    }
    
    // B orthogonalizes a vector x to a list of other xs.  
    template <
        typename Real,
        template <typename> class XX
    >
    void Borthogonalize(
        std::list <typename XX <Real>::Vector> const & vs,
        std::list <typename XX <Real>::Vector> const & Bvs,
        typename XX <Real>::Vector & x,
        typename XX <Real>::Vector & Bx,
        std::list <Real> & R
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Orthogonalize the vectors
        R.clear();
        typename std::list <Real>::iterator beta=R.begin();
        for(typename std::list <X_Vector>::const_iterator
                v=vs.begin(),
                Bv=Bvs.begin();
            v!=vs.end();
            v++,Bv++,beta++
        ) {
            Real beta=X::innr(*Bv,x);
            X::axpy(Real(-1.)*beta,*v,x);
            X::axpy(Real(-1.)*beta,*Bv,Bx);
            R.emplace_back(beta);
        }
    }
    
    // Computes the truncated MINRES algorithm in order to solve A(x)=b.
    // (input) A : Operator that computes A(x)
    // (input) b : Right hand side
    // (input) B: Operator that computes the symmetric positive definite
    //    preconditioner
    // (input) C : Operator that modifies the shape of the trust-region.
    // (input) eps : Relative stopping tolerance.  We check the relative 
    //    difference between the current and original preconditioned
    //    norm of the residual.
    // (input) iter_max : Maximum number of iterations
    // (input) delta : Trust region radius.  
    // (input) x_cntr : Center of the trust-region. 
    // (output) x : Final solution.
    // (output) x_cp : The Cauchy-Point, which is defined as the solution x
    //     after a single iteration.
    // (output) Bnorm_r : The B-norm of the residual.  In the case that we
    //     truncate, this is the B-norm of the residual on the previous
    //     iteration.
    // (output) iter : Number of iterations computed
    template <
        typename Real,
        template <typename> class XX
    >
    void truncated_minres(
        Operator <Real,XX,XX> const & A,
        typename XX <Real>::Vector const & b,
        Operator <Real,XX,XX> const & B,
        Operator <Real,XX,XX> const & C,
        Real const & eps,
        Natural const & iter_max,
        Natural orthog_max,
        Real const & delta,
        typename XX <Real>::Vector const & x_cntr,
        typename XX <Real>::Vector & x,
        typename XX <Real>::Vector & x_cp,
        Real & Bnorm_r0,
        Real & Bnorm_r,
        Natural & iter,
        KrylovStop::t & krylov_stop
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Adjust orthog_max if it's too big
        orthog_max = orthog_max > iter_max ? iter_max : orthog_max;

        // Initialize x to zero. 
        X::zero(x);

        // Allocate memory for the iterate update 
        X_Vector dx(X::init(x));

        // Allocate memory for a few more temps 
        X_Vector ABv_last(X::init(x));
        X_Vector x_tmp1(X::init(x));
        X_Vector x_tmp2(X::init(x));

        // Allocate memory for the final column of the R matrix in the
        // QR factorization of T where
        // A V = V T + e_m' w_m
        std::list <Real> R;

        // Allocate memory for the list of Krylov vectors
        std::list <X_Vector> vs;
        std::list <X_Vector> Bvs;
                    
        // Allocate memory for the vectors that compose B V inv(R)            
        std::list <X_Vector> B_V_Rinvs;

        // Allocoate memory for the Givens rotations
        std::list <std::pair<Real,Real> > Qts;

        // Allocate memory for the last two elements of the right hand side of
        // the linear system, the vector  Q' norm(w1) e1.  
        std::vector <Real> Qt_e1(2);

        // Find the preconditioned residual
        B.eval(b,x_tmp1);
        Bnorm_r0 = sqrt(X::innr(x_tmp1,b));
        Bnorm_r = Bnorm_r0;

        // Find the initial Krylov vector. 
        X::copy(b,x_tmp1);
        X::scal(Real(1.)/Bnorm_r,x_tmp1);

        // Insert the first Krylov vector
        vs.emplace_back(std::move(X::init(x_tmp1)));
        X::copy(x_tmp1,vs.back());
        
        Bvs.emplace_back(std::move(X::init(x_tmp1)));
        B.eval(x_tmp1,Bvs.back());

        // Find the initial right hand side for the vector Q' norm(w1) e1.  
        Qt_e1[0] = Bnorm_r;
        Qt_e1[1] = Real(0.);

        // Allocate memory for the quantity x-x_cntr
        X_Vector x_m_xcntr(X::init(x));
        X::copy(x_cntr,x_m_xcntr);
        X::scal(Real(-1.),x_m_xcntr);

        // Verify that ||C(x-x_cntr)|| <= delta.  This insures that our initial
        // iterate lies inside the trust-region.  If it does not, we exit.
        C.eval(x_m_xcntr,x_tmp1);
        Real norm_C_x_m_xcntr = sqrt(X::innr(x_tmp1,x_tmp1));
        if(norm_C_x_m_xcntr > delta) {
            X::zero(x_cp);
            iter=0;
            krylov_stop = KrylovStop::InvalidTrustRegionCenter;
            return;
        }

        // Iterate until the maximum iteration
        for(iter = 1; iter <= iter_max;iter++) {

            // Find the next Krylov vector,
            // x_tmp1 <- ABv_last,
            // x_tmp2 <- B(ABv_last)
            A.eval(Bvs.back(),ABv_last);
            X::copy(ABv_last,x_tmp1);
            B.eval(ABv_last,x_tmp2);
            
            // Orthogonalize this Krylov vector with respect to the rest
            Borthogonalize <Real,XX> (vs,Bvs,x_tmp1,x_tmp2,R); 
            
            // Determine the normalizing factor, || v ||_B.
            Real Bnorm_v = sqrt(X::innr(x_tmp2,x_tmp1));

            // Note, it is possible that the normalizing factor is NaN if
            // either A or B returned a NaN or if B is really indefinite. In
            // this case our new Krylov vector is invalid and we proceed 
            // to the truncated part.
            if(Bnorm_v==Bnorm_v) {
            
                // Check if we need to eliminate any vectors for
                // orthogonalization.
                if(vs.size()==orthog_max+1) {
                    vs.pop_front();
                    Bvs.pop_front();
                }

                // Store the Krylov vector 
                vs.emplace_back(X::init(x));
                X::copy(x_tmp1,vs.back());
                X::scal(Real(1.)/Bnorm_v,vs.back());
                
                Bvs.emplace_back(X::init(x));
                X::copy(x_tmp2,Bvs.back());
                X::scal(Real(1.)/Bnorm_v,Bvs.back());
                
                // At this point, R only contains the Gram-Schmidt coefficients
                // from the orthogonalization.  In theory, we only need two
                // coefficients.  However, if we over orthogonalize, we may
                // have more.  Depending on the situation, we may have a
                // column of R that looks like
                // [ * * * * * ].  Alternatively, we may have a column that
                // looks like [ 0 0 * * * ].  In the latter case, we need
                // to pad the vector R with a zero in the front.
                if(iter - R.size() > 0) 
                    R.push_front(Real(0.));
                
                // Apply the existing Givens rotations to the new column of R
                typename std::list <Real>::iterator beta=R.begin();
                typename std::list <Real>::iterator beta_next=beta; beta_next++;
                for(typename std::list <std::pair<Real,Real> >::iterator Qt
                        =Qts.begin();
                    Qt!=Qts.end();
                    Qt++,beta++,beta_next++
                )
                    rot <Real> (1,&(*beta),1,&(*beta_next),1,
                        Qt->first,Qt->second);
               
                // Remove unneeded Givens rotations
                if(Qts.size()==orthog_max+1)
                    Qts.pop_front();

                // Form the new Givens rotation
                Qts.emplace_back(std::pair <Real,Real> ());
                rotg <Real> (
                    R.back(),Bnorm_v,Qts.back().first,Qts.back().second);

                // Apply this new Givens rotation to the last element of R and 
                // norm(w).  This fixes our system R.
                rot <Real> (1,&(R.back()),1,
                    &(Bnorm_v),1,Qts.back().first,Qts.back().second);

                // Apply the new givens rotation to the RHS.  This also 
                // determines the new B-norm of the residual.
                rot <Real> (1,&(Qt_e1[0]),1,&(Qt_e1[1]),
                    1,Qts.back().first,Qts.back().second);

                // Determine B V inv(R)
                typename std::list <X_Vector>::const_reverse_iterator
                    Bv_iter_m_1=Bvs.rbegin(); Bv_iter_m_1++;
                X::copy(*(Bv_iter_m_1),x_tmp1);
                beta=R.begin();
                for(typename std::list <X_Vector>::const_iterator
                        B_V_Rinv=B_V_Rinvs.begin();
                    B_V_Rinv!=B_V_Rinvs.end();
                    B_V_Rinv++,beta++
                ) 
                    X::axpy(-(*beta),*B_V_Rinv,x_tmp1);
                X::scal(Real(1.)/R.back(),x_tmp1);
               
                // Remove unneeded vectors in B V inv(R).
                if(B_V_Rinvs.size()==orthog_max+1) 
                    B_V_Rinvs.pop_front();

                // Add in the new B V inv(R) vector.
                B_V_Rinvs.emplace_back(std::move(X::init(x_tmp1)));
                X::copy(x_tmp1,B_V_Rinvs.back());

                // Solve for the new iterate update
                X::copy(B_V_Rinvs.back(),dx);
                X::scal(Qt_e1[0],dx);

                // Find || C( (x-x_cntr) + dx) ||
                X::copy(x_m_xcntr,x_tmp1);
                X::axpy(Real(1.),dx,x_tmp1);
                C.eval(x_tmp1,x_tmp2);
                norm_C_x_m_xcntr = sqrt(X::innr(x_tmp2,x_tmp2));
                
            // Sometimes the second Krylov vector lies in the nullspace
            // of B.  In this case, we don't have a search direction to
            // fall back on.  Hence, we need to find a valid search
            // direction. In order to do this, we solve for alpha in 
            // min_{alpha} .5 || A B v1 alpha - b ||^2_B, which is 
            // alpha = <B A B v1, b> / || B A B v1 ||^2.
            // Then, we set dx=alpha B v1.
            } else if(iter==1){
                // x_tmp1 <- B A B v1
                B.eval(ABv_last,x_tmp1);

                // Solve for alpha
                Real Bnorm_ABv1_2 = X::innr(x_tmp1,ABv_last);
                Real Binnr_ABv1_b = X::innr(x_tmp1,b);
                Real alpha=Binnr_ABv1_b / Bnorm_ABv1_2; 

                // Find dx = alpha B v1
                X::copy(Bvs.back(),dx);
                X::scal(alpha,dx);

                // Since we start with a zero initial guess, Bnorm_r0=||b||_B.
                // Hence,
                //   || A(x+dx)-b ||_B
                // = || A(dx) - b ||_B
                // = || alpha A B v1 - b ||_B
                // = < B A B v1, A B v1>alpha^2-2< B A B v1, b >alpha+< B b, b>
                // = Bnorm_Abv1_2 alpha^2 - 2 Binnr_ABv1_b alpha + Bnorm_r0^2.
                // Note, if this is near zero, it's possible for the number
                // to go negative, which messes up the output.
                Bnorm_r = sqrt(
                    Bnorm_ABv1_2*alpha*alpha
                    - Real(2.)*Binnr_ABv1_b*alpha
                    + Bnorm_r0*Bnorm_r0); 

                // Find || C( (x-x_cntr) + dx) ||
                X::copy(x_m_xcntr,x_tmp1);
                X::axpy(Real(1.),dx,x_tmp1);
                C.eval(x_tmp1,x_tmp2);
                norm_C_x_m_xcntr = sqrt(X::innr(x_tmp2,x_tmp2));
            
                // Determine if we should exit since the norm of the
                // preconditioned residual is small.  Note, since we're
                // technically in this block of code since we're detected
                // instability, we need the exit statement here otherwise
                // we'll scale the current direction to the trust-region.
                if(!(norm_C_x_m_xcntr >= delta) && Bnorm_r < eps*Bnorm_r0) {
                    // Move to the new iterate, x <- x + dx
                    X::axpy(Real(1.),dx,x);
                    X::axpy(Real(1.),dx,x_m_xcntr);

                    // Save the Cauchy-Point
                    if(iter==1) 
                        X::copy(x,x_cp);

                    // Set the stopping condition
                    krylov_stop = KrylovStop::RelativeErrorSmall;
                    break;	
                }
            }

            // Determine if we're unstable.  Right now, this is just if
            // we've calculate Bnorm_v = || v ||_B to be a NaN.
            bool instability = (Bnorm_v != Bnorm_v); 

            // Exit early if we have either exceeded our trust-region or
            // if we have detected a NaN.
            if(norm_C_x_m_xcntr >= delta || instability) {
                // If we're paying attention to the trust-region, scale the
                // step appropriately.  We do not do this if we're unstable.
                if(delta < std::numeric_limits <Real>::infinity()
                    && !instability
                ) {

                    // Find sigma so that
                    //
                    // || C( (x-x_cntr) + sigma dx) ||=delta.
                    //
                    // This can be found by finding the positive root of the 
                    // quadratic
                    //
                    // || C( (x-x_cntr) + sigma dx) ||^2 = delta^2.
                    //
                    // Specifically, we want the positive root of
                    // 
                    // sigma^2 < C(dx),C(dx) >
                    //     + sigma (2 < C(dx),C(x-x_cntr) >)
                    //     + (< C(x-x_cntr),C(x-x_cntr) > - delta^2).

                    // x_tmp1 <- C(dx)
                    C.eval(dx,x_tmp1);

                    // x_tmp2 <- C(x-x_cntr);
                    C.eval(x_m_xcntr,x_tmp2);

                    // Solve the quadratic equation for the positive root 
                    Real aa = X::innr(x_tmp1,x_tmp1);
                    Real bb = Real(2.)*X::innr(x_tmp1,x_tmp2);
                    Real cc = X::innr(x_tmp2,x_tmp2)-delta*delta;
                    Natural nroots;
                    Real r1;
                    Real r2;
                    quad_equation(aa,bb,cc,nroots,r1,r2);
                    Real sigma = r1 > r2 ? r1 : r2;

                    // Take the step, x <- x + sigma dx
                    X::axpy(sigma,dx,x);
                    X::axpy(sigma,dx,x_m_xcntr);
                }

                // Determine why we stopped
                if(iter>1 && instability)
                    krylov_stop = KrylovStop::Instability;
                else
                    krylov_stop = KrylovStop::TrustRegionViolated;
 
                // If this is the first iteration, save the Cauchy-Point
                if(iter==1) 
                    X::copy(x,x_cp);
                break;
            }

            // Shift the elements of Qt_e1
            Qt_e1[0]=Qt_e1[1];
            Qt_e1[1]=Real(0.);

            // Move to the new iterate, x <- x + dx
            X::axpy(Real(1.),dx,x);
            X::axpy(Real(1.),dx,x_m_xcntr);
               
            // Determine the B-norm of the residual
            Bnorm_r = fabs(Qt_e1[0]);

            // If this is the first iteration, save the Cauchy-Point
            if(iter==1) 
                X::copy(x,x_cp);
            
            // Determine if we should exit since the norm of the preconditioned
            // residual is small
            if(Bnorm_r <= eps*Bnorm_r0) {
                krylov_stop = KrylovStop::RelativeErrorSmall;
                break;	
            }
        }
        
        // If we've exceeded the maximum iteration, make sure to denote this
        if(iter > iter_max) {
            krylov_stop=KrylovStop::MaxItersExceeded;
            iter = iter_max;
        }
    }

    // Determines the relative error between two vectors where the second vector
    // may or may not have been initialized.  This is typically used for
    // determining the relative error between a vector and some cached value.
    template <typename Real,template <typename> class XX>
    Real rel_err_cached(
        typename XX <Real>::Vector const & x,
        std::pair <bool,typename XX <Real>::Vector> const & x_cached
    ) {
        // Create a type shortcut
        typedef XX <Real> X;

        // Create some workspace
        typename X::Vector x_tmp1(X::init(x));

        // If we've not been cached yet, return infinity
        if(!x_cached.first)
            return std::numeric_limits <Real>::infinity();

        // Otherwise, figure out the relative error
        else {
            // Figure out the residual between x_cached and x 
            X::copy(x_cached.second,x_tmp1);
            X::axpy(Real(-1.),x,x_tmp1);

            // Figure out the relative error between x and x_cached 
            Real rel_err = sqrt(X::innr(x_tmp1,x_tmp1)) /
                (std::numeric_limits <Real>::epsilon() + sqrt(X::innr(x,x)));

            // Return the relative error 
            return rel_err;
        }
    } 
//---Optizelle2---
}
//---Optizelle3---

#endif
