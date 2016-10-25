#pragma once

#include <vector>
#include <list>
#include <deque>
#include <cmath>
#include <limits>
#include <utility>
#include <string>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <random>
#include <functional>

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
            r = scale*std::sqrt(sq <Real> (a/scale) + sq <Real> (b/scale));
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
    
    // A safeguard search used primarily for inequality constraints
    template <typename Real,template <typename> class XX>
    using SafeguardSimplified = std::function <
        Real(
            typename XX <Real>::Vector const & dx_base,
            typename XX <Real>::Vector const & dx_dir
        )>;

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
        norm_v=std::sqrt(dot<Real>(
            m,&(V[ijtok(1,1,m)]),1,&(V[ijtok(1,1,m)]),1));
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
                norm_v=std::sqrt(
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
            norm_v=std::sqrt (dot <Real> (m,&(v2[0]),1,&(v2[0]),1));

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
    // (output) Roots, if they exist 
    template <typename Real>
    std::vector <Real> quad_equation(
        Real const & a,
        Real const & b,
        Real const & c
    ) {

        // It's sort of hard to tell if we have a quadratic or a linear since
        // we don't have a good way with the information provided to tell if
        // the quadratic coefficient is small.  As such, we do a bad, hard
        // check if the leading coefficient is zero.  If it is not the case,
        // we assume that we have a quadratic and we use the most stable
        // equation that we can for the root.
        if( a != Real(0.) ) { 
            auto r1 = b < Real(0.) ?
                (-b + std::sqrt(b*b-Real(4.)*a*c)) / (Real(2.)*a) :
                (Real(2.)*c) / (-b - std::sqrt(b*b-Real(4.)*a*c));
            auto r2 = b < Real(0.) ?
                (Real(2.)*c) / (-b + std::sqrt(b*b-Real(4.)*a*c)) :
                (-b - std::sqrt(b*b-Real(4.)*a*c)) / (Real(2.)*a);
            return std::vector <Real> {r1,r2};

        // Now, in the case that a is zero, but b is not, we have a linear
        // function and we can solve for the root.
        } else if( b != Real(0.)) {
            Real r1 = -c/b;
            return std::vector <Real> {r1};

        // Here, we have a constant function.  Now, we could have no roots
        // if c is zero.  Alternatively, we could have an infinity number of
        // roots if c is zero.  There's not really a good way to denote all
        // of these cases, so we just assume that c is not zero and return
        // zero roots.
        } else 
            return std::vector <Real> ();
    }

    // Reasons we stop truncated CG 
    namespace TruncatedStop{
        enum t{
            //---TruncatedStop0---
            NotConverged,              // Algorithm has not converged
            NegativeCurvature,         // Negative curvature detected
            RelativeErrorSmall,        // Relative error is small
            MaxItersExceeded,          // Maximum number of iterations exceeded
            TrustRegionViolated,       // Trust-region radius violated
            NanOperator,               // NaN detected in the operator
            NanPreconditioner,         // NaN detected in the preconditioner
            NonProjectorPreconditioner,// Detected a nonprojecting
                                       // preconditioner when one is required.
                                       // Too much inexactness in the
                                       // composite-step SQP method can trigger
                                       // this.
            NonSymmetricPreconditioner,// Detected a nonsymmetric preconditioner
            NonSymmetricOperator,      // Detected a nonsymmetric operator 
            LossOfOrthogonality,       // Loss of orthogonality between the
                                       // Krylov vectors detected
            OffsetViolatesTrustRegion, // Offset is chosen such that
                                       // || x_offset || > delta where
                                       // delta is the trust-region radius
            OffsetViolatesSafeguard,   // Offset violates the safeguard
            TooManyFailedSafeguard,    // Too many safeguarded steps have failed
            ObjectiveIncrease          // CG objective, 0.5 <ABx,Bx> - <b,Bx>
                                       // increased between iterations, which
                                       // shouldn't happen.
            //---TruncatedStop1---
        };

        // Converts the Krylov stopping condition to a string 
        std::string to_string(t const & trunc_stop);
        
        // Converts a string to a Krylov stopping condition
        t from_string(std::string const & trunc_stop);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name);
    }

    // Rotate a container, so that the first element is the last
    template <typename Xs>
    void rotate(Xs & xs) {
        xs.emplace_back(std::move(xs.front()));
        xs.pop_front();
    }

    // Store a normalized vector as well as its norm
    template <
        typename Real,
        template <typename> class XX
    >
    auto archive(
        Natural const & maxsize,
        Real const & normalization,
        typename XX <Real>::Vector const & x,
        std::deque <typename XX <Real>::Vector> & xs,
        std::deque <Real> & norm_xs
    ) -> std::function<void()> {
        // Create some type shortcuts
        typedef XX <Real> X;

        return [maxsize,&normalization,&xs,&norm_xs,&x]() {
            // If we're not storing anything, exit
            if(maxsize <= 0) return;

            // If we've too many elements, rotate the containers
            if(norm_xs.size()==maxsize) {
                rotate(xs);
                rotate(norm_xs);

            // Otherwise, allocate more memory
            } else {
                xs.emplace_back(X::init(x));
                norm_xs.emplace_back(Real(0.));
            }
            
            // Copy the element into place and normalize the vector
            X::copy(x,xs.back());
            X::scal(1/normalization,xs.back());
            norm_xs.back() = std::sqrt(X::innr(xs.back(),xs.back()));
        };
    };

    // A orthogonalizes a vector x to a list of other xs.  
    template <
        typename Real,
        template <typename> class XX
    >
    void Aorthogonalize(
        std::deque <typename XX <Real>::Vector> const & vs,
        std::deque <typename XX <Real>::Vector> const & Avs,
        Natural const & iter_max,
        typename XX <Real>::Vector & x,
        typename XX <Real>::Vector & Ax
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Orthogonalize the vectors
        for(auto iter = Natural(1);iter <= iter_max;iter++) {
            auto Av = Avs.cbegin();
            for(auto const & v : vs) {
                auto beta=X::innr(*Av,x);
                X::axpy(Real(-1.)*beta,v,x);
                X::axpy(Real(-1.)*beta,*Av,Ax);
                Av++;
            }
        }
    }

    // Grow a matrix of reals in both the rows and columns until we hit a
    // certain size
    template <typename Real,typename X>
    auto grow_matrix(
        Natural const & maxsize,
        X & x
    ) -> std::function <void()> {
        return [maxsize,&x]() {
            // If we're not storing anything, exit
            if(maxsize <= 0) return;

            // If we've hit the max number of rows, rotate items
            if(x.size()==maxsize) {

                // Rotate the first row to the back
                rotate(x);

                // Rotate the elements of each row 
                for(auto & xi : x)
                    rotate(xi);

            // Otherwise, allocate more memory
            } else {
                // Allocate a new row that's the same size as the first
                if(x.size()>=1)
                    x.push_back(x.front());

                // If we have no rows, allocate an empty container
                else
                    x.emplace_back(typename std::remove_reference<decltype(x)>
                        ::type::value_type());

                // Allocate a new colum to each row 
                for(auto & xj : x)
                    xj.emplace_back(Real(0.));
            }
        };
    }

    // Grow a vector of reals until we hit a certain size
    template <typename Real,typename X>
    auto grow_vector(
        Natural const & maxsize,
        X & x
    ) -> std::function <void()> {
        return [maxsize,&x]() {
            // If we're not storing anything, exit
            if(maxsize <= 0) return;

            // If we've hit the max number of columns, rotate items
            if(x.size()==maxsize) 
                rotate(x);

            // Otherwise, allocate more memory
            else 
                x.emplace_back(Real(0.));
        };
    }

    // 2-norm of a vector
    template <typename Real,typename X>
    auto norm_2(X const & x) -> Real {
        auto sum = Real(0.);
        for(auto & ele : x)
            sum += sq(ele);
        return std::sqrt(sum);
    };

    // Frobenius-norm of a matrix 
    template <typename Real,typename X>
    auto norm_F(X const & x) -> Real{
        auto sum = Real(0.);
        for(auto const & xi : x)
            for(auto const & xij : xi)
                sum += sq(xij);
        return std::sqrt(sum);
    }

    // Frobenius-norm of the antisymmetric part of a matrix 
    template <typename Real,typename X>
    auto norm_anti_F(X const & x) -> Real{
        auto sum = Real(0.);
        for(auto i=Natural(1);i<x.size();i++)
            for(auto j=Natural(0);j<i;j++)
                sum += sq(x[i][j]-x[j][i]);
        return std::sqrt(sum);
    }

    // Returns whether or not a truncated-CG direction is salvagable based on
    // the current stopping condition.  Basically, salvagable directions are
    // directions that exceed the trust-region, have negative curvature, 
    // violated the safeguard, or have nothing wrong with them.  Non-salvagable
    // situations are where Bdx has something like a NaN or the operators we
    // used to calculate it have something bad going on.
    auto is_Bdx_salvagable(TruncatedStop::t const & stop) -> bool;

    // Returns whether or not the exit condition is related to the truncated-CG
    // direction, Bdx
    auto is_Bdx_related(TruncatedStop::t const & stop) -> bool;

    // Computes the truncated projected conjugate gradient algorithm in order
    // to solve Ax=b where we restrict x to be in the range of B and that
    // || x + x_offset || <= delta.  The parameters are as follows.
    // 
    // (input) A : Operator in the system A B x = b
    // (input) b : Right hand side in the system A B x = b
    // (input) B : Projection in the system A B x = b
    // (input) C : Operator that modifies the shape of the trust-region
    // (input) eps : Stopping tolerance
    // (input) iter_max :  Maximum number of iterations
    // (input) orthog_storage_max : Maximum number of orthgonalizations.  If
    //     this number is 1, then we do the conjugate gradient algorithm.
    // (input) delta : Trust region radius.  If this number is infinity, we
    //     do not scale the final step if we detect negative curvature.
    // (input) x_offset : Offset for checking the TR radius 
    // (input) safeguard_failed_max : Maximum number of failed safeguard steps
    //     before exiting
    // (input) safeguard : Our safeguard function
    // (input) check_B_projector: Check preconditioner that is a projector 
    // (input) check_B_properties: Check other properties of the preconditioner 
    // (input) check_A_properties: Check properties of the operator 
    // (output) x : Final solution x
    // (output) x_cp : The Cauchy-Point, which is defined as the solution x
    //     after a single iteration
    // (output) norm_Br : The norm ||B r|| of the final residual
    // (output) iter : The number of iterations required to converge.
    // (output) stop : The reason why the method was terminated
    // (output) safeguard_failed : Number of failed safeguard steps upon exiting
    // (output) alpha_safeguard : Amount we truncated the last iteration
    template <
        typename Real,
        template <typename> class XX
    >
    void truncated_cg(
        Operator <Real,XX,XX> const & A,
        typename XX <Real>::Vector const & b,
        Operator <Real,XX,XX> const & B,
        Real const & eps,
        Natural const & iter_max,
        Natural const & orthog_storage_max,
        Natural const & orthog_iter_max,
        Real const & delta,
        typename XX <Real>::Vector const & x_offset,
        Natural const & safeguard_failed_max,
        SafeguardSimplified <Real,XX> const & safeguard,
	bool const & check_B_projector,
	bool const & check_B_properties,
	bool const & check_A_properties,
        typename XX <Real>::Vector & x,
        typename XX <Real>::Vector & x_cp,
        Real & norm_Br0,
        Real & norm_Br,
        Natural & iter,
        TruncatedStop::t & stop,
        Natural & safeguard_failed,
        Real & alpha_safeguard
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Initialize x to zero
        X::zero(x);

        // Initialize the initial Cauch point to zero
        X::zero(x_cp);

        // Initialize our iteration to zero.  If we exit during setup, we
        // report the number of iterations as zero.  Otherwise, we report
        // the number of iterations the algorithm performs.
        iter = 0;

        // At the start, we haven't converged
        stop = TruncatedStop::NotConverged;

        // We also haven't truncated anything at the start
        alpha_safeguard = Real(1.0);

        // Set a tolerance for our diagnostic checks.  For this, we pick the
        // magic number 0.5.  There's probably a better number as this is sort
        // of one of those, "Well, that looks good," numbers.
        auto const eps_diag = Real(0.5);

        // Normalization and curvature quantity
        auto Anorm_Bdx_2 = Real(0.);
        auto Anorm_Bdx = Real(0.);

        // Constant normalization, when we don't need one
        auto const one = Real(1.);
        
        // Residual for the sytem 
        auto r = X::init(x); 
        auto norm_r = Real(0.);
        auto rs = std::deque <X_Vector> ();
        auto norm_rs = std::deque <Real> ();

        // Preconditioned residual for the sytem
        auto Br = X::init(x);
        // norm_Br returned from function
        auto Brs = std::deque <X_Vector> ();
        auto norm_Brs = std::deque <Real> ();
       
        // Preconditioned directions
        auto Bdx = X::init(x);
        auto norm_Bdx = Real(0.);
        auto Bdxs = std::deque <X_Vector> ();
        auto norm_Bdxs = std::deque <Real> ();
        
        // Operator applied to the preconditioned directions
        auto ABdx = X::init(x);
        auto norm_ABdx = Real(0.);
        auto ABdxs = std::deque <X_Vector> ();
        auto norm_ABdxs = std::deque <Real> ();

        // Setup a bunch of functions to store elements
        auto archive_r = archive <Real,XX> (
            check_B_projector || check_B_properties ? orthog_storage_max : 0,
            one,
            r,
            rs,
            norm_rs);
        auto archive_Br = archive <Real,XX> (
            check_B_projector || check_B_properties ? orthog_storage_max : 0,
            one,
            Br,
            Brs,
            norm_Brs);
        auto archive_Bdx = archive <Real,XX> (
            orthog_storage_max,
            Anorm_Bdx, 
            Bdx,
            Bdxs,
            norm_Bdxs);
        auto archive_ABdx = archive <Real,XX> (
            orthog_storage_max,
            Anorm_Bdx, 
            ABdx,
            ABdxs,
            norm_ABdxs);
                
        // Allocate memory for a vector where
        //
        // X_i = <B ri, ri> / || Bri ||^2 - 1
        //
        // From this vector, we check the following
        //
        // 1. || X ||_2 > epsilon 
        //
        //    B is not a projector, which we need for the composite step SQP
        //    method.  Note, when B is a projector we have || Bri ||^2 =
        //    <Bri,Bri> = <B*Bri,ri> = <B^2ri,ri> = <Bri,ri>, so we really
        //    should get 0.
        auto B_projector = std::deque <Real> ();
        auto allocate_B_projector = grow_vector <Real>(
            check_B_projector ? orthog_storage_max : 0,
            B_projector);
        auto update_B_projector = [&]() {
            // Allocate additional memory when required
            allocate_B_projector();

            // Exit if we're not doing the check 
            if(!check_B_projector) return;
            
            // Update the elements
            auto const & Bri = Brs.back();
            auto const norm_Bri = norm_Brs.back();
            auto const & ri = rs.back();
            B_projector.back() = X::innr(Bri,ri) / sq(norm_Bri) - Real(1.);
        };
        auto is_B_projector = [&]() {
            return norm_2 <Real>(B_projector) <= eps_diag;
        };

        // Allocate memory for a matrix where
        //
        // X_ij = <B ri, rj> / (|| B ri || || rj ||)    i!=j
        //      = 0                                     i==j
        //
        // From this matrix, we check the following in order
        //
        // 1. || (X - X*)/2 ||_F > epsilon
        //
        //    B is not symmetric.  In truth, this should always be true, but
        //    users give strange things.  Basically, <Bri,rj> = <ri,B*rj> =
        //    <ri,Brj> when B is symmetric.  Hence, the norm of the
        //    anti-symmetric part should be zero.
        //
        // 2. || X ||_F > epsilon
        //
        //    Residuals are not orthogonal.  In theory, <B ri, rj> = 0 for
        //    i!=j.  As such, the matrix should be zero.  Nevertheless, since
        //    we're using a variant of Gram-Schmidt, we accumulate error here
        //    overtime.  We can mitigate this somewhat by increasing
        //    orthog_iter_max. 
        //
        // For reference, we store X in row-major format.
        auto B_properties = std::deque <std::deque <Real>> ();
        auto allocate_B_properties = grow_matrix <Real> (
            check_B_properties ? orthog_storage_max : 0,
            B_properties);
        auto update_B_properties = [&]() {
            // Allocate additional memory when required
            allocate_B_properties();

            // Exit if we're not doing the check 
            if(!check_B_properties) return;

            // Update the elements
            auto i = rs.size()-1;
            B_properties[i][i] = Real(0.);
            auto ele = [&](auto const & i,auto const & j) {
                return X::innr(Brs[i],rs[j]) / (norm_Brs[i] * norm_rs[j]);
            };
            for(auto j = Natural(0);j<i;j++) {
                B_properties[i][j]=ele(i,j);
                B_properties[j][i]=ele(j,i);
            }
        };
        auto is_B_symmetric = [&]() {
            return norm_anti_F <Real> (B_properties) <= eps_diag;
        };
        auto is_rs_orthogonal = [&]() {
            return norm_F <Real> (B_properties) <= eps_diag;
        };

        // Allocate memory for a matrix where
        //
        // X_ij = <A B dxi, B dxj> / (|| A B dxi || || B dxj ||)   i!=j
        //      = 0                                                i==j
        //
        // From this matrix, we check the following in order
        //
        // 1. || (X - X*)/2 ||_F > epsilon
        //
        //    A is not symmetric.  In truth, this should always be true, but
        //    users give strange things.  Basically, <ABdxi,ABdxj> =
        //    <Bdxi,A*Bdxj> = <Bdxi,ABdxj> when A is symmetric.  Hence, the
        //    norm of the anti-symmetric part should be zero. 
        //
        // 2. || X ||_F > epsilon
        //
        //    Preconditioned directions are not A-orthogonal.  In theory, 
        //    <A B dxi, B dxj> = 0 for i!=j.  As such, the matrix should be
        //    zero.  Nevertheless, since we're using a variant of Gram-Schmidt,
        //    we accumulate error here overtime.  We can mitigate this somewhat
        //    by increasing orthog_iter_max. 
        //
        // For reference, we store X in row-major format.
        auto A_properties = std::deque <std::deque <Real>> ();
        auto allocate_A_properties = grow_matrix <Real> (
            check_A_properties ? orthog_storage_max : 0,
            A_properties);
        auto update_A_properties = [&]() {
            // Allocate additional memory when required
            allocate_A_properties();

            // Exit if we're not doing the check 
            if(!check_A_properties) return;
            
            // Update the elements
            auto i = Bdxs.size()-1;
            A_properties[i][i] = Real(0.);
            auto ele = [&](auto const & i,auto const & j) {
                return X::innr(ABdxs[i],Bdxs[j])/(norm_ABdxs[i]*norm_Bdxs[j]);
            };
            for(auto j = Natural(0);j<i;j++) {
                A_properties[i][j]=ele(i,j);
                A_properties[j][i]=ele(j,i);
            }
        };
        auto is_A_symmetric = [&]() {
            return norm_anti_F <Real> (A_properties) <= eps_diag;
        };
        auto is_Bdxs_Aorthogonal = [&]() {
            return norm_F <Real> (A_properties) <= eps_diag;
        };
        
        // Allocate memory for the shifted iterate, x + x_offset.  Generally,
        // we care if this quantity violates the safeguard or the trust-region,
        // not whether x does directly
        auto shifted_iterate = X::init(x);
        X::copy(x_offset,shifted_iterate);
        auto norm_shifted_iterate =
            std::sqrt(X::innr(shifted_iterate,shifted_iterate));

        // Verify that || x0 + x_offset || = || x_offset || <= delta.  This
        // insures that our initial iterate lies inside the trust-region.  If
        // it does not, we exit.
        if(norm_shifted_iterate > delta) {
            stop = TruncatedStop::OffsetViolatesTrustRegion;
            return;
        }

        // Verify that x_offset obeys the safeguard.  This insures that our
        // initial iterate obeys the safeguard, which we need in order to exit
        // with a safe step later.  If it does not, we exit.
        auto zero = X::init(x);
        X::zero(zero);
        if(safeguard(zero,x_offset)<Real(1.)) {
            stop = TruncatedStop::OffsetViolatesSafeguard;
            return;
        }

        // Find the initial residual, r = A*x-b = -b, and projected residual, Br
        X::copy(b,r);
        X::scal(Real(-1.),r);
        B.eval(r,Br);
        norm_r = std::sqrt(X::innr(r,r));
        norm_Br0 = std::sqrt(X::innr(Br,Br));
        norm_Br = norm_Br0; 

        // Check for NaNs in the preconditioner.  If there's a problem here, we
        // have a difficult time recovering later, so just quit. 
        if(norm_Br!=norm_Br) {
            stop = TruncatedStop::NanPreconditioner;
            return;
        }

        // Find the projected search direction 
        X::copy(Br,Bdx);
        X::scal(Real(-1.),Bdx);
       
        // Allocate memory for the shifted trial step
        //
        // || (x + x_offset) + alpha Bdx ||
        //
        // and its norm
        auto shifted_trial = X::init(x);
        auto norm_shifted_trial = std::numeric_limits <Real>::quiet_NaN();
        
        // Track the number of iterations in a row where we violated the
        // safeguard 
        safeguard_failed = 0;

        // Given an iterate feasible with respect to the safeguard, save copies
        // of everything we need from that iteration of truncated-CG.  Well,
        // not everything needed to continue the algorithm.  Really, we want to
        // save everything we need in order to:
        //
        // 1. Calculate the error in the safeguarded solution
        //
        // 2. Acutally be able to calculate a point between this safe point
        //    and whereever the algorithm currently is
        auto x_safe = X::init(x);                    // Last safe iterate
        auto r_safe = X::init(x);                    // Last safe residual
        auto shifted_iterate_safe = X::init(x); // For finding a new safe step 
        auto Bdx_safe = X::init(x);    // For new iterate, x = x + alpha Bdx
        auto ABdx_safe = X::init(x);   // For new residual, r = r + alpha ABdx

        // Archives a set of safe iterate information 
        auto archive_iterate = [&]() {
            // Save the iterate information
            X::copy(x,x_safe);
            X::copy(r,r_safe);
            X::copy(shifted_iterate,shifted_iterate_safe);

            // Don't store the direction information.  We don't know if
            // it's good yet.
            X::zero(Bdx_safe);
            X::zero(ABdx_safe);
        };

        // Archives a set of safe direction information
        auto archive_direction = [&](auto const & alpha) {
            // If the direction looks good, keep it
            if(stop == TruncatedStop::NotConverged) {
                X::copy(Bdx,Bdx_safe);
                X::scal(alpha,Bdx_safe);
                X::copy(ABdx,ABdx_safe);
                X::scal(alpha,ABdx_safe);

            // If not, dump it
            } else {
                X::zero(Bdx_safe);
                X::zero(ABdx_safe);
            }
        };

        // At this point, we know we have a good starting iterate.  As such,
        // save everything off.
        archive_iterate();

        // Track the amount that we're about to reduce the CG objective
        //
        // 0.5 <ABx,Bx> - <b,Bx>
        //
        // when we take the step x+alpha Bdx, which is
        //
        // alpha^2/2 <ABdx,Bdx> + alpha <ABdx,x> - alpha <b,Bdx>
        //
        // or
        //
        // alpha ( <ABdx,x+alpha/2 Bdx> - <b,Bdx>)
        //
        // Note, there's probably a better way to group terms for numerical
        // error.  Anyway, CG should monotonically decrease its objective.
        // However, there's a number of things that could go wrong such as
        //
        // 1. If we exceed our precision, the algorithm sometimes tries
        //    to add new Krylov vectors beyond what should constitute a basis
        //    for the entire space.  This confuses the algorithm and the
        //    objective goes up.  If the user told us how many iterations to
        //    run, we could possibly avoid this.
        //
        // 2. Breakdown in the algorithm, by finding a solution early in the
        //    iterations, but one that numerically can't exceed our stopping
        //    tolerance. 
        // 
        // Anyway, this reduction is really what we care about
        // for optimization since as long as the CG objective goes down, we
        // know we'll get a positive predicted reduction or a descent
        // direction.
        auto x_p_ao2Bdx = X::init(x);
        auto obj_red = [&](auto const & alpha, bool const & cp=false) {
            // In general, we want this term 
            if(!cp)
                X::copy(x,x_p_ao2Bdx);

            // When we do the calculation for the Cauchy point, we've
            // already updated x, so this is a hack to undue that update
            else
                X::zero(x_p_ao2Bdx);

            // Finish the calculation
            X::axpy(Real(0.5)*alpha,Bdx,x_p_ao2Bdx);
            auto red1 = X::innr(ABdx,x_p_ao2Bdx);
            auto red2 = X::innr(b,Bdx);
            auto red3 = alpha*(red1-red2);
            return red3;
        };

        // Given a step length, take a step as long as the CG objective
        // function gets reduced.  Note, this does not consider the safeguard.
        auto step_if_obj_red = [&](auto const & alpha) {
            if(obj_red(alpha) <= Real(0.)) {
                X::axpy(alpha,Bdx,x);
                X::axpy(alpha,Bdx,shifted_iterate);
                X::axpy(alpha,ABdx,r);
                B.eval(r,Br);
                norm_r=std::sqrt(X::innr(r,r));
                norm_Br=std::sqrt(X::innr(Br,Br));
                norm_shifted_iterate =
                    std::sqrt(X::innr(shifted_iterate,shifted_iterate));
            }
        };

        // Loop until we converge (or don't)
        iter = 1;
        while(stop == TruncatedStop::NotConverged) {

            // Find the ABdx application
            A.eval(Bdx,ABdx);

            // Orthogonalize this direction to the previous directions
            Aorthogonalize <Real,XX> (Bdxs,ABdxs,orthog_iter_max,Bdx,ABdx); 

            // Check if this direction is a descent direction.  If it is not,
            // flip it so that it is.  In truth, this really shouldn't ever
            // happen, but bad operators for A and B do stranger things.
            if(X::innr(Bdx,r) > Real(0.)) {
                X::scal(Real(-1.),Bdx);
                X::scal(Real(-1.),ABdx);
            }

            // Find || Bdx ||_A^2
            Anorm_Bdx_2 = X::innr(Bdx,ABdx);

            // Check for NaNs in the operator
            if(Anorm_Bdx_2!=Anorm_Bdx_2)
                stop = TruncatedStop::NanOperator;

            // Check for negative curvature, when || Bdx ||_A^2 <= 0.  Note,
            // this also encapsulates zero curvature, which is also bad.
            if(Anorm_Bdx_2 <= Real(0.)&& stop ==TruncatedStop::NotConverged)
                stop = TruncatedStop::NegativeCurvature;

            // Allocate memory for the line-search for the optimal step
            // in the direction Bdx.
            auto alpha = std::numeric_limits<Real>::quiet_NaN();

            // We only compute the following when we have not detected some
            // kind of exiting condition
            if(stop == TruncatedStop::NotConverged) {
                // Figure out the normalization for the directions 
                // orthogonalization.
                Anorm_Bdx = std::sqrt(Anorm_Bdx_2);

                // Store our history
                archive_Bdx();
                archive_ABdx();
                archive_r();
                archive_Br();
                
                // Update our numerical error checks
                update_B_projector();
                update_B_properties();
                update_A_properties();

                // Look for numerical errors
                if(check_B_projector && !is_B_projector())
                    stop = TruncatedStop::NonProjectorPreconditioner;
                else if(check_B_properties && !is_B_symmetric())
                    stop = TruncatedStop::NonSymmetricPreconditioner;
                else if(check_B_properties && !is_rs_orthogonal())
                    stop = TruncatedStop::LossOfOrthogonality;
                else if(check_A_properties && !is_A_symmetric())
                    stop = TruncatedStop::NonSymmetricOperator;
                else if(check_A_properties && !is_Bdxs_Aorthogonal())
                    stop = TruncatedStop::LossOfOrthogonality;

                // Exit now if we have any numerical errors
                if(stop != TruncatedStop::NotConverged)
                    break;

                // Do an exact linesearch in the computed direction
                alpha = -X::innr(r,Bdx) / Anorm_Bdx_2;

                // Check that our proposed direction is going to give us
                // decrease in the CG objective
                if(obj_red(alpha)>Real(0.)&&stop==TruncatedStop::NotConverged){
                    stop=TruncatedStop::ObjectiveIncrease;
                    break;
                }

                // Determine the norm of the shifted trial step
                //
                // || (x + x_offset) + alpha Bdx ||
                //
                // We use this to determine if we've stepped outside the
                // trust-region radius.
                X::copy(shifted_iterate,shifted_trial);
                X::axpy(alpha,Bdx,shifted_trial);
                norm_shifted_trial = std::sqrt(
                    X::innr(shifted_trial,shifted_trial));

                // Check if we've met or exceeded the trust-region radius
                if(norm_shifted_trial >= delta)
                    stop = TruncatedStop::TrustRegionViolated;
            }

            // If our last iterate was safe and the direction isn't corrupted
            // for any reason, store the search direction adjusted by alpha
            if(safeguard_failed == 0)
                archive_direction(alpha);

            // If we detect any kind of exit condition at this point, resolve
            // it here.  In theory, these should be exit conditions related
            // to the search direction, Bdx.
            if( stop != TruncatedStop::NotConverged ) {
                // If our stopping condition isn't related to Bdx, then we
                // shouldn't be here and we have a bug.  Exit poorly enough to
                // cause attention.
                if(!is_Bdx_related(stop))
                    throw;

                // If our direction is corrupted for any reason, don't try
                // to recover
                if(!is_Bdx_salvagable(stop)) {

                // If the situation is salvagable, then we have negative
                // curvature or have hit the trust-region bound.  To fix
                // things, we extend the step to the trust-region when one
                // exists.  Otherwise, we retreat to the last step.  In the
                // case we're on iteration 1, we use the steepest descent
                // direction.
                } else {
                    // Amount that we cut back the step
                    auto sigma = Real(0.);

                    // When we have a trust-region, run in the current
                    // direction, Bdx, until we hit the boundary
                    if(delta < std::numeric_limits <Real>::infinity()) {
                        // Find sigma so that
                        // 
                        // || (x+x_offset) + sigma Bdx || = delta.
                        //
                        // This can be found by finding the positive root of the
                        // quadratic
                        //
                        // || (x+x_offset) + sigma Bdx ||^2 = delta^2.
                        //
                        // Specifically, we want the positive root of
                        //
                        // sigma^2 || Bdx ||^2
                        //     + sigma (2 <Bdx,x+x_offset>)
                        //     + (|| x+x_offset||^2-delta^2).

                        // Solve the quadratic equation for the positive root 
                        auto aa = X::innr(Bdx,Bdx);
                        auto bb = Real(2.)*X::innr(Bdx,shifted_iterate);
                        auto cc = norm_shifted_iterate*norm_shifted_iterate
                            - delta*delta;
                        auto roots = quad_equation(aa,bb,cc);
                        for(auto const & root : roots)
                            if(root > sigma)
                                sigma = root;

                    // When we don't have a trust region, but we're on
                    // iteration 1, move in the preconditioned steepest-descent
                    // direction
                    } else if(iter == 1) {
                        // Note, Bdx contains -Br, since we're on the first
                        // iteration, so we just take a unit step in order to
                        // move in the steepest descent direction
                        sigma = Real(1.);
                    }

                    // Determine if this new iterate is feasible with respect
                    // to our safeguard.  We calculate this from x_offset,
                    // which we assume to be a safe starting place.  In any
                    // case, if the new iterate is safe, set safeguard_failed
                    // to zero and let the code take the step down below.
                    auto trial = X::init(x);
                    X::copy(x,trial);
                    X::axpy(sigma,Bdx,trial);
                    alpha_safeguard =
                        std::min(safeguard(x_offset,trial),Real(1.0));
                    if(alpha_safeguard >= Real(1.)) {
                        safeguard_failed = 0;

                    // If this new iterate isn't safe, but the last one was,
                    // see how far we can go in the current direction.  If this
                    // amount truncates us more than sigma, then we reduce the
                    // size of sigma.
                    } else if(safeguard_failed==0) { 
                        auto sigma_Bdx = X::init(x);
                        X::copy(Bdx,sigma_Bdx);
                        X::scal(sigma,sigma_Bdx);
                        alpha_safeguard = std::min(
                            safeguard(shifted_iterate,sigma_Bdx),Real(1.0));

                    // Finally, if neither the current iterate or the last
                    // iterate are safe, then exit.  Our final exit code
                    // should resolve things.  If we're on the first iteration,
                    // the last iterate should be safe, so we'll capture
                    // the Cauchy point below.
                    } else
                        break;

                    // At this point, we should know that
                    //
                    // (x + x_offset) + sigma * alpha_safeguard Bdx
                    //
                    // is safe.  Therefore, we take the step as long as it
                    // decreases the CG objective.
                    step_if_obj_red(alpha_safeguard*sigma);

                    // If we're on the first iteration, save the Cauchy-Point 
                    if(iter==1) X::copy(x,x_cp);
                }

                // We've hit an exit criteria, so exit from the routine 
                break;
            }

            // As long as we reduce the CG objective, and we should at this
            // point, take a step
            step_if_obj_red(alpha);

            // Determine if this new iterate is feasible with respect to our
            // safeguard.  We calculate this from x_offset, which we checked to 
            // be a safe starting place.  In any case, if the new iterate is
            // safe, save it for potential use later
            alpha_safeguard = std::min(safeguard(x_offset,x),Real(1.0));
            if(alpha_safeguard < Real(1.))
                safeguard_failed += 1;
            else {
                safeguard_failed = 0;
                archive_iterate();
            }

            // If this is the first iteration and the Cauchy-point reduces the
            // CG objective, save it.
            if(iter==1 && obj_red(alpha*alpha_safeguard,true)) {
                X::copy(x,x_cp);
                X::scal(alpha_safeguard,x_cp);
            }

            // Find the projected steepest descent direction
            X::copy(Br,Bdx);
            X::scal(Real(-1.),Bdx);	

            // If we have a NaN in the preconditioner, exit 
            if(norm_Br!=norm_Br)
                stop = TruncatedStop::NanPreconditioner;

            // If we have too many failed safeguard steps, exit
            else if(safeguard_failed >= safeguard_failed_max)
                stop = TruncatedStop::TooManyFailedSafeguard;
        
            // If the norm of the residual is small relative to the starting
            // residual, exit
            else if(norm_Br <= eps*norm_Br0)
                stop = TruncatedStop::RelativeErrorSmall;

            // If we've exceeded the maximum number of iterations, also exit
            else if(iter>=iter_max)
                stop = TruncatedStop::MaxItersExceeded;

            // Otherwise, increment our iteration and keep computing
            else
                iter++;
        }

        // If our last x was safe, we assume at this point that we have
        // truncated things appropriately and that our solution in x is
        // feasible with respect to our safeguard.  If our last x was not
        // safe, then go back to our last safe x and grab its Bdx.  Then,
        // we truncate Bdx until we meet our safeguard.
        if(safeguard_failed>0) {
            // Grab our old iterate, steps, and residual
            X::copy(x_safe,x);
            X::copy(r_safe,r);
            X::copy(shifted_iterate_safe,shifted_iterate);
            X::copy(Bdx_safe,Bdx);
            X::copy(ABdx_safe,ABdx);

            // Determine how far we can go safely in the saved direction
            alpha_safeguard= std::min(safeguard(shifted_iterate,Bdx),Real(1.0));

            // Take the safeguarded step and update our residuals as long as
            // this step decreases our CG objective
            step_if_obj_red(alpha_safeguard);
        }
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
        Operator <Real,XX,XX> const & B_right,
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
        B_right.eval(V_y,dx);
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
        Operator <Real,XX,XX> const & B_left,
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

        // Apply the left preconditioner to the true residual.  This
        // completes #1
        B_left.eval(rtrue,r);

        // Store the norm of the preconditioned residual.  This completes #2.
        norm_r = std::sqrt(X::innr(r,r));

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
        NO_COPY_ASSIGNMENT(GMRESManipulator)
       
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
        NO_COPY_ASSIGNMENT(EmptyGMRESManipulator)
       
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
    // (input) B_left : Operator that computes the left preconditioner
    // (input) B_right : Operator that computes the right preconditioner
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
        Operator <Real,XX,XX> const & B_left,
        Operator <Real,XX,XX> const & B_right,
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

        // Allocate memory for the Givens rotations
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
        norm_rtrue = std::sqrt(X::innr(rtrue,rtrue));

        // Initialize the GMRES algorithm
        resetGMRES<Real,XX> (rtrue,B_left,rst_freq,v,vs,r,norm_r,
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
            B_right.eval(v,w);
            A.eval(w,A_Mrinv_v);
            B_left.eval(A_Mrinv_v,w);

            // Orthogonalize this Krylov vector with respect to the rest
            orthogonalize <Real,XX> (vs,w,&(R[(i-1)*i/2]));

            // Find the norm of the remaining, orthogonalized vector
            Real norm_w = std::sqrt(X::innr(w,w));

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
                
            // Do the follow steps twice just in case we detect a NaN
            bool nan_detected = false;
            for(Natural ii = 0;ii <= 1;ii++) { 
                // Solve for the new iterate update
                solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,B_right,x,dx);

                // Find the current iterate, its residual, the residual's norm
                X::copy(x,x_p_dx);
                X::axpy(Real(1.),dx,x_p_dx);
                A.eval(x_p_dx,rtrue);
                X::scal(Real(-1.),rtrue);
                X::axpy(Real(1.),b,rtrue);
                norm_rtrue = std::sqrt(X::innr(rtrue,rtrue));

                // If our residual is real, quit
                if(norm_rtrue==norm_rtrue)
                    break;

                // If we detect a NaN, then something has gone terribly wrong
                // during the last iteration, so eliminate the last vector and
                // quit
                else {
                    vs.pop_back();
                    iter--;
                    i--;
                    nan_detected=true;
                    // Note, norm_r will not be correct on the next iteration
                }
            }

            // If we find a NaN, quit
            if(nan_detected) break;

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
                resetGMRES<Real,XX> (rtrue,B_left,rst_freq,v,vs,r,norm_r,
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

        // As long as we didn't just solve for our new iterate, go ahead and
        // solve for it now.
        if(i > 0){ 
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,B_right,x,dx);
            X::axpy(Real(1.),dx,x);
        }

        // Return the norm and the residual
        return std::pair <Real,Natural> (norm_rtrue,iter);
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
            Real rel_err = std::sqrt(X::innr(x_tmp1,x_tmp1)) /
                (std::numeric_limits <Real>::epsilon()+std::sqrt(X::innr(x,x)));

            // Return the relative error 
            return rel_err;
        }
    } 
//---Optizelle2---
}
//---Optizelle3---
