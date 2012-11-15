#ifndef LINALG_H 
#define LINALG_H 

#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <utility>
#include <string>
#include <cstddef>

namespace peopt {
    typedef size_t Natural;
    typedef ptrdiff_t Integer;

    template <typename Real>
    void copy(Integer n,const Real* x,Integer incx,Real* y,Integer incy);
    
    template <typename Real>
    void axpy(Integer n,Real alpha,const Real* x,Integer incx,
        Real* y,Integer incy);
    
    template <typename Real>
    void scal(Integer n,const Real alpha,Real* x,Integer incx);
    
    template <typename Real>
    Real dot(Integer n,const Real* x,Integer incx,const Real* y,Integer incy);
    
    template <typename Real>
    void syr2k(char uplo,char trans,Integer n,Integer k,Real alpha,
        const Real* A,Integer lda,const Real* B,Integer ldb,
        Real beta,Real* C,Integer ldc);

    template <typename Real>
    void syevr(char jobz,char range,char uplo,Integer n,Real *A,Integer lda,
        Real vl,Real vu,Integer il,Integer iu,Real abstol,Integer& m,
        Real* w,Real* z,Integer ldz,Integer* isuppz,Real* work,Integer lwork,
        Integer* iwork,Integer liwork,Integer& info);
    
    template <typename Real>
    void stemr(char jobz,char range,Integer n,Real *D,Real *E,Real vl,Real vu,
        Integer il,Integer iu,Integer& m,Real* w,Real* z,Integer ldz,
        Integer nzc,Integer* isuppz,Integer& tryrac,Real* work,
        Integer lwork,Integer* iwork,Integer liwork,Integer& info);
    
    template <typename Real>
    void stevr(char jobz,char range,Integer n,Real *D,Real *E,Real vl,Real vu,
        Integer il,Integer iu,Real abstol, Integer& m,Real* w,Real* z,
        Integer ldz,Integer* isuppz,Real* work,Integer lwork,Integer* iwork,
        Integer liwork,Integer& info);

    template <typename Real>
    Real lamch(char cmach);

    template <typename Real>
    void gemm(char transa,char transb,Integer m,Integer n,Integer k,Real alpha,
        const Real* A,Integer lda,const Real* B,Integer ldb,Real beta,
        Real* C,Integer ldc);
    
    template <typename Real>
    void symm(char side,char uplo,Integer m,Integer n,Real alpha,const Real* A,
        Integer lda,const Real* B,Integer ldb,Real beta,Real* C,Integer ldc);
    
    template <typename Real>
    void symv(char uplo,Integer n,Real alpha,const Real* A,Integer lda,
        const Real* x,Integer incx,Real beta,Real* y,Integer incy);

    template <typename Real>
    void potrf(char uplo,Integer n,Real* A,Integer lda,Integer& info);

    template <typename Real>
    void trtri(char uplo,char diag,Integer n,Real* A,Integer lda,Integer& info);

    template <typename Real>
    void rotg(Real a,Real b,Real& c,Real& s);

    template <typename Real>
    void rot(Integer n,const Real* x,Integer incx,Real* y,Integer incy,
        Real c,Real s);

    template <typename Real>
    void tpsv(char uplo,char trans,char diag,Integer n,const Real* Ap,Real* x,
        Integer incx);

    // Indexing function for matrices
    Natural ijtok(Natural i,Natural j,Natural m);

    // Indexing for packed storage
    Natural ijtokp(Natural i,Natural j); 

    // Indexing for vectors 
    Natural itok(Natural i);

    // A simple operator specification, A : X->Y
    template <
        typename Real,
        template <typename> class X,
        template <typename> class Y
    >
    struct Operator {
    private:
        // Create some type shortcuts
        typedef typename X <Real>::Vector X_Vector;
        typedef typename Y <Real>::Vector Y_Vector;
    public:
        // Basic application
        virtual void operator () (const X_Vector& x,Y_Vector &y) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~Operator() {}
    };

    /* Given a Schur decomposition of A, A=V D V', solve the Sylvester equation
    
       A X + X A = B

    */
    template <typename Real>
    void sylvester(
        const Natural m,
        const Real* V,
        const Real* D,
        const Real* B,
        Real* X
    ) {

        // Find V' B V
        std::vector <Real> tmp(m*m);
        std::vector <Real> VtBV(m*m);
        // tmp <- B V
        symm <Real> ('L','U',Integer(m),Integer(m),Real(1.),&(B[0]),Integer(m),
            &(V[0]),Integer(m),Real(0.),&(tmp[0]),Integer(m)); 
        // VtBV <- V' B V
        gemm <Real> ('T','N',Integer(m),Integer(m),Integer(m),Real(1.),&(V[0]),
            Integer(m),&(tmp[0]),Integer(m),Real(0.),&(VtBV[0]),Integer(m));

        // Solve for each column of X.  In theory, we only need half of these
        // elements since X is symmetric.
        #pragma omp parallel for schedule(static)
        for(Natural j=1;j<=m;j++) {
            for(Natural i=1;i<=j;i++) 
                X[ijtok(i,j,m)]=VtBV[ijtok(i,j,m)]/(D[i-1]+D[j-1]);
        }

        // Realransform the solution back, X = V X V'
        // tmp <- V X
        symm <Real> ('R','U',Integer(m),Integer(m),Real(1.),&(X[0]),Integer(m),
            &(V[0]),Integer(m),Real(0.),&(tmp[0]),Integer(m));
        // X <- V X V'
        gemm <Real> ('N','T',Integer(m),Integer(m),Integer(m),Real(1.),
            &(tmp[0]),Integer(m),&(V[0]),Integer(m),Real(0.),&(X[0]),
            Integer(m));
    }

    // Find a bound on the smallest eigenvalue of the given matrix A such
    // that lambda_min(A) < alpha where alpha is returned from this function.
    template <typename Real>
    Real lanczos(
        const Natural m,
        const Real* A,
        const Natural max_iter,
        const Real tol
    ) {
        // Create the initial Krylov vector
        std::vector <Real> v(m,Real(1./std::sqrt(Real(m))));

        // Get the next Krylov vector and orthgonalize it
        std::vector <Real> w(m);
        // w <- A v
        symv <Real> ('U',Integer(m),Real(1.),&(A[0]),Integer(m),&(v[0]),
            Integer(1),Real(0.),&(w[0]),Integer(1));
        // alpha[0] <- <Av,v>
        std::vector <Real> alpha;
        alpha.push_back(dot <Real> (Integer(m),&(w[0]),Integer(1),&(v[0]),
            Integer(1)));
        // w <- Av - <Av,v> v
        axpy <Real> (Integer(m),-alpha[0],&(v[0]),Integer(1),&(w[0]),
            Integer(1));

        // Store the norm of the Arnoldi vector w in the off diagonal part of T.
        // By T, we mean the triagonal matrix such that A = Q T Q'.
        std::vector <Real> beta;
        beta.push_back(std::sqrt(dot <Real> (m,&(w[0]),1,&(w[0]),1)));

        // Allocate memory for solving an eigenvalue problem for the Ritz
        // values and vectors later.
        std::vector <Integer> isuppz;
        std::vector <Real> work(1);
        std::vector <Integer> iwork(1);
        Integer lwork=-1;
        Integer liwork=-1;
        Integer info;
        Integer nevals;
        //Integer nzc=0;
        std::vector <Real> W;
        std::vector <Real> Z;
        std::vector <Real> D;
        std::vector <Real> E;

        // Start Lanczos
        std::vector <Real> v_old(m);
        for(Natural i=0;i<max_iter;i++) {
            // Save the current Arnoldi vector
            copy <Real> (Integer(m),&(v[0]),Integer(1),&(v_old[0]),Integer(1));

            // Copy the candidate Arnoldi vector to the current Arnoldi vector
            copy <Real> (Integer(m),&(w[0]),Integer(1),&(v[0]),Integer(1));

            // Get the normalized version of this vector.  Realhis is now a real
            // Arnoldi vector.
            scal <Real> (Integer(m),Real(1.)/beta[i],&(v[0]),Integer(1));

            // Get the new Arnoldi vector, w <- A v
            symv <Real> ('U',Integer(m),Real(1.),&(A[0]),Integer(m),&(v[0]),
                Integer(1),Real(0.),&(w[0]),Integer(1));

            // Orthogonalize against v_old and v using modified Gram-Schdmit.

            // First, we orthogonalize against v_old
            // w <- Av - <Av,v_old> v_old.  Due to symmetry, <Av,v_old>=beta.
            axpy <Real> (Integer(m),-beta[i],&(v_old[0]),Integer(1),&(w[0]),
                Integer(1));

            // Now, we orthogonalize against v
            // Find the Gram-Schmidt coefficient
            alpha.push_back(dot <Real> (Integer(m),&(w[0]),Integer(1),&(v[0]),
                Integer(1)));
            // Orthogonlize w to v
            axpy <Real> (Integer(m),-alpha[i+1],&(v[0]),Integer(1),&(w[0]),
                Integer(1));

            // Store the norm of the Arnoldi vector w in the off diagonal part
            // of Real.
            beta.push_back(std::sqrt(dot <Real> (Integer(m),&(w[0]),Integer(1),
                &(w[0]),Integer(1))));
   
#if 0
            // Figure out the workspaces for the eigenvalues and eigenvectors
            Integer k=alpha.size();  // Size of the eigenvalue subproblem
            D.resize(alpha.size());
            copy <Real> (k,&(alpha[0]),1,&(D[0]),1);
            E.resize(beta.size());
            copy <Real> (k,&(beta[0]),1,&(E[0]),1);
            isuppz.resize(2*k);
            lwork=-1;
            liwork=-1;
            W.resize(k);
            Z.resize(k*k);
            peopt::stemr <double> ('V','A',k,&(D[0]),&(E[0]),double(0.),
                double(0.),0,0,nevals,&(W[0]),&(Z[0]),k,k,&(isuppz[0]),
                nzc,&(work[0]),lwork,&(iwork[0]),liwork,info);

            // Resize the workspace 
            lwork = Integer(work[0])+Integer(1);
            work.resize(lwork);
            liwork = iwork[0];
            iwork.resize(liwork);

            // Find the eigenvalues and vectors 
            peopt::stemr <double> ('V','A',k,&(D[0]),&(E[0]),double(0.),
                double(0.),0,0,nevals,&(W[0]),&(Z[0]),k,k,&(isuppz[0]),
                nzc,&(work[0]),lwork,&(iwork[0]),liwork,info);
#else
            // Figure out the workspaces for the eigenvalues and eigenvectors
            Natural k=alpha.size();  // Size of the eigenvalue subproblem
            D.resize(alpha.size());
            copy <Real> (Integer(k),&(alpha[0]),Integer(1),&(D[0]),Integer(1));
            E.resize(beta.size());
            copy <Real> (Integer(k),&(beta[0]),Integer(1),&(E[0]),Integer(1));
            isuppz.resize(Natural(2)*k);
            lwork=-1;
            liwork=-1;
            W.resize(k);
            Z.resize(k*k);
            peopt::stevr <Real> ('V','A',Integer(k),&(D[0]),&(E[0]),Real(0.),
                Real(0.),Integer(0),Integer(0),peopt::lamch <Real> ('S'),
                nevals,&(W[0]),&(Z[0]),Integer(k),&(isuppz[0]),&(work[0]),
                lwork,&(iwork[0]),liwork,info);

            // Resize the workspace 
            lwork = Integer(work[0])+Integer(1);
            work.resize(Natural(lwork));
            liwork = iwork[0];
            iwork.resize(liwork);

            // Find the eigenvalues and vectors 
            peopt::stevr <Real> ('V','A',Integer(k),&(D[0]),&(E[0]),Real(0.),
                Real(0.),Integer(0),Integer(0),peopt::lamch <Real> ('S'),
                nevals,&(W[0]),&(Z[0]),Integer(k),&(isuppz[0]),&(work[0]),
                lwork,&(iwork[0]),liwork,info);
#endif

            // Find beta_i |s_{ik}| where s_{ik} is the ith (last) element
            // of the kth Ritz vector where k corresponds to the largest
            // and smallest Ritz values.  Basically, we don't know which is
            // going to converge first, but they'll both be the first two.
            // Hence, we converge until these errors estimates are accurate
            // enough.
            Real err_est_min = fabs(Z[ijtok(k,1,k)])*beta[i+1];
            Real err_est_max = fabs(Z[ijtok(k,k,k)])*beta[i+1];

            // Stop of the error estimates are small
            if(    err_est_min < fabs(W[0]) * tol
                && err_est_max < fabs(W[i]) * tol
            )
                break;
        }

        // Return the smallest Ritz value
        return W[0];
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
        const Real& a,
        const Real& b,
        const Real& c,
        Natural& nroots,
        Real& r1,
        Real& r2
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
            nroots = Natural(2);

        // Now, in the case that a is zero, but b is not, we have a linear
        // function and we can solve for the root.
        } else if( b != Real(0.)) {
            r1 = -c/b;
            nroots = Natural(1);

        // Here, we have a constant function.  Now, we could have no roots
        // if c is zero.  Alternatively, we could have an infinity number of
        // roots if c is zero.  There's not really a good way to denote all
        // of these cases, so we just assume that c is not zero and return
        // zero roots.
        } else {
            nroots = Natural(0);
        }
    }

    // Reasons we stop the Krylov method
    struct KrylovStop{
        enum t{
            NegativeCurvature,        // Negative curvature detected
            RelativeErrorSmall,       // Relative error is small
            MaxItersExceeded,         // Maximum number of iterations exceeded
            TrustRegionViolated       // Trust-region radius violated
        };

        // Converts the Krylov stopping condition to a string 
        static std::string to_string(t krylov_stop){
            switch(krylov_stop){
            case NegativeCurvature:
                return "NegativeCurvature";
            case RelativeErrorSmall:
                return "RelativeErrorSmall";
            case MaxItersExceeded:
                return "MaxItersExceeded";
            case TrustRegionViolated:
                return "TrustRegionViolated";
            default:
                throw;
            }
        }
        
        // Converts a string to a Krylov stopping condition
        static t from_string(std::string krylov_stop){
            if(krylov_stop=="NegativeCurvature")
                return NegativeCurvature;
            else if(krylov_stop=="RelativeErrorSmall")
                return RelativeErrorSmall;
            else if(krylov_stop=="MaxItersExceeded")
                return MaxItersExceeded;
            else if(krylov_stop=="TrustRegionViolated")
                return TrustRegionViolated;
            else
                throw;
        }

        // Checks whether or not a string is valid
        struct is_valid : public std::unary_function<std::string, bool> {
            bool operator () (const std::string& name) const {
                if( name=="NegativeCurvature" ||
                    name=="RelativeErrorSmall" ||
                    name=="MaxItersExceeded" ||
                    name=="TrustRegionViolated" 
                )
                    return true;
                else
                    return false;
            }
        };
    };

    // A B orthogonalizes a vector x to a list of other xs.  
    template <
        typename Real,
        template <typename> class XX
    >
    void ABorthogonalize(
        const std::list <typename XX <Real>::Vector>& vs,
        const std::list <typename XX <Real>::Vector>& Bvs,
        const std::list <typename XX <Real>::Vector>& ABvs,
        typename XX <Real>::Vector& x,
        typename XX <Real>::Vector& Bx,
        typename XX <Real>::Vector& ABx
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Orthogonalize the vectors
        for(typename std::list <X_Vector>::const_iterator
                v=vs.begin(),
                Bv=Bvs.begin(),
                ABv=ABvs.begin();
            v!=vs.end();
            v++,Bv++,ABv++
        ) {
            Real beta=X::innr(*ABv,Bx);
            X::axpy(Real(-1.)*beta,*v,x);
            X::axpy(Real(-1.)*beta,*Bv,Bx);
            X::axpy(Real(-1.)*beta,*ABv,ABx);
        }
    }

    // Computes the truncated projected conjugate direction algorithm in order
    // to solve Ax=b where we restrict x to be in the range of B and that
    // || C x || <= delta.  The parameters are as follows.
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
    // (output) x : Final solution x.
    // (output) x_cp : The Cauchy-Point, which is defined as the solution x
    //     after a single iteration.
    // (output) norm_r : The norm of the final residual.
    // (output) iter : The number of iterations required to converge. 
    // (output) krylov_stop : The reason why the Krylov method was terminated.
    template <
        typename Real,
        template <typename> class XX
    >
    void truncated_pcd(
        const Operator <Real,XX,XX>& A,
        const typename XX <Real>::Vector& b,
        const Operator <Real,XX,XX>& B,
        const Operator <Real,XX,XX>& C,
        const Real eps,
        const Natural iter_max,
        const Natural orthog_max,
        const Real delta,
        typename XX <Real>::Vector& x,
        typename XX <Real>::Vector& x_cp,
        Real& norm_r,
        Natural& iter,
        KrylovStop::t& krylov_stop
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Allocate memory for the search direction, its projection, and the
        // the operator applied to the projection
        X_Vector p; X::init(x,p);
        X_Vector Bp; X::init(x,Bp);
        X_Vector ABp; X::init(x,ABp);

        // Allocate memory for the previous search directions
        std::list <X_Vector> ps;
        std::list <X_Vector> Bps;
        std::list <X_Vector> ABps;

        // Initialize x to zero. 
        X::zero(x);
        Real norm_Cx = Real(0.);

        // Allocate memory for and find the initial residual, A*x-b = -b
        X_Vector r; X::init(x,r);
        X::copy(b,r);
        X::scal(Real(-1.),r);

        // Find the norm of the residual and save the original residual
        norm_r = sqrt(X::innr(r,r));
        Real norm_r0 = norm_r;

        // Allocate memory for two additional work vectors.  We use these
        // when calcuating the norm of the operator C applied to the
        // trial step.  In theory, we could get rid of at least one of these
        // work vectors.  However, by doing so, we have to calculate this
        // norm implicitely, which involves additional inner products and adding
        // the results together.  This can have some numerical difficulties,
        // so we just sacrifice the additional memory.
        X_Vector x_tmp1; X::init(x,x_tmp1); 
        X_Vector x_tmp2; X::init(x,x_tmp2); 

        // Loop until the maximum iteration
        for(iter=Natural(1);iter<=iter_max;iter++){
        
            // If the norm of the residual is small relative to the starting
            // residual, exit
            if(norm_r < eps*norm_r0) {
                iter--;
                krylov_stop = KrylovStop::RelativeErrorSmall;
                break;
            }

            // Find the steepest descent search direction
            X::copy(r,p);
            X::scal(Real(-1.),p);	

            // Find the Bp and ABp applications 
            B(p,Bp);
            A(Bp,ABp);

            // Orthogonalize this direction to the previous directions
            ABorthogonalize <Real,XX> (ps,Bps,ABps,p,Bp,ABp); 

            // Check if this direction is a descent direction.  If it is not,
            // flip it so that it is.
            if(X::innr(Bp,r) > Real(0.)) {
                X::scal(Real(-1.),p);
                X::scal(Real(-1.),Bp);
                X::scal(Real(-1.),ABp);
            }

            // Check if we need to eliminate any vectors for orthogonalization.
            if(ps.size()==orthog_max) {
                ps.pop_front();
                Bps.pop_front();
                ABps.pop_front();
            }

            // Store the previous directions
            Real innr_Bp_ABp = X::innr(Bp,ABp);
            Real one_over_sqrt_innr_Bp_ABp = Real(1.)/sqrt(innr_Bp_ABp);

            ps.push_back(X_Vector()); X::init(x,ps.back());
            X::copy(p,ps.back());
            X::scal(one_over_sqrt_innr_Bp_ABp,ps.back());
            
            Bps.push_back(X_Vector()); X::init(x,Bps.back());
            X::copy(Bp,Bps.back());
            X::scal(one_over_sqrt_innr_Bp_ABp,Bps.back());
            
            ABps.push_back(X_Vector()); X::init(x,ABps.back());
            X::copy(ABp,ABps.back());
            X::scal(one_over_sqrt_innr_Bp_ABp,ABps.back());

            // Do an exact linesearch in the computed direction
            Real alpha = -X::innr(r,Bp) / innr_Bp_ABp;

            // Find the trial step and the norm || C(x+alpha Bp) || 

            // x_tmp1 <- x + alpha Bp
            X::copy(x,x_tmp1);
            X::axpy(alpha,Bp,x_tmp1);

            // x_tmp2 <- C(x+alpha Bp)
            C(x_tmp1,x_tmp2);

            // norm_CxpaBp <- || C(x+alpha Bp) ||
            Real norm_CxpaBp = sqrt(X::innr(x_tmp2,x_tmp2));

            // If we have negative curvature or our trial point is outside the
            // trust-region radius, terminate truncated-PCD and find our final
            // step.  We have the <Bp,ABp> != <Bp,ABp> check in order to trap
            // NaNs.
            if( innr_Bp_ABp <= Real(0.) ||
                norm_CxpaBp >= delta ||
                innr_Bp_ABp != innr_Bp_ABp 
            ) {
                // If we're paying attention to the trust-region, scale the
                // step appropriately.
                if(delta < std::numeric_limits <Real>::infinity()) {
                    // Find sigma so that || C(x + sigma p) || =delta.  This can
                    // be found by finding the positive root of the quadratic
                    // || C(x) + sigma C(Bp) ||^2 = delta^2.  Specifically, we
                    // want the positive root of
                    // sigma^2<C(Bp),C(Bp)>
                    //     + sigma(2 <C(Bp),C(x)>)
                    //     + (<C(x),C(x)>-delta^2).

                    // x_tmp1 <- C(Bp)
                    C(Bp,x_tmp1);

                    // x_tmp2 <- C(x)
                    C(x,x_tmp2);

                    // Solve the quadratic equation for the positive root 
                    Real aa = X::innr(x_tmp1,x_tmp1);
                    Real bb = Real(2.)*X::innr(x_tmp1,x_tmp2);
                    Real cc = norm_Cx*norm_Cx-delta*delta;
                    Natural nroots;
                    Real r1;
                    Real r2;
                    quad_equation(aa,bb,cc,nroots,r1,r2);
                    Real sigma = r1 > r2 ? r1 : r2;

                    // Take the step, find its residual, and compute the 
                    // residual's norm
                    X::axpy(sigma,Bp,x);
                    X::axpy(sigma,ABp,r);
                    norm_r=sqrt(X::innr(r,r));

                // Otherwise, just take a step with a unit scale
                } else {
                    X::axpy(Real(1.),Bp,x);
                    X::axpy(Real(1.),ABp,r);
                    norm_r=sqrt(X::innr(r,r));
                }

                // Determine why we stopped
                if(innr_Bp_ABp <= Natural(0) || innr_Bp_ABp != innr_Bp_ABp)
                    krylov_stop = KrylovStop::NegativeCurvature;
                else
                    krylov_stop = KrylovStop::TrustRegionViolated;
 

                // If this is the first iteration, save the Cauchy-Point
                if(iter==Natural(1)) X::copy(x,x_cp);
                break;
            }

            // Take a step in this direction
            X::axpy(alpha,Bp,x);

            // If this is the first iteration, save the Cauchy-Point
            if(iter==Natural(1)) X::copy(x,x_cp);

            // Update the norm of x
            norm_Cx = norm_CxpaBp;

            // Find the new residual
            X::axpy(alpha,ABp,r);

            // Compute the norm of the residual
            norm_r=sqrt(X::innr(r,r));
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
        std::vector <Real>& x
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
            p[1]=Natural(1); p[2]=Natural(2);
            q[1]=Natural(1); q[2]=Natural(2);
        } else if(i==1) {
            p[1]=Natural(2); p[2]=Natural(1);
            q[1]=Natural(1); q[2]=Natural(2);
        } else {
            p[1]=Natural(2); p[2]=Natural(1);
            q[1]=Natural(2); q[2]=Natural(1);
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
        const std::vector <Real>& A,
        const std::vector <Real>& a,
        const std::vector <Real>& x
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
        const std::vector <Real>& A,
        const std::vector <Real>& a,
        const std::vector <Real>& lb,
        const std::vector <Real>& ub,
        std::vector <Real>& x
    ) {

        // Store the best objective function up until this point
        Real f_x = std::numeric_limits <Real>::infinity();
        
        // List all of the points to check for feasibility and optimality
        std::list <std::vector <Real> > zs;
        std::vector <Real> z(2);

        // Unconstrained minimum
        std::vector <Real> minus_a(2); minus_a[0]=-a[0]; minus_a[1]=-a[1];
        solve2x2 <Natural,Real> (A,minus_a,z);
        zs.push_back(z);  

        // z1 to the lower bound
        z[0] = lb[0];
        z[1] = -(a[1]+Real(2.)*A[0]*A[1]*z[0])/(Real(2.)*A[2]);
        zs.push_back(z);
       
        // z2 to the lower bound
        z[1] = lb[1];
        z[0] = -(a[0]+Real(2.)*A[0]*A[1]*z[1])/(2*A[0]); 
        zs.push_back(z);
        
        // z1 to the upper bound 
        z[0] = ub[0];
        z[1] = -(a[1]+Real(2.)*A[0]*A[1]*z[0])/(2*A[2]);
        zs.push_back(z);
        
        // z2 to the upper bound
        z[1] = ub[1];
        z[0] = -(a[0]+Real(2.)*A[0]*A[1]*z[1])/(2*A[0]); 
        zs.push_back(z);
       
        // Lower left corner
        z[0] = lb[0];
        z[1] = lb[1]; 
        zs.push_back(z);
        
        // Lower right corner
        z[0] = ub[0];
        z[1] = lb[1]; 
        zs.push_back(z);

        // Upper right corner
        z[0] = ub[0];
        z[1] = ub[1]; 
        zs.push_back(z);

        // Upper left corner
        z[0] = lb[0];
        z[1] = ub[1]; 
        zs.push_back(z);

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
        const std::list <typename XX <Real>::Vector>& vs,
        typename XX <Real>::Vector& x,
        Real* R
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
        const Natural& m,
        const Real* R,
        const Real* Qt_e1,
        const std::list <typename XX <Real>::Vector>& vs,
        const Operator <Real,XX,XX>& Mr_inv,
        const typename XX <Real>::Vector& x,
        typename XX <Real>::Vector& dx
    ) {
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;
        
        // Allocate memory for the solution of the triangular solve 
        std::vector <Real> y(m);

        // Create one temporary element required to solve for the iterate
        X_Vector V_y; X::init(x,V_y);

        // Solve the system for y
        copy <Real> (Integer(m),&(Qt_e1[0]),Integer(1),&(y[0]),Integer(1));
        tpsv <Real> ('U','N','N',Integer(m),&(R[0]),&(y[0]),Integer(1));

        // Compute tmp = V y
        X::zero(V_y);
        typename std::list <X_Vector>::const_iterator vv=vs.begin();
        for(Natural j=0;j<m;j++) {
            X::axpy(Real(y[j]),*vv,V_y);
            vv++;
        }

        // Right recondition the above linear combination
        Mr_inv(V_y,dx);
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
        const typename XX <Real>::Vector& rtrue,
        const Operator <Real,XX,XX>& Ml_inv,
        const Natural& rst_freq,
        typename XX <Real>::Vector& v,
        std::list <typename XX <Real>::Vector>& vs,
        typename XX <Real>::Vector& r,
        Real& norm_r,
        std::vector <Real>& Qt_e1,
        std::list <std::pair<Real,Real> >& Qts
    ){
        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Apply the left preconditioner to the true residual.  This
        // completes #1
        Ml_inv(rtrue,r);

        // Store the norm of the preconditioned residual.  This completes #2.
        norm_r = sqrt(X::innr(r,r));

        // Find the initial Krylov vector.  This completes #3.
        X::copy(r,v);
        X::scal(Real(1.)/norm_r,v);

        // Clear memory for the list of Krylov vectors and insert the first
        // vector.  This completes #4.
        vs.clear();
        vs.push_back(X_Vector());
        X::init(rtrue,vs.back());
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
        // Application
        virtual void operator () (
            const typename XX <Real>::Vector& b,
            const typename XX <Real>::Vector& x,
            Real& eps
        ) const {}

        // Allow the derived class to deallocate memory
        virtual ~GMRESManipulator() {}
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
        const Operator <Real,XX,XX>& A,
        const typename XX <Real>::Vector& b,
        Real eps,
        Natural iter_max,
        Natural rst_freq,
        const Operator <Real,XX,XX>& Ml_inv,
        const Operator <Real,XX,XX>& Mr_inv,
        const GMRESManipulator <Real,XX>& gmanip,
        typename XX <Real>::Vector& x
    ){

        // Create some type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Adjust the restart frequency if it is too big
        rst_freq = rst_freq > iter_max ? iter_max : rst_freq;

        // Adjust the restart frequency if none is desired.
        rst_freq = rst_freq == Natural(0) ? iter_max : rst_freq;

        // Allocate memory for the residual
        X_Vector r; X::init(x,r);
        
        // Allocate memory for the iterate update 
        X_Vector dx; X::init(x,dx);
        
        // Allocate memory for x + dx 
        X_Vector x_p_dx; X::init(x,x_p_dx);
        
        // Allocate memory for the true residual
        X_Vector rtrue; X::init(x,rtrue);
        
        // Allocate memory for the norm of the true, preconditioned, and
        // original true norm of the residual
        Real norm_rtrue;
        Real norm_r;

        // Allocate memory for the R matrix in the QR factorization of H where
        // A V = V H + e_m' w_m
        // Note, this size is restricted to be no larger than the restart
        // frequency
        std::vector <Real> R(rst_freq*(rst_freq+Natural(1))/Natural(2));

        // Allocate memory for the normalized Krylov vector
        X_Vector v; X::init(x,v);

        // Allocate memory for w, the orthogonalized, but not normalized vector
        X_Vector w; X::init(x,w);

        // Allocate memory for the list of Krylov vectors
        std::list <X_Vector> vs;

        // Allocate memory for right hand side of the linear system, the vector
        // Q' norm(w1) e1.  Since we have a problem overdetermined by a single
        // index at each step, the size of this vector is the restart frequency
        // plus 1.
        std::vector <Real> Qt_e1(rst_freq+Natural(1));

        // Allocoate memory for the Givens rotations
        std::list <std::pair<Real,Real> > Qts;

        // Allocate a temporary work element
        X_Vector A_Mrinv_v; X::init(x,A_Mrinv_v);

        // Allocate memory for the subiteration number of GMRES taking into
        // account restarting
        Natural i;

        // Find the true residual and its norm
        A(x,rtrue);
        X::scal(Real(-1.),rtrue);
        X::axpy(Real(1.),b,rtrue);
        norm_rtrue = sqrt(X::innr(rtrue,rtrue));

        // Initialize the GMRES algorithm
        resetGMRES<Real,XX> (rtrue,Ml_inv,rst_freq,v,vs,r,norm_r,
            Qt_e1,Qts);
            
        // If for some bizarre reason, we're already optimal, don't do any work 
        gmanip(x,b,eps);
        if(norm_rtrue < eps) iter_max=0;	

        // Iterate until the maximum iteration
        Natural iter;
        for(iter = Natural(1); iter <= iter_max;iter++) {

            // Find the current iterate taking into account restarting
            i = iter % rst_freq;

            // We the above remainder is zero, we're on our final iteration
            // before restarting.  However, the iterate in this case is equal to
            // the restart frequency and not zero since our factorization has
            // size rst_freq x rst_freq.
            if(i == Natural(0)) i = rst_freq;

            // Find the next Krylov vector
            Mr_inv(v,w);
            A(w,A_Mrinv_v);
            Ml_inv(A_Mrinv_v,w);

            // Orthogonalize this Krylov vector with respect to the rest
            orthogonalize <Real,XX> (vs,w,&(R[(i-1)*i/2]));

            // Find the norm of the remaining, orthogonalized vector
            Real norm_w = sqrt(X::innr(w,w));

            // Normalize the orthogonalized Krylov vector and insert it into the
            // list of Krylov vectros
            X::copy(w,v);
            X::scal(Real(1.)/norm_w,v);
            vs.push_back(X_Vector()); X::init(x,vs.back());
            X::copy(v,vs.back());

            // Apply the existing Givens rotations to the new column of R
            Natural j=1;
            for(typename std::list <std::pair<Real,Real> >::iterator
                    Qt=Qts.begin();
                Qt!=Qts.end();
                Qt++
            ) { 
                rot <Real> (Integer(1),&(R[(j-1)+(i-1)*i/2]),
                    Integer(1),&(R[j+(i-1)*i/2]),Integer(1),
                    Qt->first,Qt->second);
                j++;
            }

            // Form the new Givens rotation
            Qts.push_back(std::pair <Real,Real> ());
            rotg <Real> (R[(i-1)+i*(i-1)/Natural(2)],norm_w,
                Qts.back().first,Qts.back().second);

            // Apply this new Givens rotation to the last element of R and 
            // norm(w).  This fixes our system R.
            rot <Real> (Integer(1),&(R[(i-1)+i*(i-1)/2]),Integer(1),
                &(norm_w),Integer(1),Qts.back().first,Qts.back().second);

            // Apply the new givens rotation to the RHS.  This also determines
            // the new norm of the preconditioned residual.
            rot <Real> (Integer(1),&(Qt_e1[i-1]),Integer(1),&(Qt_e1[i]),
                Integer(1),Qts.back().first,Qts.back().second);
            norm_r = fabs(Qt_e1[i]);
                
            // Solve for the new iterate update 
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,Mr_inv,x,dx);

            // Find the current iterate, its residual, the residual's norm
            X::copy(x,x_p_dx);
            X::axpy(Real(1.),dx,x_p_dx);
            A(x_p_dx,rtrue);
            X::scal(Real(-1.),rtrue);
            X::axpy(Real(1.),b,rtrue);
            norm_rtrue = sqrt(X::innr(rtrue,rtrue));

            // Adjust the stopping tolerance
            gmanip(x_p_dx,b,eps);

            // Determine if we should exit since the norm of the true residual
            // is small
            if(norm_rtrue < eps) break;	

            // If we've hit the restart frequency, reset the Krylov spaces and
            // factorizations
            if(i%rst_freq==Natural(0)) {

                // Move to the new iterate
                X::copy(x_p_dx,x);

                // Reset the GMRES algorithm
                resetGMRES<Real,XX> (rtrue,Ml_inv,rst_freq,v,vs,r,norm_r,
                    Qt_e1,Qts);

                // Make sure to correctly indicate that we're now working on
                // iteration 0 of the next round of GMRES.  If we exit
                // immediately thereafter, we use this check to make sure we
                // don't do any additional solves for x.
                i = Natural(0);
            }
        }

        // Adjust the iteration number if we ran out of iterations
        iter = iter > iter_max ? iter_max : iter;

        // As long as we didn't just solve for our new ierate, go ahead and
        // solve for it now.
        if(i > Natural(0)){ 
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,Mr_inv,x,dx);
            X::axpy(Real(1.),dx,x);
        }

        // Return the norm and the residual
        return std::pair <Real,Natural> (norm_rtrue,iter);
    }
}

#endif
