#ifndef LINALG_H 
#define LINALG_H 

#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <utility>
#include <string>

namespace peopt {
    template <typename T>
    void copy(int n,const T* x,int incx,T* y,int incy);
    
    template <typename T>
    void axpy(int n,T alpha,const T* x,int incx,T* y,int incy);
    
    template <typename T>
    void scal(int n,const T alpha,T* x,int incx);
    
    template <typename T>
    T dot(int n,const T* x,int incx,const T* y,int incy);
    
    template <typename T>
    void syr2k(char uplo,char trans,int n,int k,T alpha,const T* A,int lda,
        const T* B,int ldb,T beta,T* C,int ldc);

    template <typename T>
    void syevr(char jobz,char range,char uplo,int n,T *A,int lda,
        T vl,T vu,int il,int iu,T abstol,int& m,
        T* w,T* z,int ldz,int* isuppz,T* work,int lwork,
        int* iwork,int liwork,int& info);
    
    template <typename T>
    void stemr(char jobz,char range,int n,T *D,T *E,T vl,T vu,int il,int iu,
        int& m,T* w,T* z,int ldz,int nzc,int* isuppz,int& tryrac,T* work,
        int lwork,int* iwork,int liwork,int& info);
    
    template <typename T>
    void stevr(char jobz,char range,int n,T *D,T *E,T vl,T vu,int il,int iu,
        T abstol, int& m,T* w,T* z,int ldz,int* isuppz,T* work,
        int lwork,int* iwork,int liwork,int& info);

    template <typename T>
    T lamch(char cmach);

    template <typename T>
    void gemm(char transa,char transb,int m,int n,int k,T alpha,
        const T* A,int lda,const T* B,int ldb,T beta,
        T* C,int ldc);
    
    template <typename T>
    void symm(char side,char uplo,int m,int n,T alpha,const T* A,
        int lda,const T* B,int ldb,T beta,T* C,int ldc);
    
    template <typename T>
    void symv(char uplo,int n,T alpha,const T* A,int lda,const T* x,int incx,
        T beta,T* y,int incy);

    template <typename T>
    void potrf(char uplo,int n,T* A,int lda,int& info);

    template <typename T>
    void trtri(char uplo,char diag,int n,T* A,int lda,int& info);

    // Indexing function for matrices
    unsigned int ijtok(unsigned int i,unsigned int j,unsigned int m);

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
    template <typename T>
    void sylvester(
        const unsigned int m,
        const T* V,
        const T* D,
        const T* B,
        T* X
    ) {

        // Find V' B V
        std::vector <T> tmp(m*m);
        std::vector <T> VtBV(m*m);
        // tmp <- B V
        symm <T> ('L','U',m,m,T(1.),&(B[0]),m,&(V[0]),m,T(0.),&(tmp[0]),m); 
        // VtBV <- V' B V
        gemm <T> ('T','N',m,m,m,T(1.),&(V[0]),m,&(tmp[0]),m,T(0.),&(VtBV[0]),m);

        // Solve for each column of X.  In theory, we only need half of these
        // elements since X is symmetric.
        #pragma omp parallel for schedule(static)
        for(unsigned int j=1;j<=m;j++) {
            for(unsigned int i=1;i<=j;i++) 
                X[ijtok(i,j,m)]=VtBV[ijtok(i,j,m)]/(D[i-1]+D[j-1]);
        }

        // Transform the solution back, X = V X V'
        // tmp <- V X
        symm <T> ('R','U',m,m,T(1.),&(X[0]),m,&(V[0]),m,T(0.),&(tmp[0]),m);
        // X <- V X V'
        gemm <T> ('N','T',m,m,m,T(1.),&(tmp[0]),m,&(V[0]),m,T(0.),&(X[0]),m);
    }

    // Find a bound on the smallest eigenvalue of the given matrix A such
    // that lambda_min(A) < alpha where alpha is returned from this function.
    template <typename T>
    T lanczos(
        const unsigned int m,
        const T* A,
        const unsigned int max_iter,
        const T tol
    ) {
        // Create the initial Krylov vector
        std::vector <T> v(m,T(1./std::sqrt(m)));

        // Get the next Krylov vector and orthgonalize it
        std::vector <T> w(m);
        // w <- A v
        symv <T> ('U',m,T(1.),&(A[0]),m,&(v[0]),1,T(0.),&(w[0]),1);
        // alpha[0] <- <Av,v>
        std::vector <T> alpha;
        alpha.push_back(dot <T> (m,&(w[0]),1,&(v[0]),1));
        // w <- Av - <Av,v> v
        axpy <T> (m,-alpha[0],&(v[0]),1,&(w[0]),1);

        // Store the norm of the Arnoldi vector w in the off diagonal part of T.
        // By T, we mean the triagonal matrix such that A = Q T Q'.
        std::vector <T> beta;
        beta.push_back(std::sqrt(dot <T> (m,&(w[0]),1,&(w[0]),1)));

        // Allocate memory for solving an eigenvalue problem for the Ritz
        // values and vectors later.
        std::vector <int> isuppz;
        std::vector <T> work(1);
        std::vector <int> iwork(1);
        int lwork=-1;
        int liwork=-1;
        int info;
        int nevals;
        //int nzc=0;
        std::vector <T> W;
        std::vector <T> Z;
        std::vector <T> D;
        std::vector <T> E;

        // Start Lanczos
        std::vector <T> v_old(m);
        for(unsigned int i=0;i<max_iter;i++) {
            // Save the current Arnoldi vector
            copy <T> (m,&(v[0]),1,&(v_old[0]),1);

            // Copy the candidate Arnoldi vector to the current Arnoldi vector
            copy <T> (m,&(w[0]),1,&(v[0]),1);

            // Get the normalized version of this vector.  This is now a real
            // Arnoldi vector.
            scal <T> (m,T(1.)/beta[i],&(v[0]),1);

            // Get the new Arnoldi vector, w <- A v
            symv <T> ('U',m,T(1.),&(A[0]),m,&(v[0]),1,T(0.),&(w[0]),1);

            // Orthogonalize against v_old and v using modified Gram-Schdmit.

            // First, we orthogonalize against v_old
            // w <- Av - <Av,v_old> v_old.  Due to symmetry, <Av,v_old>=beta.
            axpy <T> (m,-beta[i],&(v_old[0]),1,&(w[0]),1);

            // Now, we orthogonalize against v
            // Find the Gram-Schmidt coefficient
            alpha.push_back(dot <T> (m,&(w[0]),1,&(v[0]),1));
            // Orthogonlize w to v
            axpy <T> (m,-alpha[i+1],&(v[0]),1,&(w[0]),1);

            // Store the norm of the Arnoldi vector w in the off diagonal part
            // of T.
            beta.push_back(std::sqrt(dot <T> (m,&(w[0]),1,&(w[0]),1)));
   
   #if 0
            // Figure out the workspaces for the eigenvalues and eigenvectors
            int k=alpha.size();  // Size of the eigenvalue subproblem
            D.resize(alpha.size());
            copy <T> (k,&(alpha[0]),1,&(D[0]),1);
            E.resize(beta.size());
            copy <T> (k,&(beta[0]),1,&(E[0]),1);
            isuppz.resize(2*k);
            lwork=-1;
            liwork=-1;
            W.resize(k);
            Z.resize(k*k);
            peopt::stemr <double> ('V','A',k,&(D[0]),&(E[0]),double(0.),
                double(0.),0,0,nevals,&(W[0]),&(Z[0]),k,k,&(isuppz[0]),
                nzc,&(work[0]),lwork,&(iwork[0]),liwork,info);

            // Resize the workspace 
            lwork = int(work[0])+1;
            work.resize(lwork);
            liwork = iwork[0];
            iwork.resize(liwork);

            // Find the eigenvalues and vectors 
            peopt::stemr <double> ('V','A',k,&(D[0]),&(E[0]),double(0.),
                double(0.),0,0,nevals,&(W[0]),&(Z[0]),k,k,&(isuppz[0]),
                nzc,&(work[0]),lwork,&(iwork[0]),liwork,info);
#else
            // Figure out the workspaces for the eigenvalues and eigenvectors
            int k=alpha.size();  // Size of the eigenvalue subproblem
            D.resize(alpha.size());
            copy <T> (k,&(alpha[0]),1,&(D[0]),1);
            E.resize(beta.size());
            copy <T> (k,&(beta[0]),1,&(E[0]),1);
            isuppz.resize(2*k);
            lwork=-1;
            liwork=-1;
            W.resize(k);
            Z.resize(k*k);
            peopt::stevr <T> ('V','A',k,&(D[0]),&(E[0]),T(0.),
                T(0.),0,0,peopt::lamch <T> ('S'),nevals,&(W[0]),&(Z[0]),k,
                &(isuppz[0]),&(work[0]),lwork,&(iwork[0]),liwork,info);

            // Resize the workspace 
            lwork = int(work[0])+1;
            work.resize(lwork);
            liwork = iwork[0];
            iwork.resize(liwork);

            // Find the eigenvalues and vectors 
            peopt::stevr <T> ('V','A',k,&(D[0]),&(E[0]),T(0.),
                T(0.),0,0,peopt::lamch <T> ('S'),nevals,&(W[0]),&(Z[0]),k,
                &(isuppz[0]),&(work[0]),lwork,&(iwork[0]),liwork,info);
#endif

            // Find beta_i |s_{ik}| where s_{ik} is the ith (last) element
            // of the kth Ritz vector where k corresponds to the largest
            // and smallest Ritz values.  Basically, we don't know which is
            // going to converge first, but they'll both be the first two.
            // Hence, we converge until these errors estimates are accurate
            // enough.
            T err_est_min = fabs(Z[ijtok(k,1,k)])*beta[i+1];
            T err_est_max = fabs(Z[ijtok(k,k,k)])*beta[i+1];

            // Stop of the error estimates are small
            if(    err_est_min < fabs(W[0]) * tol
                && err_est_max < fabs(W[i]) * tol
            )
                break;
        }

        // Return the smallest Ritz value
        return W[0];
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
        int i=0;
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
    // || x || <= delta.  The parameters are as follows.
    // 
    // (input) A : Operator in the system A B x = b
    // (input) b : Right hand side in the system A B x = b
    // (input) B : Projection in the system A B x = b
    // (input) eps : Stopping tolerance
    // (input) iter_max :  Maximum number of iterations
    // (input) orthog_max : Maximum number of orthgonalizations.  If this
    //     number is 1, then we do the conjugate gradient algorithm.
    // (input) delta : Trust region radius.  If this number is infinity, we
    //     do not scale the final step if we detect negative curvature.
    // (input/output) x : Initial guess and final solution x.  We assume that
    //     the initial guess satisfies norm(x) <= delta.
    // (output) x_cp : The Cauchy-Point, which is defined as the solution x
    //     after a single iteration.
    // (output) norm_r : The norm of the final residual.
    // (output) iter : The number of iterations required to converge. 
    // (output) krylov_stop : The reason why the Krylov method was terminated.
    template <
        typename Real,
        template <typename> class XX
    >
    std::pair <Real,int> truncated_pcd(
        const Operator <Real,XX,XX>& A,
        const typename XX <Real>::Vector& b,
        const Operator <Real,XX,XX>& B,
        const Real eps,
        const unsigned int iter_max,
        const unsigned int orthog_max,
        const Real delta,
        typename XX <Real>::Vector& x,
        typename XX <Real>::Vector& x_cp,
        Real& norm_r,
        unsigned int& iter,
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

        // Find the norm of x
        Real norm_x = sqrt(X::innr(x,x));

        // Allocate memory for and find the initial residual, A*x-b
        X_Vector r; X::init(x,r);
        if(norm_x!=Real(0.)) {
            A(x,r);
            X::axpy(Real(-1.),b,r);

        // Frequently, we start with an initial guess of zero.  In this case,
        // let us save some work.
        } else {
            X::copy(b,r);
            X::scal(Real(-1.),r);
        }

        // Find the norm of the residual and save the original residual
        norm_r = sqrt(X::innr(r,r));
        Real norm_r0 = norm_r;

        // Loop until the maximum iteration
        for(iter=1;iter<=iter_max;iter++){
        
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
            if(X::innr(Bp,r) > 0) {
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

            // Find the norm at the trial step
            Real innr_Bp_x = X::innr(Bp,x);
            Real innr_Bp_Bp = X::innr(Bp,Bp);
            Real norm_xpaBp = sqrt(norm_x*norm_x + Real(2.)*alpha*innr_Bp_x
                + alpha*alpha*innr_Bp_Bp);

            // If we have negative curvature or our trial point is outside the
            // trust-region radius, terminate truncated-PCD and find our final
            // step.  We have the <Bp,ABp> != <Bp,ABp> check in order to trap
            // NaNs.
            if( innr_Bp_ABp <= Real(0.) ||
                norm_xpaBp >= delta ||
                innr_Bp_ABp != innr_Bp_ABp 
            ) {
                // If we're paying attention to the trust-region, scale the
                // step appropriately.
                if(delta < std::numeric_limits <Real>::infinity()) {
                    // Find sigma so that || x + sigma p || = delta.  This can
                    // be found by finding the positive root of the quadratic
                    // || x + sigma Bp ||^2 = delta^2.  Specifically, we want
                    // the positive root of
                    // sigma^2<Bp,Bp> + sigma(2 <Bp,x>) + (<x,x>-delta^2).
                    Real aa = innr_Bp_Bp;
                    Real bb = Real(2.)*innr_Bp_x; 
                    Real cc = norm_x*norm_x-delta*delta;
                    Real sigma = (-bb + sqrt(bb*bb-Real(4.)*aa*cc))/(2*aa);

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
                if(innr_Bp_ABp<=0 || innr_Bp_ABp != innr_Bp_ABp)
                    krylov_stop = KrylovStop::NegativeCurvature;
                else
                    krylov_stop = KrylovStop::TrustRegionViolated;
 

                // If this is the first iteration, save the Cauchy-Point
                if(iter==1) X::copy(x,x_cp);
                break;
            }

            // Take a step in this direction
            X::axpy(alpha,Bp,x);

            // If this is the first iteration, save the Cauchy-Point
            if(iter==1) X::copy(x,x_cp);

            // Update the norm of x
            norm_x = norm_xpaBp;

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

        // Return the norm of the residual and the iteration number we
        // completed at
        return std::pair <Real,unsigned int> (norm_r,iter);
    }
}

#endif
