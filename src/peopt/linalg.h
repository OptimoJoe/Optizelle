#ifndef LINALG_H 
#define LINALG_H 

#include <vector>
#include <cmath>
#include <list>

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

    template <typename T>
    void rotg(T a,T b,T& c,T& s);

    template <typename T>
    void rot(int n,const T* x,int incx,T* y,int incy,T c,T s);

    template <typename Real>
    void tpsv(char uplo,char trans,char diag,int n,const Real* Ap,Real* x,
        int incx);
    
    // Indexing function for matrices
    unsigned int ijtok(unsigned int i,unsigned int j,unsigned int m);
    
    // A simple operator specification, A : X->Y
    template <
        typename Real,
        template <typename> class XX,
        template <typename> class YY
    >
    struct Operator {
    private:
        // Create some type shortcuts
        typedef typename XX <Real>::Vector X_Vector;
        typedef typename YY <Real>::Vector Y_Vector;
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
        int i=0;
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
        const unsigned int m,
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
        copy <Real> (m,&(Qt_e1[0]),1,&(y[0]),1);
        tpsv <Real> ('U','N','N',m,&(R[0]),&(y[0]),1);

        // Compute tmp = V y
        X::zero(V_y);
        typename std::list <X_Vector>::const_iterator vv=vs.begin();
        for(int j=0;j<m;j++) {
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
        const unsigned int rst_freq,
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
    class GMRESManipulator {
    public:
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
    std::pair <Real,int> gmres(
        const Operator <Real,XX,XX>& A,
        const typename XX <Real>::Vector& b,
        Real eps,
        unsigned int iter_max,
        unsigned int rst_freq,
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
        rst_freq = rst_freq == 0 ? iter_max : rst_freq;

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
        // Note, this size is restricted to be no larger than the restart frequency
        std::vector <Real> R(rst_freq*(rst_freq+1)/2);

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
        std::vector <Real> Qt_e1(rst_freq+1);

        // Allocoate memory for the Givens rotations
        std::list <std::pair<Real,Real> > Qts;

        // Allocate a temporary work element
        X_Vector A_Mrinv_v; X::init(x,A_Mrinv_v);

        // Allocate memory for the subiteration number of GMRES taking into
        // account restarting
        unsigned int i;

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
        unsigned int iter;
        for(iter = 1; iter <= iter_max;iter++) {

            // Find the current iterate taking into account restarting
            i = iter % rst_freq;

            // We the above remainder is zero, we're on our final iteration
            // before restarting.  However, the iterate in this case is equal to
            // the restart frequency and not zero since our factorization has
            // size rst_freq x rst_freq.
            if(i == 0) i = rst_freq;

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
            int j=1;
            for(typename std::list <std::pair<Real,Real> >::iterator Qt=Qts.begin();
                Qt!=Qts.end();
                Qt++
            ) { 
                rot <Real> (1,&(R[(j-1)+(i-1)*i/2]),1,&(R[j+(i-1)*i/2]),1,
                    Qt->first,Qt->second);
                j++;
            }

            // Form the new Givens rotation
            Qts.push_back(std::pair <Real,Real> ());
            rotg <Real> (R[(i-1)+i*(i-1)/2],norm_w,Qts.back().first,Qts.back().second);

            // Apply this new Givens rotation to the last element of R and norm(w).
            // This fixes our system R.
            rot <Real> (1,&(R[(i-1)+i*(i-1)/2]),1,&(norm_w),1,
                Qts.back().first,Qts.back().second);

            // Apply the new givens rotation to the RHS.  This also determines the new norm
            // of the preconditioned residual.
            rot <Real> (1,&(Qt_e1[i-1]),1,&(Qt_e1[i]),1,
                Qts.back().first,Qts.back().second);
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

        // As long as we didn't just solve for our new ierate, go ahead and solve
        // for it now.
        if(i > 0){ 
            solveInKrylov <Real,XX> (i,&(R[0]),&(Qt_e1[0]),vs,Mr_inv,x,dx);
            X::axpy(Real(1.),dx,x);
        }

        // Return the norm and the residual
        return std::pair <Real,unsigned int> (norm_rtrue,iter);
    }
}

#endif
