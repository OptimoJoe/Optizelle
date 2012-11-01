#ifndef VSPACES_H
#define VSPACES_H
#include <cmath>
#include "linalg.h"
#include "peopt.h"

namespace peopt {

    // Vector space for the nonnegative orthant.  For basic vectors
    // in R^m, use this.
    template <typename Real>
    struct Rm { 
        typedef std::vector <Real> Vector;

        // Memory allocation and size setting
        static void init(const Vector& x, Vector& y) {
            y.resize(x.size());
        }
        
        // y <- x (Shallow.  No memory allocation.)
        static void copy(const Vector& x, Vector& y) {
            peopt::copy <Real> (x.size(),&(x[0]),1,&(y[0]),1);
        }

        // x <- alpha * x
        static void scal(const Real& alpha, Vector& x) {
            peopt::scal <Real> (x.size(),alpha,&(x[0]),1);
        }

        // y <- alpha * x + y
        static void axpy(const Real& alpha, const Vector& x, Vector& y) {
            peopt::axpy <Real> (x.size(),alpha,&(x[0]),1,&(y[0]),1);
        }

        // innr <- <x,y>
        static Real innr(const Vector& x,const Vector& y) {
            return peopt::dot <Real> (x.size(),&(x[0]),1,&(y[0]),1);
        }

        // x <- 0 
        static void zero(Vector& x) {
            #pragma omp parallel for schedule(static)
            for(unsigned int i=0;i<x.size();i++) 
                x[i]=Real(0.);
        }

        // Jordan product, z <- x o y
        static void prod(const Vector& x, const Vector& y, Vector& z) {
            #pragma omp parallel for schedule(static)
            for(unsigned int i=0;i<x.size();i++) 
                z[i]=x[i]*y[i];
        }

        // Identity element, x <- e such that x o e = x
        static void id(Vector& x) {
            #pragma omp parallel for schedule(static)
            for(unsigned int i=0;i<x.size();i++) 
                x[i]=Real(1.);
        }
        
        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
        static void linv(const Vector& x,const Vector& y,Vector& z) {
            #pragma omp parallel for schedule(static)
            for(unsigned int i=0;i<x.size();i++) 
                z[i]=y[i]/x[i];
        }

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e
        static Real barr(const Vector& x) {
            Real z=0;
            #pragma omp parallel for reduction(+:z) schedule(static)
            for(unsigned int i=0;i<x.size();i++)
                z+=log(x[i]);
            return z;
        }

        // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
        // where y > 0.  If the argmax is infinity, then return Real(-1.).
        static Real srch(const Vector& x,const Vector& y) {
            // Line search parameter
            Real alpha=Real(-1.);

            #pragma omp parallel
            {
                // Create a local version of alpha
                Real alpha_loc=Real(-1.);

                // Search for the optimal linesearch parameter
                #pragma omp parallel for schedule(static)
                for(unsigned int i=0;i<x.size();i++) {
                    if(x[i] < 0) {
                        Real alpha0 = -y[i]/x[i];
                        if(alpha_loc==Real(-1.) || alpha0 < alpha_loc)
                            alpha_loc=alpha0;
                    }
                }

                // After we're through with the local search, accumulate the
                // result
                #pragma omp critical
                {
                   if(alpha==Real(-1.) || alpha_loc < alpha)
                       alpha=alpha_loc;
                }
            }
            return alpha;
        }
    };

    // Different cones used in SQL problems
    struct Cone {
        enum t {
            Linear,             // Nonnegative orthant
            Quadratic,          // Second order cone
            Semidefinite        // Cone of positive semidefinite matrices 
        };
    };

    // A vector spaces consisting of a finite product of semidefinite,
    // quadratic, and linear cones.
    template <typename Real>
    struct SQL {
        struct Vector {
            // The overall variable data
            std::vector <Real> data;

            // Offsets of each cone stored in the data
            std::vector <unsigned int> offsets;

            // The type of cones stored in the data
            std::vector <Cone::t> types;
            
            // The size of the cones stored in the data
            std::vector <unsigned int> sizes;

            // The stored Schur decompositions.  Once we have the offset, we
            // store the diagonal and then the matrix.
            mutable std::vector <Real> schur;

            // The offsets of the cached Schur decompositions
            std::vector <unsigned int> schur_offsets;

            // The point where we last took the Schur decomposition
            mutable std::vector <Real> schur_base;

            // The offsets for the bases stored for the Schur decompositions
            std::vector <unsigned int> schur_base_offsets;

            // This creates an empty, unitialized vector
            Vector () {}

            // We require a vector of cone types and their sizes
            Vector (const Messaging& msg,
                const std::vector <Cone::t>& types_,
                const std::vector <unsigned int>& sizes_
            ) : types(types_), sizes(sizes_) {

                // Insure that the type of cones and their sizes lines up
                if(types.size()!=sizes.size())
                    msg.error("The vector containing the type of cones must "
                        "be the same size as the vector with the cone sizes.");

                // Make sure we have at least one cone
                if(types.size() == 0)
                    msg.error("A SQL vector requires at least one cone.");

                // Initialize the offsets.  The last element has the total
                // number of variables
                offsets.resize(sizes.size()+1);
                offsets[0]=0;
                for(unsigned int i=1;i<offsets.size()+1;i++)
                    offsets[i] = types[i-1]==Cone::Linear ||
                                 types[i-1]==Cone::Quadratic
                                     ? offsets[i-1]+sizes[i-1]
                                     : offsets[i-1]+sizes[i-1]*sizes[i-1];

                // Create the data
                data.resize(offsets.back());

                // Calculate offsets for the Schur decompositions.  Basically,
                // the way it works is that we calculate offsets for every cone
                // even though we're not ever going to use them.  In the case
                // we don't have an SDP block, we simply use the last offset.
                // This makes it easy to index to the correct place where the
                // cached information is stored.
                schur_offsets.resize(sizes.size()+1);
                schur_base_offsets.resize(sizes.size()+1);
                for(unsigned int i=1;i<types.size()+1;i++) {
                    schur_offsets[i] =
                        types[i-1]==Cone::Linear ||
                        types[i-1]==Cone::Quadratic
                            ? schur_offsets[i-1]
                            : schur_offsets[i-1]+(sizes[i-1]+1)*sizes[i-1];
                    schur_base_offsets[i] =
                        types[i-1]==Cone::Linear ||
                        types[i-1]==Cone::Quadratic
                            ? schur_base_offsets[i-1]
                            : schur_base_offsets[i-1]+sizes[i-1]*sizes[i-1];
                }

                // Create the memory required for the cached Schur
                // decompositions
                schur.resize(schur_offsets.back());
                schur_base.resize(schur_base_offsets.back());
            }

            // Simple indexing
            Real& operator () (unsigned int i) {
                return data[i-1];
            }
            const Real& operator () (unsigned int i) const {
                return data[i-1];
            }

            // Indexing with multiple cones
            Real& operator () (unsigned int k,unsigned int i) {
                return data[offsets[k-1]+(i-1)];
            }
            const Real& operator () (unsigned int k,unsigned int i) const {
                return data[offsets[k-1]+(i-1)];
            }

            // Indexing a matrix with multiple cones
            Real& operator () (unsigned int k,unsigned int i,unsigned int j) {
                return data[offsets[k-1]+ijtok(i,j,sizes[k-1])];
            }
            const Real& operator ()(unsigned int k,unsigned int i,unsigned int j
            ) const {
                return data[offsets[k-1]+ijtok(i,j,sizes[k-1])];
            }

            // Size of the block
            unsigned int blkSize(unsigned int blk) const {
                return sizes[blk-1];
            }

            // Type of the block
            Cone::t blkType(unsigned int blk) const {
                return types[blk-1];
            }

            // Number of blocks
            unsigned int numBlocks() const {
                return types.size();
            }
        };

        // Gets the Schur decomposition of a block of the SQL vector
        static void get_schur(
            const Vector& X,
            const unsigned int blk,
            std::vector <Real>& V,
            std::vector <Real>& D
        ) {
            // Fix the block index to C indexing
            unsigned int k=blk-1;

            // Get the size of the block
            const int m=X.sizes[k];

            // Next, check if we've already calculated the Schur decomposition.

            // Copy out the the base of the last decomposition
            std::vector <Real> tmp(m*m);
            peopt::copy <Real>
                (m*m,&(X.schur_base[X.schur_base_offsets[k]]),1,&(tmp[0]),1);

            // tmp <- Base_k - X_k
            peopt::axpy <Real>
                (m*m,Real(-1.),&(X.data[X.offsets[k]]),1,&(tmp[0]),1);

            // Find the relative error between the current iterate
            // and the base
            Real norm_xk = sqrt(dot <Real>
                (m,&(X.data[X.offsets[k]]),1,&(X.data[X.offsets[k]]),1));
            Real rel_err = sqrt(dot <Real> (m,&(tmp[0]),1,&(tmp[0]),1))
                / (Real(1e-16)+norm_xk);

            // If the relative error is large, refresh the cached decomposition.
            // To be sure, I'm not really sure what this tolerance should be
            // given that we have both floats and doubles.  In reality,
            // the iterates will probably change rapidly, so I don't think
            // we have to worry too much, but be careful.
            if(rel_err > Real(1e-7)) {

                // Store X_k as the new base
                peopt::copy <Real> (m*m,&(X.data[X.offsets[k]]),1,
                    &(X.schur_base[X.schur_base_offsets[k]]),1);

                // Find the Schur decomposition of X_k 

                // Since the eigenvalue routine destroys the diagonal and
                // lower triangular portion of the matrix, make a copy for
                // computation
                std::vector <Real> Xk(m*m);
                peopt::copy <Real> (m*m,&(X.data[X.offsets[k]]),1,&(Xk[0]),1);

                // Query the size of the workplaces
                std::vector <int> isuppz(2*m);
                std::vector <Real> work(1);
                std::vector <int> iwork(1);
                int lwork=-1;
                int liwork=-1;
                int info;
                int nevals;
                peopt::syevr <Real> ('V','A','U',m,&(Xk[0]),m,Real(0.),
                    Real(0.),0,0,Real(-1.),nevals,
                    &(X.schur[X.schur_offsets[k]]),
                    &(X.schur[X.schur_offsets[k]+m]),m,&(isuppz[0]),
                    &(work[0]),lwork,&(iwork[0]),liwork,info);

                // Resize the workplaces
                lwork = int(work[0])+1;
                work.resize(lwork);
                liwork = iwork[0];
                iwork.resize(liwork);

                // Find the decomposition
                peopt::syevr <Real> ('V','A','U',m,&(Xk[0]),m,Real(0.),
                    Real(0.),0,0,Real(-1.),nevals,
                    &(X.schur[X.schur_offsets[k]]),
                    &(X.schur[X.schur_offsets[k]+m]),m,&(isuppz[0]),
                    &(work[0]),lwork,&(iwork[0]),liwork,info);
            }

            // Copy out the decomposition from the cached copy
            D.resize(m);
            peopt::copy <Real> (m,&(X.schur[X.schur_offsets[k]]),1,&(D[0]),1);
            V.resize(m*m);
            peopt::copy <Real>
                (m*m,&(X.schur[X.schur_offsets[k]+m]),1,&(V[0]),1);
        }
        
        // Memory allocation and size setting
        static void init(const Vector& x, Vector& y) {
            y.data.resize(x.data.size());
            y.offsets=x.offsets;
            y.types=x.types;
            y.sizes=x.sizes;
            y.schur.resize(x.schur.size());
            y.schur_offsets=x.schur_offsets;
            y.schur_base.resize(x.schur_base.size());
            y.schur_base_offsets=x.schur_base_offsets;
        }
        
        // y <- x (Shallow.  No memory allocation.)
        static void copy(const Vector& x, Vector& y) {
            peopt::copy <Real> (x.data.size(),&(x.data[0]),1,&(y.data[0]),1);
        }

        // x <- alpha * x
        static void scal(const Real& alpha, Vector& x) {
            peopt::scal <Real> (x.data.size(),alpha,&(x.data[0]),1);
        }

        // y <- alpha * x + y
        static void axpy(const Real& alpha, const Vector& x, Vector& y) {
            peopt::axpy <Real> (x.data.size(),alpha,&(x.data[0]),1,
                &(y.data[0]),1);
        }

        // innr <- <x,y>
        static Real innr(const Vector& x,const Vector& y) {
            return peopt::dot <Real> (x.data.size(),&(x.data[0]),1,
                &(y.data[0]),1);
        }

        // x <- 0 
        static void zero(Vector& x) {
            #pragma omp parallel for schedule(static)
            for(unsigned int i=0;i<x.data.size();i++) 
                x.data[i]=Real(0.);
        }

        // Jordan product, z <- x o y
        static void prod(const Vector& x, const Vector& y, Vector& z) {
            /* It's hard to tell apriori how to parallelize this
               computuation.  Sometimes, it helps to parallelize across
               the cones, but if the cones are large and few, it helps
               to parallelize the computation.  The hardest case to determine
               is if the cones are radically different in size.  At this point
               we use a simple strategy.  Basically, as far as I can tell,
               ATLAS BLAS won't allow us to change the number of threads
               being used by its routines.  Since I don't really want to
               recode these routines, we're stuck making sure that BLAS
               controls the parallelism, which means doing parallel
               computation on each cone one after another. 
            */
            // Loop over all the blocks
            for(unsigned int blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                unsigned int m=x.blkSize(blk);

                // Depending on the block, compute a different jordan product
                switch(x.blkType(blk)) {

                // z = diag(x) y 
                case Cone::Linear:
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=1;i<=m;i++)
                        z(blk,i)=x(blk,i)*y(blk,i);
                    break;

                // z = [x'y ; x0 ybar + y0 xbar]
                case Cone::Quadratic:
                    z(blk,1)=peopt::dot <Real> (m,&(x(blk,1)),1,&(y(blk,1)),1);
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=2;i<=m;i++)
                        z(blk,i)=x(blk,1)*y(blk,i)+x(blk,i)*y(blk,1);
                    break;

                // z = (xy + yx ) /2
                case Cone::Semidefinite:
                    peopt::syr2k <Real> ('U','N',int(m),int(m),Real(0.5),
                        &(x(blk,1,1)),int(m),&(y(blk,1,1)),int(m),
                        Real(0.),&(z(blk,1,1)),int(m));

                    // Fill in the bottom half of the matrix
                    #pragma omp parallel for schedule(guided)
                    for(unsigned int j=1;j<=m;j++)
                        for(unsigned int i=j+1;i<=m;i++)
                            z(blk,i,j)=z(blk,j,i);

                    break;
                }
            }
        }

        // Identity element, x <- e such that x o e = x
        static void id(Vector& x) {
            // Loop over all the blocks
            for(unsigned int blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                unsigned int m=x.blkSize(blk);

                // Depending on the block, compute a different identity element
                switch(x.blkType(blk)) {

                // x = vector of all 1s
                case Cone::Linear:
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=1;i<=m;i++) 
                        x(blk,i)=Real(1.);
                    break;
                // x = (1,0,...,0)
                case Cone::Quadratic:
                    x(blk,1)=Real(1.);
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=2;i<=m;i++)
                        x(blk,i)=Real(0.);
                    break;
                // x = I
                case Cone::Semidefinite:
                    // We write the diagonal elements twice to avoid the
                    // conditional.
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=1;i<=m;i++) 
                        for(unsigned int j=1;j<=m;j++) 
                            x(blk,i,j)=Real(0.);

                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=1;i<=m;i++) 
                        x(blk,i,i)=Real(1.);
                    break;
                }
            }
        }

        // This applies the inverse of the Schur complement of the Arw
        // operator to a vector.  Note, y has size one less than x.
        // Hence, m is the length of y and m+1 is the length of x.
        static void invSchur(const int m,const Real* x,Real* y) {

            // innr_xbar_y <- <xbar,y>
            Real innr_xbar_y = dot <Real> (m,&(x[1]),1,&(y[0]),1);

            // denom <- x0 ( x0^2 - <xbar,xbar> )
            Real denom = x[0]*( x[0]*x[0] - dot <Real> (m,&(x[1]),1,&(x[1]),1));

            // y <- 1/x0 y
            peopt::scal <Real> (m,Real(1.)/x[0],&(y[0]),1);

            // y <- 1/x0 y + <xbar,y> / (x0 ( x0^2 - <xbar,xbar> )) xbar
            peopt::axpy <Real> (m,innr_xbar_y/denom,&(x[1]),1,&(y[0]),1);
        }
        
        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
        static void linv(const Vector& x,const Vector& y,Vector& z) {
            // We have these two vectors in case we have a SDP block
            std::vector <double> D;
            std::vector <double> V;

            // Loop over all the blocks
            for(unsigned int blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                unsigned int m=x.blkSize(blk);

                // Depending on the block, compute a different operator
                switch(x.blkType(blk)) {

                // z = inv(Diag(x)) y
                case Cone::Linear:
                    #pragma omp parallel for schedule(static)
                    for(unsigned int i=1;i<=m;i++) 
                        z(blk,i)=y(blk,i)/x(blk,i);
                    break;

                // z = inv(Arw(x)) y
                case Cone::Quadratic: {
                    // invSchur_ybar <- invSchur(x)(y_bar) 
                    std::vector <Real> invSchur_ybar(m-1);
                    peopt::copy <Real>(m-1,&(y(blk,2)),1,&(invSchur_ybar[0]),1);
                    invSchur(m-1,&(x(blk,1)),&(invSchur_ybar[0]));

                    // a <- 1 / (x0 - (1/x0) <x_bar,x_bar>) * y0
                    Real a = y(blk,1) / (x(blk,1) - (Real(1.)/x(blk,1)) *
                        peopt::dot <Real> (m-1,&(x(blk,2)),1,&(x(blk,2)),1));

                    // b <- - (1/x0) <xbar,invSchur(x)(y_bar)>
                    Real b = Real(-1.)/x(blk,1) * peopt::dot <Real>
                        (m-1,&(x(blk,2)),1,&(invSchur_ybar[0]),1);
                    
                    // z0 <- 1 / (x0 - (1/x0) <x_bar,x_bar>) y0
                    //       - (1/x0) <x_bar,invSchur(x)(y_bar)> 
                    z(blk,1) = a + b;
                    
                    // z_bar <- invSchur(x)(x_bar)
                    peopt::copy <Real> (m-1,&(x(blk,2)),1,&(z(blk,2)),1);
                    invSchur(m-1,&(x(blk,1)),&(z(blk,2)));

                    // zbar <- (-y0/x0) invSchur(x)(x_bar)
                    peopt::scal <Real> (m-1,-y(blk,1)/x(blk,1),&(z(blk,2)),1);

                    // z_bar <- (-y0/x0) invSchur(x)(x_bar) + invSchur(x)(y_bar)
                    peopt::axpy <Real> (m-1,Real(1.),&(invSchur_ybar[0]),1,
                        &(z(blk,2)),1);
                    break;

                // (xz + zx)/2 = y.  This is an implicit definition of z.
                // Basically, we solve the Sylvester equation and then
                // multiply the result by 2.
                } case Cone::Semidefinite: {
                    // Get the Schur complement of the block.  With any luck
                    // these are cached.
                    peopt::SQL <double>::get_schur(x,blk,V,D);

                    // Solve the Sylvester equation
                    sylvester(m,&(V[0]),&(D[0]),&(y(blk,1,1)),&(z(blk,1,1)));

                    // Scale the result by 2
                    peopt::scal <Real> (m*m,Real(2.),&(z(blk,1,1)),1);
                    break;
                }}
            }
        }

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e
        static Real barr(const Vector& x) {
            // This accumulates the barrier's value
            Real z=0;

            // Loop over all the blocks
            for(unsigned int blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                unsigned int m=x.blkSize(blk);

                // In case we need to take a Choleski factorization for the
                // SDP blocks
                std::vector <Real> U;

                // Depending on the block, compute a different barrier
                switch(x.blkType(blk)) {

                // z += sum_i log(x_i)
                case Cone::Linear:
                    #pragma omp parallel for reduction(+:z) schedule(static)
                    for(unsigned int i=1;i<=m;i++)
                        z+=log(x(blk,i));
                    break;

                // z += 0.5 * log(x0^2-<xbar,xbar))
                case Cone::Quadratic:
                    z+=Real(0.5) * log(x(blk,1)*x(blk,1)
                        -dot <Real> (m-1,&(x(blk,2)),1,&(x(blk,2)),1));
                    break;

                // z += log(det(x)).  We compute this by noting that
                // log(det(x)) = log(det(u'u)) = log(det(u')det(u))
                //             = log(det(u)^2) = 2 log(det(u))
                case Cone::Semidefinite: {

                    // Find the Choleski factorization of X
                    U.resize(m*m);
                    int info;
                    peopt::copy <Real> (m*m,&(x(blk,1,1)),1,&(U[0]),1);
                    peopt::potrf <Real> ('U',m,&(U[0]),m,info);

                    // Find the deterimant of the Choleski factor
                    Real det(1.);
                    #pragma omp parallel for reduction(*:det) schedule(static)
                    for(unsigned int i=1;i<=m;i++)
                        det*=U[peopt::ijtok(i,i,m)];
                    
                    // Complete the barrier computation by taking the log
                    z+= Real(2.) * log(det);
                    
                    break;
                } }
            }

            // Return the accumulated barrier value
            return z;
        }

        // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
        // where y > 0.  If the argmax is infinity, then return Real(-1.).
        static Real srch(const Vector& x,const Vector& y) {
            // Line search parameter
            Real alpha=Real(-1.);
                   
            // Variables required for the linesearch on SDP blocks 
            std::vector <Real> invU;
            std::vector <Real> tmp;
            std::vector <Real> invUtXinvU;

            // Loop over all the blocks
            for(unsigned int blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                unsigned int m=x.blkSize(blk);

                // Depending on the block, compute a different barrier
                switch(x.blkType(blk)) {

                // Pointwise, alpha_i = -y_i / x_i.  If this number is positive,
                // then we need to restrict how far we travel.
                case Cone::Linear:

                    #pragma omp parallel
                    {
                        // Create a local version of alpha
                        Real alpha_loc=Real(-1.);

                        // Search for the optimal linesearch parameter
                        #pragma omp parallel for schedule(static)
                        for(unsigned int i=1;i<=m;i++) {
                            if(x(blk,i) < 0) {
                                Real alpha0 = -y(blk,i)/x(blk,i);
                                if(alpha_loc==Real(-1.) || alpha0 < alpha_loc)
                                    alpha_loc=alpha0;
                            }
                        }

                        // After we're through with the local search,
                        // accumulate the result
                        #pragma omp critical
                        {
                           if(alpha==Real(-1.) || alpha_loc < alpha)
                               alpha=alpha_loc;
                        }
                    }

                    break;

                // We choose the smallest positive number between:
                // -y0/x0, and the roots of alpha^2 a + alpha b + c
                //
                // where
                //
                // a = x0^2 - ||xbar||^2
                // b = 2x0y0 - 2 <xbar,ybar>
                // c = y0^2 - ||ybar||^2 
                // 
                // Technically, if a is zero, the quadratic formula doesn't
                // apply and we use -c/b instead of the roots.  If b is zero
                // and a is zero, then there's no limit to the line search
                case Cone::Quadratic: {

                    // Now, first we have to insure that the leading coefficient
                    // of the second order cone problem remains nonnegative.
                    // This number tells us how far we can step before this
                    // is not true.
                    Real alpha0 = -y(blk,1)/x(blk,1);

                    // Next, assuming that the leading coefficient is fine,
                    // figure out how far we can step before we violate the
                    // rest of the SOCP constraint.  This involves solving
                    // the quadratic equation from above.
                    Real a = x(blk,1)*x(blk,1)
                        - dot <Real> (m-1,&(x(blk,2)),1,&(x(blk,2)),1);
                    Real b = Real(2.)*(x(blk,1)*y(blk,1)
                        - dot <Real> (m-1,&(x(blk,2)),1,&(y(blk,2)),1));
                    Real c = y(blk,1)*y(blk,1)
                        - dot <Real> (m-1,&(y(blk,2)),1,&(y(blk,2)),1);
                    unsigned int nroots(0);
                    Real alpha1(-1.);
                    Real alpha2(-1.);
                    quad_equation(a,b,c,nroots,alpha1,alpha2);

                    // Now, determine the step length.

                    // If we have a restriction on alpha based on the leading
                    // coefficient.
                    if(alpha0 >= Real(0.) && (alpha==Real(-1.) || alpha0<alpha))
                        alpha=alpha0;

                    // Next, if we have two roots, determine the restriction 
                    if(nroots==2) { 
                        if(alpha1>=Real(0.) && (alpha==Real(-1.)||alpha1<alpha))
                            alpha=alpha1;
                        if(alpha2>=Real(0.) && (alpha==Real(-1.)||alpha2<alpha))
                            alpha=alpha2;

                    // If we have a single root 
                    } else if(nroots==1) {
                        if(alpha1>=Real(0.) && (alpha==Real(-1.)||alpha1<alpha))
                            alpha=alpha1;
                    }

                    // If we no roots, there's no additional restriction.
                    // This can't happen since we assume that y is strictly
                    // feasible.
                }
                break;

                // Here, if the solution to the generalized eigenvalue problem
                //
                // Y v = lambda Y v
                //
                // has a negative eigenvalue, we use the negative of that.  Now,
                // sometimes we have cached a Schur decomposition of Y.  In this
                // case we can look for the smallest eigenvalue of V' X V.  
                // However, since we're doing a line search both for X and for 
                // Z, I don't think we're always going to have a cached version.
                // That means that we may be required to do a Schur 
                // decomposition for the line search on the dual variable 
                // which is expensive.  Alternatively, since Y is positive 
                // definite, we can do a Choleski factorization.  Then, we can 
                // find the smallest eigenvalue of inv(U)' X inv(U).  From
                // this, we have that the line search parameter is -1/lambda_min
                // whenever lambda_min <0.  That's what we do here.
                case Cone::Semidefinite: {

                    // First, make a copy of Y and store it in invU
                    invU.resize(m*m);
                    peopt::copy <Real> (m*m,&(y(blk,1,1)),1,&(invU[0]),1);

                    // invU <- chol(Y)
                    int info;
                    peopt::potrf <Real> ('U',m,&(invU[0]),m,info);

                    // invU <- inv(chol(Y))
                    peopt::trtri <Real> ('U','N',m,&(invU[0]),m,info);
                    
                    // Fill in the bottom half of invU with zeros 
                    #pragma omp parallel for schedule(guided)
                    for(unsigned int j=1;j<=m;j++)
                        for(unsigned int i=j+1;i<=m;i++)
                            invU[ijtok(i,j,m)]=Real(0.);

                    // Find inv(chol(Y))' X inv(chol(Y))
                    tmp.resize(m*m);
                    invUtXinvU.resize(m*m);
        
                    // tmp <- X inv(chol(Y)) 
                    peopt::gemm <Real> ('N','N',m,m,m,Real(1.),&(x(blk,1,1)),m,
                        &(invU[0]),m,Real(0.),&(tmp[0]),m); 

                    // invUtXinvU <- inv(chol(Y))' X inv(chol(Y))
                    peopt::gemm <Real> ('T','N',m,m,m,Real(1.),&(invU[0]),m,
                        &(tmp[0]),m,Real(0.),&(invUtXinvU[0]),m);

                    // Find the smallest eigenvalue of invUtXinvU.  The
                    // negative of this is the line search parameter.
                    Real alpha0=-Real(1.)/peopt::lanczos<Real>
                        (m,&(invUtXinvU[0]),m,Real(1e-1));

                    // Adjust the line search step if necessary
                    if(alpha0 >= Real(0.) && (alpha==Real(-1.) || alpha0<alpha))
                        alpha=alpha0;
                } }
            }
            return alpha;
        }
    };
}
#endif
