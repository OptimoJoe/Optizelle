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

#ifndef VSPACES_H
#define VSPACES_H
#include <cmath>
#include <random>
#include "optizelle/linalg.h"
#include "optizelle/optizelle.h"
#include "optizelle/json.h"

//---Optizelle0---
namespace Optizelle {
//---Optizelle1---

    using namespace Optizelle;

    //---Rm0---
    // Vector space for the nonnegative orthant.  For basic vectors
    // in R^m, use this.
    template <typename Real>
    struct Rm { 
        // Disallow constructors
        NO_CONSTRUCTORS(Rm);

        // Use std::vector as our vector storage
        typedef std::vector <Real> Vector;

        // Memory allocation and size setting.
        static Vector init(Vector const & x) {
            return std::move(Vector(x.size()));
        }
        
        // y <- x (Shallow.  No memory allocation.)
        static void copy(Vector const & x, Vector & y) {
            Optizelle::copy <Real> (x.size(),&(x.front()),1,&(y.front()),1);
        }

        // x <- alpha * x.
        static void scal(Real const & alpha, Vector & x) {
            Optizelle::scal <Real> (x.size(),alpha,&(x.front()),1);
        }

        // y <- alpha * x + y.
        static void axpy(Real const & alpha, Vector const & x, Vector & y) {
            Optizelle::axpy<Real>(x.size(),alpha,&(x.front()),1,&(y.front()),1);
        }

        // innr <- <x,y>.
        static Real innr(Vector const & x,Vector const & y) {
            return Optizelle::dot<Real>(x.size(),&(x.front()),1,&(y.front()),1);
        }

        // x <- 0.
        static void zero(Vector & x) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for(Natural i=0;i<x.size();i++) 
                x[i]=Real(0.);
        }

        // x <- random
        static void rand(Vector & x){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<Real> dis(Real(0.),Real(1.));
            
            // This is not parallel since it doesn't appear that our generator
            // works properly when parallel.
            for(Natural i=0;i<x.size();i++) 
                x[i]=Real(dis(gen));
        }

        // Jordan product, z <- x o y.
        static void prod(Vector const & x, Vector const & y, Vector & z) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for(Natural i=0;i<x.size();i++) 
                z[i]=x[i]*y[i];
        }

        // Identity element, x <- e such that x o e = x.
        static void id(Vector & x) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for(Natural i=0;i<x.size();i++) 
                x[i]=Real(1.);
        }
        
        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y.
        static void linv(Vector const & x,Vector const & y,Vector & z) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for(Natural i=0;i<x.size();i++) 
                z[i]=y[i]/x[i];
        }

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e.
        static Real barr(Vector const & x) {
            Real z=Real(0.);
            #ifdef _OPENMP
            #pragma omp parallel for reduction(+:z) schedule(static)
            #endif
            for(Natural i=0;i<x.size();i++)
                z+=log(x[i]);
            return z;
        }

        // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
        // where y > 0. 
        static Real srch(Vector const & x,Vector const & y) {
            // Line search parameter
            Real alpha=std::numeric_limits <Real>::infinity();

            #ifdef _OPENMP
            #pragma omp parallel
            #endif
            {
                // Create a local version of alpha.
                Real alpha_loc=std::numeric_limits <Real>::infinity();

                // Search for the optimal linesearch parameter.
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static)
                #endif
                for(Natural i=0;i<x.size();i++) {
                    if(x[i] < Real(0.)) {
                        Real alpha0 = -y[i]/x[i];
                        alpha_loc = alpha0 < alpha_loc ? alpha0 : alpha_loc;
                    }
                }

                // After we're through with the local search, accumulate the
                // result
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    alpha = alpha_loc < alpha ? alpha_loc : alpha;
                }
            }
            return alpha;
        }

        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.
        static void symm(Vector & x) { }
    };
    //---Rm1---

    namespace json {
        // Serialization utility for the Rm vector space
        template <typename Real>
        struct Serialization <Real,Rm> {
            static std::string serialize (
                typename Rm <Real>::Vector const & x
            ) {
                // Create a jsoncpp object to copy into
                Json::Value x_json;  

                // Copy the information
                for(Natural i=0;i<x.size();i++)
                    x_json[Json::ArrayIndex(i)]=x[i];

                // Return a string of the result
                Json::StyledWriter writer;
                return writer.write(x_json);
            }
            static typename Rm <Real>::Vector deserialize (
                typename Rm <Real>::Vector const & x_,
                std::string const & x_json_
            ) {
                // Create a json tree from the input string
                Json::Value x_json;
                Json::Reader reader;
                reader.parse(x_json_,x_json,true);

                // Create a vector from the json tree 
                std::vector <Real> x(x_json.size());
                for(Natural i=0;i<x.size();i++)
                    x[i]=Real(x_json[Json::ArrayIndex(i)].asDouble());
                return std::move(x);
            }
        };
    }

    // Different cones used in SQL problems
    namespace Cone {
        enum t {
            Linear,             // Nonnegative orthant
            Quadratic,          // Second order cone
            Semidefinite        // Cone of positive semidefinite matrices 
        };

        // Converts the cone to a string
        std::string to_string(t const & cone);
        
        // Converts a string to a cone 
        t from_string(std::string const & cone);

        // Checks whether or not a string is valid
        bool is_valid(std::string const & name); 
    }

    // A vector spaces consisting of a finite product of semidefinite,
    // quadratic, and linear cones.  This uses the nonsymmetric product
    // for the SDP blocks where x o y = xy.  This is not a true Euclidean-Jordan
    // algebra, but is sufficient for our purposes.
    template <typename Real>
    struct SQL {
        // Disallow constructors
        NO_CONSTRUCTORS(SQL);

        struct Vector {
            // Overall variable data.
            std::vector <Real> data;

            // Offsets of each cone stored in the data.
            std::vector <Natural> offsets;

            // Type of cones stored in the data.
            std::vector <Cone::t> types;
            
            // Size of the cones stored in the data.
            std::vector <Natural> sizes;

            // Cached matrix inverses.  Once we have the offset, we store the
            // matrix.
            mutable std::vector <Real> inverse;

            // Offsets of the cached matrix inverses 
            std::vector <Natural> inverse_offsets;

            // Point where we last took the matrix inverse 
            mutable std::vector <Real> inverse_base;

            // Offsets for the bases stored for the matrix inverses 
            std::vector <Natural> inverse_base_offsets;

            // Eliminate constructors 
            NO_DEFAULT_COPY_ASSIGNMENT(Vector);

            // We require a vector of cone types and their sizes.
            Vector (
                std::vector <Cone::t> const & types_,
                std::vector <Natural> const & sizes_,
                Messaging const msg = Optizelle::Messaging()
            ) : data(), offsets(), types(types_), sizes(sizes_),
                inverse(), inverse_offsets(), inverse_base(),
                inverse_base_offsets()
            {

                // Insure that the type of cones and their sizes lines up.
                if(types.size()!=sizes.size())
                    msg.error("The vector containing the type of cones must "
                        "be the same size as the vector with the cone sizes.");

                // Make sure we have at least one cone.
                if(types.size() == 0)
                    msg.error("A SQL vector requires at least one cone.");

                // Initialize the offsets.  The last element has the total
                // number of variables.
                offsets.resize(sizes.size()+1);
                offsets.front()=0;
                for(Natural i=2;i<=offsets.size();i++)
                    offsets[itok(i)] = types[itok(i-1)]==Cone::Linear ||
                                       types[itok(i-1)]==Cone::Quadratic
                                     ? offsets[itok(i-1)]+sizes[itok(i-1)]
                                     : offsets[itok(i-1)]+sizes[itok(i-1)]
                                         *sizes[itok(i-1)];

                // Create the data.
                data.resize(offsets.back());

                // Calculate offsets for the matrix inverses.  Basically,
                // the way it works is that we calculate offsets for every cone
                // even though we're not ever going to use them.  In the case
                // we don't have an SDP block, we simply use the last offset.
                // This makes it easy to index to the correct place where the
                // cached information is stored.
                inverse_offsets.resize(sizes.size()+1);
                inverse_offsets.front()=0;
                inverse_base_offsets.resize(sizes.size()+1);
                inverse_base_offsets.front()=0;
                for(Natural i=1;i<types.size()+1;i++) {
                    inverse_offsets[i] =
                        types[itok(i)]==Cone::Linear ||
                        types[itok(i)]==Cone::Quadratic
                            ? inverse_offsets[itok(i)]
                            : inverse_offsets[itok(i)]
                                +sizes[itok(i)]*sizes[itok(i)];
                    inverse_base_offsets[i] =
                        types[itok(i)]==Cone::Linear ||
                        types[itok(i)]==Cone::Quadratic
                            ? inverse_base_offsets[itok(i)]
                            : inverse_base_offsets[itok(i)]
                                +sizes[itok(i)]*sizes[itok(i)];
                }

                // Create the memory required for the cached Schur
                // decompositions.
                inverse.resize(inverse_offsets.back());
                inverse_base.resize(inverse_base_offsets.back());
            }
            
            // Move constructor 
            Vector(Vector&& x) noexcept : 
                data(std::move(x.data)),
                offsets(std::move(x.offsets)),
                types(std::move(x.types)),
                sizes(std::move(x.sizes)),
                inverse(std::move(x.inverse)),
                inverse_offsets(std::move(x.inverse_offsets)),
                inverse_base(std::move(x.inverse_base)),
                inverse_base_offsets(std::move(x.inverse_base_offsets))
            {}

            // Move assignment operator
            Vector const & operator = (Vector&& x) noexcept {
                data=std::move(x.data);
                offsets=std::move(x.offsets);
                types=std::move(x.types);
                sizes=std::move(x.sizes);
                inverse=std::move(x.inverse);
                inverse_offsets=std::move(inverse_offsets);
                inverse_base=std::move(inverse_base);
                inverse_base_offsets=std::move(inverse_base_offsets);
                return *this;
            }

            // Simple indexing.
            Real & operator () (Natural const & i) {
                return data[itok(i)];
            }
            Real const & operator () (Natural const & i) const {
                return data[itok(i)];
            }

            // Indexing with multiple cones.
            Real & operator () (Natural const & k,Natural const & i) {
                return data[offsets[itok(k)]+itok(i)];
            }
            Real const & operator () (Natural const & k,Natural const & i)const{
                return data[offsets[itok(k)]+itok(i)];
            }

            // Indexing a matrix with multiple cones.
            Real & operator () (
                Natural const & k,Natural const & i,Natural const & j
            ) {
                return data[offsets[itok(k)]+ijtok(i,j,sizes[itok(k)])];
            }
            Real const & operator ()(
                Natural const & k,Natural const & i,Natural const & j
            ) const {
                return data[offsets[itok(k)]+ijtok(i,j,sizes[itok(k)])];
            }

            // First element of the block
            Real const & front(Natural const & blk) const {
                return (*this)(blk,1,1);
            }
            Real & front(Natural const & blk) {
                return (*this)(blk,1,1);
            }

            // These are really shortcuts for second-order cone blocks, which
            // are typically structured as x=[x0;xbar].  The naught function
            // gives the first element whereas the bar function gives the
            // first element of the xbar portion.
            Real const & naught(Natural const & blk) const {
                return this->front(blk);
            }
            Real & naught(Natural const & blk) {
                return this->front(blk);
            }
            Real const & bar(Natural const & blk) const {
                return (*this)(blk,2);
            }
            Real & bar(Natural const & blk) {
                return (*this)(blk,2);
            }

            // Size of the block.
            Natural blkSize(Natural const & blk) const {
                return sizes[itok(blk)];
            }

            // Type of the block.
            Cone::t blkType(Natural const & blk) const {
                return types[itok(blk)];
            }

            // Number of blocks.
            Natural numBlocks() const {
                return types.size();
            }
        };

        // Gets the matrix inverse of a block of the SQL vector.
        static void get_inverse(
            Vector const & X,
            Natural const & blk,
            std::vector <Real> & Xinv 
        ) {
            // Get the size of the block
            const Natural m=X.sizes[itok(blk)];

            // Next, check if we've already calculated the matrix inverse.

            // Copy out the the base of the last inverse 
            std::vector <Real> tmp(m*m);
            Optizelle::copy <Real>
                (m*m,&(X.inverse_base[X.inverse_base_offsets[itok(blk)]]),
                1,&(tmp.front()),1);

            // tmp <- Base_k - X_k
            Optizelle::axpy <Real>
                (m*m,Real(-1.),&(X.data[X.offsets[itok(blk)]]),
                1,&(tmp.front()),1);

            // Find the relative error between the current iterate
            // and the base
            Real norm_xk = sqrt(dot <Real>
                (m*m,&(X.data[X.offsets[itok(blk)]]),1,
                &(X.data[X.offsets[itok(blk)]]),1));
            Real rel_err = sqrt(dot<Real> (m,&(tmp.front()),1,&(tmp.front()),1))
                / (std::numeric_limits <Real>::epsilon()+norm_xk);

            // If the relative error is large, refresh the cached decomposition.
            // To be sure, I'm not really sure what this tolerance should be
            // given that we have both floats and doubles.  In reality,
            // the iterates will probably change rapidly, so I don't think
            // we have to worry too much, but be careful.
            if(rel_err > std::numeric_limits <Real>::epsilon()*1e2) {

                // Store X_k as the new base
                Optizelle::copy<Real> (m*m,&(X.data[X.offsets[itok(blk)]]),1,
                    &(X.inverse_base[X.inverse_base_offsets[itok(blk)]]),1);

                // Find the matrix inverse of X_k 

                // Find the matrix inverse.  This assumes the input is
                // symmetric positive definite.
                Integer info(0);
                Optizelle::copy<Real> (m*m,&(X.data[X.offsets[itok(blk)]]),1,
                    &(X.inverse[X.inverse_offsets[itok(blk)]]),1);
                Optizelle::potrf <Real> ('U',m,
                    &(X.inverse[X.inverse_offsets[itok(blk)]]),m,info);
                Optizelle::potri <Real> ('U',m,
                    &(X.inverse[X.inverse_offsets[itok(blk)]]),m,info);
               
                // Copy the upper triangular portion to the lower.
                for(Natural i=1;i<=m;i++)
                    Optizelle::copy <Real> (m-i,
                        &(X.inverse[X.inverse_offsets[itok(blk)]
                            +ijtok(i,i+1,m)]),m,
                        &(X.inverse[X.inverse_offsets[itok(blk)]
                            +ijtok(i+1,i,m)]),1);

            }

            // Copy out the inverse from the cached copy
            Xinv.resize(m*m);
            Optizelle::copy <Real>
                (m*m,&(X.inverse[X.inverse_offsets[itok(blk)]]),1,
                &(Xinv.front()),1);
        }
        
        // Memory allocation and size setting
        static Vector init(Vector const & x) {
            return std::move(Vector(x.types,x.sizes));
        }
        
        // y <- x (Shallow.  No memory allocation.)
        static void copy(Vector const & x, Vector & y) {
            Optizelle::copy <Real> (x.data.size(),&(x.data.front()),1,
                &(y.data.front()),1);
        }

        // x <- alpha * x
        static void scal(Real const & alpha, Vector & x) {
            Optizelle::scal <Real> (x.data.size(),alpha,&(x.data.front()),1);
        }

        // y <- alpha * x + y
        static void axpy(Real const & alpha, Vector const & x, Vector & y) {
            Optizelle::axpy <Real> (x.data.size(),alpha,&(x.data.front()),1,
                &(y.data.front()),1);
        }

        // innr <- <x,y>
        static Real innr(Vector const & x,Vector const & y) {
            return Optizelle::dot<Real> (x.data.size(),&(x.data.front()),1,
                &(y.data.front()),1);
        }

        // x <- 0 
        static void zero(Vector & x) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for(Natural i=0;i<x.data.size();i++) 
                x.data[i]=Real(0.);
        }

        // x <- random
        static void rand(Vector & x){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<Real> dis(Real(0.),Real(1.));

            // This is not parallel since it doesn't appear that our generator
            // works properly when parallel.
            for(Natural i=0;i<x.data.size();i++) 
                x.data[i]=Real(dis(gen));
        }

        // Jordan product, z <- x o y
        static void prod(Vector const & x, Vector const & y, Vector & z) {
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
            // Loop over all the blocks.
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block.
                Natural m=x.blkSize(blk);

                // Depending on the block, compute a different jordan product.
                switch(x.blkType(blk)) {

                // z = diag(x) y.
                case Cone::Linear:
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++)
                        z(blk,i)=x(blk,i)*y(blk,i);
                    break;

                // z = [x'y ; x0 ybar + y0 xbar].
                case Cone::Quadratic: {
                    // Get the size of the bar section.
                    Natural mbar=m-1;

                    // Find the first element
                    z.naught(blk)=
                        Optizelle::dot<Real>(
                            m,&(x.front(blk)),1,&(y.front(blk)),1);

                    // zbar = ybar
                    Optizelle::copy <Real> (
                        mbar,&(y.bar(blk)),1,&(z.bar(blk)),1);
                        
                    // zbar = x0 ybar
                    Optizelle::scal <Real> (mbar,x.naught(blk),&(z.bar(blk)),1);

                    // zbar = x0 ybar + y0 xbar
                    Optizelle::axpy <Real> (mbar,y.naught(blk),&(x.bar(blk)),1,
                        &(z.bar(blk)),1);
                    break;
                }

                // z = xy 
                case Cone::Semidefinite:
                    Optizelle::symm <Real> ('L','U',m,m,Real(1.),
                        &(x.front(blk)),m,&(y.front(blk)),m,Real(0.),
                        &(z.front(blk)),m);
                    break;
                }
            }
        }

        // Identity element, x <- e such that x o e = x
        static void id(Vector & x) {

            // Loop over all the blocks
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                Natural m=x.blkSize(blk);

                // Depending on the block, compute a different identity element
                switch(x.blkType(blk)) {

                // x = vector of all 1s
                case Cone::Linear:
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++) 
                        x(blk,i)=Real(1.);
                    break;
                // x = (1,0,...,0)
                case Cone::Quadratic:
                    x(blk,1)=Real(1.);
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=2;i<=m;i++)
                        x(blk,i)=Real(0.);
                    break;
                // x = I
                case Cone::Semidefinite:
                    // We write the diagonal elements twice to avoid the
                    // conditional.
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural j=1;j<=m;j++) 
                        for(Natural i=1;i<=m;i++) 
                            x(blk,i,j)=Real(0.);

                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++) 
                        x(blk,i,i)=Real(1.);
                    break;
                }
            }
        }

        // This applies the inverse of the Schur complement of the Arw
        // operator to a vector.  Note, y has size one less than x.
        // Hence, m is the length of y and m+1 is the length of x.
        static void invSchur(Natural const & m,Real const * const x,Real * y) {

            // innr_xbar_y <- <xbar,y>
            Real innr_xbar_y = dot <Real> (m,&(x[1]),1,&(y[0]),1);

            // denom <- x0 ( x0^2 - <xbar,xbar> )
            Real denom = x[0]*( x[0]*x[0] - dot <Real> (m,&(x[1]),1,&(x[1]),1));

            // y <- 1/x0 y
            Optizelle::scal <Real> (m,Real(1.)/x[0],&(y[0]),1);

            // y <- 1/x0 y + <xbar,y> / (x0 ( x0^2 - <xbar,xbar> )) xbar
            Optizelle::axpy <Real> (m,innr_xbar_y/denom,&(x[1]),1,&(y[0]),1);
        }
        
        // Jordan product inverse, z <- inv(L(x)) y where L(x) y = x o y
        static void linv(Vector const & x,Vector const & y,Vector & z) {
            // We have this vector in case we have a SDP block
            std::vector <Real> Xinv;

            // Loop over all the blocks
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                Natural m=x.blkSize(blk);

                // Depending on the block, compute a different operator
                switch(x.blkType(blk)) {

                // z = inv(Diag(x)) y
                case Cone::Linear:
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++) 
                        z(blk,i)=y(blk,i)/x(blk,i);
                    break;

                // z = inv(Arw(x)) y
                case Cone::Quadratic: {
                    // Get the size of the bar piece
                    Natural mbar=m-1;

                    // invSchur_ybar <- invSchur(x)(y_bar) 
                    std::vector <Real> invSchur_ybar(mbar);
                    Optizelle::copy <Real>(mbar,&(y.bar(blk)),1,
                        &(invSchur_ybar.front()),1);
                    invSchur(mbar,&(x.front(blk)),&(invSchur_ybar[0]));

                    // a <- 1 / (x0 - (1/x0) <x_bar,x_bar>) * y0
                    Real a = y.naught(blk) / (
                        x.naught(blk) - (Real(1.)/x.naught(blk)) *
                        Optizelle::dot<Real>(
                            mbar,&(x.bar(blk)),1,&(x.bar(blk)),1));

                    // b <- - (1/x0) <xbar,invSchur(x)(y_bar)>
                    Real b = -Optizelle::dot <Real> (mbar,
                            &(x.bar(blk)),1,&(invSchur_ybar.front()),1)
                        / x.naught(blk);
                    
                    // z0 <- 1 / (x0 - (1/x0) <x_bar,x_bar>) y0
                    //       - (1/x0) <x_bar,invSchur(x)(y_bar)> 
                    z(blk,1) = a + b;
                    
                    // z_bar <- invSchur(x)(x_bar)
                    Optizelle::copy <Real> (
                        mbar,&(x.bar(blk)),1,&(z.bar(blk)),1);
                    invSchur(mbar,&(x.front(blk)),&(z.bar(blk)));

                    // zbar <- (-y0/x0) invSchur(x)(x_bar)
                    Optizelle::scal <Real> (mbar,
                        -y.naught(blk)/x.naught(blk),&(z.bar(blk)),1);

                    // z_bar <- (-y0/x0) invSchur(x)(x_bar) + invSchur(x)(y_bar)
                    Optizelle::axpy <Real> (mbar,Real(1.),
                        &(invSchur_ybar.front()),1,&(z.bar(blk)),1);
                    break;

                // Z=inv(X) Y
                } case Cone::Semidefinite: {
                    // Get the Schur complement of the block.  With any luck
                    // these are cached.
                    Optizelle::SQL <Real>::get_inverse(x,blk,Xinv);

                    // Multiply out the result
                    Optizelle::symm <Real> ('L','U',m,m,Real(1.),
                        &(Xinv.front()),m,&(y.front(blk)),m,Real(0.),
                        &(z.front(blk)),m);
                    break;
                }}
            }
        }

        // Barrier function, barr <- barr(x) where x o grad barr(x) = e
        static Real barr(Vector const & x) {
            // This accumulates the barrier's value
            Real z(0.);

            // Loop over all the blocks
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                Natural m=x.blkSize(blk);

                // In case we need to take a Choleski factorization for the
                // SDP blocks
                std::vector <Real> U;

                // Depending on the block, compute a different barrier
                switch(x.blkType(blk)) {

                // z += sum_i log(x_i)
                case Cone::Linear:
                    #ifdef _OPENMP
                    #pragma omp parallel for reduction(+:z) schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++)
                        z+=log(x(blk,i));
                    break;

                // z += 0.5 * log(x0^2-<xbar,xbar))
                case Cone::Quadratic: {
                    // Get the size of the bar part.
                    Natural mbar=m-1;

                    z+=Real(0.5) * log(x.naught(blk)*x.naught(blk)
                        -dot <Real> (mbar,&(x.bar(blk)),1,&(x.bar(blk)),1));
                    break;
                }

                // z += log(det(x)).  We compute this by noting that
                // log(det(x)) = log(det(u'u)) = log(det(u')det(u))
                //             = log(det(u)^2) = 2 log(det(u))
                case Cone::Semidefinite: {

                    // Find the Choleski factorization of X
                    U.resize(m*m);
                    Integer info;
                    Optizelle::copy <Real> (
                        m*m,&(x.front(blk)),1,&(U.front()),1);
                    Optizelle::potrf <Real> (
                        'U',m,&(U.front()),m,info);

                    Real log_det(0.);
                    #ifdef _OPENMP
                    #pragma omp parallel for reduction(+:log_det) schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++)
                        log_det += log(U[Optizelle::ijtok(i,i,m)]);
                    
                    // Complete the barrier computation by taking the log
                    z+= Real(2.) * log_det;
                    break;
                } }
            }

            // Return the accumulated barrier value
            return z;
        }

        // Line search, srch <- argmax {alpha \in Real >= 0 : alpha x + y >= 0}
        // where y > 0.
        static Real srch(Vector const & x,Vector const & y) {
            // Line search parameter
            Real alpha=std::numeric_limits <Real>::infinity();
                   
            // Variables required for the linesearch on SDP blocks 
            Integer info(0);
            std::vector <Real> Xrf;
            std::vector <Real> Yrf;
            std::vector <Real> Zrf;

            // Loop over all the blocks
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                Natural m=x.blkSize(blk);

                // Depending on the block, do a different line search 
                switch(x.blkType(blk)) {

                // Pointwise, alpha_i = -y_i / x_i.  If this number is positive,
                // then we need to restrict how far we travel.
                case Cone::Linear:

                    #ifdef _OPENMP
                    #pragma omp parallel
                    #endif
                    {
                        // Create a local version of alpha
                        Real alpha_loc=std::numeric_limits <Real>::infinity();

                        // Search for the optimal linesearch parameter
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(static)
                        #endif
                        for(Natural i=1;i<=m;i++) {
                            if(x(blk,i) < Real(0.)) {
                                Real alpha0 = -y(blk,i)/x(blk,i);
                                alpha_loc = alpha0 < alpha_loc ?
                                    alpha0 : alpha_loc;
                            }
                        }

                        // After we're through with the local search,
                        // accumulate the result
                        #ifdef _OPENMP
                        #pragma omp critical
                        #endif
                        {
                            alpha = alpha_loc < alpha ? alpha_loc : alpha;
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

                    // Get the size of the bar block.
                    Natural mbar=m-1;

                    // Now, first we have to insure that the leading coefficient
                    // of the second order cone problem remains nonnegative.
                    // This number tells us how far we can step before this
                    // is not true.
                    Real alpha0 = -y.naught(blk)/x.naught(blk);

                    // Next, assuming that the leading coefficient is fine,
                    // figure out how far we can step before we violate the
                    // rest of the SOCP constraint.  This involves solving
                    // the quadratic equation from above.
                    Real a = x.naught(blk)*x.naught(blk)
                        - dot <Real> (mbar,&(x.bar(blk)),1,&(x.bar(blk)),1);
                    Real b = Real(2.)*(x.naught(blk)*y.naught(blk)
                        - dot <Real> (mbar,&(x.bar(blk)),1,&(y.bar(blk)),1));
                    Real c = y.naught(blk)*y.naught(blk)
                        - dot <Real> (mbar,&(y.bar(blk)),1,&(y.bar(blk)),1);
                    Natural nroots(0);
                    Real alpha1(-1.);
                    Real alpha2(-1.);
                    quad_equation(a,b,c,nroots,alpha1,alpha2);

                    // Now, determine the step length.

                    // If we have a restriction on alpha based on the leading
                    // coefficient.
                    alpha = x.naught(blk) < Real(0.) && alpha0<alpha ?
                        alpha0 : alpha;

                    // Next, if we have two roots, determine the restriction 
                    if(nroots==2) { 
                        alpha = alpha1>=Real(0.)&&alpha1<alpha ? alpha1 : alpha;
                        alpha = alpha2>=Real(0.)&&alpha2<alpha ? alpha2 : alpha;

                    // If we have a single root 
                    } else if(nroots==1)
                        alpha = alpha1>=Real(0.)&&alpha1<alpha ? alpha1 : alpha;

                    // If we no roots, there's no additional restriction.
                    // This can't happen since we assume that y is strictly
                    // feasible.
                }
                break;

                // We need to find the solution of the generalized eigenvalue
                // problem alpha X v + Y v = 0.  Since Y is positive definite,
                // we want to divide by alpha to get a standard form
                // generalized eigenvalue problem X v = (-1/alpha) Y v.  This
                // means that we solve the problem X v = lambda Y v and then
                // set alpha = -1/lambda as long as lambda is negative.  Note,
                // our Krylov method will converge to lambda from the right,
                // which is going to give an upper bound on alpha.  This is
                // not good for our line search, since we want a lower bound.
                // However, since we get an absolute estimate of the error
                // in lambda, we can just back off of it by a small amount.
                case Cone::Semidefinite: {

                    // Convert X and Y to rectangular packed storage
                    Xrf.resize(m*(m+1)/2);
                    Optizelle::trttf <Real>('N','U',m,&(x(blk,1,1)),m,&(Xrf[0]),
                        info);

                    Yrf.resize(m*(m+1)/2);
                    Optizelle::trttf <Real>('N','U',m,&(y(blk,1,1)),m,&(Yrf[0]),
                        info);

                    // Solve the generalized eigenvalue problem X v = lambda Y v
                    Real abs_tol=1e-2;
                    std::pair <Real,Real> lambda_err=Optizelle::gsyiram <Real> (
                        m,&(Xrf[0]),&(Yrf[0]),20,20,abs_tol);

                    // IRAM converges from the right, but we really need a lower
                    // bound on the eigenvalue.  Hence, modify the result
                    // so that we have a lower bound
                    Real lambda=lambda_err.first-abs_tol;

                    // Now, find the line-search parameter
                    Real alpha0=-Real(1.)/lambda;

                    // Do a safeguard step because sometimes the eigenvalue
                    // solver converges to the wrong eigenvalue.  Now, if
                    // alpha0 is negative, ostensibly we can take as big
                    // as step as we want.  However, if we converged to
                    // the wrong eigenvalue, this may not be true.  Hence, 
                    // if alpha0 is negative, we do the line-search with
                    // alpha0 = 2.  If this value doesn't move, we assume
                    // that our eigenvalue estimate was fine and this
                    // direction is feasible for all alpha.
                    //
                    // Also, note that the Choleski check is not full-proof.
                    // It's possible that the Choleski check passes and yet
                    // we have an indefinite matrix.  This is sort of hard
                    // to check.  Basically, that means that the next iteration
                    // will have an infeasible solution, which is going to
                    // cause issues with this routine.  In theory, we should
                    // continue to cut alpha0 until it becomes a hard 0 and
                    // then this routine will exit.  Hopefully, the other
                    // pieces in the code will pick up on the interior point
                    // instability and exit.
                    bool completely_feasible_dir= alpha0<=Real(0.);
                    alpha0= alpha0>0 ? alpha0 : Real(2.);
                    Zrf.resize(m*(m+1)/2);
                    do {
                        // Basically, we find X+alpha0 Y and try to take
                        // the Choleski factorization.  If that fails, we're
                        // infeasible and we do a backtracking line search.
                        Optizelle::copy <Real> (
                            m*(m+1)/2,&(Yrf[0]),1,&(Zrf[0]),1);
                        Optizelle::axpy <Real> (m*(m+1)/2,alpha0,&(Xrf[0]),1,
                            &(Zrf[0]),1);
                        pftrf('N','U',m,&(Zrf[0]),info);

                        // Check if the Choleski failed
                        if(info!=0) {
                            alpha0 /= Real(2.); 
                            completely_feasible_dir=false;
                        }
                    
                    // If alpha0 ever becomes 0, then something wrong has
                    // gone on and we really ought to exit.
                    } while(info!=0 && alpha0>Real(0.));

                    // If we still have a completely feasible direction,
                    // fix alpha0 so that we don't update our line search.
                    alpha0 = completely_feasible_dir ?
                        std::numeric_limits <Real>::infinity(): alpha0;

                    // Adjust the line search step if necessary
                    alpha = alpha0<alpha ? alpha0 : alpha;
                } }
            }
            return alpha;
        }
        // Symmetrization, x <- symm(x) such that L(symm(x)) is a symmetric
        // operator.
        static void symm(Vector & x) { 
            // Allocate vectors to help with the SDP blocks
            std::vector <Real> I;
            std::vector <Real> Xk;

            // Loop over all the blocks
            for(Natural blk=1;blk<=x.numBlocks();blk++) {

                // Get the size of the block
                Natural m=x.blkSize(blk);

                // Depending on the block, do a different symmetrization 
                switch(x.blkType(blk)) {

                // Linear and quadratic cones don't have this issue
                case Cone::Linear:
                case Cone::Quadratic:
                    break;

                // Find the symmetric part of X, (X+X')/2
                case Cone::Semidefinite: {
                    // Create the identity matrix
                    I.resize(m*m);
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++) 
                        for(Natural j=1;j<=m;j++) 
                            I[ijtok(i,j,m)]=Real(0.);

                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
                    #endif
                    for(Natural i=1;i<=m;i++) 
                        I[ijtok(i,i,m)]=Real(1.);

                    // Create a copy of X
                    std::vector <Real> Xk(m*m);
                    Optizelle::copy<Real> (
                        m*m,&(x.front(blk)),1,&(Xk.front()),1);

                    // X <- (X+X')/2
                    syr2k <Real> ('U','N',m,m,Real(0.5),&(Xk.front()),m,
                        &(I.front()),m,Real(0.),&(x.front(blk)),m);

                    // Copy the upper part of X to the lower
                    for(Natural i=1;i<=m;i++)
                        Optizelle::copy <Real> (
                            m-i,&(x(blk,i,i+1)),m,&(x(blk,i+1,i)),1);
                    break;
                } }
            }
        }
    };

    namespace json {
        // Serialization utility for the SQL vector space
        template <typename Real>
        struct Serialization <Real,SQL> {
            static std::string serialize (
                typename SQL <Real>::Vector const & x
            ) {
                // Create a jsoncpp object to copy into
                Json::Value x_json;  

                // Copy the information
                for(Natural i=0;i<x.data.size();i++)
                    x_json["data"][Json::ArrayIndex(i)]=x.data[i];

                for(Natural i=0;i<x.offsets.size();i++)
                    x_json["offsets"][Json::ArrayIndex(i)]
                        =Json::Value::UInt64(x.offsets[i]);

                for(Natural i=0;i<x.types.size();i++)
                    x_json["types"][Json::ArrayIndex(i)]
                        =Cone::to_string(x.types[i]);

                for(Natural i=0;i<x.sizes.size();i++)
                    x_json["sizes"][Json::ArrayIndex(i)]
                        =Json::Value::UInt64(x.sizes[i]);

                for(Natural i=0;i<x.inverse.size();i++)
                    x_json["inverse"][Json::ArrayIndex(i)]=x.inverse[i];

                for(Natural i=0;i<x.inverse_offsets.size();i++)
                    x_json["inverse_offsets"][Json::ArrayIndex(i)]
                        =Json::Value::UInt64(x.inverse_offsets[i]);

                for(Natural i=0;i<x.inverse_base.size();i++)
                    x_json["inverse_base"][Json::ArrayIndex(i)]
                        =x.inverse_base[i];

                for(Natural i=0;i<x.inverse_base_offsets.size();i++)
                    x_json["inverse_base_offsets"][Json::ArrayIndex(i)]
                        =Json::Value::UInt64(x.inverse_base_offsets[i]);
                
                // Return a string of the result
                Json::StyledWriter writer;
                return writer.write(x_json);
            }
            static typename SQL <Real>::Vector deserialize (
                typename SQL <Real>::Vector const & x_,
                std::string const & x_json_
            ) {
                // Create a json tree from the input string
                Json::Value x_json;
                Json::Reader reader;
                reader.parse(x_json_,x_json,true);

                // Grab the types of the cones 
                std::vector <Cone::t> types;
                types.resize(x_json["types"].size());
                for(Natural i=0;i<types.size();i++)
                    types[i]=Cone::from_string(x_json["types"]
                        [Json::ArrayIndex(i)].asString());

                // Grab the sizes of the cones
                std::vector <Natural> sizes;
                sizes.resize(x_json["sizes"].size());
                for(Natural i=0;i<sizes.size();i++)
                    sizes[i]=x_json["sizes"][Json::ArrayIndex(i)]
                        .asUInt64();

                // Allocate a new SQL vector
                typename SQL <Real>::Vector x(types,sizes);

                // Read in the data
                for(Natural i=0;i<x.data.size();i++)
                    x.data[i]=Real(x_json["data"][Json::ArrayIndex(i)]
                        .asDouble());

                for(Natural i=0;i<x.offsets.size();i++)
                    x.offsets[i]=x_json["offsets"][Json::ArrayIndex(i)]
                        .asUInt64();

                for(Natural i=0;i<x.inverse.size();i++)
                    x.inverse[i]=Real(x_json["inverse"]
                        [Json::ArrayIndex(i)].asDouble());
                
                for(Natural i=0;i<x.inverse_offsets.size();i++)
                    x.inverse_offsets[i]=x_json["inverse_offsets"]
                        [Json::ArrayIndex(i)].asUInt64();
                
                for(Natural i=0;i<x.inverse_base.size();i++)
                    x.inverse_base[i]=Real(x_json["inverse_base"]
                        [Json::ArrayIndex(i)].asDouble());
                
                for(Natural i=0;i<x.inverse_base_offsets.size();i++)
                    x.inverse_base_offsets[i]=x_json["inverse_base_offsets"]
                        [Json::ArrayIndex(i)].asUInt64();

                // Return the newly constructed vector
                return std::move(x);
            }
        };
    }

    // Optimization problems instantiated on these vector spaces.  In theory,
    // this should help our compilation times.
    extern template struct Unconstrained<double,Rm>;
    extern template struct Unconstrained<float,Rm>;
    extern template struct EqualityConstrained<double,Rm,Rm>;
    extern template struct EqualityConstrained<float,Rm,Rm>;
    extern template struct InequalityConstrained<double,Rm,Rm>;
    extern template struct InequalityConstrained<float,Rm,Rm>;
    extern template struct InequalityConstrained<double,Rm,SQL>;
    extern template struct InequalityConstrained<float,Rm,SQL>;
    extern template struct Constrained<double,Rm,Rm,Rm>;
    extern template struct Constrained<float,Rm,Rm,Rm>;
    extern template struct Constrained<double,Rm,Rm,SQL>;
    extern template struct Constrained<float,Rm,Rm,SQL>;
//---Optizelle2---
}
//---Optizelle3---
#endif
