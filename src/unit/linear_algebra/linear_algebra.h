#pragma once

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "unit.h"

// For the time being, just lump everything in the unit namespace 
namespace Unit {

    // Grab the natural type from Optizelle
    typedef Optizelle::Natural Natural;

    // Different kinds of matrices 
    template <typename Real>
    struct Matrix {
        // Matrix type.  This needs to ultimately be an Optizelle operator
        struct t :
            public Optizelle::Operator <Real,Optizelle::Rm,Optizelle::Rm>
        {
        private:
            // Create some type shortcuts
            typedef std::vector <Real> X_Vector;
            
            // Size of the matrix
            Natural size;
        public:
            // Storage for the matrix
            std::vector <Real> data;

            // Create an empty matrix filled with some arbitrary data 
            t(Natural const & size_)
                : size(size_), data(size_*size_)
            {}

            // Apply the matrix to the vector 
            void eval(X_Vector const & x,X_Vector & y) const {
                for(auto i=0;i<size;i++) {
                    y[i]=0;
                    for(auto j=0;j<size;j++)
                        y[i]+=data[i+size*j]*x[j];
                }
            }
        };
        
        // Create a matrix with the scheme where the (i,j)th element is
        //
        // cos(I^2)
        //
        // where
        //
        // I = i+(j-1)*m
        // m the matrix size
        //
        // Basically, we just need something arbitrary, but fixed for testing.
        static t nonsymmetric(Natural const & size,Natural const & offset) {
            auto A = t(size);
            for(auto j=1;j<=size;j++)
                for(auto i=1;i<=size;i++) {
                    auto I = j+(i-1)*size;
                    A.data[I-1] = cos(std::pow(I+offset,2)); 
                }
            return A;
        }

        // Creates a matrix with the scheme where the (i,j)th element is
        //
        // cos(I^2))                  when i>j
        // abs(cos(I^2))) + m + i^2   when i==j
        // (i,j)th element            when j<i
        //
        // where
        //
        // I = i+(j-1)*m
        // m the matrix size
        //
        // Again, we just need something arbitrary, but fixed for testing.
        // Since we need something for CG, the m factor is there to force
        // things to be diagonally dominant and hence positive definite wheras
        // the element i^2 is there to make sure the spectrum is spread out, so
        // we don't converge too quickly.
        static t symmetric(Natural const & size,Natural const & offset) {
            auto A = t(size);
            for(auto j=1;j<=size;j++)
                for(auto i=1;i<=j;i++) {
                    auto I = i+(j-1)*size;
                    auto J = j+(i-1)*size;
                    if(i!=j) {
                        auto ele0 = std::pow(I+offset,2);
                        auto ele = cos(std::pow(I+offset,2));
                        A.data[I-1] = ele;
                        A.data[J-1] = ele;
                    } else if(i==j)
                        A.data[I-1]=std::fabs(cos(std::pow(I+offset,2)))+size
                            + i*i;
                }
            return A;
        }

        // Create a symmetric negative definite matrix
        static t symm_nd(Natural const & size,Natural const & offset) {
            auto A = symmetric(size,offset);
            for(auto k=0;k<A.data.size();k++) {
                A.data[k] = -A.data[k];
            }
            return A;
        }

        // Fills the matrix with diagonal matrix where the (i,i)th element is
        //
        // (i+1)/2     when i is odd
        // 0                i is even
        static t diagonal(Natural const & size) {
            auto A = t(size);
            A.data.assign(size*size,Real(0.));
            for(auto i=1;i<=size;i++) {
                auto I = i+(i-1)*size;
                if(i % 2 == 1) 
                    A.data[I-1] = (i+1)/Real(2.);
            }
            return A;
        }

        // Fills the matrix with diagonal matrix where the (i,i)th element is
        //
        // 2/(i+1)     when i is odd
        // 0                i is even
        //
        // Basically, it's the inverse of the matrix above.
        static t diagonal_inv(Natural const & size) {
            auto A = t(size);
            A.data.assign(size*size,Real(0.));
            for(auto i=1;i<=size;i++) {
                auto I = i+(i-1)*size;
                if(i % 2 == 1) 
                    A.data[I-1] = Real(2.)/(i+1);
            }
            return A;
        }

        // Fills the matrix with diagonal matrix where the (i,i)th element is
        //
        // 1       when i \in {1,2} 
        // 0            0 otherwise 
        //
        // Essentially, we project out the first two elements and stay in the
        // nullspace of the rest.
        static t project_2(Natural const & size) {
            auto A = t(size);
            A.data.assign(size*size,Real(0.));
            for(auto i=1;i<=size;i++) {
                auto I = i+(i-1)*size;
                if(i<=2) 
                    A.data[I-1] = Real(1.); 
            }
            return A;
        }

        // Creates a mostly diagonally dominant matrix, with a single negative
        // eigenvalue.  We set the (i,j) element to
        //
        // 2    when i=j, i>1 
        // -2   when i=j, i=1
        // 1    when i-1=j
        // 1    when i=j-1
        static t mostly_dd_indef(Natural const & size) {
            auto A = t(size);
            A.data.assign(size*size,Real(0.));
            for(auto i=1;i<=size;i++) {
                auto I = i-1+(i-1)*size;
                auto up = i-2+(i-1)*size;
                auto left = i-1+(i-2)*size;
                if(i == 1) 
                    A.data[I] = Real(-2.);
                else {
                    A.data[I] = Real(2.); 
                    A.data[up] = Real(1.); 
                    A.data[left] = Real(1.); 
                }
            }
            return A;
        }
    };

    // Operators
    template <typename Real,template <typename> class XX>
    struct Operator {
        // Create some type shortcuts
        typedef Optizelle::Rm<Real> X;
        typedef typename X::Vector X_Vector;
        typedef Optizelle::Operator <Real,XX,XX> t;

        // Throw a NaN after a specified number of calls
        struct Nan: public t { 
        private:
            // Track the number of times the operator is called
            mutable Natural calls;

            // Track the maximum number of calls
            Natural const calls_max;

            // Track the internal matrix
            typename Matrix <Real>::t A;
        public:
            // On creation, set the number of calls before NaNing
            Nan(typename Matrix <Real>::t const& A_,Natural const& calls_max_) :
                calls(0),
                calls_max(calls_max_),
                A(A_)
            {}
                
            // Apply the matrix until we NaN 
            void eval(X_Vector const & x,X_Vector & y) const {
                A.eval(x,y);
                if(calls >= calls_max)
                    X::scal(std::numeric_limits <Real>::quiet_NaN(),y);
                calls++;
            }
        };

        // Identity operator
        typedef typename Optizelle::Unconstrained <Real,XX>::Functions::Identity
            Identity;
    };

    // A bunch of structured vectors
    template <typename Real>
    struct Vector {

        // Create a vector where the the ith element is cos(i+25)
        static std::vector <Real> basic(Natural & size) {
            auto b = std::vector <double>(size);
            for(auto i=1;i<=size;i++)
                b[i-1] = cos(i+25); 
            return b;
        }

        // Create a zero vector
        static std::vector <Real> zero(Natural & size) {
            auto b = std::vector <double>(size,Real(0.));
            return b;
        }

        // Create a rhs where the odd elements are 1 and the even are 0 
        static std::vector <Real> alternate(Natural const & size) {
            auto b = std::vector <double>(size);
            for(auto i=1;i<=size;i++) {
                if(i%2==1)
                    b[i-1] = Real(1.);
                else
                    b[i-1] = Real(0.);
            }
            return b;
        }

        // Create a rhs that's just the sum of the first two columns of our
        // symmetric matrix
        static std::vector <Real> sum_2(
            Natural const & size,
            Natural const & offset
        ) {
            auto A = Matrix <Real>::symmetric(size,offset);
            auto x = std::vector <double>(size);
            for(auto i=1;i<=size;i++) {
                if(i<=2)
                    x[i-1] = Real(1.);
                else
                    x[i-1] = Real(0.);
            }
            auto b = std::vector <double>(size);
            A.eval(x,b);
            return b;
        }

        // Create a vector of all ones 
        static std::vector <Real> ones(Natural & size) {
            auto b = std::vector <double>(size,Real(1.));
            return b;
        }

        // Create a vector where the last element is 1 and the others are 0 
        static std::vector <Real> elast(Natural & size) {
            auto b = std::vector <double>(size,Real(0.));
            b[size-1]=Real(1.);
            return b;
        }
    };

    // Turns off safeguarding
    template <typename Real,template <typename> class XX>
    Real no_safeguard(
        typename XX <Real>::Vector const & dx_base,
        typename XX <Real>::Vector const & dx_dir
    ) {
        return Real(1.0);
    }

    // Calculates a error between two vectors.  The second number returned 
    // is the norm of the second vector, which we use to calculate the relative
    // error.
    template <typename Real,template <typename> class XX>
    std::tuple <Real,Real> error(
        typename XX <Real>::Vector const & x,
        typename XX <Real>::Vector const & y
    ) {
        typedef XX <Real> X;
        auto r = X::init(x);
        X::copy(x,r);
        X::axpy(Real(-1.),y,r);
        auto norm_r = std::sqrt(X::innr(r,r));
        auto norm_y = std::sqrt(X::innr(y,y));
        return std::tuple <Real,Real>(norm_r,norm_y);
    }

    // Calculates the norm of the residual of a linear system.  This second
    // number returned is the norm of the rhs, which we use to calculate the
    // relative norm of the residual. 
    template <typename Real,template <typename> class XX>
    std::tuple<Real,Real> residual(
        Optizelle::Operator <Real,XX,XX> const & A,
        typename XX <Real>::Vector const & x,
        typename XX <Real>::Vector const & b,
        Optizelle::Operator <Real,XX,XX> const & B =
            typename Optizelle::Unconstrained <Real,XX>::Functions::Identity()
    ) {
        typedef XX <Real> X;
        auto A_x = X::init(x);
        A.eval(x,A_x);
        auto BA_x = X::init(x);
        B.eval(A_x,BA_x);
        auto B_b = X::init(x);
        B.eval(b,B_b);
        return error <Real,XX> (BA_x,B_b);
    }

    // Simple linear algebra solver setup
    template <typename Real,template <typename> class XX>
    struct Solver {
        // Type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Size 
        Natural m;

        // Stopping tolerance 
        Real eps;

        // Maximum number of iterations
        Natural iter_max;

        // Operator
        std::unique_ptr <Optizelle::Operator <Real,XX,XX>> A;

        // Right hand side
        std::unique_ptr <X_Vector> b;

        // Final solution
        std::unique_ptr <X_Vector> x;

        // Final number of iterations
        Natural iter;

        // Desired solution
        std::unique_ptr <X_Vector> x_star;

        // Desired number of iterations
        Natural iter_star;

        // Desired solution accuracy
        Real eps_sol;

        // Desired residual accuracy
        Real eps_res;

        // Iteration check
        bool check_iter;

        // Solution check
        bool check_sol;

        // Residual check 
        bool check_res;

        // Setup some simple parameters
        Solver() :
            m(5),
            eps(std::pow(
                Real(10.),
                std::log10(std::numeric_limits <Real>::epsilon())
                    *Real(0.75))),
            iter_max(500),
            A(),
            b(),
            x(),
            iter(0),
            x_star(),
            iter_star(0),
            eps_sol(std::pow(
                Real(10.),
                std::log10(std::numeric_limits <Real>::epsilon())
                    *Real(0.70))),
            eps_res(eps),
            check_iter(false),
            check_sol(false),
            check_res(false)
        {}
    };

    // Base GMRES setup 
    template <typename Real,template <typename> class XX>
    struct gmres : public Solver <Real,XX> {
        // Type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Restart frequency
        Natural rst_freq;

        // Preconditioners
        std::unique_ptr <Optizelle::Operator <Real,XX,XX>> B_left;
        std::unique_ptr <Optizelle::Operator <Real,XX,XX>> B_right;

        // GMRES manipulator
        std::unique_ptr <Optizelle::GMRESManipulator <Real,XX>> gmanip;

        // Setup some simple parameters
        gmres():
            Solver <Real,XX> (),
            rst_freq(0),
            B_left(new typename Optizelle::Unconstrained <Real,XX>::
                Functions::Identity()),
            B_right(new typename Optizelle::Unconstrained <Real,XX>::
                Functions::Identity()),
            gmanip(new Optizelle::EmptyGMRESManipulator <Real,XX>())
        {}
    };

    // Base TCG setup 
    template <typename Real,template <typename> class XX>
    struct tcg : public Solver <Real,XX> {
        // Type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Trust-region radius 
        Real delta; 

        // Preconditioners
        std::unique_ptr <Optizelle::Operator <Real,XX,XX>> B;

        // Offset 
        std::unique_ptr <X_Vector> x_offset;

        // Safeguard
        std::unique_ptr <Optizelle::SafeguardSimplified <Real,XX>> safeguard;

        // Maximum number of failed safeguard iterations
        Natural failed_max;

        // Number of vectors we orthogonalize against 
        Natural orthog_storage_max; 

        // Number of orthogonalization iterations
        Natural orthog_iter_max; 

        // Whether we do the projector check
        bool proj_check;

        // Cauchy-Point check
        bool check_cp;

        // Check that the length of the solution is the TR
        bool check_tr;

        // Tolerance for the TR check
        Real eps_tr;

        // Stopping condition check
        bool check_stop;

        // Desired stopping condition
        Optizelle::TruncatedStop::t stop_star;

        // Setup some simple parameters
        tcg():
            Solver <Real,XX> (),
            delta(std::numeric_limits <Real>::infinity()),
            B(new typename Optizelle::Unconstrained <Real,XX>::
                Functions::Identity()),
            x_offset(),
            safeguard(std::make_unique<Optizelle::SafeguardSimplified<Real,XX>>(
                no_safeguard <Real,XX>)),
            failed_max(std::numeric_limits <Natural>::max()),
            orthog_storage_max(1),
            orthog_iter_max(1),
            proj_check(false),
            check_cp(false),
            check_tr(false),
            eps_tr(this->eps),
            check_stop(false),
            stop_star(Optizelle::TruncatedStop::NotConverged)
        {}
    };

    // Run and verify the problem setup 
    template <typename Real,template <typename> class XX>
    void run_and_verify(gmres <Real,XX> & setup) {
        // Type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Fill in some pieces if they're missing
        if(!setup.x) {
            setup.x = std::make_unique <X_Vector> (X::init(*setup.b));
            X::zero(*setup.x);
        }

        // Solve this linear system
        auto err = Real(0.);
        std::tie(err,setup.iter) = Optizelle::gmres <Real,XX> (
            *setup.A,
            *setup.b,
            setup.eps,
            setup.iter_max,
            setup.rst_freq,
            *setup.B_left,
            *setup.B_right,
            *setup.gmanip,
            *setup.x);

        // Check that the number of iterations matches 
        if(setup.check_iter)
            CHECK(setup.iter == setup.iter_star);

        // Check the residual is less than our tolerance.  We do this two
        // different ways in case our solver is lying to us.
        if(setup.check_res) {
            CHECK(err < setup.eps_res);
            auto norm_r = Real(0.);
            auto norm_b = Real(0.);
            std::tie(norm_r,norm_b) =
                residual <Real,XX>(*setup.A,*setup.x,*setup.b);
            CHECK(norm_r <= setup.eps_res*norm_b);
            CHECK(std::fabs(norm_r-err) <= setup.eps_res);
        }

        // Check that the solution is correct
        if(setup.check_sol) {
            auto norm_r = Real(0.);
            auto norm_xstar = Real(0.);
            std::tie(norm_r,norm_xstar)=error <Real,XX>(*setup.x,*setup.x_star);
            CHECK(norm_r <= setup.eps_sol*norm_xstar);
        }
    }

    // Run and verify the problem setup 
    template <typename Real,template <typename> class XX>
    void run_and_verify(tcg <Real,XX> & setup) {
        // Type shortcuts
        typedef XX <Real> X;
        typedef typename X::Vector X_Vector;

        // Fill in some pieces if they're missing
        if(!setup.x) {
            setup.x = std::make_unique <X_Vector> (X::init(*setup.b));
            X::zero(*setup.x);
        }
        if(!setup.x_offset) {
            setup.x_offset = std::make_unique <X_Vector> (X::init(*setup.b));
            X::zero(*setup.x_offset);
        }

        // Allocate some necessary memory

        // Solve this linear system
        auto x_cp = X::init(*setup.b);
        auto err = Real(0.);
        auto norm_Br0 = Real(0.);
        auto norm_Br = Real(0.);
        auto stop = Optizelle::TruncatedStop::NotConverged;
        auto failed = Natural(0);
        auto alpha = Real(0.);
        Optizelle::truncated_cg <Real,XX> (
            *setup.A,
            *setup.b,
            *setup.B,
            setup.eps,
            setup.iter_max,
            setup.orthog_storage_max,
            setup.orthog_iter_max,
            setup.delta,
            *setup.x_offset,
            setup.proj_check,
            setup.failed_max,
            *setup.safeguard,
            *setup.x,
            x_cp,
            norm_Br0,
            norm_Br,
            setup.iter,
            stop,
            failed,
            alpha);

        // Check that the number of iterations matches 
        if(setup.check_iter)
            CHECK(setup.iter == setup.iter_star);

        // Check the preconditioned residual is less than our tolerance.  We do
        // this two different ways in case our solver is lying to us.
        if(setup.check_res) {
            CHECK(err < setup.eps_res);
            auto Bnorm_r = Real(0.);
            auto Bnorm_b = Real(0.);
            std::tie(Bnorm_r,Bnorm_b) =
                residual <Real,XX>(*setup.A,*setup.x,*setup.b,*setup.B);
            CHECK(Bnorm_r <= setup.eps_res*Bnorm_b);
            CHECK(std::fabs(Bnorm_r-err) <= setup.eps_res);
        }

        // Check that the solution is correct
        if(setup.check_sol) {
            auto norm_r = Real(0.);
            auto norm_xstar = Real(0.);
            std::tie(norm_r,norm_xstar)=error <Real,XX>(*setup.x,*setup.x_star);
            CHECK(norm_r <= setup.eps_sol*norm_xstar);
        }

        // Check that the solution is the Cauchy point
        if(setup.check_cp) {
            auto norm_r = Real(0.);
            auto norm_xcp = Real(0.);
            std::tie(norm_r,norm_xcp)=error <Real,XX>(*setup.x,x_cp);
            CHECK(norm_r <= setup.eps_sol*norm_xcp);
        }

        // Check that the length of the solution plus the offset is the size
        // of the TR
        if(setup.check_tr) {
            auto x_p_offset = X::init(*setup.x);
            X::copy(*setup.x,x_p_offset);
            X::axpy(Real(1.),*setup.x_offset,x_p_offset);
            auto norm_xpoffset = std::sqrt(X::innr(x_p_offset,x_p_offset));
            CHECK(std::abs(norm_xpoffset-setup.delta) <=
                setup.eps_tr * setup.delta);
        }

        // Check that the stopping condition matches what we expect
        if(setup.check_stop) {
            CHECK(setup.stop_star == stop);
        }
    }

    // Reset the different checks 
    template <typename Real,template <typename> class XX>
    void reset_checks(Solver <Real,XX> & setup) {
        setup.check_iter = false;
        setup.check_sol = false;
        setup.check_res = false;
    }

    template <typename Real,template <typename> class XX>
    void reset_checks(gmres <Real,XX> & setup) {
        reset_checks(dynamic_cast <Solver <Real,XX> &>(setup));
    }

    template <typename Real,template <typename> class XX>
    void reset_checks(tcg <Real,XX> & setup) {
        reset_checks(dynamic_cast <Solver <Real,XX> &>(setup));
        setup.check_cp = false;
        setup.check_tr = false;
        setup.check_stop = false;
    }
}
