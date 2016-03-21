#pragma once
// Supporting functions for testing the augmented systems solves in their
// various forms

#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"
#include "unit.h"

// Grab the squaring function
using Optizelle::sq;

// JustA lump everything in the unit structure 
template <typename Real>
struct Unit {
    // Set some type shortcuts
    template <typename Real_>
    using XX = Optizelle::Rm <Real_>;
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector;

    template <typename Real_>
    using YY = Optizelle::Rm <Real_>;
    typedef YY <Real> Y;
    typedef typename Y::Vector Y_Vector;

    template <typename Real_>
    using ZZ = Optizelle::Rm <Real_>;
    typedef ZZ <Real> Z;
    typedef typename Z::Vector Z_Vector;

    typedef Optizelle::Natural Natural;

    // Various objective functions
    struct Objective {
        // Zero objective function
        struct Zero : public Optizelle::ScalarValuedFunction <Real,XX> {
            Real eval(X_Vector const & x) const {
                return Real(0.);
            }
            void grad(
                X_Vector const & x,
                X_Vector & grad
            ) const {
                X::zero(grad);
            }
            void hessvec(
                X_Vector const & x,
                X_Vector const & dx,
                X_Vector & H_dx
            ) const {
                X::zero(H_dx);
            }
        };
    };

    // Various constraints
    struct Constraint {
        // Define a linear constraint
        //
        // g(x,y)= [ x + y = 1 ] 
        //
        struct Linear : public Optizelle::VectorValuedFunction<Real,XX,YY> {
            // y=g(x) 
            void eval(
                X_Vector const & x,
                Y_Vector & y
            ) const {
                y[0] = x[0] + x[1] - Real(1.); 
            }

            // y=g'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector & y
            ) const {
                y[0] = dx[0] + dx[1];
            }

            // x_hat=g'(x)*dy
            void ps(
                X_Vector const & x,
                Y_Vector const & dy,
                X_Vector & x_hat
            ) const {
                x_hat[0]=dy[0];
                x_hat[1]=dy[0];
            }

            // x_hat=(g''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector const & dy,
                X_Vector & x_hat 
            ) const {
                X::zero(x_hat);
            }
        };
        
        // Define a quadratic constraint
        //
        // g(x,y)= [ a x^2 + b y^2 = 1 ] 
        //
        struct Quadratic : public Optizelle::VectorValuedFunction<Real,XX,YY> {
            // Determines the width of the ellipse
            Real a;

            // Determines the height of the ellpise
            Real b;

            // Grabs the size of the quadratic
            Quadratic(Real const & a_, Real const & b_) : a(a_), b(b_) {}

            // y=g(x) 
            void eval(
                X_Vector const & x,
                Y_Vector & y
            ) const {
                y[0] = a * sq(x[0]) + b * sq(x[1]) - Real(1.); 
            }

            // y=g'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector & y
            ) const {
                y[0] = Real(2.)*a*x[0]*dx[0] + Real(2.)*b*x[1]*dx[1];
            }

            // x_hat=g'(x)*dy
            void ps(
                X_Vector const & x,
                Y_Vector const & dy,
                X_Vector & x_hat
            ) const {
                x_hat[0]=Real(2.)*a*x[0]*dy[0];
                x_hat[1]=Real(2.)*b*x[1]*dy[0];
            }

            // x_hat=(g''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector const & dy,
                X_Vector & x_hat 
            ) const {
                x_hat[0]=Real(2.)*a*dx[0]*dy[0];
                x_hat[1]=Real(2.)*b*dx[1]*dy[0];
            }
        };
        
        // Define a constraint where we intersect two circles 
        //
        // g(x,y)= [ (x-a)^2 + (y-b)^2 = 1 ] 
        //         [ (x-c)^2 + (y-d)^2 = 1 ] 
        //
        struct CircleIntersection :
            public Optizelle::VectorValuedFunction<Real,XX,YY>
        {
            // Location of the first circle 
            Real a;
            Real b;

            // Location of the second circle 
            Real c;
            Real d;

            // Grabs the size of the quadratic
            CircleIntersection(
                Real const & a_, 
                Real const & b_,
                Real const & c_,
                Real const & d_
            ) : a(a_), b(b_), c(c_), d(d_) {}

            // y=g(x) 
            void eval(
                X_Vector const & x,
                Y_Vector & y
            ) const {
                y[0] = sq(x[0]-a) + sq(x[1]-b) - Real(1.); 
                y[1] = sq(x[0]-c) + sq(x[1]-d) - Real(1.); 
            }

            // y=g'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector & y
            ) const {
                y[0] = Real(2.)*(x[0]-a)*dx[0] + Real(2.)*(x[1]-b)*dx[1];
                y[1] = Real(2.)*(x[0]-c)*dx[0] + Real(2.)*(x[1]-d)*dx[1];
            }

            // x_hat=g'(x)*dy
            void ps(
                X_Vector const & x,
                Y_Vector const & dy,
                X_Vector & x_hat
            ) const {
                x_hat[0] = Real(2.)*(x[0]-a)*dy[0] + Real(2.)*(x[0]-c)*dy[1];
                x_hat[1] = Real(2.)*(x[1]-b)*dy[0] + Real(2.)*(x[1]-d)*dy[1];
            }

            // x_hat=(g''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector const & dy,
                X_Vector & x_hat 
            ) const {
                x_hat[0] = Real(2.)*dx[0]*dy[0] + Real(2.)*dx[0]*dy[1];
                x_hat[1] = Real(2.)*dx[1]*dy[0] + Real(2.)*dx[1]*dy[1];
            }
        };
        
        // Define a trigonometric constraint 
        //
        // g(x) = cos(x)
        //
        struct Cos : public Optizelle::VectorValuedFunction<Real,XX,YY> {
            // y=g(x) 
            void eval(
                X_Vector const & x,
                Y_Vector & y
            ) const {
                y[0] = std::cos(x[0]); 
            }

            // y=g'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector & y
            ) const {
                y[0] = -std::sin(x[0])*dx[0];
            }

            // x_hat=g'(x)*dy
            void ps(
                X_Vector const & x,
                Y_Vector const & dy,
                X_Vector & x_hat
            ) const {
                x_hat[0] = -std::sin(x[0])*dy[0];
            }

            // x_hat=(g''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector const & dy,
                X_Vector & x_hat 
            ) const {
                x_hat[0] = -std::cos(x[0])*dx[0]*dy[0];
            }
        };
        
        // Define a polynomial constraint 
        //
        // g(x) = (x-1)*(x-2)*(x-3)*(x-4) 
        //
        struct Poly : public Optizelle::VectorValuedFunction<Real,XX,YY> {
            // y=g(x) 
            void eval(
                X_Vector const & x,
                Y_Vector & y
            ) const {
                y[0] = (x[0]-Real(1.))*(x[0]-Real(2.))*
                       (x[0]-Real(3.))*(x[0]-Real(4.));
            }

            // y=g'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector & y
            ) const {
                y[0] = (Real(2.)*(Real(2.)*x[0]-Real(5.))*
                           (sq(x[0])-Real(5.)*x[0]+Real(5.)))*dx[0]; 
            }

            // x_hat=g'(x)*dy
            void ps(
                X_Vector const & x,
                Y_Vector const & dy,
                X_Vector & x_hat
            ) const {
                x_hat[0] = (Real(2.)*(Real(2.)*x[0]-Real(5.))*
                               (sq(x[0])-Real(5.)*x[0]+Real(5.)))*dy[0]; 
            }

            // x_hat=(g''(x)dx)*dy
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Y_Vector const & dy,
                X_Vector & x_hat 
            ) const {
                x_hat[0] = Real(2.)*(Real(6.)*sq(x[0])-Real(30.)*x[0]+Real(35.))
                    *dy[0]*dx[0]; 
            }
        };
        
        // Define a box constraint, which we use for inequality constraints
        //
        // h(x) = [ x - lb ] 
        //      = [ ub - x ]
        //
        struct Box: public Optizelle::VectorValuedFunction<Real,XX,YY> {
            // Lower bound
            X_Vector lb;
            
            // Upper bound
            X_Vector ub;

            // Number of variables
            Natural const m;

            // Grab the bounds
            Box(X_Vector const & lb_,X_Vector const & ub_) :
                lb(X::init(lb_)),
                ub(X::init(ub_)),
                m(lb_.size())
            {
                X::copy(lb_,lb);
                X::copy(ub_,ub);
            }

            // z=h(x) 
            void eval(
                X_Vector const & x,
                Z_Vector & z
            ) const {
                for(Natural i=0;i<m;i++) {
                    z[i]=x[i]-lb[i];
                    z[i+m]=ub[i]-x[i];
                }
            }

            // z=h'(x)dx
            void p(
                X_Vector const & x,
                X_Vector const & dx,
                Z_Vector & z
            ) const {
                for(Natural i=0;i<m;i++) {
                    z[i]=dx[i];
                    z[i+m]=-dx[i];
                }
            }

            // x_hat=h'(x)*dz
            void ps(
                X_Vector const & x,
                Z_Vector const & dz,
                X_Vector & x_hat
            ) const {
                for(Natural i=0;i<m;i++)
                    x_hat[i]=dz[i]-dz[i+m];
            }

            // x_hat=(h''(x)dx)*dz
            void pps(
                X_Vector const & x,
                X_Vector const & dx,
                Z_Vector const & dz,
                X_Vector & x_hat 
            ) const {
                X::zero(x_hat);
            }
        };
    };

    // Generic elements for augmented system solves
    struct Augsys {
        // Tolerance for all our tests
        Real eps;

        // Points where we run the test 
        X_Vector x;
        Y_Vector y;

        // Equality constraint 
        std::unique_ptr <Optizelle::VectorValuedFunction <Real,XX,YY>> g;

        // Trust-region radius
        Real delta;

        // Do diagonstic checks instead
        bool do_diagnostics;

        // Check that we computed augmented system solve iterations 
        bool check_augsys; 

        // Check that we exited the augmented system solve early
        bool check_augsys_exit;
        
        // Setup some simple parameters
        Augsys(X_Vector const & x_,Y_Vector const & y_) :
            eps(std::pow(
                Real(10.),
                std::log10(std::numeric_limits <Real>::epsilon())
                    *Real(0.75))),
            x(X::init(x_)),
            y(Y::init(y_)),
            g(),
            delta(1e16),
            do_diagnostics(false),
            check_augsys(false),
            check_augsys_exit(false)
        {
            X::copy(x_,x); 
            Y::copy(y_,y);
        }
    };

    // Simple quasinormal subproblem setup 
    struct QN : public Augsys {

        // Points where we run the test 
        Z_Vector z;

        // Inequality constraint
        std::unique_ptr <Optizelle::VectorValuedFunction <Real,XX,YY>> h;

        // Desired solution
        std::unique_ptr <X_Vector> dx_ncp_star;
        std::unique_ptr <X_Vector> dx_n_star;

        // Desired stopping condition 
        Optizelle::QuasinormalStop::t qn_stop_star;

        // Hard coded zeta for reference
        Real const zeta;

        // Hard coded gamma for reference
        Real const gamma;

        // Check the stopping condition
        bool check_stop;

        // Check that the trust-region is active 
        bool check_tr;

        // Check the quasinormal step
        bool check_dx_n;

        // Check the quasinormal cauchy point
        bool check_dx_ncp;

        // Check the feasibility
        bool check_feas;

        // Check that we cut back the step from the safeguard
        bool check_safe;

        // Copy of the Cauchy point for some more complicated tests 
        X_Vector cp;

        // Setup some simple parameters
        QN(X_Vector const & x_,Y_Vector const & y_) :
            Augsys(x_,y_),
            z(),
            h(),
            dx_ncp_star(),
            dx_n_star(),
            qn_stop_star(),
            zeta(0.8),
            gamma(0.99),
            check_stop(false),
            check_tr(false),
            check_dx_n(false),
            check_dx_ncp(false),
            check_feas(false),
            check_safe(false),
            cp(X::init(x_))
        {
            z.resize(x_.size()*2);
        }
    };

    // Calculates a error between two vectors.  The second number returned 
    // is the norm of the second vector, which we use to calculate the relative
    // error.
    static std::tuple <Real,Real> error(
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

    // Run and verify the problem setup
    static void run_and_verify(QN & setup) {
        // Create an optimization state
        typename Optizelle::Constrained <Real,XX,YY,ZZ>::State::t
            state(setup.x,setup.y,setup.z);

        // Set some appropriate state information
        state.delta = setup.delta;

        // Create a bundle of functions
        typename Optizelle::Constrained <Real,XX,YY,ZZ>::Functions::t fns;
        fns.f.reset(new typename Objective::Zero);
        fns.g = std::move(setup.g);
        if(setup.h)
            fns.h = std::move(setup.h);
        else {
            auto lb = X::init(setup.x);
            auto ub = X::init(setup.x);
            X::id(lb);
            X::id(ub);
            X::scal(Real(-1e6),lb);
            X::scal(Real(1e6),ub);
            fns.h.reset(new typename Constraint::Box(lb,ub));
        }

        // Create the messaging object
        auto msg = Optizelle::Messaging();

        // Fill out the bundle of functions
        Optizelle::Constrained <Real,XX,YY,ZZ>::Functions::init(
            msg,state,fns);

        // Evaluate the function and cache information about it
        fns.g->eval(state.x,state.g_x);
        state.norm_gxtyp = std::sqrt(Y::innr(state.g_x,state.g_x));
        auto gps_g = X::init(state.x);
        fns.g->ps(state.x,state.g_x,gps_g);
        state.norm_gpsgxtyp = std::sqrt(X::innr(gps_g,gps_g));
        fns.g->eval(state.x,state.g_x);
        fns.h->eval(state.x,state.h_x);

        // If we're doing diagnostics, run them and then quit.  Really, we
        // should just be using this for debugging our problem setups.
        if(setup.do_diagnostics) {
            state.f_diag = Optizelle::FunctionDiagnostics::SecondOrder;
            state.g_diag = Optizelle::FunctionDiagnostics::SecondOrder;
            state.h_diag = Optizelle::FunctionDiagnostics::SecondOrder;
            Optizelle::Constrained<Real,XX,YY,ZZ>::Diagnostics
                ::checkFunctions(msg,fns,state);
            return;
        }

        // Compute the quasinormal step
        Optizelle::EqualityConstrained<Real,XX,YY>::Algorithms::quasinormalStep(
            fns,state);

        // Check the stopping condition 
        if(setup.check_stop)
            CHECK(state.qn_stop == setup.qn_stop_star);
        
        // Check that we hit the trust-region modified by zeta 
        if(setup.check_tr) {
            auto norm_dxn = std::sqrt(X::innr(state.dx_n,state.dx_n));
            CHECK(std::abs(norm_dxn-state.delta*state.zeta) <= setup.eps);
        }

        // Check that the quasinormal point is correct
        if(setup.check_dx_n) {
            auto norm_r = Real(0.);
            auto norm_dxn_star = Real(0.);
            std::tie(norm_r,norm_dxn_star)=error(state.dx_n,*setup.dx_n_star);
            CHECK(norm_r <= setup.eps);
        }

        // Check that the cauchy point is correct
        if(setup.check_dx_ncp) {
            auto norm_r = Real(0.);
            auto norm_dxncp_star = Real(0.);
            std::tie(norm_r,norm_dxncp_star)=error(
                state.dx_ncp,*setup.dx_ncp_star);
            CHECK(norm_r <= setup.eps);
        }

        // Check that our step takes us to feasibility 
        if(setup.check_feas) {
            auto x_p_dxn = X::init(state.x);
            X::copy(state.x,x_p_dxn);
            X::axpy(Real(1.),state.dx_n,x_p_dxn);
            auto g_x = Y::init(state.y);
            fns.g->eval(x_p_dxn,g_x);
            auto norm_gx = std::sqrt(X::innr(g_x,g_x));
            CHECK(norm_gx <= setup.eps);
        }

        // Check that the safeguard truncated the step 
        if(setup.check_safe) {
            CHECK(state.alpha_x_qn < Real(1.));
        } else {
            CHECK(state.alpha_x_qn == Real(1.));
        }

        // Check that we're right on the fraction to the boundary rule
        {
            // Find h(x+dx_n)
            auto x_p_dx = X::init(state.x);
            X::copy(state.x,x_p_dx);
            X::axpy(Real(1.),state.dx_n,x_p_dx);
            auto h_xpdx = Z::init(state.z);
            fns.h->eval(x_p_dx,h_xpdx);

            // Find -h_x
            auto m_h_x = Z::init(state.h_x);
            Z::copy(state.h_x,m_h_x);
            Z::scal(Real(-1.),m_h_x);

            // Search from h(x + dx_n) to -h(x).  The fraction to the
            // boundary rule states that h(x+dx) >= (1-gamma*zeta) h(x).
            auto alpha = Z::srch(m_h_x,h_xpdx);

            // When we safeguard, we should hit the fraction to the boundary
            // exactly
            if(setup.check_safe) {
                CHECK(std::abs(alpha-(Real(1.)-state.gamma*state.zeta))
                    <= setup.eps);

            // When we don't safeguard, we should not be on the fraction to the
            // boundary
            } else {
                CHECK(std::abs(alpha-(Real(1.)-state.gamma*state.zeta)) > 0.1);
            }
        }

        // Check that we computed augmented system solves 
        if(setup.check_augsys) {
            CHECK(state.augsys_qn_iter > 0);

        // Otherwise, make sure that we didn't
        } else {
            CHECK(state.augsys_qn_iter == 0);
        }

        // Move the function out the bundle of functions and back into setup in
        // case we need it later
        setup.g = std::move(fns.g);
        setup.h = std::move(fns.h);

        // Copy the Cauchy point
        X::copy(state.dx_ncp,setup.cp); 
    }

    // Simple nullspace projection setup 
    struct Proj : public Augsys {
        // Direction to project
        std::unique_ptr <X_Vector> dx;

        // Desired solution
        std::unique_ptr <X_Vector> P_dx_star;

        // Check that we're in the nullspace of the constraint
        bool check_null;

        // Check that we solution that we obtained
        bool check_sol;

        // Setup some simple parameters
        Proj(X_Vector const & x_,Y_Vector const & y_) :
            Augsys(x_,y_),
            dx(),
            P_dx_star(),
            check_null(false),
            check_sol(false)
        {}
    };

    // Run and verify the problem setup
    static void run_and_verify(Proj & setup) {
        // Create an optimization state
        typename Optizelle::EqualityConstrained <Real,XX,YY>::State::t
            state(setup.x,setup.y);

        // Create a bundle of functions
        typename Optizelle::Constrained <Real,XX,YY,ZZ>::Functions::t fns;
        fns.f.reset(new typename Objective::Zero);
        fns.g = std::move(setup.g);

        // Create the messaging object
        auto msg = Optizelle::Messaging();

        // Fill out the bundle of functions
        Optizelle::EqualityConstrained <Real,XX,YY>::Functions::init(
            msg,state,fns);

        // Evaluate the function and cache information about it
        fns.g->eval(state.x,state.g_x);
        state.norm_gxtyp = std::sqrt(Y::innr(state.g_x,state.g_x));
        fns.g->eval(state.x,state.g_x);
        state.norm_dxtyp = Real(1.);

        // If we're doing diagnostics, run them and then quit.  Really, we
        // should just be using this for debugging our problem setups.
        if(setup.do_diagnostics) {
            state.f_diag = Optizelle::FunctionDiagnostics::SecondOrder;
            state.g_diag = Optizelle::FunctionDiagnostics::SecondOrder;
            Optizelle::EqualityConstrained<Real,XX,YY>::Diagnostics
                ::checkFunctions(msg,fns,state);
            return;
        }

        // Compute the nullspace projection
        auto P_dx = X::init(setup.x);
        typename Optizelle::EqualityConstrained<Real,XX,YY>::Algorithms
            ::NullspaceProjForTrunc(state,fns).eval(*setup.dx,P_dx);

        // Check that the quasinormal point is correct
        if(setup.check_sol) {
            auto norm_r = Real(0.);
            auto norm_dx = Real(0.);
            std::tie(norm_r,norm_dx)=error(P_dx,*setup.P_dx_star);
            CHECK(norm_r <= setup.eps);
        }

        // Check that our step takes us into the nullspace 
        if(setup.check_null) {
            auto gp_x_Pdx = Y::init(setup.y);
            fns.g->p(setup.x,P_dx,gp_x_Pdx);
            auto norm_gpxPdx = std::sqrt(Y::innr(gp_x_Pdx,gp_x_Pdx));
            CHECK(norm_gpxPdx <= setup.eps);
        }

        // Check that we computed augmented system solves 
        if(setup.check_augsys) {
            CHECK(state.augsys_proj_iter > 0);

        // Otherwise, make sure that we didn't
        } else {
            CHECK(state.augsys_proj_iter == 0);
        }
        
        // Check that we exited the augmented system solve early
        if(setup.check_augsys_exit) {
            CHECK(state.augsys_proj_err_target == Real(1.));
        }

        // Move the function out the bundle of functions and back into setup in
        // case we need it later
        setup.g = std::move(fns.g);
    }
};
