// A unit test with three different vector spaces.  Mostly, we're using this
// unit to catch errors where we mix vector spaces t##hat the compiler didn't
// catch.  Our strategy for this is to make sure t##hat each vector space has
// incompatible functions defined on the vector itself and then to call these
// functions within each of the functions.  This will not run, but it should
// compile if we're type clean.

#include "optizelle/optizelle.h"
#include "optizelle/json.h"

// Base part of a vector space including the vector definition
#define BASE(v) \
    struct Vector { \
        void check_##v() const {} \
        Vector() = default; \
        Vector(Vector && x) = default; \
        Vector & operator = (Vector && x) = default; \
        Vector(Vector const & x) = delete; \
        Vector & operator = (Vector const & x) = delete; \
    }; \
    static Vector init(Vector const & x) { \
        x.check_##v(); \
        throw; \
    } \
    static void copy(Vector const & x, Vector & y) { \
        x.check_##v(); \
        y.check_##v(); \
        throw; \
    } \
    static void scal(Real const & alpha, Vector & x) { \
        x.check_##v(); \
        throw; \
    } \
    static void zero(Vector & x) { \
        x.check_##v(); \
        throw; \
    } \
    static void axpy(Real const & alpha, Vector const & x, Vector & y) { \
        x.check_##v(); \
        y.check_##v(); \
        throw; \
    } \
    static Real innr(Vector const & x,Vector const & y) { \
        x.check_##v(); \
        y.check_##v(); \
        throw; \
    } \
    static void rand(Vector & x){ \
        x.check_##v(); \
        throw; \
    }

// Inequality algebra pieces of a vector space 
#define INEQ(v) \
    static void prod(Vector const & x, Vector const & y, Vector & z){ \
        x.check_##v(); \
        y.check_##v(); \
        z.check_##v(); \
        throw; \
    } \
    static void id(Vector & x) { \
        x.check_##v(); \
        throw; \
    } \
    static void linv(Vector const & x,Vector const & y,Vector & z) { \
        x.check_##v(); \
        y.check_##v(); \
        z.check_##v(); \
        throw; \
    } \
    static Real barr(Vector const & x) { \
        x.check_##v(); \
        throw; \
    } \
    static Real srch(Vector const & x,Vector const & y) { \
        x.check_##v(); \
        y.check_##v(); \
        throw; \
    } \
    static void symm(Vector & x) { \
        x.check_##v(); \
        throw; \
    }

// Definitions for a scalar valued function
#define SCALAR_FN(V,v) \
    typedef VS_##V <Real> V; \
    typedef typename V::Vector V##_Vector; \
    Real eval(V##_Vector const & x) const { \
        v.check_##v(); \
        throw; \
    } \
    void grad( \
        V##_Vector const & v, \
        V##_Vector & grad \
    ) const { \
        v.check_##v(); \
        grad.check_##v(); \
        throw; \
    } \
    void hessvec( \
        V##_Vector const & v, \
        V##_Vector const & d##v, \
        V##_Vector & H_d##v \
    ) const { \
        v.check_##v(); \
        d##v.check_##v(); \
        H_d##v.check_##v(); \
        throw; \
    }

// Definitions for a vector valued function
#define VECTOR_FN(V,v,W,w) \
    typedef VS_##V <Real> V; \
    typedef typename V::Vector V##_Vector; \
    typedef VS_##W <Real> W; \
    typedef typename W::Vector W##_Vector; \
    void eval( \
        V##_Vector const & v, \
        W##_Vector & w \
    ) const { \
        v.check_##v(); \
        w.check_##w(); \
        throw; \
    } \
    void p( \
        V##_Vector const & v, \
        V##_Vector const & d##v, \
        W##_Vector & w \
    ) const { \
        v.check_##v(); \
        d##v.check_##v(); \
        w.check_##w(); \
        throw; \
    } \
    void ps( \
        V##_Vector const & v, \
        W##_Vector const & d##w, \
        V##_Vector & v##hat  \
    ) const { \
        v.check_##v(); \
        d##w.check_##w(); \
        v##hat.check_##v(); \
        throw; \
    } \
    void pps( \
        V##_Vector const & v, \
        V##_Vector const & d##v, \
        W##_Vector const & d##w, \
        V##_Vector & v##hat  \
    ) const { \
        v.check_##v(); \
        d##v.check_##v(); \
        d##w.check_##w(); \
        v##hat.check_##v(); \
        throw; \
    }
#define OP_FN(V,v) \
    typedef VS_##V <Real> V; \
    typedef typename V::Vector V##_Vector; \
    void eval(V##_Vector const & d##v,V##_Vector & v##hat) const { \
        d##v.check_##v(); \
        v##hat.check_##v(); \
    }

// Three different vector spaces, which should be incompatible 
template <typename Real>
struct VS_X { 
    BASE(x)
};
template <typename Real>
struct VS_Y { 
    BASE(y)
};
template <typename Real>
struct VS_Z { 
    BASE(z)
    INEQ(z)
};

// Define functions on each vector space
template <typename Real>
struct F : public Optizelle::ScalarValuedFunction <Real,VS_X> {
    SCALAR_FN(X,x)
};
template <typename Real>
struct G : public Optizelle::VectorValuedFunction<Real,VS_X,VS_Y> {
    VECTOR_FN(X,x,Y,y);
};
template <typename Real>
struct H : public Optizelle::VectorValuedFunction<Real,VS_X,VS_Z> {
    VECTOR_FN(X,x,Z,z);
};

// Define all the preconditioners 
template <typename Real>
struct PH : public Optizelle::Operator <Real,VS_X,VS_X> {
    OP_FN(X,x)
};
template <typename Real>
struct PSchurLeft : public Optizelle::Operator <Real,VS_Y,VS_Y> {
    OP_FN(Y,y)
};
template <typename Real>
struct PSchurRight : public Optizelle::Operator <Real,VS_Y,VS_Y> {
    OP_FN(Y,y)
};

// Put together an optimization problem that uses the above vector spaces and
// functions
int main(int argc,char * argv[]) {
    // Set our spaces 
    typedef double Real;
    typedef VS_X <Real>::Vector X_Vector;
    typedef VS_Y <Real>::Vector Y_Vector;
    typedef VS_Z <Real>::Vector Z_Vector;

    // Don't actually run anything
    if(argc==1) return EXIT_SUCCESS;

    // Create some initial vectors 
    auto x = X_Vector(); 
    auto y = Y_Vector(); 
    auto z = Z_Vector(); 

    // Unfortunately, we need to do this twice since there are some different
    // code paths for unconstrained and constrained problems.  For example, on
    // how they deal with preconditioners.  Inequalities are somewhat bolted
    // on, so we really just need to check with equality constraints and
    // without.

    // Unconstrained
    {
        // Create a state
        Optizelle::Unconstrained <Real,VS_X>::State::t state(x);
      
        // Read a restart file
        Optizelle::json::Unconstrained <Real,VS_X>::read_restart(
            "junk.json",x,state);

        // Read additional parameters from file
        Optizelle::json::Unconstrained <Real,VS_X>::read("junk.json",state);

        // Create the bundle of functions 
        Optizelle::Constrained <Real,VS_X,VS_Y,VS_Z>::Functions::t fns;
        fns.f.reset(new F <Real>);
        fns.PH.reset(new PH <Real>);
        
        // Solve the optimization problem
        Optizelle::Unconstrained <Real,VS_X>::Algorithms
            ::getMin(Optizelle::Messaging::stdout,fns,state);

        // Write out the final answer to file
        Optizelle::json::Unconstrained <Real,VS_X>::write_restart(
            "junk.json",state);

    }

    // Constrained
    {
        // Create a state
        Optizelle::Constrained <Real,VS_X,VS_Y,VS_Z>::State::t state(x,y,z);
      
        // Read a restart file
        Optizelle::json::Constrained <Real,VS_X,VS_Y,VS_Z>::read_restart(
            "junk.json",x,y,z,state);

        // Read additional parameters from file
        Optizelle::json::Constrained <Real,VS_X,VS_Y,VS_Z>
            ::read("junk.json",state);

        // Create the bundle of functions 
        Optizelle::Constrained <Real,VS_X,VS_Y,VS_Z>::Functions::t fns;
        fns.f.reset(new F <Real>);
        fns.g.reset(new G <Real>);
        fns.h.reset(new H <Real>);
        fns.PH.reset(new PH <Real>);
        fns.PSchur_left.reset(new PSchurLeft <Real>);
        fns.PSchur_right.reset(new PSchurRight <Real>);
        
        // Solve the optimization problem
        Optizelle::Constrained <Real,VS_X,VS_Y,VS_Z>::Algorithms
            ::getMin(Optizelle::Messaging::stdout,fns,state);

        // Write out the final answer to file
        Optizelle::json::Constrained <Real,VS_X,VS_Y,VS_Z>::write_restart(
            "junk.json",state);
    }

    // Return success
    return EXIT_SUCCESS;
}
