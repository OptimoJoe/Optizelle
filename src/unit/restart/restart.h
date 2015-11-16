#ifndef RESTART_H
#define RESTART_H

#include "optizelle/optizelle.h"
#include "optizelle/json.h"

// Create some type shortcuts
typedef double Real;
template <typename Real> using XX = Optizelle::Rm <Real>;
template <typename Real> using YY = Optizelle::Rm <Real>;
template <typename Real> using ZZ = Optizelle::Rm <Real>;

// Create a blank state mainpulator
template <typename ProblemClass>
struct BlankManipulator : Optizelle::StateManipulator <ProblemClass> {
    void eval(
        typename ProblemClass::Functions::t const & fns,
        typename ProblemClass::State::t & state,
        Optizelle::OptimizationLocation::t const & loc
    ) const { } 
};

// Create some empty functions
struct F : public Optizelle::ScalarValuedFunction <Real,XX> {
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector Vector;

    // <- f(x) 
    Real eval(Vector const & x) const { throw; }

    // grad = grad f(x) 
    void grad(Vector const & x,Vector & grad) const {}

    // H_dx = hess f(x) dx 
    void hessvec(Vector const & x,Vector const & dx,Vector & H_dx) const {}
};

struct G : public Optizelle::VectorValuedFunction <Real,XX,YY> {
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector;
    typedef YY <Real> Y;
    typedef typename Y::Vector Y_Vector;

    // y=f(x)
    void eval(X_Vector const & x,Y_Vector & y) const {}

     // y=f'(x)dx 
     void p(
         X_Vector const & x,
         X_Vector const & dx,
         Y_Vector & y
     ) const {}

     // z=f'(x)*dy
     void ps(
         X_Vector const & x,
         Y_Vector const & dy,
         X_Vector & z
     ) const {}

     // z=(f''(x)dx)*dy
     void pps(
         X_Vector const & x,
         X_Vector const & dx,
         Y_Vector const & dy,
         X_Vector & z
     ) const {} 
};

struct H : public Optizelle::VectorValuedFunction <Real,XX,ZZ> {
    // Create some type shortcuts
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector;
    typedef ZZ <Real> Z;
    typedef typename Z::Vector Z_Vector;

    // y=f(x)
    void eval(X_Vector const & x,Z_Vector & z) const {}

     // z=f'(x)dx
     void p(
         X_Vector const & x,
         X_Vector const & dx,
         Z_Vector & z
     ) const {}

     // z=f'(x)*dz
     void ps(
         X_Vector const & x,
         Z_Vector const & dz,
         X_Vector & z
     ) const {}

     // z=(f''(x)dx)*dz
     void pps(
         X_Vector const & x,
         X_Vector const & dx,
         Z_Vector const & dz,
         X_Vector & z
     ) const {} 
};

#endif
