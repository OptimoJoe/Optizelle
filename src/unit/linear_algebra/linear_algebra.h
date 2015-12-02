#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include "optizelle/optizelle.h"
#include "optizelle/vspaces.h"

// Functions and definitions used by all the linear algebra tests 
typedef Optizelle::Natural Natural;

template <typename Real>
struct BasicOperator :
    public Optizelle::Operator <Real,Optizelle::Rm,Optizelle::Rm>
{
private:
    // Create some type shortcuts
    typedef std::vector <Real> X_Vector;
    typedef std::vector <Real> Y_Vector;
    
    // Size of the matrix
    Natural m;
public:
    // Storage for the matrix
    std::vector <Real> A;

    // Create an empty matrix, which must be filled by the user 
    BasicOperator(const Natural m_) : m(m_) {
        A.resize(m*m); 
    }
    
    // Apply the random matrix to the vector 
    void eval(const X_Vector& x,Y_Vector &y) const  {
        for(Natural i=0;i<m;i++){
            y[i]=0;
            for(Natural j=0;j<m;j++)
                y[i]+=A[i+m*j]*x[j];
        }
    }
};

// Create the identity operator
template <typename Real>
struct IdentityOperator :
    public Optizelle::Operator <Real,Optizelle::Rm,Optizelle::Rm>
{
private:
    // Create some type shortcuts
    typedef std::vector <Real> X_Vector;
    typedef std::vector <Real> Y_Vector;

public:
    // We don't need to do anything on creation
    IdentityOperator() {};

    // Just copy the input to the output 
    void eval(const X_Vector& x,Y_Vector &y) const  {
        Optizelle::Rm <Real>::copy(x,y);
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

#endif
