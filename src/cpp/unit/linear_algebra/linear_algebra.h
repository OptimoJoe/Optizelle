#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include "peopt/peopt.h"
#include "peopt/vspaces.h"

// Functions and definitions used by all the linear algebra tests 
typedef peopt::Natural Natural;

template <typename Real>
struct BasicOperator : public peopt::Operator <Real,peopt::Rm,peopt::Rm> {
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
    void operator () (const X_Vector& x,Y_Vector &y) const  {
        for(Natural i=0;i<m;i++){
            y[i]=0;
            for(Natural j=0;j<m;j++)
                y[i]+=A[i+m*j]*x[j];
        }
    }
};

// Create the identity operator
template <typename Real>
struct IdentityOperator : public peopt::Operator <Real,peopt::Rm,peopt::Rm> {
private:
    // Create some type shortcuts
    typedef std::vector <Real> X_Vector;
    typedef std::vector <Real> Y_Vector;

public:
    // We don't need to do anything on creation
    IdentityOperator() {};

    // Just copy the input to the output 
    void operator () (const X_Vector& x,Y_Vector &y) const  {
        peopt::Rm <Real>::copy(x,y);
    }
};

#endif
