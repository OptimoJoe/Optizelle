#ifndef DIFF_OPERATOR_H
#define DIFF_OPERATOR_H

#include "core.h"

// Operator possessing two derivatives and adjoints of those derivatives 
template <class Domain, class Codomain>
class DiffOperator {
public:
    // Basic application
    virtual void operator () (const int i,const Domain& x,Codomain &y) const=0;

    // Derivative of the operator
    virtual void p(const int i,const Domain& x,const Domain& eta,Codomain& y)
    	const=0;

    // Derivative adjoint
    virtual void ps(const int i,const Domain& x,const Codomain& xi,Domain& y)
    	const=0;

    // Second derivative in the direction eta, adjoint, in the direction xi
    virtual void pps(const int i,const Domain& x, const Domain& eta,
	const Codomain& xi, Domain& y) const=0;

    // Get's the maximum index for each of the above operators
    virtual int max_index() const=0;
};

/* Functions that compute Hessian approximations */
namespace Hessians{

    // The Gauss-Newton Hessian approximation
    template <class Y,class U>
    class GaussNewton : public Operator <U,U> {
    private:
	const DiffOperator <Y,Y>& f;
	const DiffOperator <U,Y>& h;
	const U& u;
	list <Y>& workY;
	list <U>& workU;
    public:
	GaussNewton(
	    const DiffOperator <Y,Y>& f_,
	    const DiffOperator <U,Y>& h_,
	    const U& u_,
	    list <Y>& workY_,
	    list <U>& workU_) :
	    f(f_), h(h_), u(u_), workY(workY_), workU(workU_)
	{
	    // Check that we have enough work space
	    if(workY.size() < 3)
		pe_error("The Gauss-Newton Hessian approximation requires "
		    "at least three work elements in the state space.");
	    if(workU.size() < 1)
		pe_error("The Gauss-Newton Hessian approximation requires "
		    "at least one work element in the control space.");
	}
	void operator () (const U& p,U& result) const{
	    // Create shortcuts for each of the elements that we need
	    typename list <Y>::iterator y_iter=workY.begin();
	    Y& workY_0=*y_iter;
	    Y& workY_1=*(++y_iter);
	    Y& workY_2=*(++y_iter);
	    typename list <U>::iterator u_iter=workU.begin();
	    U& workU_0=*u_iter;

	    // Zero out the result 
	    Operations::zero(result);

	    // Accumulate the affect of the GN approximation on p one piece
	    // at a time
	    for(int i=0;i<f.max_index();i++){

		// Find the solution, y.  Store in workY_0.
		h(i,u,workY_0);

		// Find h'(u)p.  Store in workY_1.
		h.p(i,u,p,workY_1);

		// Find f'(y) workY_1.  Store in workY_2.
		f.p(i,workY_0,workY_1,workY_2);	

		// Find f'(y)* workY_2.  Store in workY_1.
		f.ps(i,workY_0,workY_2,workY_1);

		// Find h'(u)* workY_1.  Store in workU_0.
		h.ps(i,u,workY_1,workU_0);

		// Accumulate the result
		Operations::axpy(1.,workU_0,result);
	    }
	}
    };
    
    // The full-Newton Hessian (full Lagrangian)
    template <class Y,class U>
    class Newton : public Operator <U,U> {
    private:
	const DiffOperator <Y,Y>& f;
	const DiffOperator <U,Y>& h;
	const U& u;
	list <Y>& workY;
	list <U>& workU;
    public:
	Newton(
	    const DiffOperator <Y,Y>& f_,
	    const DiffOperator <U,Y>& h_,
	    const U& u_,
	    list <Y>& workY_,
	    list <U>& workU_) :
	    f(f_), h(h_), u(u_), workY(workY_), workU(workU_)
	{
	    // Check that we have enough work space
	    if(workY.size() < 4)
		pe_error("The full-Newton Hessian approximation requires "
		    "at least four work elements in the state space.");
	    if(workU.size() < 1)
		pe_error("The full-Newton Hessian approximation requires "
		    "at least one work element in the control space.");
	}
	void operator () (const U& p,U& result) const{
	    // Create shortcuts for each of the elements that we need
	    typename list <Y>::iterator y_iter=workY.begin();
	    Y& workY_0=*y_iter;
	    Y& workY_1=*(++y_iter);
	    Y& workY_2=*(++y_iter);
	    Y& workY_3=*(++y_iter);
	    typename list <U>::iterator u_iter=workU.begin();
	    U& workU_0=*u_iter;

	    // Zero out the result 
	    Operations::zero(result);

	    // Accumulate the affect of the Newton Hessian on p one piece
	    // at a time
	    for(int i=0;i<f.max_index();i++){

		// Find the solution, y.  Store in workY_0.
		h(i,u,workY_0);

		// Accumulate the Gauss-Newton part of the Hessian

		// Find h'(u)p.  Store in workY_1.
		h.p(i,u,p,workY_1);

		// Find f'(y) workY_1.  Store in workY_2.
		f.p(i,workY_0,workY_1,workY_2);	

		// Find f'(y)* workY_2.  Store in workY_1.
		f.ps(i,workY_0,workY_2,workY_1);

		// Find h'(u)* workY_1.  Store in workU_0.
		h.ps(i,u,workY_1,workU_0);

		// Accumulate the result
		Operations::axpy(1.,workU_0,result);

		// Accumlate the second order terms for the solution operator

		// Find f(y).  Store in workY_1.
		f(i,workY_0,workY_1);

		// Find f'(y)* workY_1.  Store in workY_2.
		f.ps(i,workY_0,workY_1,workY_2);

		// Find (h''(u)p)*workY_2.  Store in workU_0.
		h.pps(i,u,p,workY_2,workU_0);

		// Accumulate the result
		Operations::axpy(1.,workU_0,result);

		// Accumulate the second order terms for the projection operator

		// Assume that workY_1 still contains f(y)

		// Find h'(u)p.  Store in workY_2.
		h.p(i,u,p,workY_2);

		// Find (f''(y) workY_2)* workY_1.  Store in workY_3.
		f.pps(i,workY_0,workY_2,workY_1,workY_3);

		// Find h'(u)* workY_3.  Store in workU_0.
		h.ps(i,u,workY_3,workU_0);

		// Accumulate the result
		Operations::axpy(1.,workU_0,result);
	    }
	}
    };
}

namespace General{
    // Finds the graident in the parameter estimation problem
    template <class Y,class U>
    class getGradient: public Operator <U,U> {
    private:
	const DiffOperator <Y,Y>& f;
	const DiffOperator <U,Y>& h;
	list <Y>& workY;
	list <U>& workU;
    public:
	getGradient(
	    const DiffOperator <Y,Y>& f_,
	    const DiffOperator <U,Y>& h_,
	    list <Y>& workY_,
	    list <U>& workU_) :
	    f(f_), h(h_), workY(workY_), workU(workU_)
	{
	    // Check that we have enough work space
	    if(workY.size() < 3)
		pe_error("The gradient computation requires "
		    "at least three work elements in the state space.");
	    if(workU.size() < 1)
		pe_error("The gradient computation requires "
		    "at least one work element in the control space.");
	
	    // Check that the sizes of f and h are compatible
	    if(f.max_index()!=h.max_index())
		pe_error("The gradient computation requires that the "
		    "solution operator and the matching operator index over "
		    "the same number of elements.");
	}
	
	void operator () (const U& u,U& g) const{
	    // Create shortcuts for each of the elements that we need
	    typename list <Y>::iterator y_iter=workY.begin();
	    Y& workY_0=*y_iter;
	    Y& workY_1=*(++y_iter);
	    Y& workY_2=*(++y_iter);
	    typename list <U>::iterator u_iter=workU.begin();
	    U& workU_0=*u_iter;

	    // Zero out the gradient
	    Operations::zero(g);

	    // Accumulate the gradient one piece at a time
	    for(int i=0;i<f.max_index();i++){

		// Find the solution y, store in workY_0
		h(i,u,workY_0);

		// Determine how well this matches the experimental responses
		// (y-d).  Then, store in workY_1
		f(i,workY_0,workY_1);

		// Next, find (f'(y))^* workY_1.  Store in workY_2.
		f.ps(i,workY_0,workY_1,workY_2);

		// Finally, find h'(u)^* workY_2.  Store in workU_0.
		h.ps(i,u,workY_2,workU_0);

		// Add this to the current gradient accumulation
		Operations::axpy(1.,workU_0,g);
	    }
	}
    };

    // Computes the objective value of the problem given a control
    template <class Y,class U>
    class getObjValue : public Functional <U> {
    private:
	const DiffOperator <Y,Y>& f;
	const DiffOperator <U,Y>& h;
	list <Y>& workY;
    public:
	getObjValue(
	    const DiffOperator <Y,Y>& f_,
	    const DiffOperator <U,Y>& h_,
	    list <Y>& workY_) :
	    f(f_), h(h_), workY(workY_)
	{
	    // Check that we have enough work space
	    if(workY.size() < 2)
		pe_error("In order to find the objective value, we require"
		    " at least two work elements in the state space.");
	}

	double operator () (const U& u) const{
	    // Grab temp storage for h(u) and f(h(u)) 
	    typename list <Y>::iterator y_iter=workY.begin();
	    Y& workY_0 = *y_iter;
	    Y& workY_1 = *(++y_iter);

	    // Determine the objective function evaluated at h(u)
	    double obj_u=0.;
	    for(int i=0;i < f.max_index();i++){
		// Determine h_i(u).  Store it in workY_0.
		h(i,u,workY_0);

		// Determine f_i(h_i(u)).  Store it in workY_1.
		f(i,workY_0,workY_1);

		// Accumulate .5*|| f_i(h_i(u)) ||^2
		obj_u += .5*Operations::innr(workY_1,workY_1);
	    }

	    // Return the accumulated answer
	    return obj_u;
	}
    };
}
#endif
