#ifndef PEOPT_H
#define PEOPT_H

/** \file peopt.h
    \brief Parameter Estimation Optimization interface and algorithms
**/

#include<list>
#include<cmath>

using std::list;

// Error reporting for all of these functions
void pe_error(char message[]);

/// A simple operator specification 
template <class Domain, class Codomain>
class Operator {
public:
    /// Basic application
    virtual void operator () (const Domain& x,Codomain &y) const = 0;
};

/// A simple functional interface
template <class Domain>
class Functional {
public:
    /// Basic application
    virtual double operator () (const Domain& x) const = 0;
};

/// Basic linear algebra operations 
namespace Operations{
    /// Scalar multiple
    template <class Domain>
    void scal(const double alpha,Domain& x);

    /// Copy
    template <class Domain>
    void copy(const Domain& from, Domain& to); 

    /// alpha*x+y
    template <class Domain>
    void axpy(const double alpha,const Domain& x,Domain& y);

    /// Takes some element and sets it to all zeros
    template <class Domain>
    void zero(Domain& x);

    /// Inner product
    template <class Domain>
    double innr(const Domain& x,const Domain& y);
}

/// Functions that act has Hessian approximations
namespace Hessians{
    /// The identity Hessian approximation 
    template <class U>
    class Identity : public Operator <U,U> {
    public:
	void operator () (const U& p,U& result) const{
	    Operations::copy(p,result);
	}
    };

    /// The BFGS Hessian approximation.  
    /** Note, the formula we normally see for BFGS denotes the inverse
        Hessian approximation.  This is not the inverse, but the true
        Hessian approximation.  The oldY and oldS lists have the same
        structure as the BFGS preconditioner. **/
    template <class U>
    class BFGS : public Operator <U,U> {
    private:
    	list <U>& oldY;
	list <U>& oldS;
	list <U>& work;
    public:
    	BFGS(list <U>& oldY_,list <U>& oldS_,list <U>& work_)
	    : oldY(oldY_) , oldS(oldS_), work(work_) {};
	
        /// Operator interface
        /** It's not entirely clear to me what the best implementation for
            this method really is.  In the following implementation, we
            require an additional k work elements where k is the number of
            stored gradient and position differences.  It's possible to
            reduce this to 1 or 2, but we need to compute redundant
            information.  It's also possible to implementation the compact
            representation, see "Representations of quasi-Newton matrices
            and their use in limited memory methods" from Byrd, Nocedal,
            and Schnabel.  The problem with that algorithm is that is
            requires machinery such as linear system solves that we don't
            current have.  It also works much better with matrices or
            multivectors of data and we don't require the user to provide
            these abstractions. **/
	void operator () (const U& p,U& result) const{
	    // Check that the number of stored gradient and trial step
	    // differences is the same.
	    if(oldY.size() != oldS.size())
	    	pe_error("In the BFGS Hessian approximation, the number "
		    "of stored gradient differences must equal the number of "
		    "stored trial step differences.");

	    // Check that we have enough work space
	    if(work.size() != oldY.size())
	    	pe_error("In the BFGS Hessian approximation, we require a "
		    "number of work elements equal to the number of stored "
		    "gradient differences.");

	    // If we have no vectors in our history, we return the direction
	    Operations::copy(p,result);
	    if(oldY.size() == 0)
		return;

	    // As a safety check, insure that the inner product between all
	    // the (s,y) pairs is positive
	    typename list <U>::iterator y0=oldY.begin();
	    typename list <U>::iterator s0=oldS.begin();
	    while(y0!=oldY.end()){
		double inner_y_s=Operations::innr(*y0++,*s0++);
		if(inner_y_s<0)
		    pe_error("Detected a (s,y) pair in BFGS that possesed a "
			"nonpositive inner product");
	    }

	    // Othwerwise, we copy all of the trial step differences into the
	    // work space
	    typename list <U>::iterator Bisj_iter=work.begin();
	    typename list <U>::iterator sk_iter=oldS.begin();
	    while(Bisj_iter!=work.end())
	    	Operations::copy((*sk_iter++),(*Bisj_iter++));

	    // Keep track of the element Bisi
	    typename list <U>::iterator Bisi_iter=work.end(); Bisi_iter--;

	    // Keep iterating until Bisi equals the first element in the work
	    // list.  This means we have computed B1s1, B2s2, ..., Bksk.
	    Bisj_iter=work.begin();
	    typename list<U>::iterator si_iter=oldS.end(); si_iter--;
	    typename list<U>::iterator yi_iter=oldY.end(); yi_iter--;
	    typename list<U>::iterator sj_iter=oldS.begin();
	    while(1){

		// Create some reference to our iterators that are easier to
		// work with
		U& si=*si_iter;
		U& yi=*yi_iter;
		U& Bisi=*Bisi_iter;

		// Determine <Bisi,si>
		double inner_Bisi_si=Operations::innr(Bisi,si);

		// Determine <yi,si>
		double inner_yi_si=Operations::innr(yi,si);

		// Determine <si,Bip>
		double inner_si_Bip=Operations::innr(si,result);

		// Determine <yi,p>
		double inner_yi_p=Operations::innr(yi,p);

		// Determine -<si,Bip>/<Bisi,si> Bisi + Bip.  Store in Bip.
		// This will become B_{i+1}p.
		Operations::axpy(-inner_si_Bip/inner_Bisi_si,Bisi,result);

		// Determine <yi,p>/<yi,si> yi + w where we calculated w
		// in the line above.  This completes the calculation of
		// B_{i+1}p
		Operations::axpy(inner_yi_p/inner_yi_si,yi,result);

		// Check whether or not we've calculated B_{i+1}p for the
		// last time
		if(Bisi_iter==work.begin()) break;

		// Begin the calculation of B_{i+1}sj
	    	while(si_iter!=sj_iter){
		    // Add some additional references to the iterators 
		    U& sj=*sj_iter;
		    U& Bisj=*Bisj_iter;

		    // Determine <si,Bisj>
		    double inner_si_Bisj=Operations::innr(si,Bisj);

		    // Determine <yi,sj>
		    double inner_yi_sj=Operations::innr(yi,sj);

		    // Determine -<si,Bisj>/<Bisi,si> Bisi + Bisj
		    // Store in Bisj.  This will become B_{i+1}sj.
		    Operations::axpy(-inner_si_Bisj/inner_Bisi_si,Bisi,Bisj);

		    // Determine <yi,sj>/<yi,si> yi + w where we calculated w
		    // in the line above.  This completes the computation of
		    // B_{i+1}sj.
		    Operations::axpy(inner_yi_sj/inner_yi_si,yi,Bisj);

		    // Change j to be j-1 and adjust Bisj and sj accordingly
		    sj_iter++;
		    Bisj_iter++;
		}

		// At this point, we've computed all Bisj entries on the current
		// row.  As a result, we increment i and set j to be k.  This
		// requires us to modify si, yi, sj, Bisj, and Bisi accordingly.
		
		// Increment i and adjust si
		si_iter--;

		// Increment i and adjust yi
		yi_iter--;

		// Set j=k and adjust sj
		sj_iter=oldS.begin();

		// Set j=k, increment i, and adjust Bisj
		Bisj_iter=work.begin();

		// Increment i and adjust Bisi
		Bisi_iter--;
	    }
	}
    };
    
    /// The SR1 Hessian approximation.  
    /** The oldY and oldS lists have the same structure as the BFGS
        preconditioner. **/
    template <class U>
    class SR1 : public Operator <U,U> {
    private:
    	list <U>& oldY;
	list <U>& oldS;
	list <U>& work;
    public:
        /// Constructor
    	SR1(list <U>& oldY_,list <U>& oldS_,list <U>& work_)
	    : oldY(oldY_) , oldS(oldS_), work(work_) {};
        /// Operator interface
	void operator () (const U& p,U& result) const{
	    // Check that the number of stored gradient and trial step
	    // differences is the same.
	    if(oldY.size() != oldS.size())
	    	pe_error("In the SR1 Hessian approximation, the number "
		    "of stored gradient differences must equal the number of "
		    "stored trial step differences.");

	    // Check that we have enough work space
	    if(work.size() != oldY.size())
	    	pe_error("In the SR1 Hessian approximation, we require a "
		    "number of work elements equal to the number of stored "
		    "gradient differences.");

	    // If we have no vectors in our history, we return the direction
	    Operations::copy(p,result);
	    if(oldY.size() == 0)
		return;

	    // Othwerwise, we copy all of the trial step differences into the
	    // work space
	    typename list <U>::iterator Bisj_iter=work.begin();
	    typename list <U>::iterator sk_iter=oldS.begin();
	    while(Bisj_iter!=work.end())
	    	Operations::copy((*sk_iter++),(*Bisj_iter++));

	    // Keep track of the element Bisi
	    typename list <U>::iterator Bisi_iter=work.end(); Bisi_iter--;

	    // Keep iterating until Bisi equals the first element in the work
	    // list.  This means we have computed B1s1, B2s2, ..., Bksk.
	    Bisj_iter=work.begin();
	    typename list<U>::iterator si_iter=oldS.end(); si_iter--;
	    typename list<U>::iterator yi_iter=oldY.end(); yi_iter--;
	    typename list<U>::iterator sj_iter=oldS.begin();
	    while(1){

		// Create some reference to our iterators that are easier to
		// work with
		U& si=*si_iter;
		U& yi=*yi_iter;
		U& Bisi=*Bisi_iter;

		// Determine <yi,p>
		double inner_yi_p=Operations::innr(yi,p);

		// Determine <Bisi,p>
		double inner_Bisi_p=Operations::innr(Bisi,p);

		// Determine <yi,si>
		double inner_yi_si=Operations::innr(yi,si);

		// Determine <Bisi,si>
		double inner_Bisi_si=Operations::innr(Bisi,si);

		// Determine (<yi,p>-<Bisi,p>) / (<y_i,s_i>-<Bisi,si>).
		// Store in alpha
		double alpha=
		    (inner_yi_p-inner_Bisi_p)/(inner_yi_si-inner_Bisi_si);

		// Determine alpha y_i + Bip.  Store in result (which
		// accumulate Bip).
		Operations::axpy(alpha,yi,result);

		// Then, add -alpha*Bisi to this result
		Operations::axpy(-alpha,Bisi,result);

		// Check whether or not we've calculated B_{i+1}p for the
		// last time
		if(Bisi_iter==work.begin()) break;

		// Begin the calculation of B_{i+1}sj
	    	while(si_iter!=sj_iter){
		    // Add some additional references to the iterators 
		    U& sj=*sj_iter;
		    U& Bisj=*Bisj_iter;

		    // Determine <yi,sj>
		    double inner_yi_sj=Operations::innr(yi,sj);

		    // Determine <Bisi,sj>
		    double inner_Bisi_sj=Operations::innr(Bisi,sj);

		    // Determine (<yi,p>-<Bisi,p>) / (<y_i,s_i>-<Bisi,si>).
		    // Store in beta 
		    double beta=
			(inner_yi_sj-inner_Bisi_sj)/(inner_yi_si-inner_Bisi_si);
		
		    // Determine beta y_i + Bisj.  Store in Bisj. 
		    Operations::axpy(beta,yi,Bisj);

		    // Add -beta*Bisi to this result
		    Operations::axpy(-beta,Bisi,Bisj);

		    // Change j to be j-1 and adjust Bisj and sj accordingly
		    sj_iter++;
		    Bisj_iter++;
		}

		// At this point, we've computed all Bisj entries on the current
		// row.  As a result, we increment i and set j to be k.  This
		// requires us to modify si, yi, sj, Bisj, and Bisi accordingly.
		
		// Increment i and adjust si
		si_iter--;

		// Increment i and adjust yi
		yi_iter--;

		// Set j=k and adjust sj
		sj_iter=oldS.begin();

		// Set j=k, increment i, and adjust Bisj
		Bisj_iter=work.begin();

		// Increment i and adjust Bisi
		Bisi_iter--;
	    }
	}
    };
}


/// Functions that precondition truncated CG
namespace Preconditioners{
    
    /// The identity preconditioner
    template <class U>
    class Identity : public Operator <U,U> {
    public:
	void operator () (const U& p,U& result) const{
	    Operations::copy(p,result);
	}
    };

    /// The BFGS preconditioner
    /** The oldY list has the following structure
        oldY[0] = y_k = grad f(h(u_k)) - grad f(h(u_{k-1}))
        oldY[1] = y_{k-1} = grad f(h(u_{k-1})) - grad f(h(u_{k-2}))
        The oldS list has the following structure
        oldS[0] = s_k = u_k - u_k{-1}
        oldS[1] = s_{k-1} = u_{k-1} - u_k{k-2} **/
    template <class U>
    class BFGS : public Operator <U,U> {
    private:
    	list <U>& oldY;
	list <U>& oldS;
    public:
        /// Constructor 
    	BFGS(list <U>& oldY_,list <U>& oldS_) : oldY(oldY_) , oldS(oldS_) {}
        /// Operator interface
	void operator () (const U& p,U& result) const{
	    // Check that we have an even number of elements in the info space
	    if(oldY.size() != oldS.size())
	    	pe_error("In the BFGS preconditioner, the number"
		" of stored gradients must equal the number of stored"
		" trial steps.");
#if 0
	    // As a safety check, insure that the inner product between all
	    // the (s,y) pairs is positive
	    typename list <U>::iterator y0=oldY.begin();
	    typename list <U>::iterator s0=oldS.begin();
	    while(y0!=oldY.end()){
		double inner_y_s=Operations::innr(*y0++,*s0++);
		if(inner_y_s<0)
		    pe_error("Detected a (s,y) pair in BFGS that possesed a "
			"nonpositive inner product");
	    }
#endif

	    // Create two vectors to hold some intermediate calculations
	    double alpha[oldY.size()];
	    double rho[oldY.size()];

	    // Before we begin computing, copy p to our result 
	    Operations::copy(p,result);

	    // In order to compute, we first iterate over all the stored
	    // element in the forward direction.  Then, we iterate over them
	    // backward.
	    typename list <U>::iterator y_iter=oldY.begin();
	    typename list <U>::iterator s_iter=oldS.begin();
	    int i=0;
	    while(y_iter != oldY.end()){
		// Find y_k, s_k, and their inner product
		U& y_k=*(y_iter++);
		U& s_k=*(s_iter++);
		rho[i]=1./Operations::innr(y_k,s_k);

		// Find rho_i <s_i,result>.  Store in alpha_i
		alpha[i]=rho[i]*Operations::innr(s_k,result);

		// result = - alpha_i y_i + result 
		Operations::axpy(-alpha[i],y_k,result);

		// Make sure we don't overwrite alpha and rho
		i++;
	    }

	    // Assume that H_0 is the identity operator (which may or may not
	    // work in Hilbert space)

	    // Now, let us iterate backward over our elements to complete the
	    // computation
	    while(y_iter != oldY.begin()){
		// Find y_k and s_k
		U& s_k=*(--s_iter);
		U& y_k=*(--y_iter);

		// beta=rho_i <y_i,result>
		double beta= rho[--i] * Operations::innr(y_k,result);

		// result=  (alpha_i-beta) s_i + result
		Operations::axpy(alpha[i]-beta,s_k,result);
	    }
	}
    };
    
    /// The SR1 preconditioner.  
    /** In this definition, we take a short cut and simply use the SR1
        Hessian approximation where we swap Y and S.  The oldY and oldS lists
        have the same structure as the BFGS preconditioner. **/
    template <class U>
    class SR1 : public Operator <U,U> {
    private:
	Hessians::SR1 <U> sr1;
    public:
    	SR1(list <U>& oldY_,list <U>& oldS_,list <U>& work_)
	    : sr1(oldS_,oldY_,work_) {};
	
	void operator () (const U& p,U& result) const{
	    sr1(p,result);
	}
    };
}

/// Functions that support the trust region algorithm
namespace TrustRegion{
    	
    /// Computes the Steihaug-Toint trial step
    template <class U>
    void getStep(
	Operator<U,U>& Minv,
	Operator<U,U>& H,
	const U& u,
	const U& g,
	const double delta,
	const int max_iter,
	const double eps_cg,
	list <U>& workU,
	U& s){

	// Check that we have enough work space
	if(workU.size() < 4)
	    pe_error("In order to compute the truncated CG algorithm, we "
		"require at least four work elements in the control space.");

	// Create shortcuts for each of the elements that we need
	typename list <U>::iterator u_iter=workU.begin();
	U& s_k=s;
	U& g_k=*u_iter;
	U& v_k=*(++u_iter);
	U& p_k=*(++u_iter);
	U& H_pk=*(++u_iter);

	// Allocate memory for a few constants that we need to track 
	double kappa;
	double sigma;
	double alpha;
	double beta;
	double norm_sk_M2,norm_skp1_M2,norm_pk_M2,norm_g2;
	double inner_sk_M_pk,inner_gk_vk,inner_gkp1_vkp1;

	// Initialize our variables
	Operations::zero(s_k);			// s_0=0
	Operations::copy(g,g_k);		// g_0=g
	Minv(g_k,v_k);				// v_0=inv(M)*g_0
	Operations::copy(v_k,p_k);		// p_0=-v_0
	Operations::scal(-1.,p_k);
	norm_sk_M2=0.;				// || s_0 ||_M^2 = 0
	norm_pk_M2=Operations::innr(g_k,v_k);	// || p_0 ||_M^2 = <g_0,v_0>	
	inner_sk_M_pk=0.;			// <s_0,M p_0>=0
	inner_gk_vk=norm_pk_M2;			// <g_0,v_0> = || p_0 ||_M^2
	norm_g2=Operations::innr(g,g);		// || -g ||^2

	// Run truncated CG until we hit our max iteration or we converge
	for(int i=0;i<max_iter;i++){
	    // H_pk=H p_k
	    H(p_k,H_pk);

	    // Compute the curvature for this direction.  kappa=<p_k,H p_k>
	    kappa=Operations::innr(p_k,H_pk);

	    // If we have negative curvature, don't bother with the next two "
	    // steps since we're going to exit and we won't need them.  
	    if(kappa > 0){
		// Determine a trial point
		alpha = Operations::innr(g_k,v_k)/kappa;

		// || s_k+alpha_k p_k ||
		norm_skp1_M2=norm_sk_M2+2*alpha*inner_sk_M_pk
		    +alpha*alpha*norm_pk_M2;
	    }

	    // If we have negative curvature or our trial point is outside the
	    // trust region radius, terminate truncated-CG and find our final
	    // step.
	    if(kappa <= 0 || norm_skp1_M2 >= delta*delta){
		// sigma = positive root of || s_k + sigma p_k ||_M = delta
		sigma= (-inner_sk_M_pk + sqrt(inner_sk_M_pk*inner_sk_M_pk
		    + norm_pk_M2*(delta*delta-norm_sk_M2)))/norm_pk_M2;

		// s_kp1=s_k+sigma p_k
		Operations::axpy(sigma,p_k,s_k);

		// Exit the loop
		break;
	    }

	    // Take a step in the computed direction. s_k=s_k+alpha p_k
	    Operations::axpy(alpha,p_k,s_k);

	    // Update the norm of sk
	    norm_sk_M2=norm_skp1_M2;
	    
	    // g_k=g_k+alpha H p_k
	    Operations::axpy(alpha,H_pk,g_k);

	    // Test whether we've converged CG
	    if(Operations::innr(g_k,g_k) / norm_g2 <= eps_cg*eps_cg)
		break;

	    // v_k = Minv g_k
	    Minv(g_k,v_k);

	    // Compute the new <g_kp1,v_kp1>
	    inner_gkp1_vkp1=Operations::innr(g_k,v_k);

	    // beta = <g_kp1,v_kp1> / <g_k,v_k>
	    beta= inner_gkp1_vkp1 / inner_gk_vk;

	    // Store the new inner product between g_k and p_k
	    inner_gk_vk=inner_gkp1_vkp1;
	    
	    // Find the new search direction.  p_k=-v_k + beta p_k
	    Operations::scal(beta,p_k);
	    Operations::axpy(-1.,v_k,p_k);

	    // Update the inner product between s_k and M p_k
	    inner_sk_M_pk=beta*(inner_sk_M_pk+alpha*norm_pk_M2);

	    // Update the norm of p_k
	    norm_pk_M2=inner_gk_vk+beta*beta*norm_pk_M2; 
	}
    }

    /// Checks whether we accept or reject a step
    template <class U>
    void checkStep(
	const U& u,
	const U& s,
	const Functional<U>& obj_fn,
	const Operator<U,U>& H,
	const U& g,
	const double obj_u,
	const double eta1,
	const double eta2,
	list <U>& workU,
	double& delta,
	bool& accept,
	double& obj_ups,
	double& rho
    ){
	// Check that we have enough work space
	if(workU.size() < 2)
	    pe_error("In order to check the step from the TR subproblem, we "
		"require at least two work elements in the control space.");

	// Grab temp storage for u+s, and H(u)(u+s)
	typename list <U>::iterator u_iter=workU.begin();
	U& ups = *u_iter;
	U& Hu_s = *(++u_iter);
	
	// Determine u+s 
	Operations::copy(s,ups);
	Operations::axpy(1.0,u,ups);

	// Determine the objective function evaluated at h(u)
	// We now pass this in as a parameter.

	// Determine the objective function evaluated at h(u+s)
	obj_ups=obj_fn(ups);
	
	// Determine H(u)s
	H(s,Hu_s);

	// Determine alpha+<g,s>+.5*<H(u)s,s>
	double model_s=obj_u+Operations::innr(g,s)+.5*Operations::innr(Hu_s,s);

	// As a safety check, insure that we receive reduction in the objective.
	// Sometimes, we have some odd looking Hessian approximations due to
	// the quasi-Newton methods.  Further, we don't exactly solve the TR
	// subproblem.  Hence, it's possible that the model minizer is greater
	// than the current objective value and the objective value at the
	// trial step is greater which makes it appear as though we have a
	// reduction in value which we don't.
	// This could also be a bug in our code.  It's worthwhile to
	// investigate.
	if(obj_ups >= obj_u){
	  accept=false;
	  delta /= 2.;
	  return;
	}

	// Determine the ratio of reductions
	rho = (obj_u - obj_ups) / (obj_u - model_s);

	// Determine if we should keep the step
	accept = rho >= eta1 ? true : false;

	// Update the trust region radius
	if(rho >= eta2)
	    delta *= 2.;
	else if(rho >= eta1 && rho < eta2)
	    delta = delta;
	else
	    delta /= 2.;
    }

    /// Reasons why we stop the algorithm
    enum StoppingCondition{
      NotConverged,               ///< Algorithm did not converge
      GradientSmall,              ///< Gradient was sufficiently small
      RelativeGradientSmall,      ///< Relative gradient was sufficiently small
      TrustRegionSmall,           ///< Trust region became sufficiently small
      RelativeTrustRegionSmall,   ///< Relative small trust region
      RelativeObjectiveSmall      ///< Relative value of objective is small
    };

    /// Checks a set of stopping conditions
    template <class U>
    StoppingCondition checkStop(
    	Operator <U,U>& Minv,
	const U& u,
	const U& g,
	const U& g_typ,
	const double obj_u,
	const double obj_ums,
	const double delta,
	const double eps_g,
	const double eps_d,
	const double eps_f,
	list <U>& workU
    ){
	// Check that we have enough work space
	if(workU.size() < 1)
	    pe_error("In order to check the stopping conditions, we require"
		" at least one work element in the control space.");

	// Grab temp storage for Minv(g) 
	typename list <U>::iterator u_iter=workU.begin();
	U& Minv_g = *u_iter;

	// Determine the norm of the gradient squared for debugging purposes
	double norm_g2 = Operations::innr(g,g);

	// Determine the norm of a typical gradient squared
	double norm_gtyp2 = Operations::innr(g_typ,g_typ);

	// Determine the norm of the gradient squared using the M^{-1} norm
	Minv(g,Minv_g);
	double Minv_norm_g2 = Operations::innr(g,Minv_g);

	// Check whether the norm of the gradient has become small.  If so,
	// terminate.
	if(norm_g2 < eps_g*eps_g) return GradientSmall;

	// Check whether the norm is small relative to some typical gradient
	if(norm_g2 < eps_g*eps_g*norm_gtyp2) return RelativeGradientSmall;

	// Determine the norm of the current iterate
	double norm_u2=Operations::innr(u,u);

	// Check whether the relative size of trust region radius has become
	// too small.  If so, terminate
	if(delta*delta <= eps_d*eps_d*norm_u2) return RelativeTrustRegionSmall;

	// Check whether the absolute size of the trust region radius has
	// become to small.  If so, terminate
	if(delta*delta <= eps_d*eps_d) return TrustRegionSmall;

	// If the relative difference between objective values becomes too
	// small, terminate
	if((obj_ums-obj_u)/(1.+fabs(obj_ums))<eps_f)
	    return RelativeObjectiveSmall;

	// Otherwise, return that we're not converged 
	return NotConverged;
    }
}

#endif  // PEOPT_H
