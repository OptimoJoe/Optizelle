#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include "peopt/peopt.h"
#include "peopt/vspaces.h"
#include "peopt/json.h"

typedef size_t Natural;
typedef ptrdiff_t Integer;

// Index vectors starting from 1
inline Natural itok(Natural i) {
    return i-1;
}
// Index packed matrices starting from 1
inline Natural ijtokp(Natural i,Natural j) {
    return i+j*(j-1)/2-1;
}

// A simple structure that contains either a digonal matrix or a sparse matrix
// in triple format.  If we have a diagonal matrix, the js are empty.
template <typename Real>
struct SparseMat{
    std::vector <Natural> is;
    std::vector <Natural> js;
    std::vector <Real> data;
};

// Stores a sparse SDPA format problem of the form
//
// min b1*x1 + ... + bm*xm
// st  A1*x1 + ... + Am*xm - A0 >= 0
//
// where each A has a block structure with sizes blksizes_1,...,blksizes_nblocks
template <typename Real>
struct SparseSDP {
    // Block sizes.  Negative means a diagonal block.
    std::vector <Integer> blk_sizes;

    // Objective function
    std::vector <Real> b;

    // Constraint matrices
    std::vector < std::vector <SparseMat<Real> > > A;
};

// Clears out whitespace
void eat_whitespace(std::stringstream& sin) {
    while(!sin.eof() && (sin.peek() == ' ' || sin.peek()=='\t'))
        sin.get();
}

// Clears out the formating characters
void eat_formatting(std::stringstream& sin) {
    while(!sin.eof() && (sin.peek() == ' ' || sin.peek()=='\t' ||
        sin.peek()==',' || sin.peek() =='(' || sin.peek() == ')' ||
        sin.peek()=='{' || sin.peek() == '}')
    )
        sin.get();
}

// Reads in the sparse SDPA format
template <typename Real>
void parse_sdpa(const std::string fname,SparseSDP <Real>& prob) {
    // Open the file
    std::ifstream fin(fname.c_str());

    // Check if we were able to open the file
    if(fin.fail()) {
        std::cerr << "Unable to open the file: " << fname << '.' << std::endl;
        exit(EXIT_FAILURE);
    }

    // Get rid of all the lines with comments
    std::string line;
    do {
        std::getline (fin,line);
        if(fin.fail()) {
            std::cerr << "Error while parsing initial comments." << std::endl;
            exit(EXIT_FAILURE);
        } else if(fin.eof()) {
            std::cerr << "Found end of file while parsing initial comments."
                << std::endl;
            exit(EXIT_FAILURE);
        }
    } while(line[0]=='"' || line[0]=='*');

    // Next, read in the number of constraint matrices
    std::stringstream sin(line);
    Natural m(0);
    sin >> m;
    if(sin.fail()) {
        std::cerr << "Error while parsing the number of constraint matrices."
            << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in the number of blocks
    std::getline (fin,line);
    if(fin.fail()) {
        std::cerr << "Error while reading the number of blocks." << std::endl;
        exit(EXIT_FAILURE);
    } else if(fin.eof()) {
        std::cerr << "Found end of file while reading the number of blocks."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    sin.str(line);
    sin.clear();
    Natural nblocks(0);
    sin >> nblocks;
    if(sin.fail()) {
        std::cerr << "Error while parsing the number of blocks." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read in the sizes of the blocks
    prob.blk_sizes.clear();
    std::getline (fin,line);
    if(fin.fail()) {
        std::cerr << "Error while reading the sizes of blocks." << std::endl;
        exit(EXIT_FAILURE);
    } else if(fin.eof()) {
        std::cerr << "Found end of file while reading the sizes of blocks."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    sin.str(line);
    sin.clear();
    while(!sin.eof()) {
        // Clear out the formating characters
        eat_formatting(sin);
        if(sin.eof()) break;

        // Parse the block size
        Integer blk_size(0);
        sin >> blk_size;
        if(sin.fail()) {
            std::cerr << "Error while parsing the block size." << std::endl;
            exit(EXIT_FAILURE);
        }
        prob.blk_sizes.push_back(blk_size);
    }

    // Check that the number of blocks corresponds to what we parsed
    if(nblocks != prob.blk_sizes.size()) 
        std::cerr << "The number of parsed block sizes differs from the "
            "specified number." << std::endl;

    // Read in the objective function
    prob.b.clear();
    std::getline (fin,line);
    if(fin.fail()) {
        std::cerr << "Error while reading the objective." << std::endl;
        exit(EXIT_FAILURE);
    } else if(fin.eof()) {
        std::cerr << "Found end of file while reading the objective."
            << std::endl;
        exit(EXIT_FAILURE);
    }
    sin.str(line);
    sin.clear();
    while(!sin.eof()) {
        // Clear out the formating characters
        eat_formatting(sin);
        if(sin.eof()) break;

        // Read the objective
        Real val(0.);
        sin >> val;
        if(sin.fail()) {
            std::cerr << "Error while parsing the objective." << std::endl;
            exit(EXIT_FAILURE);
        }
        prob.b.push_back(val);
    }

    // Create structure for the constraint matrices
    prob.A= std::vector < std::vector <SparseMat<Real> > >
        (m+1, std::vector <SparseMat<Real> > (nblocks));
    for(Natural i=0;i<m+1;i++)
        for(Natural j=0;j<prob.blk_sizes.size();j++) {
            prob.A[i][j].data.clear();
            prob.A[i][j].is.clear();
            prob.A[i][j].js.clear();
        }

    // Read constraints until we finish
    while(1) {
        // Read in the constraint matrices
        std::getline (fin,line);
        if(fin.eof())
            break;
        else if(fin.fail()) {
            std::cerr << "Error while reading the constraints." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Read in the matno, blkno, i, j and entry
        Integer matno(0), blkno(0), i(0), j(0);
        Real entry;
        sin.str(line);
        sin.clear();
        sin >> matno >> blkno >> i >> j >> entry;
        if(sin.fail()) {
            std::cerr << "Error while parsing the constraint." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Insert the element into the appropriate slot
        if(prob.blk_sizes[itok(blkno)] < 0) {
            if(i!=j) {
                std::cerr
                    << "Specified a off-diagonal element of a diagonal block."
                    << std::endl;
                exit(EXIT_FAILURE);
            }
            prob.A[matno][itok(blkno)].is.push_back(i);
            prob.A[matno][itok(blkno)].data.push_back(entry);
        } else {
            if(i>j) std::swap(i,j);
            prob.A[matno][itok(blkno)].is.push_back(i);
            prob.A[matno][itok(blkno)].js.push_back(j);
            prob.A[matno][itok(blkno)].data.push_back(entry);
        }
    }

    // Close the file
    fin.close();
}

// Used for doing a tagged sort on sparse matrices
struct MatComparison{
    const std::vector <Natural>& is;
    const std::vector <Natural>& js;
    MatComparison(
        const std::vector <Natural>& is_,
        const std::vector <Natural>& js_
    ) : is(is_), js(js_) {}
    bool operator () (Natural k,Natural l) {
        return ijtokp(is[k],js[k]) < ijtokp(is[l],js[l]);
    }
};

// Used for doing a tagged sort on digonal matrices 
struct DiagComparison{
    const std::vector <Natural>& is;
    DiagComparison(const std::vector <Natural>& is_) : is(is_) {}
    bool operator () (Natural k,Natural l) {
        return is[k] < is[l];
    }
};

// Sorts the indices used in the SDP problem first by column and then by row.
template <typename Real>
void sort_sdp (SparseSDP <Real>& prob) {
    // Loop over the constraints
    for(Natural i=0;i<prob.A.size();i++) {
        // Loop over the blocks
        for(Natural j=0;j<prob.blk_sizes.size();j++) {
            // Diagonal blocks
            if(prob.blk_sizes[j]<0) {
                // Do the tagged sort
                std::vector <Natural> tag(prob.A[i][j].is.size());
                for(Natural k=0;k<prob.A[i][j].is.size();k++)
                    tag[k]=k;
                DiagComparison comp(prob.A[i][j].is);
                std::sort(tag.begin(),tag.end(),comp);

                // Create new memory for the indices and data
                std::vector <Natural> is_new(prob.A[i][j].is.size());
                std::vector <Real> data_new(prob.A[i][j].data.size());

                // Insert the data where it needs to be
                for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                    is_new[k]=prob.A[i][j].is[tag[k]];
                    data_new[k]=prob.A[i][j].data[tag[k]];
                }

                // Reassign the indices and data
                prob.A[i][j].is=is_new;
                prob.A[i][j].data=data_new;

            // Sparse blocks
            } else {
                // Do the tagged sort
                std::vector <Natural> tag(prob.A[i][j].is.size());
                for(Natural k=0;k<prob.A[i][j].is.size();k++)
                    tag[k]=k;
                MatComparison comp(prob.A[i][j].is,prob.A[i][j].js);
                std::sort(tag.begin(),tag.end(),comp);

                // Create new memory for the indices and data
                std::vector <Natural> is_new(prob.A[i][j].is.size());
                std::vector <Natural> js_new(prob.A[i][j].js.size());
                std::vector <Real> data_new(prob.A[i][j].data.size());

                // Insert the data where it needs to be
                for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                    is_new[k]=prob.A[i][j].is[tag[k]];
                    js_new[k]=prob.A[i][j].js[tag[k]];
                    data_new[k]=prob.A[i][j].data[tag[k]];
                }

                // Reassign the indices and data
                prob.A[i][j].is=is_new;
                prob.A[i][j].js=js_new;
                prob.A[i][j].data=data_new;
            }
        }
    }
}

// Cartesian product between Rm and the SQL spaces 
template <typename Real>
struct RmxR {
private:
    // This is a templated namespace.  Do not allow construction.
    RmxR();

    // Create a shortcut for these routines
    typedef peopt::Rm <Real> Rm;

public:
    // Use a Cartesian product of two Rm spaces for our storage 
    typedef std::pair <typename Rm::Vector,Real> Vector;

    // Memory allocation and size setting.
    static void init(const Vector& x, Vector& y) {
        Rm::init(x.first,y.first);
    }
    
    // y <- x (Shallow.  No memory allocation.)
    static void copy(const Vector& x, Vector& y) {
        Rm::copy(x.first,y.first);
        y.second=x.second;
    }

    // x <- alpha * x.
    static void scal(const Real& alpha, Vector& x) {
        Rm::scal(alpha,x.first);
        x.second*=alpha;
    }

    // y <- alpha * x + y.
    static void axpy(const Real& alpha, const Vector& x, Vector& y) {
        Rm::axpy(alpha,x.first,y.first);
        y.second+=alpha*x.second;
    }

    // innr <- <x,y>.
    static Real innr(const Vector& x,const Vector& y) {
        return Rm::innr(x.first,y.first)+x.second*y.second;
    }

    // x <- 0.
    static void zero(Vector& x) {
        Rm::zero(x.first);
        x.second=Real(0.);
    }
};

// Define the SDP objective where 
// 
// f(x)=<b,x> 
//
template <typename Real>
struct SDPObj : public peopt::ScalarValuedFunction <Real,peopt::Rm> {
private:
    const SparseSDP <Real>& prob;
        
public:
    typedef peopt::Rm <Real> Rm;
    typedef typename Rm::Vector X_Vector;

    // Grab a reference to the underlying SDP problem
    SDPObj(const SparseSDP <Real>& prob_) : prob(prob_) {}

    // Evaluation 
    double operator () (const X_Vector& x) const {
        return Rm::innr(prob.b,x);
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad 
    ) const {
        Rm::copy(prob.b,grad);
    }

    // Hessian-vector product
    void hessvec(
        const X_Vector& x,
        const X_Vector& dx,
        X_Vector& H_dx
    ) const {
        Rm::zero(H_dx);
    }
};


// Define the SDP inequality where 
//
// h(x) = A1*x1 + ... + Am*xm - A0 >= 0
//
template <typename Real>
struct SDPIneq : public peopt::VectorValuedFunction <Real,peopt::Rm,peopt::SQL>{
public:
    typedef peopt::Rm <Real> X;
    typedef peopt::SQL <Real> Z;
    typedef typename X::Vector X_Vector;
    typedef typename Z::Vector Z_Vector;

private:
    const SparseSDP <Real>& prob;

    // z=h(x), except that we de start with A_start.  Mostly, this is
    // to toggle whether we start from A0 or A1.
    template <Natural start>
    void eval_from(
        const X_Vector& x,
        Z_Vector& z
    ) const {
        // Zero out the solution
        Z::zero(z);

        // Loop over the constraints
        for(Natural i=start;i<prob.A.size();i++) {
            // Loop over the blocks
            for(Natural j=0;j<prob.blk_sizes.size();j++) {
                // Diagonal blocks
                if(prob.blk_sizes[j]<0) {
                    // Loop over the indices
                    for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                        // Constant block
                        if(i==0)
                            z(j+1,prob.A[i][j].is[k]) -= prob.A[i][j].data[k];
                        // Variable block
                        else 
                            z(j+1,prob.A[i][j].is[k])
                                += prob.A[i][j].data[k]*x[itok(i)];
                    }

                // Sparse blocks
                } else {
                    // Loop over the indices
                    for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                        // Constant block
                        if(i==0) {
                            z(j+1,prob.A[i][j].is[k],prob.A[i][j].js[k])
                                -= prob.A[i][j].data[k];
                            if(prob.A[i][j].is[k] != prob.A[i][j].js[k])
                                z(j+1,prob.A[i][j].js[k],prob.A[i][j].is[k])
                                    -= prob.A[i][j].data[k];
                        // Variable block
                        } else {
                            z(j+1,prob.A[i][j].is[k],prob.A[i][j].js[k])
                                += prob.A[i][j].data[k]*x[itok(i)];
                            if(prob.A[i][j].is[k] != prob.A[i][j].js[k])
                                z(j+1,prob.A[i][j].js[k],prob.A[i][j].is[k])
                                    += prob.A[i][j].data[k]*x[itok(i)];
                        }
                    }
                }
            }
        }
    }

public:
    // Grab a reference to the underlying SDP problem
    SDPIneq(const SparseSDP <Real>& prob_) : prob(prob_) {}

    // z=h(x) 
    void operator () (
        const X_Vector& x,
        Z_Vector& z
    ) const {
        eval_from <0> (x,z);
    }

    // z=h'(x)dx
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Z_Vector& z
    ) const {
        eval_from <1> (dx,z);
    }

    // xhat=h'(x)*dz
    void ps(
        const X_Vector& x,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        // Zero out the solution
        X::zero(xhat);

        // Loop over the constraints
        for(Natural i=1;i<prob.A.size();i++) {
            // Loop over the blocks
            for(Natural j=0;j<prob.blk_sizes.size();j++) {
                // Diagonal blocks
                if(prob.blk_sizes[j]<0) {
                    // Loop over the indices
                    for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                        xhat[itok(i)]+=
                            dz(j+1,prob.A[i][j].is[k]) * prob.A[i][j].data[k];
                    }

                // Sparse blocks
                } else {
                    // Loop over the indices
                    for(Natural k=0;k<prob.A[i][j].is.size();k++) {
                        xhat[itok(i)]+=
                            dz(j+1,prob.A[i][j].is[k],prob.A[i][j].js[k])
                            * prob.A[i][j].data[k];
                        if(prob.A[i][j].is[k] != prob.A[i][j].js[k])
                            xhat[itok(i)]+=
                                dz(j+1,prob.A[i][j].js[k],prob.A[i][j].is[k])
                                * prob.A[i][j].data[k];
                    }
                }
            }
        }
    }

    // xhat=(h''(x)dx)*dz
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        X::zero(xhat);
    }
};

// Initializes an SQL vector 
template <typename Real>
void initSQL(
    const SparseSDP <Real>& prob,
    typename peopt::SQL <Real>::Vector& x
) {
    // Create a type shortcut
    typedef peopt::SQL <Real> SQL;

    // Figure out the structure of the codomain of the inequality
    // constraint h 
    std::vector <peopt::Natural> sizes(prob.blk_sizes.size());
    std::vector <peopt::Cone::t> types(prob.blk_sizes.size());
    for(Natural i=0;i<sizes.size();i++) {
        sizes[i]=abs(prob.blk_sizes[i]);
        if(prob.blk_sizes[i]<0)
            types[i]=peopt::Cone::Linear;
        else
            types[i]=peopt::Cone::Semidefinite;
    }

    // Create a new element 
    typename SQL::Vector xx(peopt::Messaging(),types,sizes);

    // Initialize the memory for the user input
    SQL::init(xx,x);
}

// Define the phase-1 objective where 
// 
// f(x,y) = 1/2 || y - epsilon ||^2 + beta/2 || h(x) ||^2
//
template <typename Real>
struct Phase1Obj : public peopt::ScalarValuedFunction <Real,RmxR> {
private:
    // Inequality constraint
    const SDPIneq <Real> h;
    
    // Workspace for working with h(x) and its derivatives
    mutable typename peopt::SQL <Real>::Vector z_tmp1;
   
    // Extent to which we push for positive definiteness 
    const Real epsilon;

    // Regularization constant on h
    const Real beta;
        
public:
    typedef peopt::Rm <Real> Rm;
    typedef peopt::SQL <Real> SQL;
    typedef typename RmxR <Real>::Vector X_Vector;

    // Set some optimization parameters, allocate memory for the workspace,
    // and grab a reference for the inequality constraint
    Phase1Obj(
        const SparseSDP <Real>& prob,
        const Real& epsilon_,
        const Real& beta_
    ) : h(prob), epsilon(epsilon_), beta(beta_) {
        // Initialize memory for the internal workspace
        initSQL(prob,z_tmp1); 
    }

    // Evaluation 
    double operator () (const X_Vector& x) const {
        // f_x <- 1/2 || y - epsilon ||^2
        Real f_x = Real(0.5) * (x.second-epsilon)*(x.second-epsilon);

        // z_tmp1 <- h(x)
        h(x.first,z_tmp1);

        // f_x <- 1/2 || y - epsilon e ||^2 + beta/2 || h(x) ||^2
        f_x += beta/Real(2.) * SQL::innr(z_tmp1,z_tmp1); 

        // <- f_x
        return f_x;
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad 
    ) const {
        // grad_2 <- y-e
        grad.second=x.second-epsilon;

        // z_tmp1 <- h(x)
        h(x.first,z_tmp1);

        // grad_1 <- h'(x)* h(x)
        h.ps(x.first,z_tmp1,grad.first);

        // grad_1 <- beta h'(x)* h(x)
        Rm::scal(beta,grad.first);
    }

    // Hessian-vector product
    void hessvec(
        const X_Vector& x,
        const X_Vector& dx,
        X_Vector& H_dx
    ) const {
        // H_dx_2 <- dy
        H_dx.second=dx.second;

        // z_tmp1 <- h'(x)dx
        h.p(x.first,dx.first,z_tmp1);

        // H_dx_1 <- h'(x)*h'(x)dx
        h.ps(x.first,z_tmp1,H_dx.first);

        // H_dx_1 <- beta h'(x)*h'(x)dx
        Rm::scal(beta,H_dx.first);
    }
};

// Define a combination of the SDP inequality as well as a piece that
// helps with feasibility
//
// hh(x,y) = [ h(x) >= y e ]
//
template <typename Real>
struct Phase1Ineq : public peopt::VectorValuedFunction <Real,RmxR,peopt::SQL>{
public:
    typedef peopt::Rm <Real> Rm;
    typedef peopt::SQL <Real> SQL;
    
    typedef typename RmxR <Real>::Vector X_Vector;
    typedef peopt::SQL <Real> Z;
    typedef typename Z::Vector Z_Vector;

private:
    const SDPIneq <Real> h;
    mutable typename peopt::SQL <Real>::Vector e;

public:
    // Grab a reference to the SDP inequality, the identity element, and
    // the amount of infeasibility we want to allow
    Phase1Ineq(
        const SparseSDP <Real>& prob
    ) : h(prob) { 
        // Initialize memory for the identity element 
        initSQL(prob,e); 
        SQL::id(e);
    }

    // z=hh(x,y) 
    void operator () (
        const X_Vector& x,
        Z_Vector& z
    ) const {
        // z <- h(x)
        h(x.first,z);

        // z <- h(x) - y e
        SQL::axpy(-x.second,e,z);
    }

    // z=hh'(x,y)(dx,dy)
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Z_Vector& z
    ) const {
        // z <- h'(x)dx
        h.p(x.first,dx.first,z);
       
        // z <- h'(x)dx - dy e 
        SQL::axpy(-dx.second,e,z);
    }

    // xhat=hh'(x,y)*dz
    void ps(
        const X_Vector& x,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        // xhat_1 <- h'(x)*dx
        h.ps(x.first,dz,xhat.first);

        // xhat_2 <- -<dz,e>
        xhat.second = -SQL::innr(dz,e);
    }

    // xhat=(hh''(x,y)(dx,dy)*dz
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        Rm::zero(xhat.first);
        xhat.second=Real(0.);
    }
};

// Creates an initial guess for x
template <typename Real>
bool initPhase1X(const SparseSDP <Real> prob,typename RmxR <Real>::Vector& x){
    // Create some type shortcuts
    typedef typename peopt::Rm <Real> Rm;
    typedef typename peopt::SQL <Real> SQL;

    // Set the size of the primary part of x
    x.first.resize(prob.A.size()-1);
    for(Natural i=0;i<x.first.size();i++)
        x.first[i]=drand48();

    // Create the identity element
    typename SQL::Vector e;
    initSQL <Real> (prob,e);
    SQL::id(e);
    
    // Determine how infeasible we are.  Basically, we find delta such that
    // e + delta h(x) >=0.  Dividing by delta, we have
    //     (1/delta) e + h(x) >= 0
    // ==> h(x) >= (-1/delta) e > (-2/delta) e
    // This gives us a couple of scenarios.  First, we don't have to worry
    // about delta being 0 since e is strictly feasible.  Second, if delta > 0,
    // we can simply set y=-2/delta then we're strictly feasible.  Third,
    // if delta < 0, we're strictly feasible and we can set y=0.
    SDPIneq <Real> h(prob);

    // xx <- x_1
    typename Rm::Vector xx;
        Rm::init(x.first,xx);
        Rm::copy(x.first,xx);

    // h_xx <- h(xx)
    typename SQL::Vector h_xx;
        SQL::init(e,h_xx);
        h(xx,h_xx);

    // Figure out the extent of our infeasibility.  Use the formula above
    // to transform delta into this value.
    Real delta = SQL::srch(h_xx,e);

    // Determine if we're feasible
    bool feasible = delta < Real(0.);

    // Set y so that we're guaranteed to be feasible
    x.second=!feasible ? -Real(2.)/delta : Real(0.);

    // Return whether or not we're feasible
    return feasible;
}

// Creates an initial guess for dx
template <typename Real>
void initPhase1DX(const SparseSDP<Real>prob,typename RmxR <Real>::Vector& dx){
    // First, initialize the perturbation just like x
    initPhase1X <Real> (prob,dx);

    // Randomize the elements in y 
    dx.second=Real(drand48());
}

// Create an initial guess for z
template <typename Real>
void initZ(
    const SparseSDP <Real> prob,
    typename peopt::SQL <Real>::Vector& z
) {
    // Allocate memory for z
    initSQL <Real> (prob,z);

    // Randomize the elements in z
    for(Natural i=0;i<z.data.size();i++)
        z.data[i]=Real(drand48());
}

// Parse the values beta and epsilon for the phase-1 problem.  In addition,
// parse whether or not we want finite difference tests.
template <typename Real>
void parseSDPSettings(
    const peopt::Messaging& msg,
    const std::string& fname,
    Real& epsilon,
    Real& beta,
    bool& fd_tests
) {
    Json::Value root=peopt::json::parse(msg,fname);
    epsilon=Real(root["sdp_settings"].get("epsilon",1.).asDouble());
    beta=Real(root["sdp_settings"].get("beta",1e-5).asDouble());
    fd_tests=root["sdp_settings"].get("fd_tests",false).asBool();
}

// Sets up and runs the problem
int main(int argc,char* argv[]) {
    // Type shortcuts
    typedef double Real;
    typedef peopt::Rm <Real> Rm;
    typedef peopt::SQL <Real> SQL;

    // Check that we have sufficient inputs
    if(argc!=2) {
        std::cerr << "Usage: sdpa_sparse_format <problem>"
            << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Grab the filename
    std::string fname(argv[1]);

    // Grab the settings for the phase-1 problem and whether or not we
    // need to do finite difference tests.
    double epsilon;
    double beta;
    bool phase1_fd_tests;
    bool phase2_fd_tests;

    // Note, we're going to ignore the values of epsilon and beta from this
    // parsing.  Mostly, it's just easier not to have two different routines.
    parseSDPSettings(peopt::Messaging(),"sdpa_sparse_format.peopt",
        epsilon,beta,phase2_fd_tests);

    parseSDPSettings(peopt::Messaging(),"sdpa_sparse_format_phase1.peopt",
        epsilon,beta,phase1_fd_tests);

    // Parse the file sparse SDPA file
    SparseSDP <Real> prob;
    parse_sdpa <Real> (fname,prob);

    // Sort the indices of the resulting problem
    sort_sdp <Real> (prob);

    // Initialize our random seed
    srand48(1);

    // Create an initial guess for the problem
    RmxR <Real>::Vector x;
    bool feasible = initPhase1X <Real> (prob,x);

    // Create the directions for the FD test
    RmxR <Real>::Vector dx;
        initPhase1DX <Real> (prob,dx);
    RmxR <Real>::Vector dxx;
        initPhase1DX <Real> (prob,dxx);

    // Create an initial guess for the inequality multiplier
    SQL::Vector z;
        initZ <Real> (prob,z);

    // Create the phase-1 state 
    peopt::InequalityConstrained <Real,RmxR,peopt::SQL>::State::t
        phase1_state(x,z);

    // Read the parameters from file
    peopt::json::InequalityConstrained <Real,RmxR,peopt::SQL>::read(
        peopt::Messaging(),"sdpa_sparse_format_phase1.peopt",phase1_state);

    // Create the bundle of phase-1 functions
    peopt::InequalityConstrained <Real,RmxR,peopt::SQL>::Functions::t
        phase1_fns;
    phase1_fns.f.reset(new Phase1Obj <Real> (prob,epsilon,beta)); 
    phase1_fns.h.reset(new Phase1Ineq <Real> (prob));

    // Solve the phase-1.  Right now, we always run this even if we start
    // feasible.  The idea is to help find a better starting position, which
    // may or may not be true.
    std::cout << std::endl <<
        "Solving the phase-1 problem for an initial solution." << std::endl;

    if(phase1_fd_tests) {
        // Run some finite difference tests on this problem 
        std::cout << "Finite difference test on the objective." << std::endl;
        peopt::Diagnostics::gradientCheck <> (
            peopt::Messaging(),*phase1_fns.f,x,dx);
        peopt::Diagnostics::hessianCheck <> (
            peopt::Messaging(),*phase1_fns.f,x,dx);
        peopt::Diagnostics::hessianSymmetryCheck <> (
            peopt::Messaging(),*phase1_fns.f,x,dx,dxx);
        
        std::cout << std::endl
            << "Finite difference test on the inequality constraint."
            << std::endl;
        peopt::Diagnostics::derivativeCheck <> (
            peopt::Messaging(),*phase1_fns.h,x,dx,z);
        peopt::Diagnostics::derivativeAdjointCheck <> (
            peopt::Messaging(),*phase1_fns.h,x,dx,z);
        peopt::Diagnostics::secondDerivativeCheck <> (
            peopt::Messaging(),*phase1_fns.h,x,dx,z);
    }

    // Solve the SDP 
    peopt::InequalityConstrained <Real,RmxR,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),phase1_fns,phase1_state);

    // Tell us why the problem converged
    std::cout << "Phase-1 problem converged due to: "
        << peopt::StoppingCondition::to_string(phase1_state.opt_stop)
        << std::endl;

    // Check if we're feasible
    if(phase1_state.x.front().second <= Real(0.)) {
        std::cout << "Phase-1 problem failed to find a feasible solution."
            << std::endl;
        exit(EXIT_FAILURE);
    }

    // Copy the solution back into the variable x
    Rm::copy(phase1_state.x.front().first,x.first);

    // Create the optimization state 
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::State::t
        state(x.first,z);

    // Read the parameters from file
    peopt::json::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::read(
        peopt::Messaging(),"sdpa_sparse_format.peopt",state);

    // Create the bundle of functions
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new SDPObj <Real> (prob));
    fns.h.reset(new SDPIneq <Real> (prob));
    
    // Keep our user informed
    std::cout << std::endl << "Solving the SDP probem: " << fname << std::endl;

    if(phase2_fd_tests) {
        // Run some finite difference tests on this problem 
        std::cout << "Finite difference test on the objective." << std::endl;
        peopt::Diagnostics::gradientCheck <> (
            peopt::Messaging(),*fns.f,x.first,dx.first);
        peopt::Diagnostics::hessianCheck <> (
            peopt::Messaging(),*fns.f,x.first,dx.first);
        peopt::Diagnostics::hessianSymmetryCheck <> (
            peopt::Messaging(),*fns.f,x.first,dx.first,dxx.first);
        
        std::cout << std::endl
            << "Finite difference test on the inequality constraint."
            << std::endl;
        peopt::Diagnostics::derivativeCheck <> (
            peopt::Messaging(),*fns.h,x.first,dx.first,z);
        peopt::Diagnostics::derivativeAdjointCheck <> (
            peopt::Messaging(),*fns.h,x.first,dx.first,z);
        peopt::Diagnostics::secondDerivativeCheck <> (
            peopt::Messaging(),*fns.h,x.first,dx.first,z);
    }

    // Solve the SDP 
    peopt::InequalityConstrained<Real,peopt::Rm,peopt::SQL>::Algorithms::getMin(
        peopt::Messaging(),fns,state);

    // Tell us why the problem converged
    std::cout << "SDP problem converged due to: "
        << peopt::StoppingCondition::to_string(state.opt_stop)
        << std::endl;

    // Return the objective function
    std::cout << "Objective value: " << std::setprecision(16)
        << std::scientific << state.f_x << std::endl;
}
