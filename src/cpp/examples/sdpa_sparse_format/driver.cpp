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

// Stores a sparse SDPA format problem of the form
//
// min b1*x1 + ... + bm*xm
// st  A1*x1 + ... + Am*xm - A0 >= 0
//
// where each A has a block structure with sizes blksizes_1,...,blksizes_nblocks

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

// Store the SDP problem in a sparse format
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
    std::getline (fin,line);
    do {
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

// Define the SDP objective where 
// 
// f(x)=<b,x>
//
template <typename Real>
struct SDPObj : public peopt::ScalarValuedFunction <Real,peopt::Rm> {
private:
    const SparseSDP <Real>& prob;
        
public:
    typedef peopt::Rm <Real> X;
    typedef typename X::Vector X_Vector;

    // Grab a reference to the underlying SDP problem
    SDPObj(SparseSDP <Real>& prob_) : prob(prob_) {}

    // Evaluation 
    double operator () (const X_Vector& x) const {
        return X::innr(prob.b,x);
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad 
    ) const {
        X::copy(prob.b,grad);
    }

    // Hessian-vector product
    void hessvec(
        const X_Vector& x,
        const X_Vector& dx,
        X_Vector& H_dx
    ) const {
        X::zero(H_dx);
    }
};

// Define the SDP inequality where 
//
// h(x) = A1*x1 + ... + Am*xm - A0 >= 0
//
template <typename Real>
struct SDPIneq 
    : public peopt::VectorValuedFunction <Real,peopt::Rm,peopt::SQL> 
{
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
    SDPIneq(SparseSDP <Real>& prob_) : prob(prob_) {}

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

// Stage-1 vector space
template <typename Real>
struct Stage1VS {
private:
    // This is a templated namespace.  Do not allow construction.
    Stage1VS();

    // Create a shortcut for these routines
    typedef peopt::Rm <Real> Rm;

public:
    // Use std::vector as our vector storage
    typedef std::pair <typename peopt::Rm <Real>::Vector,Real> Vector;

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
        x.second *= alpha;
    }

    // y <- alpha * x + y.
    static void axpy(const Real& alpha, const Vector& x, Vector& y) {
        Rm::axpy(alpha,x.first,y.first);
        y.second += alpha*x.second;
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

// Define the stage-1 objective where 
// 
// f(x)=.5*(x0-1)^2
//
template <typename Real>
struct Stage1Obj : public peopt::ScalarValuedFunction <Real,Stage1VS> {
    typedef Stage1VS <Real> X;
    typedef typename X::Vector X_Vector;

    // Evaluation 
    double operator () (const X_Vector& x) const {
        return Real(0.5)*(x.second-Real(1.))*(x.second-Real(1.));
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad 
    ) const {
        X::zero(grad);
        grad.second=x.second-Real(1.);
    }

    // Hessian-vector product
    void hessvec(
        const X_Vector& x,
        const X_Vector& dx,
        X_Vector& H_dx
    ) const {
        X::zero(H_dx);
        H_dx.second=dx.second;
    }
};

// Define the stage-1 inequality where 
//
// hh(x,x0) = h(x) - x0 e
//
// As long as we can find x0 > 0, this guarantees that h(x) >= x0 e > 0.  This
// gives us a strictly feasible solution.
template <typename Real>
struct Stage1Ineq 
    : public peopt::VectorValuedFunction <Real,Stage1VS,peopt::SQL> 
{
public:
    typedef Stage1VS <Real> X;
    typedef peopt::SQL <Real> Z;
    typedef typename X::Vector X_Vector;
    typedef typename Z::Vector Z_Vector;

private:
    const SDPIneq <Real>& h;

public:
    // Grab a reference to the SDP inequality 
    Stage1Ineq(SDPIneq <Real>& h_) : h(h_) {}

    // z=hh(x) 
    void operator () (
        const X_Vector& x,
        Z_Vector& z
    ) const {
        h(x.first,z);
        Z_Vector e;
            Z::init(z,e);
            Z::id(e);
        Z::axpy(-x.second,e,z);
    }

    // z=hh'(x)dx
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Z_Vector& z
    ) const {
        h.p(x.first,dx.first,z);
        Z_Vector e;
            Z::init(z,e);
            Z::id(e);
        Z::axpy(-dx.second,e,z);
    }

    // xhat=hh'(x)*dz
    void ps(
        const X_Vector& x,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        h.ps(x.first,dz,xhat.first);
        Z_Vector e;
            Z::init(dz,e);
            Z::id(e);
        xhat.second=-Z::innr(dz,e);
    }

    // xhat=(hh''(x)dx)*dz
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        X::zero(xhat);
    }
};

#if 0
template <typename ProblemClass>
class DebugManip : public peopt::StateManipulator<ProblemClass> {
public:
    // Application
    virtual void operator () (
        const typename ProblemClass::Functions::t& fns,
        typename ProblemClass::State::t& state,
        peopt::OptimizationLocation::t loc
    ) const {
        typedef typename ProblemClass::X X;
        typedef typename ProblemClass::Z Z;
        typedef typename X::Vector X_Vector;
        typedef typename Z::Vector Z_Vector;
        const Z_Vector& z=state.z.front();
        const Z_Vector& dz=state.dz.front();
        const Z_Vector& h_x=state.h_x.front();

        switch(loc) {
        case peopt::OptimizationLocation::BeforeGetStep: {
            {
                std::ofstream fout("foo_hx.dat");
                for(int i=1;i<=161;i++) {
                    for(int j=1;j<=161;j++) {
                        fout << h_x(1,i,j) << std::endl;
                    }
                }
                fout.close();
            }
            {
                std::ofstream fout("bar_hx.dat");
                for(int i=1;i<=174;i++) {
                    fout << h_x(2,i) << std::endl;
                }
                fout.close();
            }
            {
                std::ofstream fout("foo_z.dat");
                for(int i=1;i<=161;i++) {
                    for(int j=1;j<=161;j++) {
                        fout << z(1,i,j) << std::endl;
                    }
                }
                fout.close();
            }
            {
                std::ofstream fout("bar_z.dat");
                for(int i=1;i<=174;i++) {
                    fout << z(2,i) << std::endl;
                }
                fout.close();
            }
            break;
        } case peopt::OptimizationLocation::BeforeActualVersusPredicted: {
            {
                std::ofstream fout("foo_dz.dat");
                for(int i=1;i<=161;i++) {
                    for(int j=1;j<=161;j++) {
                        fout << dz(1,i,j) << std::endl;
                    }
                }
                fout.close();
            }
            {
                std::ofstream fout("bar_dz.dat");
                for(int i=1;i<=174;i++) {
                    fout << dz(2,i) << std::endl;
                }
                fout.close();
            }
        } default:
            break;
        }
    }
};
#endif

// Sets up and runs the problem
int main(int argc,char* argv[]) {
    // Type shortcuts
    typedef double Real;
    typedef peopt::Rm <Real> X;
    typedef peopt::SQL <Real> Z;
    typedef Stage1VS <Real> XStage1;
    typedef X::Vector X_Vector;
    typedef Z::Vector Z_Vector;
    typedef XStage1::Vector XStage1_Vector;

    // Check that we have sufficient inputs
    if(argc!=2) {
        std::cerr << "Usage: sdpa_sparse_format <problem>" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Grab the filename
    std::string fname(argv[1]);

    // Parse the file
    SparseSDP <Real> prob;
    parse_sdpa <Real> (fname,prob);

    // Sort the indices of the resulting problem
    sort_sdp <Real> (prob);

    // Initialize our random seed
    srand48(1);

    // Create an initial guess for the problem
    X_Vector x(prob.A.size()-1);
    for(Natural i=0;i<x.size();i++)
        x[i]=drand48();
    
    // Create an initial guess for the dual
    std::vector <peopt::Natural> sizes(prob.blk_sizes.size());
    std::vector <peopt::Cone::t> types(prob.blk_sizes.size());
    for(Natural i=0;i<sizes.size();i++) {
        sizes[i]=abs(prob.blk_sizes[i]);
        if(prob.blk_sizes[i]<0)
            types[i]=peopt::Cone::Linear;
        else
            types[i]=peopt::Cone::Semidefinite;
    }
    Z_Vector z(peopt::Messaging(),types,sizes);
    Z::id(z);
    
    // Create a constraint for determining the amount of infeasibility
    SDPIneq <Real> h(prob);

    // Determine how infeasible we are.  Basically, we find alpha such that
    // e + alpha h(x) >=0.  Dividing by alpha, we have
    //     (1/alpha) e + h(x) >= 0
    // ==> h(x) >= (-1/alpha) e > (-2/alpha) e
    // Hence, we look at alpha and if alpha > 0, we're infeasible, but we
    // can solve a stage-1 problem to recover infeasibility.  Otherwise, if
    // alpha < 0, we're feasible and we can proceed to solve the problem.
    // Since e is strictly feasible, alpha will not be zero.
    X::zero(x);
    Z_Vector e;
        Z::init(z,e);
        Z::id(e);
    Z_Vector h_x;
        Z::init(z,h_x);
        h(x,h_x);
    Real alpha = Z::srch(h_x,e);

    // If we're infeasible, do a stage-1 process in order to find a feasible
    // starting solution
    if(alpha > 0) {
        std::cout << std::endl <<
            "Problem is infeasible.  Solving stage-1 problem for feasibility."
            << std::endl;

        // The constraint for the stage-1 problem is h(x) - x0 e >= 0
        // or h(x) >= x0 e.  From above, we know that h(x) > (-2/alpha) e.
        // Hence, if we initialize with x0 = -2/alpha, we should be striclty
        // feasible.
        XStage1_Vector x_stage1;
            X::init(x,x_stage1.first);
            X::copy(x,x_stage1.first);
            x_stage1.second=Real(-2.)/alpha;
#if 0
        // Run some finite difference tests on the stage1 problem
    
        // Create the directions for the FD test
        XStage1_Vector dx_stage1;
            X::init(dx,dx_stage1.first);
            X::copy(dx,dx_stage1.first);
            dx_stage1.second=drand48();
        XStage1_Vector dxx_stage1;
            X::init(dxx,dxx_stage1.first);
            X::copy(dxx,dxx_stage1.first);
            dxx_stage1.second=drand48();

        // Create the functions for the FD test
        Stage1Obj <Real> f_stage1;
        Stage1Ineq <Real> h_stage1(h);

        // Do the actual tests
        peopt::Diagnostics::gradientCheck <> (
            peopt::Messaging(),f_stage1,x_stage1,dx_stage1);
        peopt::Diagnostics::hessianCheck <> (
            peopt::Messaging(),f_stage1,x_stage1,dx_stage1);
        peopt::Diagnostics::hessianSymmetryCheck <> (peopt::Messaging(),
            f_stage1,x_stage1,dx_stage1,dxx_stage1);
        peopt::Diagnostics::derivativeCheck <> (
            peopt::Messaging(),h_stage1,x_stage1,dx_stage1,z);
        peopt::Diagnostics::derivativeAdjointCheck <>(
            peopt::Messaging(),h_stage1,x_stage1,dx_stage1,z);
        peopt::Diagnostics::secondDerivativeCheck <>(
            peopt::Messaging(),h_stage1,x_stage1,dx_stage1,z);
#endif

        // Create the stage-1 state
        peopt::InequalityConstrained <Real,Stage1VS,peopt::SQL>::State::t
            state_stage1(x_stage1,z);

        // Read the parameters from file
        peopt::json::InequalityConstrained <Real,Stage1VS,peopt::SQL>::read(
            peopt::Messaging(),"sdpa_sparse_format_stage1.peopt",state_stage1);

        // Create the stage-1 bundle of functions
        peopt::InequalityConstrained <Real,Stage1VS,peopt::SQL>::Functions::t
            fns_stage1;
        fns_stage1.f.reset(new Stage1Obj <Real>); 
        fns_stage1.h.reset(new Stage1Ineq <Real>(h));

        // Solve the stage-1 problem
        peopt::InequalityConstrained <Real,Stage1VS,peopt::SQL>::Algorithms
            ::getMin(peopt::Messaging(),fns_stage1,state_stage1);

        // Check that the stage-1 problem was successful
        if(state_stage1.x.front().second <= Real(0.)) {
            std::cout << "Stage-1 problem failed to find a feasible solution."
                << std::endl;
            return(EXIT_SUCCESS);
        }

        // Copy this answer into x
        X::copy(state_stage1.x.front().first,x);
    }
   
    // Keep our user informed
    std::cout << std::endl << "Solving the SDP probem: " << fname << std::endl;

#if 0
    // Run some finite difference tests on this problem 

    // Create the directions for the FD test
    X_Vector dx(prob.A.size()-1);
    X_Vector dxx(prob.A.size()-1);
    for(Natural i=0;i<x.size();i++) {
        dx[i]=drand48();
        dxx[i]=drand48();
    }

    // Create the functions for the FD test
    SDPObj <Real> f(prob);

    // Do the actual tests
    peopt::Diagnostics::gradientCheck <> (peopt::Messaging(),f,x,dx);
    peopt::Diagnostics::hessianCheck <> (peopt::Messaging(),f,x,dx);
    peopt::Diagnostics::hessianSymmetryCheck <> (peopt::Messaging(),f,x,dx,dxx);
    peopt::Diagnostics::derivativeCheck <> (peopt::Messaging(),h,x,dx,z);
    peopt::Diagnostics::derivativeAdjointCheck <>(peopt::Messaging(),h,x,dx,z);
    peopt::Diagnostics::secondDerivativeCheck <>(peopt::Messaging(),h,x,dx,z);
#endif
    

    // Create the SDP state 
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Read the parameters from file
    peopt::json::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::read(
        peopt::Messaging(),"sdpa_sparse_format.peopt",state);

    // Create the bundle of functions
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new SDPObj <Real> (prob)); 
    fns.h.reset(new SDPIneq <Real> (prob));
        
    // Solve the SDP 
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Algorithms
        ::getMin(peopt::Messaging(),fns,state);

    // Return the objective function
    std::cout << "Objective value: " << std::setprecision(16)
        << std::scientific << state.f_x << std::endl;
}
