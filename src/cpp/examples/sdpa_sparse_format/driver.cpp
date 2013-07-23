#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <random>
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
        prob.blk_sizes.emplace_back(blk_size);
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
        prob.b.emplace_back(val);
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
            prob.A[matno][itok(blkno)].is.emplace_back(i);
            prob.A[matno][itok(blkno)].data.emplace_back(entry);
        } else {
            if(i>j) std::swap(i,j);
            prob.A[matno][itok(blkno)].is.emplace_back(i);
            prob.A[matno][itok(blkno)].js.emplace_back(j);
            prob.A[matno][itok(blkno)].data.emplace_back(entry);
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
    typename peopt::SQL <Real>::Vector& x,
    const bool phase1=false
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

    // If we're in phase-1, add the extra cone for feasibility. 
    if(phase1) {
        types.emplace_back(peopt::Cone::Linear);
        sizes.emplace_back(2);
    }

    // Create a new element 
    typename SQL::Vector xx(peopt::Messaging(),types,sizes);

    // Initialize the memory for the user input
    SQL::init(xx,x);
}

// Define the phase-1 objective where 
// 
// f(x,y) = y2 
//
template <typename Real>
struct Phase1Obj : public peopt::ScalarValuedFunction <Real,peopt::Rm> {
    typedef peopt::Rm <Real> Rm;
    typedef typename Rm::Vector X_Vector;

    // We basically have an empty constructor .
    Phase1Obj() {} 

    // Evaluation 
    double operator () (const X_Vector& x) const {
        // Just return y2
        return x.back();
    }

    // Gradient
    void grad(
        const X_Vector& x,
        X_Vector& grad 
    ) const {
        // grad <- 0
        Rm::zero(grad);

        // grad_y2 <- 1
        grad.back()=Real(1.);
    }

    // Hessian-vector product
    void hessvec(
        const X_Vector& x,
        const X_Vector& dx,
        X_Vector& H_dx
    ) const {
        // H_dx <- 0
        Rm::zero(H_dx);
    }
};

// Define a combination of the SDP inequality as well as a piece that
// helps with feasibility
//
// hh(x,y) = [ h(x) >= y e ]
//           [ y2 >= y1 - epsilon ]
//           [ y2 >= -y1 + epsilon ]
//
template <typename Real>
struct Phase1Ineq
    : public peopt::VectorValuedFunction <Real,peopt::Rm,peopt::SQL>
{
public:
    typedef peopt::Rm <Real> Rm;
    typedef peopt::SQL <Real> SQL;
    
    typedef typename Rm::Vector X_Vector;
    typedef peopt::SQL <Real> Z;
    typedef typename Z::Vector Z_Vector;

private:
    // SDP inequality constraint
    const SDPIneq <Real> h;

    // Identity vector
    mutable typename peopt::SQL <Real>::Vector e;
    
    // Extent to which we push for positive definiteness 
    const Real epsilon;

public:
    // Grab a reference to the SDP inequality, the identity element, and
    // the amount of infeasibility we want to allow
    Phase1Ineq(
        const SparseSDP <Real>& prob,
        const Real& epsilon_
    ) : h(prob), epsilon(epsilon_) { 
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
        h(x,z);

        // Get the size of x
        Natural m = x.size()-2;

        // z <- h(x) - y1 e
        SQL::axpy(-x[m],e,z);
        
        // Get the number of cones
        Natural ncones = z.types.size();

        // z_2 <- (-y1 + y2 + epsilon,y1 + y2 - epsilon)
        z(ncones,1) = -x[m] + x[m+1] + epsilon; 
        z(ncones,2) = x[m] + x[m+1] - epsilon; 
    }

    // z=hh'(x,y)(dx,dy)
    void p(
        const X_Vector& x,
        const X_Vector& dx,
        Z_Vector& z
    ) const {
        // z <- h'(x)dx
        h.p(x,dx,z);

        // Get the size of x
        Natural m = x.size()-2;
       
        // z <- h'(x)dx - dy e 
        SQL::axpy(-dx[m],e,z);
        
        // Get the number of cones
        Natural ncones = z.types.size();

        // z_2 <- (-dy1 + dy2, dy1 + dy2)
        z(ncones,1) = -dx[m] + dx[m+1];
        z(ncones,2) = dx[m] + dx[m+1];
    }

    // xhat=hh'(x,y)*dz
    void ps(
        const X_Vector& x,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        // xhat_1 <- h'(x)*dx
        h.ps(x,dz,xhat);
        
        // Get the size of x
        Natural m = x.size()-2;
        
        // Get the number of cones
        Natural ncones = dz.types.size();

        // xhat_2 <- -<e,dz>.  We need to do e first to make sure we don't
        // take the inner product with the last cone in dz.
        xhat[m] = -SQL::innr(e,dz)-dz(ncones,1)+dz(ncones,2);
        xhat[m+1] = dz(ncones,1)+dz(ncones,2);
    }

    // xhat=(hh''(x,y)(dx,dy)*dz
    void pps(
        const X_Vector& x,
        const X_Vector& dx,
        const Z_Vector& dz,
        X_Vector& xhat 
    ) const {
        Rm::zero(xhat);
    }
};

// Creates the ith cannonical vector
template <typename Real>
void create_ei(const Natural& i,std::vector <Real>& ei) {
    for(Natural j=0;j<ei.size();j++) {
        ei[j] = Real(i==j);
    }
}

// Projects X to Rm
template <typename Real,template <typename> class XX>
struct ProjectX {
    virtual typename peopt::Rm <Real>::Vector * operator () (
        typename XX <Real>::Vector& x
    ) const = 0;
    virtual ~ProjectX() {}
};

// Projects Rm to Rm
template <typename Real>
struct ProjectRm : public ProjectX <Real,peopt::Rm> {
    typename peopt::Rm <Real>::Vector * operator () (
        typename peopt::Rm <Real>::Vector& x
    ) const {
        return &x;
    }
};

template <typename Real,template <typename> class XX>
struct SDPPreconditioner : public peopt::Operator <Real,XX,XX> {
private:
    // Create some type shortcuts
    typedef peopt::Rm <Real> Rm;
    typedef XX <Real> X;
    typedef typename X::Vector X_Vector;
    typedef typename Rm::Vector Rm_Vector;

    // Projection from X to Rm
    const std::unique_ptr <ProjectX <Real,XX> > proj;

    // Function modifications used by the inequality constrained problem.  This
    // needs to be a reference to the unique_ptr since the actual function is
    // not initialized yet when we grab the reference.
    const std::unique_ptr
        <peopt::ScalarValuedFunctionModifications <Real,XX> >& f_mod;

    // Current iterate
    const X_Vector& x;

    // Workspace 
    mutable X_Vector x_tmp1;
    mutable X_Vector x_tmp2;

    // Variables used for caching.  The boolean values denote whether or not
    // we've started caching yet.
    mutable std::pair <bool,X_Vector> x_last;

    // Dense matrix that stores the interior point piece of the Hessian and
    // then its Choleski factorization
    mutable std::vector <Real> H;

    // Cannonical vector
    mutable X_Vector ei;

    // Inverse of the condition number of H
    mutable Real invCondH;
 
public:
    SDPPreconditioner(
        ProjectX <Real,XX>* proj_,
        const std::unique_ptr<peopt::ScalarValuedFunctionModifications<Real,XX> >&
            f_mod_,
        const X_Vector& x_
    ) : proj(proj_), f_mod(f_mod_), x(x_), invCondH(1.) {
        X::init(x,x_tmp1);
        X::init(x,x_tmp2);
        X::init(x,x_last.second);
            x_last.first=false;
        Rm_Vector const * const P_x=(*proj)(const_cast <X_Vector&>(x));
        H.resize(P_x->size()*P_x->size());
        X::init(x,ei);
    }

    // Basic application
    void operator () (const X_Vector& dx,X_Vector &PH_dx) const {
        // Determine the size of the projected vector
        Natural m = (*proj)(ei)->size();

        // See if we need to recalculate the preconditioner
        if( peopt::rel_err_cached <Real,XX> (x,x_last) >=
            std::numeric_limits <Real>::epsilon()*1e1
        ){
            // Cache the values
            x_last.first=true;
            X::copy(x,x_last.second);

            // Zero out H and our temp vector
            Rm::zero(H);
            X::zero(x_tmp1);

            // Form H
            for(Natural i=0;i<m;i++) {
                // Create the cannonical vector
                create_ei <Real> (i,*((*proj)(ei))); 

                // Grab the ith column of H
                f_mod->hessvec_step(x,ei,x_tmp1,x_tmp2);

                // Project out the bits in Rm
                Rm_Vector* const P_x_tmp2=(*proj)(x_tmp2);

                // Copy in the column into the correct place
                peopt::copy <Real> (m,&((*P_x_tmp2)[0]),1,&(H[i*m]),1);
            }

            // Find the condition number of H
            Integer info(0);
            std::vector <Real> work(3*m);
            std::vector <Integer> iwork(m);
            peopt::trcon('I','U','N',m,&(H[0]),m,invCondH,&(work[0]),
                &(iwork[0]),info);

            // Find the Choleski factorization of H
            peopt::potrf <Real> ('U',m,&(H[0]),m,info);
        }

        // Start by copying over the direction
        X::copy(dx,PH_dx);

        // If the matrix is well enough conditioned, do a triangular solve
        // for y
        if(invCondH >= std::numeric_limits <Real>::epsilon()*1e3) {
            // Do the triangular solve for y
            Rm_Vector* P_PHdx=(*proj)(PH_dx);
            peopt::trsv <Real> ('U','T','N',m,&(H[0]),m,&((*P_PHdx)[0]),1);
            peopt::trsv <Real> ('U','N','N',m,&(H[0]),m,&((*P_PHdx)[0]),1);
        }
    }
};

// Creates an initial guess for x
template <typename Real>
bool initPhase1X(
    const SparseSDP <Real> prob,
    typename peopt::Rm <Real>::Vector& x
){
    // Create some type shortcuts
    typedef typename peopt::Rm <Real> Rm;
    typedef typename peopt::SQL <Real> SQL;

    // Set the size of the primary part of x
    Natural m = prob.A.size()-1;
    x.resize(m+2);
    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    for(Natural i=0;i<x.size();i++)
        x[i]=Real(dis(gen));

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
    // if delta < 0, we're strictly feasible and we can set y=0.  Finally,
    // if delta=infinity, we're also feasible.
    SDPIneq <Real> h(prob);

    // xx <- x_1
    typename Rm::Vector xx;
        Rm::init(x,xx);
        Rm::copy(x,xx);
        xx.pop_back();
        xx.pop_back();

    // h_xx <- h(xx)
    typename SQL::Vector h_xx;
        SQL::init(e,h_xx);
        h(xx,h_xx);

    // Figure out the extent of our infeasibility.  Use the formula above
    // to transform delta into this value.
    Real delta = SQL::srch(h_xx,e);

    // Determine if we're feasible
    bool feasible = delta < Real(0.)
        || delta > std::numeric_limits <Real>::max();

    // Set y so that we're guaranteed to be feasible
    x[m] = !feasible ? -Real(2.)/delta : Real(0.);
    x[m+1] = Real(fabs(x[m]))*10; 

    // Return whether or not we're feasible
    return feasible;
}

// Creates an initial guess for dx
template <typename Real>
void initPhase1DX(
    const SparseSDP<Real>prob,
    typename peopt::Rm <Real>::Vector& dx
){
    // First, initialize the perturbation just like x
    initPhase1X <Real> (prob,dx);
}

// Create an initial guess for z
template <typename Real>
void initZ(
    const SparseSDP <Real> prob,
    typename peopt::SQL <Real>::Vector& z,
    const bool phase1=false 
) {
    // Allocate memory for z
    initSQL <Real> (prob,z,phase1);

    // Randomize the elements in z
    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    for(Natural i=0;i<z.data.size();i++)
        z.data[i]=Real(dis(gen));
}

// Parse the value epsilon for the phase-1 problem.  In addition,
// parse whether or not we want finite difference tests.
template <typename Real>
void parseSDPSettings(
    const peopt::Messaging& msg,
    const std::string& fname,
    Real& epsilon,
    bool& fd_tests
) {
    Json::Value root=peopt::json::parse(msg,fname);
    epsilon=Real(root["sdp_settings"].get("epsilon",1.).asDouble());
    fd_tests=root["sdp_settings"].get("fd_tests",false).asBool();
}

// Sets up and runs the problem
int main(int argc,char* argv[]) {
    // Type shortcuts
    typedef double Real;
    typedef peopt::Rm <Real> Rm;
    typedef peopt::SQL <Real> SQL;

    // Check that we have sufficient inputs
    if(argc!=4) {
        std::cerr << "Usage: sdpa_sparse_format <problem> <phase-1 parameters> "
            << "<phase-2 parameters>" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Grab the filenames
    std::string fname(argv[1]);
    std::string phase1_params(argv[2]);
    std::string phase2_params(argv[3]);

    // Grab the settings for the phase-1 problem and whether or not we
    // need to do finite difference tests.
    double epsilon;
    bool phase1_fd_tests;
    bool phase2_fd_tests;

    // Note, we're going to ignore the values of epsilon and beta from this
    // parsing.  Mostly, it's just easier not to have two different routines.
    parseSDPSettings(peopt::Messaging(),phase2_params,epsilon,phase2_fd_tests);
    parseSDPSettings(peopt::Messaging(),phase1_params,epsilon,phase1_fd_tests);

    // Parse the file sparse SDPA file
    SparseSDP <Real> prob;
    parse_sdpa <Real> (fname,prob);

    // Sort the indices of the resulting problem
    sort_sdp <Real> (prob);

    // Create an initial guess for the problem
    Rm::Vector x;
    bool feasible = initPhase1X <Real> (prob,x);

    // Create the directions for the FD test
    Rm::Vector dx;
        initPhase1DX <Real> (prob,dx);
    Rm::Vector dxx;
        initPhase1DX <Real> (prob,dxx);

    // Create an initial guess for the inequality multiplier
    SQL::Vector z_phase1;
        initZ <Real> (prob,z_phase1,true);

    // Create the phase-1 state 
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::State::t
        phase1_state(x,z_phase1);

    // Read the parameters from file
    peopt::json::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::read(
        peopt::Messaging(),phase1_params,phase1_state);

    // Create the bundle of phase-1 functions
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Functions::t
        phase1_fns;
    phase1_fns.f.reset(new Phase1Obj <Real> ()); 
    phase1_fns.h.reset(new Phase1Ineq <Real> (prob,epsilon));
    phase1_fns.PH.reset(new SDPPreconditioner <Real,peopt::Rm> (
        new ProjectRm <Real> (),phase1_fns.f_mod,phase1_state.x.front()));

    if(phase1_fd_tests) {
        // Run some finite difference tests on this problem 
        std::cout << std::endl 
            << "Running the finite difference tests on the phase-1 problem."
            << std::endl << std::endl;

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
            peopt::Messaging(),*phase1_fns.h,x,dx,z_phase1);
        peopt::Diagnostics::derivativeAdjointCheck <> (
            peopt::Messaging(),*phase1_fns.h,x,dx,z_phase1);
        peopt::Diagnostics::secondDerivativeCheck <> (
            peopt::Messaging(),*phase1_fns.h,x,dx,z_phase1);
    }

    // Solve the phase-1 problem if we're infeasible.
    if(!feasible) {
        std::cout << std::endl <<
            "Solving the phase-1 problem for an initial solution." << std::endl;

        // Solve the SDP 
        peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Algorithms
            ::getMin(peopt::Messaging(),phase1_fns,phase1_state);

        // Tell us why the problem converged
        std::cout << "Phase-1 problem converged due to: "
            << peopt::StoppingCondition::to_string(phase1_state.opt_stop)
            << std::endl;

        // Check if we're feasible
        Natural m = phase1_state.x.front().size()-2;
        if(phase1_state.x.front()[m] <= Real(0.)) {
            std::cout << "Phase-1 problem failed to find a feasible solution."
                << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Copy the solution back into the variable x
    Rm::copy(phase1_state.x.front(),x);

    // Eliminate the extra entries between phase-1 and 2 
    x.pop_back(); x.pop_back();
    dx.pop_back(); dx.pop_back();
    dxx.pop_back(); dxx.pop_back();
    
    // Create an inequality multiplier
    SQL::Vector z;
        initZ <Real> (prob,z);

    // Create the optimization state 
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::State::t
        state(x,z);

    // Read the parameters from file
    peopt::json::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::read(
        peopt::Messaging(),phase2_params,state);

    // Create the bundle of functions
    peopt::InequalityConstrained <Real,peopt::Rm,peopt::SQL>::Functions::t
        fns;
    fns.f.reset(new SDPObj <Real> (prob));
    fns.h.reset(new SDPIneq <Real> (prob));
    fns.PH.reset(new SDPPreconditioner <Real,peopt::Rm> (
        new ProjectRm <Real> (),fns.f_mod,state.x.front()));
    
    // Run some finite difference tests on this problem 
    if(phase2_fd_tests) {
        std::cout << std::endl 
            << "Running the finite difference tests on the phase-2 problem."
            << std::endl << std::endl;;

        std::cout << "Finite difference test on the objective." << std::endl;
        peopt::Diagnostics::gradientCheck <> (
            peopt::Messaging(),*fns.f,x,dx);
        peopt::Diagnostics::hessianCheck <> (
            peopt::Messaging(),*fns.f,x,dx);
        peopt::Diagnostics::hessianSymmetryCheck <> (
            peopt::Messaging(),*fns.f,x,dx,dxx);
        
        std::cout << std::endl
            << "Finite difference test on the inequality constraint."
            << std::endl;
        peopt::Diagnostics::derivativeCheck <> (
            peopt::Messaging(),*fns.h,x,dx,z);
        peopt::Diagnostics::derivativeAdjointCheck <> (
            peopt::Messaging(),*fns.h,x,dx,z);
        peopt::Diagnostics::secondDerivativeCheck <> (
            peopt::Messaging(),*fns.h,x,dx,z);
    }
    
    // Keep our user informed
    std::cout << std::endl << "Solving the SDP probem: " << fname << std::endl;

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

    // Write out the final answer to file
    peopt::json::InequalityConstrained <Real,peopt::Rm,peopt::SQL>
        ::write_restart(peopt::Messaging(),"sdpa_sparse_format.perst",state);
}
