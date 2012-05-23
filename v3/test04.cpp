#include "linalg.h"
#include "vspaces.h"
#include "peopt.h"
#include <cstdlib>
#include <fstream>

// Fills a vector with random data
template <typename T>
void randomize(typename peopt::SQL <T>::Vector& X) {
    // Loop over all the blocks in the variable
    for(unsigned int blk=1;blk<=X.numBlocks();blk++) {

        // Get the size of the block
        unsigned int m=X.blkSize(blk);

        // Randomize the block
        switch(X.blkType(blk)) {
        case peopt::Cone::Linear:
        case peopt::Cone::Quadratic:
            for(unsigned int i=1;i<=m;i++)
                X(blk,i)=drand48();
            break;
        case peopt::Cone::Semidefinite:
            for(unsigned int j=1;j<=m;j++)
                for(unsigned int i=1;i<=j;i++){
                    X(blk,i,j)=drand48();
                    X(blk,j,i)=X(blk,i,j);
                }
            break;
        }
    }
}

// Makes a vector feasible in a very simple manner
template <typename T>
void feasibilize(typename peopt::SQL <T>::Vector& X) {
    // Loop over all the blocks in the variable
    for(unsigned int blk=1;blk<=X.numBlocks();blk++) {

        // Get the size of the block
        unsigned int m=X.blkSize(blk);

        // Make the block feasible 
        switch(X.blkType(blk)) {

        // For the linear blocks, we just negate any negative elements.  For
        // zero elements, we set them to a random number between 0.5 and 1.5.
        case peopt::Cone::Linear:
            for(unsigned int i=1;i<=m;i++) {
                if(X(blk,i) < 0)
                    X(blk,i) = -X(blk,i);
                else if(X(blk,i)==0)
                    X(blk,i)=drand48()+T(0.5);
            }
            break;

        // For quadratic blocks, we determine ||xbar|| and if x0 is not larger
        // than this, we set x0 to be ||xbar||+rand where random is a random
        // number between 0.5 and 1.5.
        case peopt::Cone::Quadratic: {

            T norm_xbar= peopt::sqrt <T> (
                peopt::dot <T> (m-1,&(X(blk,2)),1,&(X(blk,2)),1));

            if(X(blk,1) <= norm_xbar)
                X(blk,1) = norm_xbar + drand48() + T(0.5);
        }

        // For semidefinite blocks, we find the smallest eigenvalue.  If this
        // is negative, then we do an identity shift by this amount plus a 
        // random number between 1 and 2.
        case peopt::Cone::Semidefinite: {

            T lambda = peopt::lanczos(m,&(X(blk,1,1)),m,1e-2);

            if(lambda <= T(0.)) {
                T extra = drand48() + T(1.);
                for(unsigned int i=1;i<=m;i++) 
                    X(blk,i,i)+=-lambda+extra;
            }
            break;
        } }
    }
}


int main() {
    // Check some simple algebra
    double x[4]={1.,2.,3.,4.};
    double y[4];
    peopt::copy <double> (4,x,1,y,1);
    peopt::axpy <double> (4,1.0,x,1,y,1);

    // Create an easy shortcut for the vector space
    typedef peopt::SQL <double> Vspace;

    // Create an SQL object
    std::vector <peopt::Cone::t> types;
    std::vector <unsigned int> sizes;
    types.push_back(peopt::Cone::Linear);
        sizes.push_back(3);
    types.push_back(peopt::Cone::Semidefinite);
        sizes.push_back(5);
    types.push_back(peopt::Cone::Quadratic);
        sizes.push_back(7);
    types.push_back(peopt::Cone::Semidefinite);
        sizes.push_back(11);
    Vspace::Vector XX(peopt::Messaging(),types,sizes);
    Vspace::Vector YY; Vspace::init(XX,YY);
    Vspace::Vector ZZ; Vspace::init(XX,ZZ);
    Vspace::Vector WW; Vspace::init(XX,WW);

    // Randomize the data
    randomize <double>(XX); 
    randomize <double>(YY); 

#if 0
    // Fill in the SDP blocks with random data
    for(unsigned int j=1;j<=X.blkSize(1);j++)
        for(unsigned int i=1;i<=j;i++) {
            X(2,i,j)=drand48();
            X(2,j,i)=X(2,i,j);
        }
#endif

    // Find the Schur decomposition of the first SDP block
    std::vector <double> D1;
    std::vector <double> V1;
    peopt::SQL <double>::get_schur(XX,2,V1,D1);

    // Do it again to make sure we're getting the cached copy
    std::vector <double> D2;
    std::vector <double> V2;
    peopt::SQL <double>::get_schur(XX,2,V2,D2);

    // Create a random tridiagonal matrix
    unsigned int m=40;
    std::vector <double> D(m);
    std::vector <double> E(m);
    for(unsigned int i=0;i<m;i++){
        D[i]=drand48(); 
        E[i]=drand48();
    }

    // Figure out the workspaces for the eigenvalues and eigenvectors
    std::vector <int> isuppz(2*m);
    std::vector <double> work(1);
    std::vector <int> iwork(1);
    int lwork=-1;
    int liwork=-1;
    int info;
    int nevals;
    int nzc=0;
    std::vector <double> W(m);
    std::vector <double> Z(m*m);
    peopt::stemr <double> ('V','A',m,&(D[0]),&(E[0]),double(0.),double(0.),0,0,
        nevals,&(W[0]),&(Z[0]),m,m,&(isuppz[0]),nzc,&(work[0]),lwork,
        &(iwork[0]),liwork,info);

    // Resize the workspace 
    lwork = int(work[0])+1;
    work.resize(lwork);
    liwork = iwork[0];
    iwork.resize(liwork);

    // Find the eigenvalues and vectors 
    peopt::stemr <double> ('V','A',m,&(D[0]),&(E[0]),double(0.),double(0.),0,0,
        nevals,&(W[0]),&(Z[0]),m,m,&(isuppz[0]),nzc,&(work[0]),lwork,
        &(iwork[0]),liwork,info);
    
    // Create a random matrix
    std::vector <double> A(m*m);
    for(unsigned int j=1;j<=m;j++)
        for(unsigned int i=1;i<=j;i++) {
            A[peopt::ijtok(i,j,m)]=drand48();
            A[peopt::ijtok(j,i,m)]=A[peopt::ijtok(i,j,m)];
        }

    // Find its smallest eigenvalue
    double eval_min=peopt::lanczos(m,&(A[0]),20,1e-1);

    // Write this matrix to file for checking
    std::ofstream fout("junk.txt");
    for(unsigned int j=1;j<=m;j++) {
        for(unsigned int i=1;i<=m;i++) {
            fout << A[peopt::ijtok(i,j,m)] << ' ';
        }
        fout << std::endl;
    }
    fout.close();

    // Find Z = X o Y
    Vspace::prod(XX,YY,ZZ);

    // Find inv(L(X)) Z, hopefully this is Y
    Vspace::linv(XX,ZZ,WW);

    // Find the relative error between W and Y
    Vspace::Vector diff;  Vspace::init(XX,diff);
    Vspace::copy(WW,diff);
    Vspace::axpy(-1.,YY,diff);
    double abs_err=std::sqrt(Vspace::innr(diff,diff));

    // Find the barrier at the identity element 
    Vspace::Vector VV;  Vspace::init(XX,VV);  Vspace::id(VV);
    double barr=Vspace::barr(VV);

    // Generate some strictly feasible vector
    Vspace::Vector F; Vspace::init(XX,F);
    randomize <double> (F);
    feasibilize <double> (F);

    // Find the barrier at this feasible point
    barr=Vspace::barr(F);

    // Do a line search with F as the base in some random direction
    double alpha = Vspace::srch(XX,F);

    // Do a step in this direction with the length alpha
    Vspace::axpy(alpha,XX,F);
    
    // Do another line search.  If everything is as it is supposed to be,
    // the parameter here should be zero.
    alpha = Vspace::srch(XX,F);

    // Create a place for a break point
    int junk=1;
}
