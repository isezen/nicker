#include "numdev.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

#define VAL(f,i,j) ( *(f+((i)*m+(j))) )

numdev::numdev(){}
numdev::~numdev(){}
void numdev::lap5p_o2(double* difu, double* t, int n, int m, double h) {
    int n1 = n-1, m1 = m-1;
    double da = h*h;
    for(int i=1; i<n1; ++i)
        for(int j=1; j<m1; ++j)
            VAL(difu, i, j) = (VAL(t,i,j+1) + VAL(t, i, j-1) + VAL(t, i+1, j) + VAL(t, i-1, j) - 4*VAL(t,i,j))/da;

    // Side boundaries
    for(int j=1; j<m1; ++j) {
        VAL(difu, 0, j)   = (2*VAL(t,1,j)    + VAL(t, 0, j+1)   + VAL(t, 0, j-1)   - 4*VAL(t,0,j))/da;
        VAL(difu, n-1, j) = (2*VAL(t,n1-1,j) + VAL(t, n-1, j+1) + VAL(t, n-1, j-1) - 4*VAL(t,n-1,j))/da;
    }

    // Top and bottom boundaries
    for(int i=1; i<n1; ++i) {
        VAL(difu, i, 0)   = (2*VAL(t,i,1)    + VAL(t, i-1, 0)   + VAL(t, i+1, 0)   - 4*VAL(t,i,0))/da;
        VAL(difu, i, m-1) = (2*VAL(t,i,m1-1) + VAL(t, i-1, m-1) + VAL(t, i+1, m-1) - 4*VAL(t,i,m-1))/da;
    }

    double db = 2/da;
    // Four Corners
    VAL(difu, 0, 0)     = db*(VAL(t,0,1)      + VAL(t,1,0)      - 2*VAL(t,0,0));
    VAL(difu, 0, m-1)   = db*(VAL(t,0,m1-1)   + VAL(t,1,m-1)    - 2*VAL(t,0,m-1));
    VAL(difu, n-1, m-1) = db*(VAL(t,n1-1,m-1) + VAL(t,n-1,m1-1) - 2*VAL(t,n-1,m-1));
    VAL(difu, n-1, 0)   = db*(VAL(t,n1-1,0)   + VAL(t,n-1,1)    - 2*VAL(t,n-1,0));

}

void numdev::jac9p_o2(double* a, double* p, double* z, int l, int m, double h) {
    j1_first_col(a, p, z, l, m, h);
    j1_last_col(a, p, z, l, m, h);
    j1(a, p, z, l, m, h);
}

void numdev::j1_first_col(double* a, double* p, double* z, int l, int m, double h) {
    for(int i=0; i<l; ++i) {
        int im1, ip1;
        if(i==0 || i==(l-1)) {
            im1 = l-2; ip1 = 1;
        }else{
            im1 = i-1; ip1 = i+1;
        }
       VAL(a,i,0) = (VAL(p,i,0)   + VAL(p,ip1,0)  - VAL(p,i,1)    - VAL(p,ip1,1)) * (VAL(z,i,0)   + VAL(z,ip1,0)) -
                    (VAL(p,im1,0) + VAL(p,i,0)    - VAL(p,im1,1)  - VAL(p,i,1))   * (VAL(z,im1,0) + VAL(z,i,0)) +
                    (VAL(p,ip1,0) + VAL(p,ip1,1)  - VAL(p,im1,0)  - VAL(p,im1,1)) * (VAL(z,i,0)   + VAL(z,i,1)) +
                    (VAL(p,ip1,0) - VAL(p,i,1))   * (VAL(z,i,0)   + VAL(z,ip1,1)) +
                    (VAL(p,i,1)   - VAL(p,im1,0)) * (VAL(z,im1,1) + VAL(z,i,0));
        VAL(a,i,0) /= (12*h*h);
    }
}

void numdev::j1_last_col(double* a, double* p, double* z, int l, int m, double h) {
    for(int i=0; i<l; ++i) {
        int im1, ip1;
        if(i==0 || i==(l-1)) {
            im1 = l-2; ip1 = 1;
        }else{
            im1 = i-1; ip1 = i+1;
        }
       VAL(a,i,m-1) = (VAL(p,i,m-2)   + VAL(p,ip1,m-2)  - VAL(p,i,m-1)      - VAL(p,ip1,m-1)) * (VAL(z,i,m-1)   + VAL(z,ip1,m-1)) -
                      (VAL(p,im1,m-2) + VAL(p,i,m-2)    - VAL(p,im1,m-1)    - VAL(p,i,m-1))   * (VAL(z,im1,m-1) + VAL(z,i,m-1)) -
                      (VAL(p,ip1,m-2) + VAL(p,ip1,m-1)    - VAL(p,im1,m-2)  - VAL(p,im1,m-1)) * (VAL(z,i,m-2) + VAL(z,i,m-1)) -
                      (VAL(p,i,m-2)   - VAL(p,im1,m-1))   * (VAL(z,im1,m-2) + VAL(z,i,m-1)) -
                      (VAL(p,ip1,m-1)   - VAL(p,i,m-2))   * (VAL(z,i,m-1)     + VAL(z,ip1,m-2));
        VAL(a,i,m-1) /= (12*h*h);
    }
}

void numdev::j1(double* a, double* p, double* z, int l, int m, double h) {
    for(int i=0; i<l; ++i) {
        int im1, ip1;
        if(i==0 || i==(l-1)) {
            im1 = l-2; ip1 = 1;
        }else{
            im1 = i-1; ip1 = i+1;
        }
        for(int j=1; j<(m-1); j++) {
            VAL(a,i,j) = (VAL(p,i,j-1)   + VAL(p,ip1,j-1)  - VAL(p,i,j+1)    - VAL(p,ip1,j+1) ) * (VAL(z,ip1,j) - VAL(z,i,j)   ) +
                         (VAL(p,im1,j-1) + VAL(p,i,j-1)    - VAL(p,im1,j+1)  - VAL(p,i,j+1)   ) * (VAL(z,i,j)   - VAL(z,im1,j) ) +
                         (VAL(p,ip1,j)   + VAL(p,ip1,j+1)  - VAL(p,im1,j)    - VAL(p,im1,j+1) ) * (VAL(z,i,j+1) - VAL(z,i,j)   ) +
                         (VAL(p,ip1,j-1) + VAL(p,ip1,j)    - VAL(p,im1,j-1)  - VAL(p,im1,j)   ) * (VAL(z,i,j)   - VAL(z,i,j-1) ) +
                         (VAL(p,ip1,j)   - VAL(p,i,j+1))   * (VAL(z,ip1,j+1) - VAL(z,i,j)     ) +
                         (VAL(p,i,j-1)   - VAL(p,im1,j))   * (VAL(z,i,j)     - VAL(z,im1,j-1) ) +
                         (VAL(p,i,j+1)   - VAL(p,im1,j))   * (VAL(z,im1,j+1) - VAL(z,i,j)     ) +
                         (VAL(p,ip1,j)   - VAL(p,i,j-1))   * (VAL(z,i,j)     - VAL(z,ip1,j-1) );
            VAL(a,i,j) /= (12*h*h);
        }
    }
}

void numdev::RELAX1(double* p, double* eta, int n, int m, double alpha) {
    double eps = 0.001;
    double dx  = 10;
    int nscan  = 0;
    int n1     = n-1;
    int m1     = m-1;
    
    for(;;) {
        nscan++;
        double rmax = 0;
        int isave = 0, jsave =0;
        for(int i=1; i<n1; ++i) {
            for(int j=1; j<m1; ++j) {
                double r1 =  0.25 * (VAL(p,i+1,j) + VAL(p,i-1,j) + VAL(p,i,j+1) + VAL(p,i,j-1));
                double r2 = -VAL(p,i,j) -  0.25 * dx * dx * VAL(eta,i,j);
                double r  = r1 + r2;
                if(rmax < abs(r)) {
                    isave = i;
                    jsave = j;
                    rmax  = abs(r);
                }
                VAL(p,i,j) = VAL(p,i,j) + alpha * r;
            }
        }
        if(nscan == 300){
            cout << "No Convergence achieved\n";
            exit(1);
        }
        double ps = abs(VAL(p,isave,jsave));
        if(((rmax/ps-eps)<0)) break;
    }
    //while((rmax/ps-eps)>=0);
}

void numdev::test(double* a, int l, int m) {
    for(int i=0;i<l;i++){
        cout << endl;
        for(int j=0;j<m;j++)
            cout << VAL(a,i,j) << " ";
    }
}
