#include "numdev.h"
#include <iostream>
#include <cmath>
using namespace std;

#define VAL(f,i,j) (*(f+(j*m+i)))

numdev::numdev(){}
numdev::~numdev(){}
double numdev::lap5p_o2(double* f, int m, double h) {
    //return (f[0]+f[2]+f[3]+f[4]-4*f[1])/(h*h);
    return (VAL(f,0,1)+VAL(f,-1,0)+VAL(f,0,-1)+VAL(f,1,0)-4*VAL(f,0,0))/(h*h);
}

double numdev::jac9p_o2(double* p, double* z, int m, double h) {

//     cout << VAL(p,-1,-1)<< " " << VAL(p,0,-1)<< " " << VAL(p,1,-1)<< "\n";
//     cout << VAL(p,-1,0) << " " << VAL(p,0,0) << " " << VAL(p,1,0) << "\n";
//     cout << VAL(p,-1,1) << " " << VAL(p,0,1) << " " << VAL(p,1,1) << "\n";

    double jj1 = ((VAL(p,1,0)-VAL(p,-1,0))*(VAL(z,0, 1)-VAL(z,0,-1))
                -((VAL(p,0,1)-VAL(p,0,-1))*(VAL(z,1,0)-VAL(z,-1,0))))/(4*h*h);
    double jj2 = -(VAL(p,1,0)*(VAL(z,1,1)-VAL(z,-1,1))
                 -VAL(p,-1,0)*(VAL(z,1,-1)-VAL(z,-1,-1))
                 -VAL(p,1,0)*(VAL(z,1,1)-VAL(z,1,-1))
                 -VAL(p,-1,0)*(VAL(z,-1,1)-VAL(z,-1,-1)))/(4*h*h);
    double jj3 = -(VAL(z,1,0)*(VAL(p,1,1)-VAL(p,1,-1))
                 -VAL(z,-1,0)*(VAL(p,-1,1)-VAL(p,-1,-1))
                 -VAL(z,0,1)*(VAL(p,1,1)-VAL(p,-1,1))
                 -VAL(z,0,-1)*(VAL(p,1,-1)-VAL(p,-1,-1)))/(4*h*h);
//    cout << jj1 << "\n";
//    cout << jj2 << "\n";
//    cout << jj3 << "\n";

    return (jj1+jj2+jj3)/3;

}

