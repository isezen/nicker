#include <iostream>
#include <cmath>
#include "numdev.h"
using namespace std;

numdev nd;

void test_lap5p_o2() {
    double arr[3][3] ={
        {1,2,3},
        {4,6,6},
        {7,8,9}
    };
    double res_lap = nd.lap5p_o2(&arr[1][1], 3, 1);

    cout << (res_lap==-4?"OK":"ERROR") << "\n";
}

void test_jac9p_o2() {
    int l=10,m=20;
    double psi[l][m], zta[l][m], dx[m], x[l], y[m], JJ[l][m];
    double pi = 4.*atan(1.0);
    double h = 200;
    double dy = h;

    for(int i=0;i<m;i++) dx[i] = dy;
    double yk = 2.*pi/1000;
    double yl = pi/1000;

    x[0] = 0; y[0] = 0;
    for(int i=1;i<l;i++){
        int im1 = i-1;
        for(int j=1;j<m;j++) {
            int jm1 = j-1;
            x[i] = x[im1] + h;
            y[j] = y[jm1] + h;
        }
    }

    double sum = 0;
    for(int i=0;i<l;i++){
        for(int j=0;j<m;j++) {
            // cout << sin(yk*x[i]) << " ";
            psi[i][j] = sin(yk*x[i]) * sin(yl*y[j]) + cos(yl*y[j]);
            zta[i][j] = -(yk*yk+yl*yl)*sin(yk*x[i])*sin(yl*y[j])- yl*yl * cos(yl*y[j]);
            JJ[0][j] = zta[0][j];
            JJ[l-1][j] = zta[l-1][j];
            JJ[i][1] = zta[i][1];
            JJ[i][m-1] = zta[i][m-1];
            sum += pow((zta[i][j]/(l*m)),2);
        }
    }

    for(int j=1;j<m;j++) {
        for(int i=1;i<l;i++) {
            //JJ[i][j] = nd.jac9p_o2(&psi[i][j],&zta[i][j], l, h);
        }
    }

    for(int j=0;j<m;j++){
        for(int i=0;i<l;i++){
            cout << i << " " << j << " " << JJ[i][j] << "\n";
        }
    }

    //double res_jac = nd.jac9p_o2(&p[1][1], &z[1][1], 3, 1);
    //cout << res_jac << "\n";
}

int main() {
    cout << "Hello world!" << endl;
    // test_lap5p_o2();
    test_jac9p_o2();
    return 0;
}
