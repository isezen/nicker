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
    //double res_lap = nd.lap5p_o2(&arr[1][1], 3, 1);

    //cout << (res_lap==-4?"OK":"ERROR") << "\n";
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

void nicker() {
    int l = 39, n = l, m = 61;
    int n1 = n-1, m1 = n-1, n2 = n-2, m2 = m-2;
    double dt = 3;
    int loopt = 40;
    double gnu = 0.5;
    double dx = 10, e = dx;
    double tneut = 300;
    double g = 9.81;
    int nold = 1;
    int nnew = 2;
    double alpha = 1.89;

    double psi[n][m], q[n][m], arak[n][m], difu[n][m], u[n][m], w[n][m], tmap[n][m], wrk[n][m];
    double eta[n][m][2], t[n][m][2];

    for(int i=0; i<n1;i++){
        for(int j=0; j<m1; j++){
            psi[i][j] = 0;
            eta[i][j][nold] = 0;
            t[i][j][nold] = 0;
            q[i][j] = 0;
            u[i][j] = 0;
            w[i][j] = 0;
        }
    }

    for(int i=0;i<17;i++)
        for(int j=10;j<31;j++)
            t[i][j][nold] = 0.5*cos(M_PI*(i-1)/32)*pow((cos(M_PI*((j-1)-10)/40)),2);

    for(int i = 0;i<17;i++) t[i][10][nold] = t[i][11][nold];

    // The non-adiabatic heating (deg/sec) is defined as:
    for(int i=0;i<17;i++)
        for(int j=8;j<13;j++)
            q[i][j] = 4.e-03*cos(M_PI*(i-1)/32)*pow((cos(M_PI*((j-1)-10)/4)),2);

    int time = 0;


    for(int i=0;i<l;i++) for(int j=0;j<m;j++) wrk[i][j] = t[i][j][nold];

    for(int i=0;i<16;i++) for(int j=9;j<30;j++) t[i][j][nold]=t[i][j][nold]/tneut;

    int loop=0;
    do{
        loop++;
        if(loop!=1){
            int nsave=nold;
            nold=nnew;
            nnew=nsave;
        }
    // Define the heating as a function of time
        for(int j=8;j<13;j++){
            for(int i=0;i<17;i++){
                double qq=4.e-03*cos(M_PI*(i-1)/32)*pow((cos(M_PI*(j-1)-10/4)),2);
                q[i][j]=2*qq/M_PI;
                if(loop>100 && loop<=200) q[i][j]=-q[i][j];
            }
        }
        nd.jac9p_o2(&arak[0][0], &psi[0][0], &eta[0][0][nold],l,m,dx);
        nd.lap5p_o2(&difu[0][0],&eta[0][0][nold],n,m,dx);

        for(int i=1;i<n1;i++)
            for(int j=1;j<m1;j++)
                eta[i][j][nnew]=eta[i][j][nold]+dt*arak[i][j]+dt*gnu*difu[i][j]-dt*g*(t[i+1][j][nold]-t[i-1][j][nold])/(2*dx);

        for(int i=0;i<n;i++) eta[i][1][nnew]=eta[i][m][nnew]=0;
        for(int j=0;j<m;j++) eta[1][j][nnew]=eta[n][j][nnew]=0;

    // Calculate initial estimate of phi from heat transfer equation

        nd.jac9p_o2(&arak[0][0], &psi[0][0], &t[0][0][nold],l,m,dx);
        nd.lap5p_o2(&difu[0][0],&t[0][0][nold],n,m,dx);

        for(int i=1;i<n;i++)
            for(int j=1;j<m;j++)
                t[i][j][nnew]=t[i][j][nold]+dt*arak[i][j]+dt*q[i][j]/tneut + dt*gnu*difu[i][j];

    // Relax eta to get psi.routine RELAX1 solves the poisson equation

        nd.RELAX1(&psi[0][0],&eta[0][0][nnew],n,m,alpha);

    // The horizontal velocity and vertical velocity are defined from the stream function

        for(int i=0;i<n;i++) for(int j=0;j<m1;j++) u[i][j]=(psi[i][j+1]-psi[i][j]/dx);
        for(int j=0;j<m;j++) for(int i=0;i<n1;i++) w[i][j]=-(psi[i+1][j]-psi[i][j]/dx);

    // Calculate the final (corrector) estimate of eta and phi for this time step
    // from the predicted vorticity and tempature fields

        nd.jac9p_o2(&arak[0][0], &psi[0][0], &eta[0][0][nnew],l,m,dx);
        nd.lap5p_o2(&difu[0][0],&eta[0][0][nnew],n,m,dx);

        for(int i=1;i<n1;i++)
            for(int j=1;j<m1;j++)
                eta[i][j][nnew]=eta[i][j][nold]+dt*arak[i][j]+dt*gnu*difu[i][j]-dt*g*(t[i+1][j][nnew]-t[i-1][j][nnew])/(2*dx);

    // Relax final estimate of eta to get final psi and consequently u and w fields.

        nd.RELAX1(&psi[0][0],&eta[0][0][nnew],n,m,alpha);

        for(int i=0;i<n;i++) for(int j=0;j<m1;j++) u[i][j]=(psi[i][j+1]-psi[i][j]/dx);
        for(int j=0;j<m;j++) for(int i=0;i<n1;i++) w[i][j]=-(psi[i+1][j]-psi[i][j]/dx);

    // Final estimate of phi for this time step.

        nd.jac9p_o2(&arak[0][0], &psi[0][0], &t[0][0][nnew],l,m,dx);
        nd.lap5p_o2(&difu[0][0],&t[0][0][nnew],n,m,dx);

        for(int i=0;i<n;i++)
            for(int j=0;j<m;j++)
                t[i][j][nnew]=t[i][j][nold]+dt*arak[i][j]+dt*q[i][j]/tneut + dt*gnu*difu[i][j];

        time+=dt;
        if(loop % loopt !=0) continue;
        time=3*loop;
        for(int i=0;i<n;i++) for(int j=0;j<m;j++) tmap[i][j]=t[i][j][nnew]*tneut;
        for(int j=0;j<m;j++) for(int i=0;i<l;i++) wrk[i][j]=tmap[i][j];

        for(int j=0;j<m;j++) w[n][j]=w[n1][j];
        for(int i=0;i<n;i++) u[i][m]=u[i][m1];






    }while(loop<200);

}

int main() {
    cout << "Hello world!" << endl;
    //test_jac9p_o2();
    nicker();
    return 0;
}
