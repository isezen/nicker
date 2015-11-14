#include <iostream>
#include <cmath>    // Needed for mathematical calculations
#include <cstdlib>  // needed for exit(1)
#include <iomanip>  // needed to use manipulators with parameters (precision, width)
#include "numdev.h"
using namespace std;

numdev nd;

#define VAL(f,i,j) ( *(f+((i)*m+(j))) )
#define SWAP(x,y) int xy;xy=x;x=y;y=xy;
#define loop(i, n, m) for(int i = n; i < m; ++i)
#define loop2(i, imin, imax, j, jmin, jmax) loop(i, imin, imax) loop(j, jmin, jmax)

void test() {
//    int arr[3][4] ={
//        {1,2,3,4},
//        {5,6,7,8},
//        {9,10,11,12}
//    };
    double arr[2][3][4] = { { {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4} },
                     { {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4} } };
    nd.test(&arr[0][0][4],2,3);
}

void test_lap5p_o2() {
//    double arr[3][3] ={{1,2,3},
//                       {4,6,6},
//                       {7,8,9}};
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

void print_arr(double* a, int l, int m, const char* arr_name, bool showColRow) {
    
    int wdth = 14;
    cout << setprecision(6);
    
    if(showColRow) {
      cout << arr_name << "[i,j]" << endl;
      cout << "    ";
      for(int i=0;i<m;++i) cout << setw(wdth) << i;
      cout << endl << "----";
      for(int i=0;i<m;++i) cout << setw(wdth) << "--------------";
      cout << endl;
        
    }
    
    for(int i=0;i<l;++i) {
        if(showColRow) cout << setw(3) << i << "|";
        for(int j=0;j<m;++j){
            double v = VAL(a,i,j);
            if(v!=0)
                cout << setw(wdth) << scientific << VAL(a,i,j);
            else
                cout << setw(wdth) << "0";
        }
        cout << endl;
    }
    cout << endl;
}

void nicker() {
    int l = 39, n = l, m = 61; // n,l: Vertical m: Horizontal dimension
    int n1 = n-1, m1 = m-1;
    double dt                = 1;     // time-step
    int output_time_interval = 60;    // output time interval (seconds)
    int LOOP_MAX             = 200;   // End time of process (dt*loop=3*200=600 sec)
    double gnu               = 0.5;   // eddy kinematic coefficient of heat and momentum
    double dx                = 10;    // grid spacing. 10 meters in horizontal and vertical
    double tneut             = 300;   // potential temperature of neutral environment
    double g                 = 9.81;  // acceleration of gravity
    double alpha             = 1.89;  // relaxation factor
    
    // Array Definitions
    double psi[n][m];    // streamfunction
    double difu[n][m];   // eddy dissipation (m**2/sec)
    double eta[2][n][m]; // vorticity (per sec)
    double t[2][n][m];   // potential temperature excess (deg celcius)
    double q[n][m];      // The non-adiabatic heating (deg/sec)
    double arak[n][m], u[n][m], w[n][m], tmap[n][m], wrk[n][m];
    
    
    // Initialize arrays
    loop2(i, 0, n, j, 0, m) {
        psi[i][j]    = 0;
        eta[0][i][j] = 0;
        eta[1][i][j] = 0;
        t[0][i][j]   = 0;
        t[1][i][j]   = 0;
        q[i][j]      = 0;
        u[i][j]      = 0;
        w[i][j]      = 0;
        arak[i][j]   = 0;
        difu[i][j]   = 0;
        tmap[i][j]   = 0;
        wrk[i][j]    = 0;
    }
    

    loop2(i, 1, 18, j, 11, 32) t[0][i-1][j-1] = 0.5 * cos(M_PI*(i-1)/32) * pow(cos(M_PI*((j-1)-10)/40), 2);
    loop(i, 1, 18)               t[0][i-1][9] = t[0][i-1][10];
    loop2(i, 0, 16, j, 9, 30)      t[0][i][j] = t[0][i][j]/tneut;
    loop2(i, 0, l, j, 0, m)         wrk[i][j] = t[0][i][j];
    
    // The non-adiabatic heating (deg/sec) is defined as:
    loop2(i,0,16,j,8,12) q[i][j] = 2*(4e-03 * cos(M_PI*i/32) * pow(cos(M_PI*(j-10)/4), 2))/M_PI;
    
    int time = 0;
    int loop = 0; // number of time steps. here max 200 = 10 minutes
    int nold = 0, nnew = 1; //swap indices for 3D arrays eta[] and t[]
    do {
        loop++;
        if(loop != 1) {SWAP(nold, nnew);}

        // Change non-adiabatic heating as a function of time
        // after 300 seconds, non-adiabatic heating will be negative.
        if(loop==101) loop2(i,1,18,j,9,14) q[i-1][j-1] = -q[i-1][j-1];

        nd.jac9p_o2(&arak[0][0], &psi[0][0], &eta[nold][0][0],l, m, dx);
        nd.lap5p_o2(&difu[0][0], &eta[nold][0][0], n, m, dx);
        
        loop2(i, 1, n1, j, 1, m1)
            eta[nnew][i][j] = eta[nold][i][j] + dt*(arak[i][j] - g*(t[nold][i+1][j]-t[nold][i-1][j])/(2*dx) + gnu*difu[i][j]);

        loop(i, 0, n) eta[nnew][i][0] = eta[nnew][i][m-1] = 0;
        loop(j, 0, m) eta[nnew][0][j] = eta[nnew][n-1][j] = 0;
        
        // Calculate initial estimate of phi from heat transfer equation
        nd.jac9p_o2(&arak[0][0], &psi[0][0], &t[nold][0][0], l, m, dx);
        nd.lap5p_o2(&difu[0][0], &t[nold][0][0], n, m, dx);
        
        loop2(i, 0, n, j, 0, m) t[nnew][i][j] = t[nold][i][j] + dt*(arak[i][j] + q[i][j]/tneut + gnu*difu[i][j]);

        // Relax eta to get psi.routine RELAX1 solves the poisson equation
        nd.RELAX1(&psi[0][0], &eta[nnew][0][0], n, m, alpha);
        
        // The horizontal velocity and vertical velocity are defined from the stream function

        loop2(i,0,n,j,0,m1) u[i][j] =  (psi[i][j+1] - psi[i][j])/dx;
        loop2(j,0,m,i,0,n1) w[i][j] = -(psi[i+1][j] - psi[i][j])/dx;
        
        // Calculate the final (corrector) estimate of eta and phi for this time step
        // from the predicted vorticity and tempature fields
        nd.jac9p_o2(&arak[0][0], &psi[0][0], &eta[nnew][0][0], l, m, dx);
        nd.lap5p_o2(&difu[0][0], &eta[nnew][0][0], n, m, dx);
        
        loop2(i,1,n1,j,1,m1)
            eta[nnew][i][j] = eta[nold][i][j] + dt*(arak[i][j] + gnu*difu[i][j] - g*(t[nnew][i+1][j]-t[nnew][i-1][j])/(2*dx));
        
        // Relax final estimate of eta to get final psi and consequently u and w fields.
        nd.RELAX1(&psi[0][0], &eta[nnew][0][0], n, m, alpha);
        
        loop2(i,0,n,j,0,m1) u[i][j] =  (psi[i][j+1]-psi[i][j])/dx;
        loop2(j,0,m,i,0,n1) w[i][j] = -(psi[i+1][j]-psi[i][j])/dx;
        
        // Final estimate of phi for this time step.
        nd.jac9p_o2(&arak[0][0], &psi[0][0], &t[nnew][0][0], l, m, dx);
        nd.lap5p_o2(&difu[0][0], &t[nnew][0][0], n, m, dx);
        
        loop2(i,0,n,j,0,m) t[nnew][i][j] = t[nold][i][j] + dt * (arak[i][j] + q[i][j]/tneut + gnu * difu[i][j]);
        
        time += dt;
        if(time % output_time_interval !=0) continue;
        
        loop2(i,0,n,j,0,m) tmap[i][j] = t[nnew][i][j] * tneut;
        loop2(j,0,m,i,0,l)  wrk[i][j] = tmap[i][j];
        
        cout << time << endl;
        print_arr(&wrk[0][0], l, m, "wrk", false);
        
        loop(j,0,m) w[n][j] = w[n1][j];
        loop(i,0,n) u[i][m] = u[i][m1];

    }while(loop<LOOP_MAX);

}

int main() {
//    test_jac9p_o2();
    nicker();
    //test();
    return 0;
}
