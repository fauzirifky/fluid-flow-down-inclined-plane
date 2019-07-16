//Persamaan Saint-Venant : Simulasi Roll Waves dengan Angin

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include <algorithm>
// # define M_PI           3.14159265358979323846  /* pi */
using namespace std;

float fun1(float a, float h0, float k) {
    return( (0 + 0.05 * sin(k * a)));
}

float fun2(float a, float h0, float k, float co, float F) {
    return( sqrt( (0 + 0.00 * sin(k * a))) );
}

int main() {
    // double L = 1, dx = 0.000833333333333, dt = 0.000083333;
    // double L = 1, dx = 0.000833333333333, dt = 0.00046667; // 0.56dx;
    double L = 1, dx = 0.002, dt = dx*0.4;
    // double L = 1, dx = 0.000833333333333, dt = dx*1.2;
    // double L = 1, dx = 0.000833333333333, dt = dx*0.78;
    int X = round(L/dx);
    int T = round(20/dt);
    double x[X-1], h[X - 1],  u[X - 1], hi[X - 1], ui[X - 1], hh[X - 1], uu[X - 1], hn[X - 1], un[X - 1], H1[X - 1], H2[X - 1], H3[X - 1], HH[X - 1], U1[X - 1], U2[X - 1], U3[X - 1], UU[X - 1], Hs[4][X - 1], Us[4][X - 1], umax;
    double ht[23][X-1],  ta, tb, tc, td;
    double theta, co, si, u0,nu,mu;
    int n,j, vv=1;
    int pembagi = 100, maks = T/100;
    int v = 1, skip = round( T / 1000 );
    
    if (skip == 0){
        skip = 1;
    }
    
    // space & time discretization
    for (j = 0; j <= X ; j++) {
        x[j + 1] = x[j] + dx;
    }
    
    
    
    // time
    int a = 3, b = round(T / 3), c = round(2 * T / 3), d = T - 1;
    
    
    // Parameters
    double k = 2 * M_PI, F = 3.0 , R = 500, Cf = 0.006, g = 9.8, h0 = 0.0042, sigma = 0.5;
    printf("Simulasi Roll Waves hingga T = %4.2f\n", T*dt);
    printf("k = %d pi, F = %f, R = %f, sigma = %f\n", int(k / M_PI), F, R, sigma);
    
    /// calculate several paramaters
    theta = atan(F * Cf);
    co = g * cos(theta);
    si = g * sin(theta);
    u0 = sqrt( h0);
    nu = 0; u0 / (R * L);
    mu = 0.0005;0.001/5; //-0.075;//0.1 ;
    
    /// initial condition
    x[0] = dt;
    for (j = 0; j <= X - 1; j++) {
        h[j] = fun1(x[j], h0,  k );
        u[j] = 0.5*fun1(x[j], h0,  k ); //fun2(x[j], h0,  k,  co,  F);
    }
    
    for (j = 0; j <= X - 1; j++) {
        hi[j] = fun1(x[j], h0,  k );
        ui[j] = 0.5*fun1(x[j], h0,  k ); //fun2(x[j], h0,  k,  co,  F);
    }
    
    for (j = 0; j <= X - 1; j++) {
        H1[j] = 0;
        U1[j] = 0;
        H2[j] = 0;
        U2[j] = 0;
    }
    // Export Plot
    ofstream data_plot, umaxx;
    data_plot.open ("data_plot.csv");
    umaxx.open ("umax.csv");
    
    for (int j = 0; j <= X - 2; j++){
        data_plot << x[j] <<"," ;
    }
    data_plot << x[X - 1] <<"\n";
    
    for (int j = 0; j <= X - 2; j++){
        data_plot << h[j] <<"," ;
    }
    data_plot << h[X - 1] <<"\n";
    

    
    /// main iteration starts
    for (n = 0; n <= T; n++) {
        
        // Predictor Operator
        for (j = 1; j <= X - 2 ; j++){
            H3[j] = - (u[j + 1] - u[j - 1]) / (2 * dx) - (h[j + 1] - h[j - 1]) / (2 * dx);
            U3[j] = - (u[j + 1] - u[j - 1]) / (2 * dx) - (1/F) * (h[j + 1] - h[j - 1]) / (2 * dx) + (h[j] - 2*u[j])/sigma + (1/R) * (u[j + 1] - 2 * u[j] + u[j - 1]) / (dx * dx) + mu * (h[j + 1] - 2 * h[j] + h[j - 1]) / (dx * dx);
        }
        
        // Periodic Boundary Condition
        /*
         H3[0] = H3[X - 2];
         H3[X - 1] = H3[1];
         U3[0] = U3[X - 2];
         U3[X - 1] = U3[1];
         */
        
        // Predictor steps
        for (j = 1; j <= X - 2; j++){
            hh[j] = h[j] + dt * (23 * H3[j] - 16 * H2[j] + 5 * H1[j])/12;
            uu[j] = u[j] + dt * (23 * U3[j] - 16 * U2[j] + 5 * U1[j])/12;
        }
        
        hh[X - 1] = hh[1];
        hh[0] = hh[X - 2];
        
        uu[X - 1] = uu[1];
        uu[0] = uu[X - 2];
        
        
        // Corrector Operator
        for (j = 1; j <= X - 2; j++) {
            HH[j] = - (uu[j + 1] - uu[j - 1]) / (2 * dx) - (hh[j + 1] - hh[j - 1]) / (2 * dx);
            UU[j] = - (uu[j + 1] - uu[j - 1]) / (2 * dx) - (1/F) * (hh[j + 1] - hh[j - 1]) / (2 * dx) +  (hh[j] - 2*uu[j] ) / sigma + (1/R) * (uu[j + 1] - 2 * uu[j] + uu[j - 1]) / (dx * dx) + mu * (hh[j + 1] - 2 * hh[j] + hh[j - 1]) / (dx * dx);
        }
        
        // Periodic Boundary Condition
        /*
         HH[0] = HH[X - 2];
         HH[X - 1] = HH[1];
         UU[0] = UU[X - 2];
         UU[X - 1] = UU[1];
        */
        
        // Corrector steps
        for (j = 1; j <= X - 2; j++) {
            hn[j] = h[j] + dt * (5 * HH[j] + 8 * H3[j] - H2[j])/12;
            un[j] = u[j] + dt * (5 * UU[j] + 8 * U3[j] - U2[j])/12;
        }
        
        hn[X - 1] = hn[1];
        hn[0] = hn[X - 2];
        
        un[X - 1] = un[1];
        un[0] = un[X - 2];
        
        
        
        //printf("hhX-1 =  %f \n", hh[X-1]);
        // update operator
        
        for (j = 0; j <= X - 1; j++) {
            H1[j] = H2[j];
            H2[j] = H3[j];
            U1[j] = U2[j];
            U2[j] = U3[j];
        }
        
        // Time Capture
        if (n == 0) {
            //printf("n = %d\n", n);
            for (int i = 0; i <= X - 1; i++) {
                Hs[0][i] = hn[i];
                Us[0][i] = un[i];
            }
        } else if (n == b) {
            //printf("n = %d\n", n);
            for (int i = 0; i <= X - 1; i++) {
                Hs[1][i] = hn[i];
                Us[1][i] = un[i];
            }
        } else if (n == c) {
            //printf("n = %d\n", n);
            for (int i = 0; i <= X - 1; i++) {
                Hs[2][i] = hn[i];
                Us[2][i] = un[i];
            }
        } else if (n == d) {
            //printf("n = %d\n", n);
            for (int i = 0; i <= X - 1; i++) {
                Hs[3][i] = hn[i];
                Us[3][i] = un[i];
            }
        }
        
        // update
        for (int i = 0; i <= X - 1; i++) {
            h[i] = hn[i];
            u[i] = un[i];
        }
        
        //plot
        if ( n == skip * v ) {
            for (j = 0; j <= X - 2; j++){
                data_plot << hn[j] <<"," ;
            }
            data_plot << hn[X - 1];
            data_plot <<"\n";
            v++;
        }
        
        
        umax = *std::max_element(un, un + X - 1);
        umaxx << umax << "\n";
    } //main iteartion stops
    cout<< "\n" <<endl;
    data_plot.close();
    umaxx.close();
    
    // Print Results
    for (int j = 0; j <= X - 1; j++){
        printf(" hi = %f, hn = %f, hh = %f, H0 = %f, H1 = %f\n", hi[j], hn[j], hh[j], Hs[0][j], Hs[1][j]);
    } // Print Results
    
    // Export Results
    ofstream data_array;
    data_array.open ("data_array.csv");
    for (int j = 0; j <= X - 1; j++){
        data_array << Hs[0][j] <<","  << Hs[1][j] <<","<< Hs[2][j] <<","<< Hs[3][j] <<"\n";
    } // Export Results
    data_array.close();
    
    
    
    
    
    cout << "Iterasi Berhasil \n";
}
