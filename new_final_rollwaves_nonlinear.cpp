//Persamaan Saint-Venant : Simulasi Roll Waves dengan Angin

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include <algorithm>
#include <complex>
#include <time.h>
#include<tuple>
#include <sstream>
// # define M_PI           3.14159265358979323846  /* pi */
using namespace std;


float fun1(float a, float h0, float k) {
    return(h0 * (1 + 0.005 * sin(k * a)));
}

float fun2(float a, float h0, float k, float co, float F) {
    return( sqrt(F * co * h0 * (1 + 0.005 * sin(k * a))) );
}

float fun3(float a, float h0, float u0, float rp, float tp, float k, float co, float F) {
    return( u0 +  0.005*rp*sin(k*a + tp) );
}

std::tuple<double, double> HU(double h[], double u[], double Cf, double co, double si, double nu, double mu, double dx){
	double a = -h[1] * (u[2] - u[0]) / (2 * dx) - u[1] * (h[2] - h[0]) / (2 * dx);
	double b = -u[1] * (u[2] - u[0]) / (2 * dx) - co * (h[2] - h[0]) / (2 * dx) + (si - Cf * (u[1] * u[1]) / h[1]) + nu * (u[2] - 2 * u[1] + u[0]) / (dx * dx) + mu * (h[2] - 2 * h[1] + h[0]) / (dx * dx);
	return make_tuple(a,b); //<h3, u3>
}

string string_round(float F, int n)//round up a float type and show one decimal place
{
		string Fr, FF = to_string(F);
		Fr =  FF.at(0);
		for (int i=1; i <= n; i++)
		{
			 Fr = Fr + FF.at(i);
		}
		return Fr;
}

void SWE(double F, double U, double nu, double k) {
//int main() {
    // double L = 1, dx = 0.000833333333333, dt = 0.000083333;
    // double L = 1, dx = 0.000833333333333, dt = 0.00046667; // 0.56dx;
    double L = 1, dx = 0.0005, dt = dx*0.1;
    // double L = 1, dx = 0.000833333333333, dt = dx*1.2;
    // double L = 1, dx = 0.000833333333333, dt = dx*0.78;
    long int X = round(L/dx);
    long int T = round(1/dt);
    double x[X-1], h[X - 1],  u[X - 1], hi[X - 1], ui[X - 1], hh[X - 1], uu[X - 1], hn[X - 1], un[X - 1], H1[X - 1], H2[X - 1], H3[X - 1], HH[X - 1], U1[X - 1], U2[X - 1], U3[X - 1], UU[X - 1], Hs[4][X - 1], Us[4][X - 1], umax;
    double ht[23][X-1],  ta, tb, tc, td;
    int n,j, vv=1;
    int pembagi = 100, maks = T/100;
    int v = 1, skip = round( T / 2000 );

    if (skip == 0){
        skip = 1;
    }

    // space & time discretization
    for (j = 0; j <= X ; j++) {
        x[j + 1] = x[j] + dx;
    }

    /////////////////////////////////////////////////// w calculation starts here

    // Parameters
    double Cf = 0.006, g = 9.8, q0 = 0.001;//k = 4 * M_PI,  nu = 0.001, F = 9.0;

    /// calculate several paramaters
    double theta = atan(F * Cf), co = g * cos(theta), si = g * sin(theta);

    double h0 = pow((q0*q0/(co * F)), 0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333);
    double u0 = q0/h0;

    /// wind-related value
    // double U = 15; //wind velocity in m/s
    double r = 1000,  ra = 1.23, s = 0.25;

    /// calculate nondimensional parameter
    double R = u0/(nu*L), sigma = h0/(L*Cf), mu = s*(ra*h0/(r*L))*abs(U)*U/(u0*u0);


    double a = pow((1/sigma + k*k/(2*R)),2) - k*k/F;
    double b = -k/sigma + mu*k*k*k;
    double wrp = (1/sigma + k*k/(2*R)) + sqrt(0.5*(a + sqrt(a*a + b*b)));
    double wrm = (1/sigma + k*k/(2*R)) - sqrt(0.5*(a + sqrt(a*a + b*b)));
    double wip = k + sqrt(0.5*( -a + sqrt(a*a + b*b) ));
    double wim = k - sqrt(0.5*( -a + sqrt(a*a + b*b) ));

    double wa, amp;
    if (wip < wim){
        wa = wip/k;
        amp = wrp;
    }
    else {
        wa = wim/k;
        amp = wrm;
    }

    /// critical Froude number for stability of uniform flow
    double aa = (2 + sigma*k*k/R)/(1-mu*sigma*k*k), Fc = pow( aa  , 2);


    //// amplitudo and phase
    double rp = sqrt((amp - u0)*(amp - u0) + wa*wa);
    double tp = atan( wa / (amp - u0) );
    printf("############################################ Physical Parameter############################## \n");
	printf("Cf = %0.3f \n", Cf);
    printf("Fluid depth h0 = %f \n", h0);
    printf("Fluida velocity u0 = %f \n", u0);
    printf("Fluid kinematic Viscosity nu = %0.3f \n", nu);
    printf("Wind velocity U = %0.2f \n", U);
    printf("Wave number k = %0.fpi \n", round(k/(M_PI)));
    printf("############################################ Nondimensional Parameter############################## \n");
    printf("Froude Number F = %0.3f \n", F);
    printf("Wind Parameter mu = %f \n", mu);
    printf("Reynold Number R = %f \n", R);
    printf("Parameter sigma = %f \n", sigma);
    printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Stability criterion F <= %f ////////////////////////////////\n", Fc);
    printf("############################################ Monocromatic Wave ############################## \n");
    printf("Initial fluid velocity amplitude rp = %f \n", rp);
    printf("Initial fluid velocity phase tp = %f \n", tp);
    printf("wave speed1 = %f \n", wip/k);
    printf("wave speed2 = %f \n", wim/k);
    printf("amplification_factor_1 = %f \n", wrp);
    printf("amplification_factor_2 = %f \n", wrm);
    printf("Please wait, I'm calculating for you ;) \n");



    /////////////////////////////////////////////////// w calculation ends here

    /// initial condition
    x[0] = dt;
    for (j = 0; j <= X - 1; j++) {
        h[j] = fun1(x[j], h0,  k );
        //u[j] = fun2(x[j], h0,  k,  co,  F);
        u[j] = fun3(x[j], h0,  u0, rp, tp, k,  co,  F);
    }
    // initial wave profile
    for (j = 0; j <= X - 1; j++) {
        hi[j] = fun1(x[j], h0,  k );
        ui[j] = fun3(x[j], h0,  u0, rp, tp, k,  co,  F);
    }

    for (j = 0; j <= X - 1; j++) {
        H1[j] = 0;
        U1[j] = 0;
        H2[j] = 0;
        U2[j] = 0;
    }
    // Export Plot
		ofstream data_plot, u_plot;
		string Fr = string_round(F, 3); // round decimal to 3
		string Ustr = string_round(U, 3); // round decimal to 3
		string nus = string_round(nu, 5); // round decimal to 3
    data_plot.open ("data_plot_F="+Fr+"_U="+Ustr+"_nu="+nus+".csv");
    u_plot.open ("u_plot_F="+Fr+"_U="+Ustr+"_nu="+nus+".csv");

    for (int j = 0; j <= X - 2; j++){
        data_plot << x[j] <<"," ;
    }
    data_plot << x[X - 1] <<"\n";

    for (int j = 0; j <= X - 2; j++){
        data_plot << hi[j] <<"," ;
    }
    data_plot << hi[X - 1] <<"\n";

    for (int j = 0; j <= X - 2; j++){
        u_plot << ui[j] <<"," ;
    }
    u_plot << ui[X - 1] <<"\n";

    /// main iteration starts
    for (n = 0; n <= T; n++) {

        // Predictor Operator
        for (j = 1; j <= X - 2 ; j++){
						double ha[] = {h[j-1], h[j], h[j+1]}, ua[] = {u[j-1], u[j], u[j+1]};
						auto v = HU(ha, ua, Cf, co, si, nu, mu, dx);
            H3[j] = std::get<0>(v);
						U3[j] = std::get<1>(v);
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
						double ha[] = {hh[j-1], hh[j], hh[j+1]}, ua[] = {uu[j-1], uu[j], uu[j+1]};
						auto v = HU(ha, ua, Cf, co, si, nu, mu, dx);
						HH[j] = std::get<0>(v);
						UU[j] = std::get<1>(v);
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


        // update
        for (int i = 0; i <= X - 1; i++) {
            h[i] = hn[i];
            u[i] = un[i];
        }

        //plot
        if ( n == skip * v ) {
        	// fluid depth data save
            for (j = 0; j <= X - 2; j++){
                data_plot << hn[j] <<"," ;
            }
            data_plot << hn[X - 1];
            data_plot <<"\n";
            // fluid velocity data save
            for (j = 0; j <= X - 2; j++){
                u_plot << un[j] <<"," ;
            }
            u_plot << un[X - 1];
            u_plot <<"\n";
            v++;
        }


        //umax = *std::max_element(un, un + X - 1);
        //umaxx << umax << "\n";
    } //main iteartion stops
    cout<< "\n" <<endl;


    // Print Results
    for (int j = 0; j <= X - 1; j++){
				if (hn[j] != hn[j]){
					printf(" blow up\n");
				}
    } // Print Results
		data_plot.close();
		u_plot.close();
    cout << "Iterasi Berhasil \n";
}

int main() {
	double F[] = {5.0, 5.5, 4.5};
	double U = 15.0;
	double nu = 0.001;
	double k = 4 * M_PI;
	for (int i = 0; i <= 2; i++){
		SWE(F[i], U, nu, k);
	}

}
