#kestabilan metode prediktor-korektor adam bashforth moulton untuk persamaan swe
#numerical stability for predictor-corrector adam for swe
##%%% Pred
##%%% 529 a b x^4 - 736 a b x^3 + 486 a b x^2 - 160 a b x + 25 a b - 529 a c x^4 + 736 a c x^3 - 486 a c x^2 + 160 a c x - 25 a c + 23 a x^5 - 39 a x^4 + 21 a x^3 - 5 a x^2 + 23 b x^5 - 39 b x^4 + 21 b x^3 - 5 b x^2 + x^6 - 2 x^5 + x^4
##%%% + x^6 + 23 a x^5 + 23 b x^5 - 2 x^5 + 529 a b x^4 - 529 a c x^4 - 39 a x^4 + x^4 - 39 b x^4 - 736 a b x^3 + 21 a x^3 + 736 a c x^3 + 21 b x^3 + 486 a b x^2  - 5 a x^2    - 5 b x^2 - 486 a c x^2 - 160 a b x + 160 a c x + 25 a b     - 25 a c
##%%%  x^6 + ( 23*a* + 23*b  - 2 ) * x^5 + ( 529*a*b - 529*a*c - 39*a + 1 - 39*b ) * x^4 + ( -736*a*b + 21*a + 736*a*c + 21*b ) *  x^3 + ( 486*a*b  - 5*a - 5*b - 486*a*c ) * x^2 + ( -160*a*b + 160*a*c ) * x + ( 25*a*b - 25*a*c )
##
##
##%%% Corr
##%%% 25 a b x^6 + 80 a b x^5 + 4 a b x^4 - 96 a b x^3 + 36 a b x^2 - 25 a c x^6 - 80 a c x^5 - 4 a c x^4 + 96 a c x^3 - 36 a c x^2 + 5 a x^6 + 3 a x^5 - 14 a x^4 + 6 a x^3 + 5 b x^6 + 3 b x^5 - 14 b x^4 + 6 b x^3 + x^6 - 2 x^5 + x^4
##%%% + x^6 + 25 a b x^6 - 25 a c x^6 + 5 a x^6 + 5 b x^6 + 80 a b x^5 - 80 a c x^5  + 3 a x^5 + 3 b x^5 - 2 x^5 + 4 a b x^4 - 4 a c x^4 - 14 a x^4 - 14 b x^4   + x^4 - 96 a b x^3 + 96 a c x^3 + 6 a x^3 + 6 b x^3 + 36 a b x^2  - 36 a c x^2
##%%% ( 1 + 25*a*b - 25*a*c + 5*a + 5*b ) * x^6 + ( 80*a*b - 80*a*c  + 3*a + 3*b - 2) * x^5 + ( 4*a*b - 4*a*c - 14*a - 14*b + 1 ) * x^4 + ( -96*a*b + 96*a*c + 6*a + 6*b ) * x^3 + ( 36*a*b - 36*a*c) * x^2
##
import numpy as np
import matplotlib.pyplot as plt

ls = ['-', '--', '-.', ':', ',']

##% Cr = [0.05 0.1 0.5 1.0 1.5];
Cr = np.linspace(0.01, 1.5, 4);
theta = np.linspace(0, 2*np.pi, 100);
F = 3;
s = 0.5;
dx = 0.01/np.sqrt(F);
dtn = Cr*dx;
R   = 500;
mu = 0.001

r1 = np.zeros((len(Cr),len(theta)))
rr1 = np.zeros((len(Cr),len(theta)))


for n, dt in enumerate( dtn ):
    a = 1j*dt*np.sin(theta)/(12*dx);
    b = dt/(12) * ( 1j*np.sin(theta)/(F*dx) - 1/(s) +  mu*(4*np.sin(theta/2)**2)/( dx**2 ) );
    c = dt/(12) * ( 1j*np.sin(theta)/(dx) +  2/(s) +  (4*np.sin(theta/2)**2)/(R * dx**2 ) );

##    %%% Pred
    for j in range(len(theta)):
        x6 = 1;
        x5 = 23*a[j] + 23*c[j] - 2;
        x4 = -529*a[j]*b[j] + 529*a[j]*c[j] - 39*a[j] - 39*c[j] + 1;
        x3 = 736*a[j]*b[j] - 736*a[j]*c[j] + 21*a[j] + 21*c[j];
        x2 = -486*a[j]*b[j] + 486*a[j]*c[j] - 5*a[j] - 5*c[j];
        x1 = 160*a[j]*b[j] - 160*a[j]*c[j];
        x0 = -25*a[j]*b[j] + 25*a[j]*c[j];
        p1 = [x6, x5, x4, x3, x2, x1, x0];
        r1[n][j] = max(abs(np.roots(p1)));


##    %%% Corr
    for j in range(len(theta)):
        x6 = 1 - 25*a[j]*b[j] + 25*a[j]*c[j] + 5*a[j] + 5*c[j];
        x5 = -80*a[j]*b[j] + 80*a[j]*c[j] + 3*a[j] + 3*c[j] - 2;
        x4 = -54*a[j]*b[j] + 54*a[j]*c[j] - 9*a[j] - 9*c[j] + 1;
        x3 = 16*a[j]*b[j] - 16*a[j]*c[j] + a[j] + c[j];
        x2 = -a[j]*b[j] + a[j]*c[j];
        x1 = 0;
        x0 = 0;
        p2 = [x6, x5, x4, x3, x2, x1, x0];
        rr1[n][j] = max(abs(np.roots(p2)));


##
plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')
leg = []
g = plt.figure(1)
for n in range(len(dtn)):
    plt.plot(theta, r1[n], ls[n]);
    leg.append('Cr = ' + str( round(Cr[n],2) ))
##plt.axis([0, max(theta), 0.7, 3.6])
plt.ylabel(r'$|\lambda_{max}|$')
plt.xlabel(r'$\theta$')
plt.legend(leg, loc = 2)
plt.savefig('pred.eps')
plt.show(block=False)

f = plt.figure(2)
for n in range(len(dtn)):
    plt.plot(theta, rr1[n], ls[n]);
#plt.axis([0, max(theta), 0.7, 3.6])
plt.ylabel(r'$|\lambda_{max}|$')
plt.xlabel(r'$\theta$')
plt.legend(leg, loc = 1)
plt.savefig('corr.eps')
plt.show(block=False)
