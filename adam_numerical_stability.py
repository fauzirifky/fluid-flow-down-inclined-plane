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


##% Cr = [0.05 0.1 0.5 1.0 1.5]; 
Cr = np.linspace(0.001, 1.5, 5); 
theta = np.linspace(0, 2*pi, 50);
F = 3.5; 
s = 0.7;
dx = 0.0025/np.sqrt(F);
dtn = Cr*dx;
R   = 1000;
a = b = c = r1 = rr1 zeros((len(Cr),len(theta)))

for n, dt in enumerate( dtn ):
    a[n] = 1j*dt*sin(theta)./(12*dx);
    b[n] = dt/(12) * ( 1j*sin(theta)/dx + 2/(s) +  (4*sin(theta/2)**2)/(R * dx**2 ) );
    c[n] = dt/(12) * ( 1j*sin(theta)./(dx * F) -  1/(s) );
end

sep = round(length(Cr)/5);

p = 1; q = 1;
for n in len( dtn ):
##    %%% Pred
    for j in len(theta):
        x6 = 1;
        x5 = 23*a[n][j] + 23*b(n,j) - 2;
        x4 = 529*a[n][j]*b[n][j] - 529*a[n][j]*c[n][j] - 39*a[n][j] - 39*b[n][j] + 1;
        x3 = -736*a[n][j]*b[n][j] + 736*a[n][j]*c[n][j] + 21*a[n][j] + 21*b[n][j];
        x2 = 486*a[n][j]*b[n][j] - 486*a[n][j]*c[n][j] - 5*a[n][j] - 5*b[n][j];
        x1 = -160*a[n][j]*b[n][j] + 160*a[n][j]*c[n][j];
        x0 = 25*a[n][j]*b[n][j] - 25*a[n][j]*c[n][j];
        p1 = [x6, x5, x4, x3, x2, x1, x0];
        r1[n][j] = max(abs(roots(p1)));  
    

##    %%% Corr
    for j in len(theta):
        x6 = 1 + 25*a[n][j]*b[n][j] - 25*a(n,j)*c(n,j) + 5*a(n,j) + 5*b(n,j);
        x5 = 80*a(n,j)*b(n,j) - 80*a(n,j)*c(n,j) + 3*a(n,j) + 3*b(n,j) - 2;
        x4 = 4*a(n,j)*b(n,j) - 4*a(n,j)*c(n,j) - 14*a(n,j) - 14*b(n,j) + 1;
        x3 = -96*a(n,j)*b(n,j) + 96*a(n,j)*c(n,j) + 6*a(n,j) + 6*b(n,j);
        x2 = 36*a(n,j)*b(n,j) - 36*a(n,j)*c(n,j);
        x1 = 0;
        x0 = 0;
        p2 = [x6 x5 x4 x3 x2 x1 x0];
        rr1(n,j) = max(abs(roots(p2)));  
    
    
    
[r1_max(n) r1_loc(n)] = max(r1(n,:));
if r1_max(n) <= 1.0001
    q = q+1;
end

[rr1_max(n) rr1_loc(n)] = max(rr1(n,:));
if rr1_max(n) <= 1.0001
    p = p+1;
end

end




w = 1 : sep : length(Cr);

T       = length(theta);
step    = 5;
figure(1)
    plot(theta(1 : 2*step : T), r1(w(1), 1 : 2*step : T),'o');
    hold on
    plot(theta(1 : 4*step : T), r1(q-1, 1 : 4*step : T),'*b');
    plot(theta(1 : step : T), r1(w(3), 1 : step : T), 'k');
    plot(theta(1 : step : T), r1(p, 1 : step : T), '--', 'linewidth',2);
    plot(theta(1 : step : T), r1(w(5), 1 : step : T), '-+', 'linewidth',2);
    ylabel('|\lambda_{max}|')
    xlabel('\theta')
    title('Norma amplifikasi maksimum untuk skema prediktor')
   legend(['Cr =' num2str(Cr(w(1)))],['Cr = ' num2str(Cr(q-1))],['Cr = ' num2str(Cr(w(3)))], ['Cr = ' num2str(Cr(p))], ['Cr = ' num2str(Cr(w(5)))])
hold off
axis([0 max(theta) 0.7 3.6])


 

figure(2)
    polar(theta(1 : 2*step : T), rr1(w(1), 1 : 2*step : T),'o');
    hold on
    polar(theta(1 : 4*step : T), rr1(q-1, 1 : 4*step : T),'*b');
    polar(theta(1 : step : T), rr1(w(3), 1 : step : T), 'k');
    polar(theta(1 : step : T), rr1(p, 1 : step : T), '--');
    polar(theta(1 : step : T), rr1(w(5), 1 : step : T), '-o');
%     ylabel('|\lambda_{max}|')
%     xlabel('\theta')
    title('Norma amplifikasi maksimum untuk skema korektor')
%     legend(['Cr =' num2str(Cr(w(1)))],['Cr = ' num2str(Cr(q))],['Cr = ' num2str(Cr(w(3)))], ['Cr = ' num2str(Cr(p))], ['Cr = ' num2str(Cr(w(5)))])
    hold off
%     rlim('auto')
axis([0 max(theta) 0.85 1.25])
