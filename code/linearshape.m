function ainit = linearshape(R,Bo,Re,ep,L,n,del)
A = -(2*pi/L)^3+(1/R^2-9*Bo*Re/40)*2*pi/L;
B = 4*Bo/R + 3*Bo/ep;
g1 = - 3*Bo/R;
g2 = 2*pi/L*(1/R^2 - (2*pi/L)^2) ;

a = 1/(A^2+B^2)*(A*g1+B*g2);
b = 1/(A^2+B^2)*(B*g1 - A* g2);

z = 0:L/n:L - L/n;
ainit = 1+ a*del*sin(2*pi/L*z)+b*del*cos(2*pi/L*z);