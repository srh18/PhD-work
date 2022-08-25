function [R,Re,Bo,eps,h] = getndparams(R1,Q,L,eqn)
L = L/(2*pi);
h = (3*Q*1e-6./(2*pi*9.81*R1)).^(1/3);
R = 1
eps = h./R1
Bo = 997*L.^2*9.81/0.072./(eqn*(eps -1)+1)
Re = 3*Q./(2*pi*R1*1e-6)
k = 2*pi*R1/L

             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h^3*g/(nu*D2)
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
              a = (5*eps.*(20*eps.*R.^2.*Bo.^2-30*Bo.^2.*R.^3+eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4-300*Bo.^2.*eps.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));
                
                theta = atan(Bo.*k.*R.^2.*(15*k.^2.*R.^2-15-4*eps.*Re.*R.*Bo)./(30*Bo.^2.*R.^3-eps.*(20*Bo.^2.*R.^2+k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi
                A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4-300*eps.*Bo.^2.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)))
                h2 = A*exp(1i*theta)/eps
                kc = kp*h^2/D2
                Ca = 2*C1*km/(kp*C2)
end