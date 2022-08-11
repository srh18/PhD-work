function [A,theta,B,phi,Be,phie] = smallA2(L,R,Bo,Re,eps)
k = 2*pi./L;
a = -(5*eps.*(20*eps.*R.^2.*Bo.^2-30*Bo.^2.*R.^3+eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4-300*Bo.^2.*eps.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));

theta = atan(Bo.*k.*R.^2.*(15*k.^2.*R.^2-15-4*eps.*Re.*R.*Bo)./(30*Bo.^2.*R.^3-eps.*(20*Bo.^2.*R.^2+k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi;
A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4-300*eps.*Bo.^2.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)));
phi = atan(R.*A.*tan(theta)./(R.*A-eps.*sec(theta)));
B = A.*R.*sin(theta)./((R+eps).*sin(phi));
phi = phi +(1-sign(B))/2*pi;
B = abs(B);
theta = theta+(1-sign(theta))*pi;
phi = phi+(1-sign(phi))*pi;
phie = atan((k-R.^2.*k.^3)./(Bo.*R));
Be = eps./(3*R).*sec(phie);
phie = phie+(1-sign(B))/2*pi;
Be = abs(Be);
    
end
