function [A,theta] = smallA(L,R,Bo)
k = 2*pi./L;
theta = atan(3*Bo.*R.^2./(k.^3.*R.^2-k));
A = sqrt((9.*Bo.^2.*R.^4.*(k-k.^3.*R.^2).^2 + (k-k.^3.*R.^2).^4)./(((k - k.^3.*R.^2).^2+9*Bo.^2.*R.^4).^2));
end