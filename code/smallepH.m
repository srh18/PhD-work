function [A,theta] = smallepH(L,R,Bo)

theta = atan((8*pi^3*R.^2 - 2*pi*L.^2)./(R.*L.^3.*Bo))+pi;
A = -1/(3*R)*sec(theta);
end