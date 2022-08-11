function [R,Re,Bo,eps] = getndparams(R1,Q,L,eqn)
L = L/(2*pi);
h = (3*Q*1e-6./(2*pi*9.81*R1)).^(1/3);
R = R1./L;
eps = h./L;
Bo = 997*L.^2*9.81/0.072./(eqn*(eps -1)+1);
Re = 3*Q./(2*pi*R1*1e-6);
end