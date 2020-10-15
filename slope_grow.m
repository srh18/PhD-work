function slope_grow(m,B)
z = 0:pi/100:2*pi;
eta = 0;
deta = 0;
dddeta = 0;
deta1 = 0*z;
a1 = m/3*(z-4/(1+m^2)*deta);
a2 = - 1/9*(2*m^2*z.^2-3*eta-4*m^2/(1+m^2)*z.*deta+2*(3+m^2)/(1+m^2)^2*deta.^2-12*m/(1+m^2)*deta1-3*B*dddeta/(1+m^2)^(3/2));
clf, hold on

plot(a1+0.1*a2,z);
plot(a1,z)
plot(a2,z)
plot(m*z,z)

ax = gca;
ax.YDir = 'reverse';
end