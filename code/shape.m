function shape(R,A,k,Bo,del)
z = -pi:pi/100:pi;

a = 1 - del*(1/(3*R) + A*cos(z)/R+1/3/Bo*(-A*k/R^2+A*k^3)*sin(z));
plot(a,z) 
end