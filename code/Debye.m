function gamma = Debye(T,I, a,z) 
if T>200
    T = T-273;
end
A = 0.4883+8.074e-4*T;
B = 0.3241+1.6e-4*T;
gamma = 10^(-A*z^2*sqrt(I)/(1+B*a*sqrt(I)));