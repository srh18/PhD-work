Q = 10:10:1000;
Qa = (10:10:1000)*1e-6/3600;
nu = 0.01*1e-4;
R = [2.5e-3 1e-2 5e-2 1e-1];

clf
for r = R
h = (3*Qa*nu/(2*pi*9.8*r)).^(1/3);
loglog(Qa,h*1e6,'DisplayName',sprintf('$R = %g$mm',r*1000)),hold on
end
legend('northwest')
title('Fluid thickness for increasing flow rates at different radii')
xlabel('Flow rate $Q$ (