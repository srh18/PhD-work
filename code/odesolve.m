Re = 0.01;
R = 1; 
Bo = 10;
Q = 1;
eps = 1; 
zspan = [-pi pi];
y0 = [1; 0 ; -1 ];
func = @(z,y) [y(2);y(3); 3/eps*Bo*Q./y(1).^3-Bo*(1/eps+(y(1)/R+ 3*cos(z)/R-9/40*Re*y(1).^3*y(2)))+(sin(z)-y(2))/R^2-sin(z)];
[t,y] = ode45(func,zspan,y0);
plot(t,y(:,1))

