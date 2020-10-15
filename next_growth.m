function next_growth(a,b)
z= 0:pi/100:10*pi;
eta0 =  sin(z);
deta0 =  cos(z);
dddeta0 = -cos(z);
eta0(1:400)= 0 ;
eta0(601:1001)= 0;
deta0(1:400)= 0 ;
deta0(601:1001)= 0;
dddeta0(1:400)= 0 ;
dddeta0(601:1001)= 0;
for t = 0:0.01:10
    eta = (-a/3 +(a*z+b))*t+eta0;
    inteta = (-a/3 +(a*z+b))*t^2/2+eta0*t;
    intdeta02 = a*t^3/3+a*deta0*t^2+deta0.^2*t;
    eta2 = 1/315*((z*a+b-a/3).*(420*inteta+105*intdeta02+35*dddeta0)+51*a*t);
    figure(1)
    clf,hold on

    plot(eta0,z)
    plot(eta+0.1*eta2,z)
    xlim([-1,20])
    ylim([0,10*pi])
    pause(0.001)
end
end