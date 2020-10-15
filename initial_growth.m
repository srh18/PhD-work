function initial_growth(a,b)
z= 0:pi/100:2*pi;
for t = 0:0.01:10
    eta = (a/3 +(a*z+b))*t;
    figure(1)
    clf
    plot(eta,z)
    xlim([0,25])
    ylim([0,2*pi])
    pause(0.1)
end
end