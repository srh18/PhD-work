function flat_growth
clf,hold on
L = 2*pi:pi/100:10*pi;
k = 2*pi./L;
for kn = 1:4
y = ((kn*k).^2-(kn*k).^4)/3;
plot(L/pi,y,'DisplayName',sprintf('k = %g',kn))
end
legend('Location','North','NumColumns',4)
ylim([0, 0.1]);
xlabel('$\frac{L}{\pi}$')
ylabel('Growth rate')
title('Growth rate of unstable mode $k$')
end
function other_growth
value = loadbigpR1(1,4*pi,1);

for delta = [0.1 0.5 1 2 3]
omega2 = [];
omega3 = [];
omega4 = [];
L = 2*pi:pi/16:10*pi;
omega1 = [];
value.del = delta;
for l = L
value.L = l;
value = value.get_h;
if value.eflag<0
    omega1 = [omega1,0];
    omega2 = [omega2,0];
    omega3 = [omega3,0];
    omega4 = [omega4,0];
else
    
    
[~,d] = value.Floquet(value.h0);
if length(d)>7
omega4 = [omega4 -real(d(end-7))];
else
omega4 = [omega4,0];
end
if length(d)>5
omega3 = [omega3 -real(d(end-5))];
else
omega3 = [omega3,0];
end
if length(d)>3
omega2 = [omega2 -real(d(end-3))];
else
omega2 = [omega2,0];
end
if length(d)>1
omega1 = [omega1 -real(d(end-1))];
else
omega1 = [omega1,0];
end
end
end
clf,hold on
plot(L/pi,omega1,'DisplayName','$k=1$')
plot(L/pi,omega2,'DisplayName','$k=2$')
plot(L/pi,omega3,'DisplayName','$k=3$')
plot(L/pi,omega4,'DisplayName','$k=4$')
str = replace(sprintf('%g',delta),'.','-');
legend('Location','North','NumColumns',4)
ylim([0, 0.1]);
xlabel('$\frac{L}{\pi}$')
ylabel('Growth rate')
title(sprintf('Growth rate of unstable mode $k$ for $\\delta =%g$',delta))
saveas(gcf,sprintf('../plots/psbigp/unstablemodegrowthd%s',str),'epsc')
end