function smallAplots

fig = gcf;
clf

hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 500];
xlabel('$\frac{L}{\pi}$')
ylabel('Amplitude')
title('Amplitude of fluid disturbance for small wall disturbance')
legend
fig2 = figure;
hold on 
fig2.WindowStyle = 'normal';
fig2.Position = [100 100 600 500];
xlabel('$\frac{L}{\pi}$')
ylabel('$\frac{\theta}{\pi}$')
title('Phase shift of fluid disturbance from small wall disturbance')
legend
L = pi/4:pi/100:8*pi;
%L(176) = 2pi
for R = [0.5,1,2]
    
[A,theta] = smallA(L,R,1);
n = (2*R - 1/4)*100 + 1;

theta = theta +pi;
theta(n) = nan;
figure(fig);
plot(L/pi,A,'DisplayName',sprintf('$R =%g$',R) )
figure(fig2);
plot(L/pi,theta/pi,'DisplayName',sprintf('$R =%g$',R))

xline(L(n)/pi,'--','Color' ,fig2.Children(2).Children(1).Color,'HandleVisibility','off' )
end
saveas(fig,'../plots/smallAA','epsc')
saveas(fig2,'../plots/smallAtheta','epsc')
end