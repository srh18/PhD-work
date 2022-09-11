function smallA2plots

fig = gcf;
clf

hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 500];
xlabel('$\frac{L}{\pi}$','FontSize',20)
ylabel('$B$')
title('Amplitude of fluid disturbance for small wall disturbance')
legend
ax = gca;
ax.YScale ='log';
fig2 = figure;
hold on 
fig2.WindowStyle = 'normal';
fig2.Position = [100 100 600 500];
xlabel('Re')
ylabel('$\frac{\theta}{\pi}$','FontSize',20)
title('Phase shift of fluid disturbance from small wall disturbance')
legend('Location','east')
L = pi/4:pi/100:8*pi;

%L(176) = 2pi
figure(fig);
ax = gca;
for R = [0.5 1 2]
    Bo = 1;
    eps = 1e-1;
    Re = 1;
[A,~,B,~,Be,~] = smallA2(L,R,Bo,1,eps);

k = 2*pi./L;



plot(L/pi,abs(B),'DisplayName',sprintf('$R =%g$',R) )
 ax.ColorOrderIndex = ax.ColorOrderIndex -1;
 plot(L/pi,Be,'--','HandleVisibility','off')
end
xlim([1/4 8])
figure(fig2);
Re = 0.1:0.01:20;
for eps = [1e-1, 1e-2]
[~,~,~,phi,~,~] = smallA2(2*pi,1,1,Re,eps);
plot(Re,phi/pi,'DisplayName',sprintf('$\\epsilon =%g$',eps))

end
saveas(fig,'../plots/smallA2A','epsc')
saveas(fig2,'../plots/smallA2theta','epsc')
end