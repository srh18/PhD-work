function epsphiplots(R,Bo,Re)
clf
L = 0.1*pi:pi/100:8*pi;
eps = logspace(-0.2,-5,100)';
[~,theta,B,phi,~,~] = smallA2(L,R,Bo,Re,eps);
contourf(L/pi,eps,phi/pi,0:0.25:2)
ax = gca;
ax.YScale = 'log';
c = colorbar;
caxis([0,2])
c.Ticks =[0:0.25:2];
set(c.Label,'String','$\frac{\varphi}{\pi}$','Interpreter','Latex', 'Rotation', 270, 'VerticalAlignment','bottom','FontSize',22)
hold on 
%contour(L/pi,eps,phi/pi,[0:0.5:2],'k','LineWidth',2)
xlabel('$\frac{L}{\pi}$','FontSize',20)
ylabel('$\epsilon$')
title('Phase shift of fluid disturbance from wall disturbance')
fig = gcf;
fig.WindowStyle = 'normal';
fig.Position = [100 100 600 500];
saveas(fig,'../plots/epsphi','epsc')



end
