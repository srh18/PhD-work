function epsphiRplots(Bo,Re)
clf
R= 0.1:1/100:10;
eps = logspace(-0.1,-5,100)';
[~,~,B,phi,~,~] = smallA2(2*pi,R,Bo,Re,eps);
contourf(R,eps,phi/pi,0:0.25:2)
ax = gca;
ax.YScale = 'log';
c = colorbar;
caxis([0,2])
c.Ticks =[0:0.25:2];
set(c.Label,'String','$\frac{\varphi}{\pi}$','Interpreter','Latex', 'Rotation', 270, 'VerticalAlignment','bottom','FontSize',22)
hold on 
%contour(L/pi,eps,phi/pi,[0:0.5:2],'k','LineWidth',2)
xlabel('$R$')
ylabel('$\epsilon$')
title('Phase shift of fluid disturbance from wall disturbance')
fig = gcf;
fig.WindowStyle = 'normal';
fig.Position = [100 100 600 500];
saveas(fig,'../plots/epsphiR','epsc')


end
