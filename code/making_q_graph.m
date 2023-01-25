% %load('/Volumes/srh18/home/eq0Re1ep-1QH2Nunstablemat2highres.mat','Q','H2','M','L','del')
% load('/Volumes/srh18/home/QH2Nunstablemat2highres.mat')
% f = figure;
% hold on
% s = pcolor(L/pi,del,Q);
% shading interp
% s.HandleVisibility = 'off';
% xlim([L(1)/pi L(end)/pi])
% Qnew = Q-(Q<1e-4).*Q;
% %contour(L/pi,del,Qnew,[0 1e-4],'k--','LineWidth',1,'HandleVisibility','off')
% view(2)
% xlabel('$\frac{L}{\pi}$')
% ylabel('$q$')
% ylabel('$\delta$')
% title({'Flow rate for the steady states','where the wall $\eta = \delta\cos\left(\frac{2\pi z}{L}\right)$'})
% text(3.15,4.35,'No steady state','FontSize',20,'Color',[1 1 1])
% 
% 
% legend
% colorbar
% plot([1 1 1],[0.1 0.5 1 ],'o','DisplayName','(b)','MarkerSize',14,'MarkerEdgeColor','r','MarkerFaceColor' , 'r')
% plot([4 4 4 4 4],[0.1 0.5 1  2 3],'s','DisplayName','(c)','MarkerSize',14,'MarkerEdgeColor','m','MarkerFaceColor' , 'm')
% plot([1,3/2, 2,4,8],[1 1 1 1 1],'x','DisplayName','(d)','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor' , 'k')
% hold on 
% make_nice_contours(0);
% fig = gca;
% for i = 2:-1:1
%     fig.Children(i).LineWidth = 2;
%     fig.Children(i).LineStyle = '--';
%     fig.Children(i).HandleVisibility = 'off';
% end
% f.WindowStyle = 'normal';
% f.Position = [100 100 900 900];
% saveas(f,'../plots/steadyqcontourhighres','epsc')
% % fig.Children(1).MarkerSize = 12
% % fig.Children(1).MarkerEdgeColor = 'm'
% % fig.Children(1).MarkerFaceColor = 'm'
% % fig.Children(2).MarkerSize = 12
% % fig.Children(2).MarkerEdgeColor = 'r'
% % fig.Children(2).MarkerFaceColor = 'r'
% % fig.Children(2).MarkerEdgeColor = 'k'
% % fig.Children(2).MarkerEdgeColor = 'r'
% % fig.Children(3).MarkerEdgeColor = 'k'
% % fig.Children(3).MarkerSize = 12
% % fig.Children
% % fig.Children(4)
% % fig.Children(4).Color = [1  1 1]
% % 
% % text(3.1,4.25,'No steady state','FontSize',20)
% % fig.Children
% % fig.Children(1).Color = [1  1 1]
% % fig.Children(5).delete
% % fig.Children
% % fig.Children = [fig.Children(4) fig.Children(2) fig.Children(3) fig.Children(1) fig.Children(5) ]
% % fig.Children(5).DisplayName = ''
% % fig.Children(5).DisplayName = '.'
% % fig.Children(5).DisplayName = '\'
% % fig.Children(5).DisplayName = '/'
% % fig.Children(5).HandleVisiblity = 'off'
% % fig.Children(5).HandleVisibility = 'off'
% % fig.Children(4).DisplayName = '(a)'
% % fig.Children(3).DisplayName = '(a)'
% % fig.Children(3).DisplayName = '(b)'
% % fig.Children(2).DisplayName = '(c)'
% % fig.Children(1).DisplayName = '(d)'
g = figure;
clf, hold on
g.WindowStyle = 'normal';
g.Position = [100 100 900/0.55*0.41 300];
value =shape_solveh;
value.equation = 1;
value.del = 1;
value.L = pi;
value.Bo = 1;
value.R = 1;
value.n = 512;
for delta = [0.1,0.5,1]
    value.del = delta;
    value = value.get_h;
    plot(value.z/value.L,value.h0,'DisplayName',sprintf('$\\delta =%g$',delta))
end
legend
xlabel('$\frac{z}{L}$')
ylabel('$h$')
title({'Fluid thickness for different amplitudes','where $\eta = \delta\cos(2z)$'})
grid
xticks(0:0.25:1)
saveas(g,'../plots/steadydelchangeLpi','epsc')
clf, hold on
value.L = 4*pi;
for delta = [0.1,0.5,1,2,3]
    value.del = delta;
    value = value.get_h;
    plot(value.z/value.L,value.h0,'DisplayName',sprintf('$\\delta =%g$',delta))
end
legend
xlabel('$\frac{z}{L}$')
ylabel('$h$')
title({'Fluid thickness for different amplitudes','where $\eta = \delta\cos(\frac{z}{2})$'})
grid
xticks(0:0.25:1)
saveas(g,'../plots/steadydelchangeL4pi','epsc')

clf, hold on
value.del = 1;
plot(value.z/value.L,value.eta,'k','LineWidth',2)
for L = [1,1.5,2,4,8]*pi
    value.L = L;
    value = value.get_h;
    plot(value.z/value.L,value.h0+value.eta)
end
legend('Wall','$L = \pi$','$L = \frac{3\pi}{2}$','$L = 2\pi$','$L = 4\pi$','$L = 8\pi$')
xlabel('$\frac{z}{L}$')
ylabel('$S$')
title({'Fluid profiles for different amplitudes','where $\eta = \delta\cos(\frac{2\pi z}{L})$'})
grid
ylim([-1,3.5])
xticks(0:0.25:1)
saveas(g,'../plots/steadyLchangedel1','epsc')
