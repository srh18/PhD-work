% %where plots are generated
% 
%waveperiod plots


fig = figure(22)
fig.WindowStyle = 'normal';

fig.Position = [100 100 900 400]
% value = loadbigpR1(1,2.25*pi);
% value.plot_period(300)
% 2/5.25
% ylim([0.6,1.6])
% ylim([0.5,1.6])
% fig = figure(22)
% fig.WindowStyle = 'normal';
% 
% fig.Position = [100 100 900 400]
% saveas(gcf,'../plots/psbigp/waveperioddel1L2-25pi','epsc')
% saveas(figure(22),'../plots/psbigp/waveperioddel1L2-25pi','epsc')
% value = loadbigpR1(1,4*pi);
% clf
% value.plot_period(300)
% ylim([0.5,1.7])
% ylim([0.5,1.75])
% saveas(figure(22),'../plots/psbigp/waveperioddel1L4pi','epsc')
% value = loadbigpR1(1,5.25*pi);
% clf
% value.plot_period(300)
% value.plot_period(300,2)
% clf
% value.plot_period(300,2)
% clf
% value.plot_period(300)
% saveas(figure(22),'../plots/psbigp/waveperioddel1L5-25pi','epsc')

%speed
fig = figure(100);
clf,hold on 
fig.WindowStyle = 'normal';
fig.Position = [100 100 700 600];
value = loadbigpR1(1,5.25*pi);
value.surfdata(0.9395,2.1)
ylim([300 330])
[m,l] = max(fig.CurrentAxes.Children.CData(:,1:128),[],2);
[mm,ll] = max(fig.CurrentAxes.Children.CData(:,129:end),[],2);
ll = ll+128;
[m2,l2] = min(fig.CurrentAxes.Children.CData(:,65:192),[],2);
[mm2,ll2] = min(fig.CurrentAxes.Children.CData(:,193:end),[],2);
l2 = l2+64;
ll2 = ll2+192;
fig.CurrentAxes.Children(1).HandleVisibility = 'off';
plot(smooth(smooth(value.z(l),7),7),value.t,'k:','LineWidth',2,'DisplayName','Peak locations')
plot(smooth(smooth(value.z(l2),7),7),value.t,'g:','LineWidth',2,'DisplayName','Trough locations')

plot(smooth(smooth(value.z(ll),7),7),value.t,'k:','LineWidth',2,'HandleVisibility','off')
plot(smooth(smooth(value.z(ll2),7),7),value.t,'g:','LineWidth',2,'HandleVisibility','off')

plot(value.z,304.9+0*value.z,'r--','LineWidth',2,'DisplayName','Peak 1 at wall peak')
plot(value.z,315+0*value.z,'r-.','LineWidth',2,'DisplayName','Peak 1 at wall trough')
plot(value.z,313.6+0*value.z,'m--','LineWidth',2,'DisplayName','Peak 2 at wall peak')
plot(value.z,323.8+0*value.z,'m-.','LineWidth',2,'DisplayName','Peak 2 at wall trough')
plot(value.z,304.9+17.6+0*value.z,'r--','LineWidth',2,'HandleVisibility','off')
plot(value.z,323.8-17.6+0*value.z,'m-.','LineWidth',2,'HandleVisibility','off')

legend('NumColumns',2,'Location','northwest')
saveas(fig,'../plots/psbigp/surfsecdel1L5-25pi','epsc')

fig = figure(100);

fig.WindowStyle ='normal';
fig.Position = [100 100 700 600];
clf,hold on 
value = loadbigpR1(1,4*pi);
value.surfdata(1.1116,3.8)
plot(value.z,300+0*value.z,'k','LineWidth',2)
plot(value.z,320+0*value.z,'k','LineWidth',2)
saveas(fig,'../plots/psbigp/surfdel1L4pi','epsc')

fig = figure(100);

fig.WindowStyle ='normal';
fig.Position = [100 100 700 600];
clf,hold on 
value = loadbigpR1(1,4*pi);
value.surfdata(1.1116,3.8)
[m,l] = max(fig.CurrentAxes.Children.CData,[],2);
[m2,l2] = min(fig.CurrentAxes.Children.CData(:,129:end-32),[],2);
l2 = l2+128;
fig.CurrentAxes.Children(1).HandleVisibility = 'off';
plot(smooth(smooth(value.z(l),7),7),value.t,'k:','LineWidth',2,'DisplayName','Peak location')
plot(smooth(smooth(value.z(l2),7),7),value.t,'g:','LineWidth',2,'DisplayName','Trough location')
ylim([300 320])
plot(value.z,302.5+0*value.z,'r--','LineWidth',2,'DisplayName','Peak at wall peak')
plot(value.z,309.1+0*value.z,'r-.','LineWidth',2,'DisplayName','Peak at wall trough')
plot(value.z,302.5+11.3+0*value.z,'r--','LineWidth',2,'HandleVisibility','off')

legend('NumColumns',2,'Location','northwest')
saveas(fig,'../plots/psbigp/surfsecdel1L4pi','epsc')


fig = figure(100);

fig.WindowStyle ='normal';
fig.Position = [100 100 700 600];
clf,hold on 
value = loadbigpR1(1,2.25*pi);
value.surfdata(0.9359,-0.4)
[m,l] = max(fig.CurrentAxes.Children.CData,[],2);
[m2,l2] = min(fig.CurrentAxes.Children.CData,[],2);
fig.CurrentAxes.Children(1).HandleVisibility = 'off';
plot(smooth(smooth(value.z(l),7),7),value.t,'k:','LineWidth',2,'DisplayName','Peak location')
plot(smooth(smooth(value.z(l2),7),7),value.t,'g:','LineWidth',2,'DisplayName','Trough location')
ylim([300 315])
plot(value.z,304.9+0*value.z,'r--','LineWidth',2,'DisplayName','Peak at wall peak')
plot(value.z,301.6+0*value.z,'r-.','LineWidth',2,'DisplayName','Peak at wall trough')
plot(value.z,301.6+7.6+0*value.z,'r-.','LineWidth',2,'HandleVisibility','off')
plot(value.z,304.9+7.6+0*value.z,'r--','LineWidth',2,'HandleVisibility','off')

legend('NumColumns',2,'Location','northwest')
saveas(fig,'../plots/psbigp/surfsecdel1L2-25pi','epsc')


h2 plots
value = loadbigpR1(0,4*pi);
value.ploth2
hold on 
value = loadbigpR1(0.1,4*pi);
value.ploth2
value = loadbigpR1(1,4*pi);
value.ploth2
title('$||h||_2$ for wall $\eta = \delta\cos(0.5z)$')
ylabel('$||h||_2$')
xlabel('$t$')
legend('$\delta = 0$','$\delta = 0.1$','$\delta = 1$','Location','northwest')
saveas(fig,'../plots/psbigp/L4pih2','epsc')

comparing h2
value = loadbigpR1(1,4*pi);
value.ploth2
xlim([300, 500])
hold on 
value = loadbigpR1(0.5,5.75*pi);
plot(value.t, (value.h2norm-mean(value.h2norm))*6 + 1.06)
value = loadbigpR1(0.5,10*pi);
plot(value.t, (value.h2norm-mean(value.h2norm))*4 + 1.1)
value = loadbigpR1(3,4*pi);
plot(value.t, (value.h2norm-mean(value.h2norm))/8 + 1.17)
ylabel('scaled $||h||_2$')
title('Examples of shapes of $||h||_2$ after initial growth period')
legend('periodic','quasi-periodic','chaotic','drop forming')
xlabel('$t$')


value = loadbigpR1(3,4*pi,2);
plot(value.z,value.mean_h-value.h0)
title({'Difference between time averaged', 'and the steady state fluid thickness'})

xlim([0, 4*pi])
xlabel('$z$')
ylabel('$\bar{h} - h_0$')
saveas(gcf,'../plots/psbigp/meandiffdrop','epsc')
%time periodic L9

load('L9cases','v1','v2','v3','v3b','v4')
clf
v1 = v1.get_h;
v2 = v2.get_h;
v3 = v3.get_h;
v3b = v3b.get_h;
v4 = v4.get_h;

plot(v1.z,v1.mean_h-v1.h0)
hold on
plot(v1.z,v4.mean_h-v1.h0)
plot(v1.z,v3.mean_h-v1.h0)
title('Difference between time averaged and the steady state fluid thickness')
xlabel('$z$')
ylabel('$\bar{h} - h_0$')
legend('(a)','(b)', '(e)')
xlim([0, 9*pi])
saveas(gcf,'../plots/psbigp/meandiff','epsc')

clf
v1.plot_fourier_modes
saveas(gcf,'../plots/psbigp/L91f','epsc')
clf
v2.plot_fourier_modes
saveas(gcf,'../plots/psbigp/L92farf','epsc')
clf
v4.plot_fourier_modes
saveas(gcf,'../plots/psbigp/L92closef','epsc')
clf
v3b.plot_fourier_modes
saveas(gcf,'../plots/psbigp/L93straighterf','epsc')
clf
v3.T = 2000;
v3.plot_fourier_modes
saveas(gcf,'../plots/psbigp/L93curvef','epsc')




[npks,time_periodic,c,T] = v2.get_peak_info(300)
v2.surfdata(c)
v2.surfdata(c,-6)
figure
saveas(gcf,'../plots/psbigp/L92close','epsc')
[npks,time_periodic,c,T] = v4.get_peak_info(300)


v4.surfdata(c)
ylim([1900 2000])
saveas(gcf,'../plots/psbigp/L92far','epsc')
figure
[npks,time_periodic,c,T] = v3.get_peak_info(1600)

v3.surfdata(9*pi/T,4)
ylim([1900 2000])
ylim([1500 2000])
[npks,time_periodic,c,T] = v3b.get_peak_info(1600)


v3b.surfdata(9*pi/T)
saveas(gcf,'../plots/psbigp/L93straighter','epsc')
clf
[npks,time_periodic,c,T] = v3.get_peak_info(1600)
v3.surfdata(9*pi/T,4)
ylim([1500 2000])
figure
v3.ploth2
hold on
v3b.ploth2
saveas(gcf,'../plots/psbigp/L92close','epsc')
saveas(figure(23),'../plots/psbigp/L92far','epsc')
saveas(figure(24),'../plots/psbigp/L93swerve','epsc')
saveas(figure(25),'../plots/psbigp/L93straighter','epsc')
saveas(figure(26),'../plots/psbigp/L93straighter','epsc')
[npks,time_periodic,c,T] = v1.get_peak_info(200)
v1.surfdata(c)
v1.surfdata(c,-10)
ylim([300,400])
saveas(gcf,'../plots/psbigp/L91','epsc')