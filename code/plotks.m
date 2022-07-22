function plotks
pH = 6:0.05:10;
clf
fig = gcf;
ax = gca;
hold on 
ax.YScale = 'log';
fig.WindowStyle = 'normal';
fig.Position = [100 100 900 400];
for T = [12,20]
    
    [r,kp,km] = C02rate(pH,T+273,1,1);
    plot(pH,kp,'DisplayName',sprintf('$k_+, T =%g^{\\circ} $C',T))
    ax.ColorOrderIndex =  ax.ColorOrderIndex -1;
    plot(pH,km,'--','DisplayName',sprintf('$k_-, T =%g^{\\circ} $C',T))
end
legend('Location','northwest')
xlabel('pH')
ylabel('Rate s$^{-1}$')

saveas(fig,'../plots/kpH','epsc')
end