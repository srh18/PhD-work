function plotsco2rate
pCO2 = logspace(-3.5,-1);
clf
fig = gcf;
ax = gca;
hold on 

ax.XScale = 'log';
fig.WindowStyle = 'normal';
fig.Position = [100 100 900 400];
for Ca = [2,5]
    
    [r,~,~] = C02rate(8,12+273,pCO2,Ca*1e-3);
    plot(pCO2,r*1e3,'DisplayName',sprintf('$k_+$, Ca $=%g$ mol/m$^3$ ',Ca))
    
end
ylim([-1e-2 1e-2])
xlim([10^-3.5,4e-2])
yline(0,'--','HandleVisibility', 'off')
legend('Location','northeast')
title('pH = 8, T = $12^{\circ}$C')
xlabel('pCO$_2$')
ylabel('Reaction rate of CO$_2$ mol m$^{-3}$ s$^{-1}$')
saveas(fig,'../plots/rCO2','epsc')
end