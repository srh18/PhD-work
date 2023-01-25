function steadystepplots
value = shape_solveh;
value.wall_shape = 1;
value.equation = 1;
fig = gcf;
clf

hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 300];
value.del = 0.1;
value.del2 = 0.01;

grid

for L = [4*pi ,6*pi]
    value.L = L;
    value = value.get_h;
    plot(value.nz,value.h0-1,'DisplayName',sprintf('$L = %g\\pi$',L/pi))
end
xticks(0:0.25:1)
legend
xlabel('$\frac{z}{L}$')
ylabel('$h-1$')
title('Fluid thickness over a step')
ylim([-0.12,0.12])
plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
saveas(gcf,'../plots/steadystep','epsc')

clf

hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 300];
value.del = 0.1;
value.del2 = 0.01;

grid
for L = [2*pi ,6*pi]
    value.L = L;
    value = value.get_h;
    plot(value.nz,value.h0+value.eta,'DisplayName',sprintf('$L = %g\\pi$',L/pi))
end
legend('Location','west')
xlabel('$\frac{z}{L}$')
ylabel('$S$')
title('Fluid profile over a step')
xticks(0:0.25:1)
plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
saveas(gcf,'../plots/steadystepS','epsc')

fig.Position = [100 100 600 300];
clf

hold on 
value.L = 2*pi;

grid

for l = [ 1/4 3/4 ]
    clf,hold on
    grid
    ax = gca;
    value.l = l;
    value = value.get_h;
    plot(value.nz,value.h0-1,'DisplayName',sprintf('$\\ell = %g\\pi$',l))
    
    plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
    %ax.ColorOrderIndex = ax.ColorOrderIndex-1;

legend
xlabel('$\frac{z}{L}$')
ylabel('$h-1$')
title('Fluid thickness over a step')
ylim([ax.YLim(1)-0.02,ax.YLim(2)+0.02])
xticks(0:0.25:1)
saveas(gcf,sprintf('..//plots//steadystepl%i',4*l),'epsc')
end
end
