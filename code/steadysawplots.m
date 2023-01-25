function steadysawplots
value = shape_solveh;
value.equation =1;
value.wall_shape = 2;
fig = gcf;
clf
value.l = 0.8;
hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 300];
value.del = 0.1;
value.del2 = 0.01;

grid
for L = [2*pi ,6*pi]
    value.L = L;
    value = value.get_h;
    plot(value.nz,value.h0-1,'DisplayName',sprintf('$L = %g\\pi$',L/pi))
end
legend
xlabel('$\frac{z}{L}$')
ylabel('$h-1$')
title('Fluid thickness over a sawtooth')
ylim([-0.12,0.12])
plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
%saveas(gcf,'../plots/steadysaw','epsc')

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
title('Fluid profile over a sawtooth')

plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
%saveas(gcf,'../plots/steadysawS','epsc')

fig.Position = [100 100 600 300];
clf

hold on 
value.L = 2*pi;

grid

for s1 = [ 1 -1 ]
  
    value.dir = s1;
    for s2 = [1 -1]
          clf,hold on
    grid
    ax = gca;
        value.dir2 = s2;
        value.l = 0.8;
    value = value.get_h;
    plot(value.nz,value.h0-1,'DisplayName','$h-1$')
    
    plot(value.z/value.L,value.eta,'k--','DisplayName','$\eta$')
    %ax.ColorOrderIndex = ax.ColorOrderIndex-1;

legend
xlabel('$\frac{z}{L}$')
ylabel('$h-1$')
title('Fluid thickness over a sawtooth')
ylim([ax.YLim(1)-0.02,ax.YLim(2)+0.02])

saveas(gcf,sprintf('..//plots//steadysaw%g%g',s1,s2),'epsc')
end
end
