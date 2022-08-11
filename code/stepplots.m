function stepplots
value = shape_solveh;
value.wall_shape = 1;
fig = gcf;
clf

hold on 


fig.WindowStyle = 'normal';
fig.Position = [100 100 600 500];
for xi = [0.1 0.05 0.01]
    
    value.del2 = xi;
    plot(value.z,value.eta,'DisplayName',sprintf('$\\xi = %g$',xi))
end
xlabel('$L$')
ylabel('$\eta$')
title('How the steepness parameter affects the wall shape')
ylim([-1.2,1.2])
xlim([0,value.L])
yline(0,'--','HandleVisibility','off')
legend('Location','south')
saveas(gcf,'../plots/stepxi','epsc')
clf, hold on 
for i1 = [0.25,0.5 0.75]
    value.l = i1;
    plot(value.z,value.eta,'DisplayName',sprintf('$\\ell = %g$',i1))
end
xlabel('$L$')
ylabel('$\eta$')
title('How $\ell$ affects the wall shape')
ylim([-1.6 1.6])
xlim([0,value.L])
yline(0,'--','HandleVisibility','off')
legend('Location','south')
saveas(gcf,'../plots/stepell','epsc')

clf, hold on 
value.wall_shape = 2;
value.dir = 1;
value.dir2 = 1;
value.del = 1;
value.l = 0.5;
for xi = [0.1 0.05 0.01]
    
    value.del2 = xi;
    plot(value.z,value.eta,'DisplayName',sprintf('$\\xi = %g$',xi))
end
xlabel('$L$')
ylabel('$\eta$')
title('How the steepness parameter affects the wall shape')
ylim([-0.6,1.6])
xlim([0,value.L])
yline(0,'--','HandleVisibility','off')
legend('Location','northwest')
saveas(gcf,'../plots/sawxi','epsc')
clf,hold on 
for i1 = [1 -1]
    value.dir = i1;
    for i2 = [1,-1]
        if i2 == 1
            str = '-';
        else
            str = '--';
        end
        value.dir2 = i2;
        value.l = 0.5;
    plot(value.z,value.eta,str,'DisplayName',sprintf('$(\\sigma_1,\\sigma_2) = (%g,%g)$',i1,i2))
    end
end
xlabel('$L$')
ylabel('$\eta$')
title('Different direction Sawtooths')
ylim([-1.5 1.85])
xlim([0,value.L])
yline(0,'-.','HandleVisibility','off')
legend('Location','north')
saveas(gcf,'../plots/sawsigma','epsc')

end