function compare_animate( speed,fluidsmallR,fluidnormR,fluidbigR)
l = length(fluidbigR.a);

% xmax = max(max(data));
% xmin = min(min(data));
for i = 1:speed:l
    clf, hold on
    
    flip_y
    subplot(1,3,1)
   
    plot(fluidsmallR.eta,fluidsmallR.nz,'DisplayName','wall')
    plot(fluidsmallR.a(i,:),fluidsmallR.nz,'DisplayName',"fluid")
    title('radius $O(1)$')
        plot_n_periods
    plot_minx,plot_maxx
    legend
    subplot(1,3,2)
    plot(fluidnormR.eta,fluidnormR.nz,'DisplayName','wall')
    plot(fluidnormR.a(i,:),fluidnormR.nz,'DisplayName',"fluid")
    title('radius $O(\frac{1}{\epsilon})$')
        plot_n_periods
    plot_minx,plot_maxx
    legend
        subplot(1,3,3)
    plot(fluidbigR.eta,fluidbigR.nz,'DisplayName','wall')
    plot(fluidbigR.a(i,:),fluidbigR.nz,'DisplayName',"fluid")
    title('radius $O(\frac{1}{\epsilon^2})$')
    plot_n_periods
    plot_minx,plot_maxx
    %xlim([xmin, xmax])

    legend
    sgtitle(sprintf('$t =%g$',fluidbigR.delt*i))
    pause(0.000001)

end
end