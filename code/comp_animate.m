function comp_animate( speed,case1,case2)
l = length(case1.a);

% xmax = max(max(data));
% xmin = min(min(data));
for i = 1:speed:l
    clf, hold on
    
    flip_y
    
   
    plot(case1.eta,case1.nz,'r--','DisplayName','wall1')
    plot(case1.S(i,:),case1.nz,'r','DisplayName',"fluid1")
    plot(case2.S(i,:),case2.nz,'b','DisplayName',"fluid2")
    plot(case2.eta,case2.nz,'b--','DisplayName','wall2')

        plot_n_periods
    plot_minx,plot_maxx
    

    
    legend
    title(sprintf('$t =%g$',case1.delt*i))
    pause(0.000001)

end
end