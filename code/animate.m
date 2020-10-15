function animate(data,z,wall, speed)
if iscell(data)
    l = length(data{1});
    m = length(data);
else 
    l = length(data);
    m= 0;
end
% 
% xmax = max(max(data));
% xmin = min(min(data));
for i = 1:speed:l
    clf, hold on
    flip_y
    if wall ~= 0
    plot(wall,z)
    end
    if m==0;
    
    plot(data(i,:),z)
    else
        for j = 1:m
            plot(data{j}(i,:),z)
        end
    end
    plot_n_periods
    plot_minx,plot_maxx
    %xlim([xmin, xmax])

    title(sprintf('$t =%g$',0.01*i)) 
    pause(0.01)

end
end