function plot_means(m)
del = [0 0.1 0.5 1 1.5, 2];
L = 2.75:0.25:10;
figure,hold on

for d = del
    meanval = [];
    for i  =1:length( L)
        value = loadn500(d,L(i));
        meanval =[meanval mean(value.h(floor(length(value.t)/2):end,m*value.n/4+1))];
    end
    plot(L,meanval)
end
legend(strcat('$\delta =$ ',string(del)))
end

        
        
