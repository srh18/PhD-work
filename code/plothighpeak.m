function plothighpeak(L) 
delta = [0 0.1 0.5 1 1.5 2];
for i = 0:3
    figure,clf,hold on 
for del = delta
    value = loadn500(del,L);
    value.peakpoints(i)
end
legend(strcat('$\delta = $',string(delta)))
end
end

    