delta = 1;
figure(22), clf, hold on
title('position of peaks')
L =4*pi;
for R = [0.9 1 1.1]
    load(erase(sprintf('outputs/L%gR%gdel%g',L,R,delta),'.'))
    value.T = 400;
    value = value.peak_speed;
    figure(21)
    plot(0:value.delt:value.T,value.peakpos)
    value.plot_features(103)
    
end