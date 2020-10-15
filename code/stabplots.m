function stabplots(Re, Bo, R)
clf, hold on 
for r = R
k = 0:0.001 :1.2;
w = 2/15*k.^2*Re+ 1/(3*Bo) *k.^2/r^2-1/(3*Bo) *k.^4;
plot(k,w,'DisplayName',sprintf('$R= %f$',r))
end
legend
end