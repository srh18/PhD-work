
function plot_amph(m,range)
figure(2*m+12)
clf, hold on
legend
value = shape_solveh;
value.ep = 0.1;

value.R = 1;
value.Re = 1;
value.L = 2*pi;
[A,lambda] = value.shift(range,m,1);
plot(range,A,'--','DisplayName','Small $\delta$ Approximation')
figure(2*m+13)
clf,hold on
legend
plot(range,lambda,'--','DisplayName','Small $\delta$ Approximation')
for del = [0.1 0.5 1 2]
    value.del = del;
    [A,lambda] = value.shift(range,m);
    figure(2*m+12)
    plot(range,A,'DisplayName',sprintf('$\\delta = %f$',del))
    figure(2*m+13)
    plot(range,lambda,'DisplayName',sprintf('$\\delta = %f$',del))

end
end