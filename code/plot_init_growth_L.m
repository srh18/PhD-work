function plot_init_growth_L(del)
hold on
value = shape_solveh;
if del == 0
    value.wall_shape = 3;
end
value.del = del;
L = 2*pi:pi/16:10*pi;
value.n = 256;
value.equation = 1;
value.ps = 1;
value.R = 1;
value.Bo = 1;
value.Re = 1;
k0 = [];
k1 = [];
for l = L
    value.L = l;
    [k0s,k1s] = value.init_growth;

    k0= [k0 max(abs(k0s))];
    k1= [k1 k1s(end)];
end
plot(L,k0)
plot(L,-k1,'--')
end