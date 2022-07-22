function plot_init_speed(L)
del = 0:0.1:2;
value = shape_solveh;
value.L = L;
value.n = 256;
value.equation = 1;
value.ps = 1;
value.R = 1;
value.Bo = 1;
value.Re = 1;
c = [];
for d = del
    value.del = d;
    cs = value.init_speed;
    c = [c cs(end)];
end
plot(del,c)
end