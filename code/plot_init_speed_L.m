function plot_init_speed_L(del)
hold on
value = shape_solveh;
if del == 0
    value.wall_shape = 3;
end
value.del = del;
L = 2*pi:pi/16:6*pi;
value.n = 256;
value.equation = 1;
value.ps = 1;
value.R = 1;
value.Bo = 1;
value.Re = 1;
c = [];
for l = L
    value.L = l;
    cs = value.init_speed;
    plot(l,cs,'bx')
    %c = [c cs(end)];
end
%plot(L,c)
end