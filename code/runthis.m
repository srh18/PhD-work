value = shape_solveh;
value.L = 4*pi;
value.ep = 0.01;
value = value.get_h;
value = value.dynamics(400, value.h+ 0.1*sin(value.z));
save('test','value')