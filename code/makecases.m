fluid = shape_solve;
for i = 1:3
for L = [ 2*pi 10] %1
for ep = [ 1e-2 1e-1]% 1e-4 1e-3
for R = [1/ep 1 ep]
for del = [0.1 1 10]
for del2 = [ 0.1 1]
fluid.wall_shape = i;
fluid.L = L;
fluid.ep = ep;
fluid.R = R;
fluid.del = del;
fluid = fluid.get_a;
try
fluid = fluid.dynamics(50,fluid.a+ del2*cos(fluid.z));
catch Me
fprintf('Wall%gL%gep%gR%gdel%gdel2%g returned error %s',i,L,ep,R,del,del2,Me.message)
continue
end
save(erase(sprintf('outputs/Wall%gL%gep%gR%gdel%gdel2%g',i,L,ep,R,del,del2),'.'),'fluid')
end
end
end
end
end
end
