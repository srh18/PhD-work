value = shape_solveh;
value.n = 200;

value.ep = 0.1;

value.T = 200;
value.reltol = 1e-6;
value.abstol = 1e-8;
value.R = 0.75;

for del = [0.1 0.5 1 2 0]
    if del == 0 
        value.wall_shape = 3;
    else
        value.wall_shape = 0;
        value.del = del;
    end
for L = [pi,2*pi,3*pi,4*pi ]
    value.L = L;
    value = value.odedyn;

    save(erase(sprintf('L%gBo1del%g',L,del),'.'),'value')
            
end
end

