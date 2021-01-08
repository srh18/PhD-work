value = shape_solveh;
value.n = 200;

value.ep = 0.1;

value.T = 600;
value.reltol = 1e-6;
value.abstol = 1e-8;

for R = [0.75, 1, 1.25]
    value.R = R;
    

for del = [0.1 0.5 1 2 0]
    if del == 0 
        value.wall_shape = 3;
    else
        value.wall_shape = 0;
        value.del = del;
    end
for L = [1.25 1.5 1.75]
    value.L = L*pi;
    value = value.odedyn;

    save(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'),'value')
            
end
end
end

