function getcases2 
value = shape_solveh;
value.n = 500;

value.ep = 0.1;

value.T = 600;
value.reltol = 1e-5;
value.abstol = 1e-8;

value.R = 0.75;
for L = 1:0.25:2.75
for del = [ 0 0.1 0.5 1 1.5 2]
    value.del = del;
    if del ==0 
        value.wall_shape = 3;
    else
    
        value.wall_shape = 0;
    end
for l = L
    value.L = l*pi;
    value = value.odedyn;

    save(erase(sprintf('outputs/n500/Casesn500del%gL%gpi',del,l),'.'),'value')
end

end
end