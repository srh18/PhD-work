function getcases 
value = shape_solveh;
value.n = 500;

value.ep = 0.1;

value.T = 600;

value.wall_shape = 1;
value.R = 0.75;
L = 2:0.5:10;
value.del = 0.1;
for l = L
    value.L = l*pi;
    value = value.odedyn;

    save(erase(sprintf('outputs/Casesn500del01FlatL%gpi',l),'.'),'value')
end
end