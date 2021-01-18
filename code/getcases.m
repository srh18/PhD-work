function getcases 
value = shape_solveh;
value.n = 200;
value.notol = 1;
value.ep = 0.1;

value.T = 400;

value.wall_shape = 3;
value.R = 0.75;
L = 1:0.1:10;
for l = L
    value.L = l*pi;
    value = value.odedyn;

    save(erase(sprintf('outputs/CasesvvlqFlatL%gpi',l),'.'),'value')
end
end