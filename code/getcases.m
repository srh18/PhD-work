function getcases 
value = shape_solveh;
value.n = 200;

value.ep = 0.1;

value.T = 600;
value.reltol = 1e-6;
value.abstol = 1e-10;
value.wall_shape = 3;
value.R = 0.75;
L = [1 1.25 1.35 3*sqrt(10)/7 1.36 1.5 1.8 2 2.2 2.5 3 3.5 4 4.5 5 5.5 6 7 8 9 10];
for l = L
    value.L = l*pi;
    value = value.odedyn;

    save(erase(sprintf('CasesFlatL%gpi',l),'.'),'value')
end
end