function getcaseslow
value = shape_solveh;
value.n = 200;

value.ep = 0.1;

value.T = 400;
value.reltol = 1e-4;
value.abstol = 1e-6;

value.wall_shape = 0;
value.R = 0.75;
l = 6;
value.L = l*pi;
del = 0.1:0.1:2;
for d = del
    value.del = d;
    value = value.odedyn;

    save(erase(sprintf('outputs/lowres/Casesn200del%gL%gpi',d,l),'.'),'value')
end
end