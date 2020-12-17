value = shape_solveh;
value.n = 400;

value.ep = 0.1;

value.T = 400;
value.reltol = 1e-7;
value.abstol = 1e-9;


for del = [0.1 1 0]
    if del == 0 
        value.wall_shape = 3;
    else
        value.wall_shape = 0;
        value.del = del;
    end
for R = [1/4,1/2,2/3,3/4 ,1, 1.25]
    value.R = R;
    value = value.odedyn;

    save(erase(sprintf('n400R%gBo1del%g',R,del),'.'),'value')
            
end
end

