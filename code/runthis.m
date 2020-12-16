value = shape_solveh;
value.n = 200;

value.ep = 0.1;

value.T = 1000;
value.reltol = 1e-8;
value.abstol = 1e-10;


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

    save(erase(sprintf('outputs/R%gBo%gdel%g',R,Bo,del),'.'),'value')
            
end
end

