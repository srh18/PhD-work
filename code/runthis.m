value = shape_solveh;
value.n = 200;

value.ep = 0.01;
for R = [0.9 1 1.1]
    for L = [2*pi]
        for delta = [0.1 1]
            value.R = R;
            value.L = L;
            value.del = delta;
            
            value = value.get_h;
            value = value.dynamics(400, value.h+ 0.1*sin(value.z));
            save(erase(sprintf('outputs/L%gR%gdel%g',L,R,delta),'.'),'value')
            
        end
    end
end