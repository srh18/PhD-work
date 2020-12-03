value = shape_solveh;
value.n = 200;
value.R = 0.75;
value.ep = 0.1;

value.T = 500;
for del = [0.5 1 2]
    value.del = del;
for Bo = [1,25/18,5*16/9,77/9*5/2]
    value.Bo = Bo;
    value = value.odedyn;

    save(erase(sprintf('R075Bo%gdel%g',Bo,del),'.'),'value')
            
end
end
