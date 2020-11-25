value = shape_solveh;
value.n = 200;

value.ep = 0.1;
value.del = 0.5;
value.T = 500;
value = value.odedyn;
save(erase(sprintf('testingcode',L,R,delta),'.'),'value')
            
