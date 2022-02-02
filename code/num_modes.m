function M = num_modes(value,del,L)
if nargin == 0
    value = shape_solveh;
    value.equation = 1;
    value.R = 1;
    
    del = 0:0.05:5;
    L = pi*(2.125:0.125:10);
end
M = zeros(length(del),length(L));
for i = 1:length(del)
    for j = 1:length(L)
        value.L = L(j);
        value.del = del(i);
        value = value.get_h;
        if value.eflag<=0 | value.q<1e-4 | sum(value.h0<0)>=1
            M(i,j) = -1;
        else
        try
            M(i,j) = value.num_unstable(value.h0);
        catch
            M(i,j) = 0;
        end
        end
    end
end
%contourf(L,del,M)
end
