function M = num_modes(value,del,L)
if nargin == 0
    value = shape_solveh;
    value.n = 256;
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
shading interp
 text(3.25,4.25,'No steady state','FontSize',20)
 text(2.75,2.9,'Stable','FontSize',20)
 text(10.05,0.15,'$5$','FontSize',15)
 text(2.75,1.5,'$k =1$','FontSize',20)
text(4.75,1.5,'$k = 2$','FontSize',20)
text(6.75,1.5,'$k = 3$','FontSize',20)
text(8.75,1.5,'$k = 4$','FontSize',20)
title('Number of unstable modes for walls $\eta = \delta\cos\left(\frac{2\pi z}{L}\right)$')
xlabel('$\frac{L}{\pi}$')
ylabel('$\delta$')
saveas(gcf,'../plots/psbigp/modeshighres','epsc')
end
