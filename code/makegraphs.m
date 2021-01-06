function makegraphs()
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot ,'defaulttextInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
for R = [1 1.25]
for L = 1:4
    for i = 1:2
        figure(i), clf, hold on;
    end
    for del = [0 0.1 0.5 1 2]
        
        load(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'))
        value.t = value.sol.x;
        figure(1)
        value.ploth2(0)
        figure(2)
        value.ploth2(1)

    end
    for i = 1:2
        figure(i)
        legend('$\delta = 0$','$\delta = 0.1$','$\delta = 0.5$','$\delta = 1$','$\delta = 2$')
        xlabel('$t$')
    end
    figure(1)
    title(sprintf('Energy for $L = %g \\pi, R = %g$ with different wall amplitudes',L,R))
    ylabel('$\int|h^2|$')
    saveas(figure(1),erase(sprintf('Energy L%gpi R%g sine amplitudes',L,R),'.'),'epsc')
    figure(2)
    title(sprintf('Energy difference from steady state for $L = %g \\pi, R = %g$ with different wall amplitudes',L,R))
    ylabel('$\int|(h-h_0)^2|$')
    saveas(figure(2),erase(sprintf('Energy diff L%gpi R%g sine amplitudes',L,R),'.'),'epsc')


end

for del = [0 0.1 1 ]
    for i = 1:2
        figure(i), clf, hold on;
    end
    for L = 1:2
        load(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'))
        value.t = value.sol.x;
        figure(1)
        value.ploth2(0)
        figure(2)
        value.ploth2(1)
       
    end
    for i = 1:2
        figure(i)
        legend('$L = \pi$','$L = 2\pi$','$L= 3\pi$','$L = 4\pi$')
        xlabel('$t$')
    end
    figure(1)
    title(sprintf('Energy for $\\delta = %g , R = %g$ with different wall amplitudes',del,R))
    ylabel('$\int|h^2|$')
    saveas(figure(1),erase(sprintf('Energy delta%g R%g sine wavelengths',del,R),'.'),'epsc')
    figure(2)
    title(sprintf('Energy difference from steady state for $\\delta = %g, R = %g $ with different wall amplitudes',del,R))
    ylabel('$\int|(h-h_0)^2|$')
    saveas(figure(2),erase(sprintf('Energy diff delta%g R%g sine amplitudes',del,R),'.'),'epsc')
    
end
end
for L = 1:4
    
    for del = [0.1 1 2]
         for i = 1:2
        figure(i), clf, hold on;
    end
        for R = [0.75 1 1.25] 
            load(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'))
        value.t = value.sol.x;
        figure(1)
        value.ploth2(0)
        figure(2)
        value.ploth2(1)
        end
            for i = 1:2
        figure(i)
        legend('$R = 0.75$','$R = 1$','$R = 1.25')
        xlabel('$t$')
        
            end
    figure(1)
    title(sprintf('Energy for $\\delta = %g , L = %g$ with different radii',del,L))
    ylabel('$\int|h^2|$')
    saveas(figure(1),erase(sprintf('Energy delta%g L%g radii',del,L),'.'),'epsc')
    figure(2)
    title(sprintf('Energy difference from steady state for $\\delta = %g, L = %g $ with different radii',del,L))
    ylabel('$\int|(h-h_0)^2|$')
    saveas(figure(2),erase(sprintf('Energy diff delta%g L%g radii',del,L),'.'),'epsc'))
    
    end
end

end