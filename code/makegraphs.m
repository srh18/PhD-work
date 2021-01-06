function makegraphs()
R = 0.75;
for L = 1:4
    for i = 1:4
        figure(i), clf, hold on;
    end
    for del = [0 0.1 0.5 1 2]
        
        load(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'))
        value.t = value.sol.x;
        figure(1)
        value.ploth2(0)
        figure(2)
        value.ploth2(1)
        figure(3)
        plot(value.t, value.mass)
        value.h = value.hdiff;
        figure(4)
        plot(value.t, value.mass);
    end
    for i = 1:4
        figure(i)
        legend('$\delta = 0$','$\delta = 0.1$','$\delta = 0.5$','$\delta = 1$','$\delta = 2$')
        xlabel('$t$')
    end
    figure(1)
    title(sprintf('Energy for $L = %g \\pi$ with different wall amplitudes',L))
    ylabel('$\int|h^2|$')
    saveas(figure(1),sprintf('Energy L%gpi sine amplitudes',L),'epsc')
    figure(2)
    title(sprintf('Energy difference from steady state for $L = %g \\pi$ with different wall amplitudes',L))
    ylabel('$\int|(h-h_0)^2|$')
    saveas(figure(2),sprintf('Energy diff L%gpi sine amplitudes',L),'epsc')
    figure(3)
    title(sprintf('Mass for $L = %g \\pi$ with different wall amplitudes',L))
    ylabel('$\int|h|$')
    saveas(figure(3),sprintf('Mass L%gpi sine amplitudes',L),'epsc')
    figure(4)
    title(sprintf('Mass difference from steady state for $L = %g \\pi$ with different wall amplitudes',L))
    ylabel('$\int|h-h_0|$')
    saveas(figure(4),sprintf('Mass diff L%gpi sine amplitudes',L),'epsc')
end

for del = [0 0.1 1 ]
    for i = 1:4
        figure(i), clf, hold on;
    end
    for L = 1:4
        load(erase(sprintf('resultR%gL%gpiBo1del%g',R,L,del),'.'))
        value.t = value.sol.x;
        figure(1)
        value.ploth2(0)
        figure(2)
        value.ploth2(1)
        figure(3)
        plot(value.t, value.mass)
        value.h = value.hdiff;
        figure(4)
        plot(value.t, value.mass);
    end
    for i = 1:4
        figure(i)
        legend('$L = \pi$','$L = 2\pi$','$L= 3\pi$','$L = 4\pi$')
        xlabel('$t$')
    end
    figure(1)
    title(sprintf('Energy for $\\delta = %g $ with different wall amplitudes',del))
    ylabel('$\int|h^2|$')
    saveas(figure(1),sprintf('Energy delta%g sine wavelengths',del),'epsc')
    figure(2)
    title(sprintf('Energy difference from steady state for $\\delta = %g $ with different wall amplitudes',del))
    ylabel('$\int|(h-h_0)^2|$')
    saveas(figure(2),sprintf('Energy diff delta%g sine amplitudes',del),'epsc')
    figure(3)
    title(sprintf('Mass for $\\delta = %g $ with different wall amplitudes',del))
    ylabel('$\int|h|$')
    saveas(figure(3),sprintf('Mass delta%g sine amplitudes',del),'epsc')
    figure(4)
    title(sprintf('Mass difference from steady state for $\\delta = %g $ with different wall amplitudes',del))
    ylabel('$\int|h-h_0|$')
    saveas(figure(4),sprintf('Mass diff delta%g sine amplitudes',del),'epsc')
    
end
end