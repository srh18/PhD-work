function makegraphs2
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot ,'defaulttextInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');


    
    figure, clf, hold on
    fig = gca;
    for l = 3.1:0.1:3.7
        load(erase(sprintf('outputs/CasesvlqFlatL%gpi',l),'.'),'value')
        value.ploth2
        fig.Children(1).DisplayName = sprintf('$L = %g$',l);
    end
    ylabel('$\int|(h)^2|$')
    xlabel('$t$')
    legend
    
end
