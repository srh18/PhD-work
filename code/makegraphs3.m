function makegraphs3
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot ,'defaulttextInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
L = 1:0.1:5.7;

figure, clf, hold on
fig = gca;
for i = 1:length(L)
    load(erase(sprintf('outputs/CasesvlqFlatL%gpi',L(i)),'.'),'value')
    value = value.find_n_peaks;
    mpks(i) = mean(value.npks);
    modpks(i) = mode(value.npks);
    
end
plot(L,mpks,'DisplayName','mean')
plot(L,modpks,'DisplayName','mode')
legend


end