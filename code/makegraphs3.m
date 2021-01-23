function makegraphs3
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot ,'defaulttextInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
L = [1 1.25 1.35 1.36 1.5 1.8 2 2.2 2.5]

figure, clf, hold on
fig = gca;
for i = 1:length(L)
    load(erase(sprintf('outputs/Caseslqdel01L%gpi',L(i)),'.'),'value')
    value = value.find_n_peaks;
    mpks(i) = mean(value.npks);
    modpks(i) = mode(value.npks);
    maxpks(i) = max(value.npks);
    
end
plot(L,mpks,'DisplayName','mean')
plot(L,modpks,'DisplayName','mode')
plot(L,maxpks,'DisplayName','max')
legend


end