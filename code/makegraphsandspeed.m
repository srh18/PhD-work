function makegraphsandspeed


n = 256;
tol = 1e-8;
l = (2.25:0.25:10)*pi;
Del = [0.1 0.5 1 1.5 2];
npk_mat = zeros(length(Del),length(l));
c_mat = zeros(length(Del),length(l));
time_periodic_mat = zeros(length(Del),length(l));
for i = 2:length(Del)
    del = Del(i);
    for j = 1:length(l)
        L = l(j);
        filestring = replace(sprintf('CasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',n,del,L,log10(tol)),'.','-');
        
        load(replace(sprintf('/Volumes/srh18/home/psbigp/R1/highres/CasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',n,del,L,log10(tol)),'.','-'),'value');
        [npks,time_periodic,c,T] = value.get_peak_info;
        npk_mat(i,j) = npks;
        time_periodic_mat(i,j) = time_periodic;
        c_mat(i,j) = c;
        value.ploth2
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/h2/%s',filestring),'epsc')
        clf
        value.plotz0
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/z0/%s',filestring),'epsc')
        clf
        value.phaseplot(0,2)
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/phase/%s',filestring),'epsc')
        clf
        value.phaseplotEt
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/phasee/%s',filestring),'epsc')
        clf
        value.plot_period(500,0.1) 
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/period/%s',filestring),'epsc')
        clf
        value.surfdata(c)
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/surfdata/%s',filestring),'epsc')
        ylim([500, floor(2.2*T)+500])
        saveas(gcf,sprintf('/Volumes/srh18/home/psbigp/R1/highres/plots/subsurfdata/%s',filestring),'epsc')
    end
    save('/Volumes/srh18/home/psbigp/R1/highres/matrix/info','npk_mat','time_periodic_mat','c_mat')
end
end