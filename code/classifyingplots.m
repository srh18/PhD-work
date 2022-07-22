function classifyingplots

%time_periodic
fig = figure(100);
fig.WindowStyle ='normal';
fig.Position = [100 100 800 300];
fig2 = figure(101);
fig2.WindowStyle ='normal';
fig2.Position = [900 100 400 300];
value = loadbigpR1(1,4*pi,1);
figure(100)
value.ploth2
xlim([400 500])
figure(101);
value.phaseplotEt(300)
saveas(figure(100),'../plots/psbigp/timeperiodich2','epsc')
saveas(figure(101),'../plots/psbigp/timeperiodicphase','epsc')
value = loadbigpR1(1.5,7.5*pi,8);
figure(100)
value.ploth2
xlim([1800,2000])
saveas(figure(100),'../plots/psbigp/4timeperiodich2','epsc')


figure(101)
value.phaseplotEt(1700)
saveas(figure(101),'../plots/psbigp/4timeperiodicphase','epsc')
value = loadbigpR1(0.5,8.5*pi,8);
figure(100)
value.ploth2
xlim([1300,2000])
saveas(figure(100),'../plots/psbigp/quasih2','epsc')


figure(101)
value.phaseplotEt(1300)
saveas(figure(101),'../plots/psbigp/quasiphase','epsc')

value = loadbigpR1(1,8.25*pi,8);
figure(100)
value.ploth2
xlim([1300,2000])
saveas(figure(100),'../plots/psbigp/chaotich2','epsc')


figure(101)
value.phaseplotEt(1300)
saveas(figure(101),'../plots/psbigp/chaoticphase','epsc')
value = loadbigpR1(3,4*pi,1);
figure(100)
value.ploth2
xlim([300 500])
figure(101);
value.phaseplotEt(300)
saveas(figure(100),'../plots/psbigp/dropleth2','epsc')
saveas(figure(101),'../plots/psbigp/dropletphase','epsc')