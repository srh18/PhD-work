
[hmax,hloc] = findpeaks(value.h2norm);
[hmin,hmloc] = findpeaks(-value.h2norm);
hloc = hloc(hloc>2000);
hm1 = hloc(end-3);
hm2 = hloc(end-1);
hl = hmloc(hmloc<hm2&hmloc>hm1);
value = value.get_h;
hold on,
for i = hm1:hm2-1
plot(value.nz,value.h(i,:),'color',[0.8,0.8,0.8],'linewidth',0.5,'HandleVisibility','off')
end
plot(value.nz,value.h(hm1,:),'linewidth',2,'DisplayName','maximum energy')
plot(value.nz,value.h(hl,:),'linewidth',2,'DisplayName','minimum energy')
%plot(value.nz,value.h(floor((hm1+hl)/2),:),'linewidth',2,'HandleVisibility','off')
%plot(value.nz,value.h(floor((hm2+hl)/2),:),'linewidth',2,'HandleVisibility','off')
plot(value.nz,sum(value.h(hloc(1):hloc(end)-1,:))/(hloc(end)-1-hloc(1)),'linewidth',2,'DisplayName','time averaged thickness')
plot(value.nz,value.h0,'linewidth',2,'DisplayName','steady state thickness')
legend

grid
xticks([0 0.25 0.5 0.75 1])
xlabel('$\frac{4}{L}$')
ylabel('$h$')