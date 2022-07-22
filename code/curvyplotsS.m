
[hmax,hloc] = findpeaks(value.h2norm);
[hmin,hmloc] = findpeaks(-value.h2norm);
hm1 = hloc(end-2);
hm2 = hloc(end-1);
hl = hmloc(hmloc<hm2&hmloc>hm1);

hold on,
for i = hm1:hm2-1
plot(value.nz,value.S(i,:),'color',[0.8,0.8,0.8],'linewidth',0.5,'HandleVisibility','off')
end
plot(value.nz,value.S(hm1,:),'linewidth',2,'DisplayName','maximum energy')
plot(value.nz,value.S(hl,:),'linewidth',2,'DisplayName','minimum energy')
plot(value.nz,value.S(floor((hm1+hl)/2),:),'linewidth',2,'HandleVisibility','off')
plot(value.nz,value.S(floor((hm2+hl)/2),:),'linewidth',2,'HandleVisibility','off')
legend
xlabel('$\frac{1}{L}$')
ylabel('$h$')