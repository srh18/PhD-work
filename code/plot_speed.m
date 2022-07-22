function c = plot_speed(del,L,spectral)
if nargin ==1
    L = (2.25:0.25:10)*pi;
    spectral = 0;
end
c = [];
for l = L
    %try
        value = loadbigpR1(del,l);
        if spectral == 1
            
            c = [c;value.get_c_spectrally];
        else
            c = [c;value.get_c(0)];
        end
%     catch
%         warning(sprintf('case del %g, L %g not found',del,l))
%         value.L = l;
%         value = value.odedyn(1+ 0.1*sin(value.nz*2*pi));
%             
%             
%         save(replace(sprintf('/Volumes/srh18/home/psbigp/R1/CasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',value.n,value.del,value.L,log10(value.reltol)),'.','-'),'value')
%         c = [c,value.get_c(0)];
%     end
    
end
if spectral == 1
    plot(L,c(:,2),'x')
else
    plot(L,c,'x')
end
end