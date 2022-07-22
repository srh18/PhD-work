for i = 30:34
    figure(i),clf,hold on 
end
for suppression = [1e-12]

    for pad = [4]

        for filter_oscillation = [0 1]

            load(sprintf('outputs/testcase/Testode15s_sup%.0fpad%ifilter%i',-log10(suppression),pad,filter_oscillation),'test_value')
           
            figure(30)
            test_value.plot_mass
            figure(31)
            test_value.ploth2
            figure(32)
            plot(test_value.t,test_value.h(:,1))
            
            Fh = fft(test_value.h,[],2);
            figure(33)
            plot(log10(abs(Fh(:,value.n/2+1))))
            figure(34)
            plot(log10(abs(Fh(end,:))))
        end
        
    end
end
