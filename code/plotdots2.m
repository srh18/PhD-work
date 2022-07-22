function plotdots2
hold on
A = readmatrix('../dataclassification2');
[j,k] = size(A);
ldg1 = 0;
ldg2 = 0;
ldg3 = 0;
ldg4 = 0;
for i = 2:j
    for l = 2:k

        addlegend = 0;
        if not(isnan(A(i,l)))
            
            
            val = floor(A(i,l));
            rem = A(i,l) - val;
            if val == 1
                if rem ==0.75
                    shape = '<';
                    c = [0.3010 0.7450 0.9330];
                    if ldg1 ==0
                        addlegend = 1;
                        d = '2 peaks to 1 peak';
                        ldg1 = ldg1+1;
                    end
                    
                else
                    shape = 'o';
                    c = [0 0.4470 0.7410];
                    if ldg2 ==0
                        addlegend = 1;
                        d = '1 peak';
                        ldg2 = ldg2+1;
                    end
                end
            elseif val == 2
                shape = '^';
                c = [0.8500 0.3250 0.0980];
                if ldg3 ==0
                    addlegend = 1;
                    d = '2 peaks';
                    ldg3 = ldg3+1;
                end
            elseif val == 3
                shape = 's';
                c = [0.4660 0.6740 0.1880];
                if ldg4 ==0
                    addlegend = 1;
                    d = '3 peaks';
                    ldg4 = ldg4+1;
                end
            end
            
            
            %         if (val >= 1) == (val <= 4)
            %             c = [0 0.4470 0.7410];
            %         elseif val == 2
            %             c = [0.8500 0.3250 0.0980];
            %         elseif val == 3
            %             c = [0.4940 0.1840 0.5560];
            %         elseif val == 4
            %             c = [0.9290 0.6940 0.1250];
            %        end
            if addlegend ==1
                
                scatter(A(1,l), A(i,1),200,c,'filled','Marker',shape,'DisplayName',d)
            else
                scatter(A(1,l), A(i,1),200,c,'filled','Marker',shape,'HandleVisibility','off')
            end
        end
    end
    
end
xlabel('Wavelength/$\pi$')
ylabel('Amplitude')
xlim([1.75,10.25])
ylim([-0.05,2.25])
title('Number of Peaks')
legend('Location','north','Orientation','horizontal')
saveas(gcf,'../plots/psbigp/classifypeaks','epsc')
end

