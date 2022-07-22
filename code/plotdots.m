function plotdots 
hold on 
 A = readmatrix('../dataclassification');
[j,k] = size(A);
for i = 2:j
    for l = 2:k
        if not(isnan(A(i,l)))
            
        
        val = floor(A(i,l));
        rem = A(i,l) - val;
        if rem ==0.5
            shape = 's';
            c = [0.4940 0.1840 0.5560];
        elseif rem == .25
            shape = 's';
            c = [0.4940 0.1840 0.5560];
        elseif rem ==.75
            shape = 's';
            c = [0.4940 0.1840 0.5560];
        elseif val == 0 
            shape = '<';
            c = [0.4660 0.6740 0.1880];
            
        elseif val == 5
            shape = 'h';
            c = [0 0 0];
        elseif val == 6
            shape = '>';
            c = [0.6350 0.0780 0.1840];
        else
            shape = 'o';
            c = [0 0.4470 0.7410];
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
        scatter(A(1,l), A(i,1),200,c,'filled','Marker',shape)
        end
    end

end
xlabel('Wavelength')
ylabel('Amplitude')
xlim([0.75,10.25])
ylim([-0.05,2.05])
title('Long term behaviour depending on amplitude and wavelength')
saveas(gcf,'../paper/plots/classifier','epsc')
end
        
            