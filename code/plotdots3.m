function plotdots3
hold on 
load('/Volumes/srh18/home/psbigp/R1/highres/matrix/info_new','npk_mat','time_periodic_mat','c_mat')
L_vec = 2.25:0.25:10;
c ={'b','g','r'};
del_vec= [ 0.1 0.5 1 1.5 2];
for i = 1:length(del_vec)
    del = del_vec(i);
    
    for j = 1:length(L_vec)
        L = L_vec(j);
        if time_periodic_mat(i,j) == 0
            shape = 's';
        
            
        else
            shape = 'o';
        end

        scatter(L, del,200,c{npk_mat(i,j)},'filled','Marker',shape)

    end

end
xlabel('$\frac{L}{\pi}$')
ylabel('$\delta$')
xlim([2,10.25])
ylim([-0.05,2.05])
title('Long term behaviour depending on amplitude and wavelength')

end
        
            