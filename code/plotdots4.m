function plotdots4(n,Msize,shift)
if nargin<3
    shift = 0;
if nargin<2
    Msize = 200;
    if nargin<1
        n = 0;
    end
end
end
hold on 
mat = load_good_mat(n,2);
shape ={'h','o','v','s'};
for i =1:length(mat)
    if mat(i,4) == 1
        c = 'b';
    elseif mat(i,4)>=0.9
        c = 'c';
    elseif mat(i,4) >= 0.5
        c = 'g';
    elseif mat(i,4) == -1
        c = 'k';
    else
        c = 'r';
    end
    
    
    scatter(mat(i,1)/pi, mat(i,2)+shift,Msize,c,'filled','Marker',shape{mat(i,3)+1})
    
end
xlabel('$\frac{L}{\pi}$')
ylabel('$\delta$')
xlim([2,10.25])
ylim([-0.05,2.05])
title('Long term behaviour depending on amplitude and wavelength')

end
        