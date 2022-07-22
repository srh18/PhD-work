function makeplotsbigp
LASTN = maxNumCompThreads(48);
files = split(ls('/Volumes/srh18/home/psbigp'));
% h1 = figure();
% set(h1, 'Visible', 'off');

for i = 1:length(files)
    file = files{i}
    if file(1:4) == 'Case'
    %load(strcat('/Volumes/srh18/home/psbigp/',file),'value')
    load(strcat('/Volumes/srh18/home/psbigp/',file),'value')
    clf
    value.ploth2
    saveas(gcf,strcat('/Volumes/srh18/home/psbigp/plots/h2/',file(1:end-4)),'epsc')
    
    clf
    value.phaseplota(0,2)
    saveas(gcf,strcat('/Volumes/srh18/home/psbigp/plots/phase/',file(1:end-4)),'epsc')
    end
end
end
