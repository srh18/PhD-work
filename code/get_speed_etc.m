function get_speed_etc()
LASTN = maxNumCompThreads(32);
for str_i =3:4
%str_i = 5;
str = {'/Volumes/srh18/home/psbigp/R1/','/Volumes/srh18/home/psbigp/R1/highres/','/Volumes/srh18/home/psbigp/R1/force1/','/Volumes/srh18/home/psbigp/R1/force2/','/Volumes/srh18/home/psbigp/R1/force3/','/Volumes/srh18/home/psbigp/R1/forcealt1/','/Volumes/srh18/home/psbigp/R1/longrun/'};
%files = split(ls('psbigp/R1'));
%files = split(ls('psbigp/R1/highres'));
files = split(ls(str{str_i}(1:end-1)));
% h1 = figure();
% set(h1, 'Visible', 'off');
mat = [];
for i = 1:length(files)
    file = files{i};
    if length(file)>4
    if file(1:4) == 'Case'
    %load(strcat('/Volumes/srh18/home/psbigp/',file),'value')
    %load(strcat('psbigp/R1/',file),'value')
%     load(strcat('psbigp/R1/highres/',file),'value')
    load(strcat(str{str_i},file),'value')
    if value.mass_err<1e-10
    which get_peak_info(value)
    [npks,time_periodic,c,T] = value.get_peak_info();
    value = value.get_h;
    if value.eflag>0
    [d,A,ps,ts,pd,td] = value.steady_mean_diff;
    end
    mat = [mat; value.L,value.del, npks,time_periodic,c,T,d,A,ps,ts,pd,td];
    end
    end
   
    end
end
%save('psbigp/R1/matrix/info','mat')
%save('psbigp/R1/highres/matrix/flatinfo_plus','mat')
save(strcat(str{str_i},'matrix/infodiff'),'mat')
end
end
