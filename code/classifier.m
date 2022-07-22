function classifier(folder,i)

for i = i:length(folder)
    txt = strcat('outputs/n500/',folder(i).name);
    load(txt)
    
    clf
    value.plotz0(0)
    
    
    clf
    value.phaseplot(0,2)
    disp(value.class)
    value.class = input('0 for steady, 1 for travelling, 2 for periodic, 3 for quasiperiodic, 4 for blows up, 5 for chaos, 6 for other');
    disp(value.notes)
    value.notes = input('Anything else');
    save(txt,'value')
    disp(i+1)
end
end
    