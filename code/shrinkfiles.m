function shrinkfiles
files = dir('*.mat');
fname = {files.name};
for i = fname
    try
        load(i{1})
    catch
        delete(i{1})
    end
    value.sol = 0;
    save(i{1},'value')
end
end

