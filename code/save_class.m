function save_class(value)
mc = metaclass(value);
proplist = {mc.PropertyList(logical(1-[mc.PropertyList(:).Dependent])).Name};
m = 1;
for i = proplist
    outdata = {i{1},value.(i{1})};
        m= 1+m;
        save(strcat('outputs/test',i{1}),'outdata')
end

end