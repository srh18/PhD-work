for del = [0 0.1 0.5 1 1.5 2] 
    for L =1:0.25:10
        try
            loadn500(del,L);
        catch
            disp([del,L])
        end
    end
end

            