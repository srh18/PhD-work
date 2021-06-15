function value = loadn500(del,L)
del = erase(string(del),'.');
L = erase(string(L), '.');
load(sprintf('//Volumes//srh18//home//n500//Casesn500del%sL%spi.mat',del,L))
end
