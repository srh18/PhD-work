function value = loadn500(del,L)
del = erase(string(del),'.');
L = erase(string(L), '.');
load(sprintf('//Users//srh18//Documents//PhD//Code//outputs//n500//Casesn500del%sL%spi.mat',del,L))
end
