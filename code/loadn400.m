function value = loadn400(del,L)
del = replace(string(del),'.','-');
L = replace(string(L), '.','-');
load(sprintf('//Users//srh18//Documents//PhD//Code//outputs//n400//Casesn400del%sL%spi.mat',del,L))
end
