function  saven500(value)
del = erase(string(value.del),'.');
L = erase(string(value.L/pi), '.');

save(sprintf('//Users//srh18//Documents//PhD//Code//outputs//n500//Casesn500del%sL%spi.mat',del,L),'value')
end
