function value = loadps(del,l)
load(replace(sprintf('/Volumes/srh18/home/PSout/CasePSdel%.2gL%.2gT600',del,l),'.','-'),'value')
end
