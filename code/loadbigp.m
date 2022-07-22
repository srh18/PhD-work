function value = loadbigp(del,L,tol,n)

if nargin<4
    n = 256;
    if nargin<3
        tol = 1e-8;
    end
end
load(replace(sprintf('/Volumes/srh18/home/psbigp/hightol/CasePSn%gdel%.2fL%.2fT500tol%ginit0-1sin',n,del,L,log10(tol)),'.','-'),'value')
end