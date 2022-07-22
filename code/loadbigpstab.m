function value = loadbigpstab(del,L,tol,n)

if nargin<4
    n = 256;
    if nargin<3
        tol = 1e-8;
    end
end
load(replace(sprintf('/Volumes/srh18/home/psbigp/stabchange/CasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',n,del,L,log10(tol)),'.','-'),'value')
end