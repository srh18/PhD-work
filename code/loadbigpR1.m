function value = loadbigpR1(del,L,folder,tol,n)
inits  = {'','steadyinit/','highres/','force1/','force2/','force3/','forcealt1/','longrun/'};
if nargin<5
    n = 256;
    if nargin<4
        tol = 1e-8;
        if nargin<3
            folder = 1;
        end
    end
 
end
if folder == 100
    load(replace(sprintf('/Volumes/srh18/home/epseq/R1/CasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',n,del,L,log10(tol)),'.','-'),'value')
else
    
load(replace(sprintf('/Volumes/srh18/home/psbigp/R1/%sCasePSn%gR1del%.2fL%.3fT500tol%ginit0-1sin',inits{folder},n,del,L,log10(tol)),'.','-'),'value')
end
end