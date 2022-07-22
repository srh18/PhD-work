function mat = load_good_mat(n,d)
if nargin<2
    d = 0;
if nargin==0
    n = 0;
end
end
d = d+1;
dstring ={'','d','diff'};
if n==0
    load(sprintf('/Volumes/srh18/home/psbigp/R1/matrix/info%s.mat',dstring{d}))
elseif n==1
    load(sprintf('/Volumes/srh18/home/psbigp/R1/highres/matrix/flat%sinfo.mat',dstring{d}))
elseif n==2
    load(sprintf('/Volumes/srh18/home/psbigp/R1/highres/matrix/info%s.mat',dstring{d}))
elseif n==3
    load(sprintf('/Volumes/srh18/home/psbigp/R1/force1/matrix/info%s.mat',dstring{d}))
elseif n==4
    load(sprintf('/Volumes/srh18/home/psbigp/R1/force2/matrix/info%s.mat',dstring{d}))

elseif n ==5
    load(sprintf('/Volumes/srh18/home/psbigp/R1/force3/matrix/info%s.mat',dstring{d}))

elseif n ==6 
    load(sprintf('/Volumes/srh18/home/psbigp/R1/forcealt1/matrix/info%s.mat',dstring{d}))
elseif n == 7
    load(sprintf('/Volumes/srh18/home/psbigp/R1/longrun/matrix/info%s.mat',dstring{d}))
end
mat = mat(mat(:,1)>2*pi,:);
[~,idx] = sort(mat(:,1));
mat = mat(idx,:);
end
