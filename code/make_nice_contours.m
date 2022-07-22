function  l = make_nice_contours(m)
if nargin <1
    m = 1;
end
load('/Volumes/srh18/home/QH2Nunstablemat2highres.mat')
l = {};
n = 5;
v = contour(L/pi,del,M, [-1 0 1 2 3 4 5 ]);
fig = gca;
fig.Children(1).delete
hold on
i = 1;
j = v(2,i);
vec = v(:,i+1:j+i);
p = polyfit(vec(1,:),vec(2,:),n);
x = 0.5:0.05:1.75;
y = polyval(p,x);
plot(x,y,'k')
l{1} = [x;y];

i = i+j+1;
j = v(2,i);

vec = v(:,i+1:j+i);
p = polyfit(vec(1,:),vec(2,:),n);
x = 2.4:0.05:6.4;
y = polyval(p,x);
plot(x,y,'k')
l{2} = [x;y];
if m==1
i = i+j+1;
j = v(2,i);
xline(2,'k')
i = i+j+1;
j = v(2,i);
vec = v(:,i+1:j+i);
vec1 = vec(:,1:72);
vec2 = vec(:,72:82);
vec3 = vec(:,82:end);
p = polyfit(vec1(1,:),vec1(2,:),n);
y1 = polyval(p,[4.11:0.01:4.19 4.2:0.1:7.5]);
%plot([4.1:0.01:4.19 4.2:0.1:7.5],y1)
p = polyfit(vec2(1,:),vec2(2,:),n);
y2 = polyval(p,4.11:0.01:4.78);
%plot(4.11:0.01:4.78,y2)
p = polyfit(vec3(1,:),vec3(2,:),n);
y3 = polyval(p,[2:0.1:4.7 4.71:0.01:4.78]);
%plot([2:0.1:4.7 4.71:0.01:4.8],y3)
x = [2:0.1:4.7 4.71:0.01:4.78 flip(4.11:0.01:4.78) 4.11:0.01:4.19 4.2:0.1:7.5 ];
y = [y3 flip(y2) y1];
l{3} = [x;y];
plot(x,y,'k')

i = i+j+1;
j = v(2,i);

vec = v(:,i+1:j+i);
p = polyfit(vec(1,:),vec(2,:),n);
x = 4:0.05:7.8;
y = polyval(p,x);
y(1) = 0;
plot(x,y,'k')
l{4} = [x;y];
i = i+j+1;
j = v(2,i);

vec = v(:,i+1:j+i);
p = polyfit(vec(1,:),vec(2,:),n);
x = 6:0.05:8.2;
y = polyval(p,x);
y(1) = 0;
plot(x,y,'k')
l{5} = [x;y];
i = i+j+1;
j = v(2,i);

vec = v(:,i+1:j+i);
p = polyfit(vec(1,:),vec(2,:),n);
x = 8:0.05:9.2;
y = polyval(p,x);
y(1) = 0;
plot(x,y,'k')
l{6} = [x;y];
i = i+j+1;
j = v(2,i);


x = [10 10 10.1];
y = [0 0.05 1.2];
plot(x,y,'k')
l{7} = [x;y];
end
ylim([0 5])

xlim([0.5,10.1])

end