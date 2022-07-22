function [y,t] = shift_ode(delta,theta,init)

if nargin<3
    n = 256;
    x = 0:pi/(n/2):2*pi-pi/(n/2);
    init = cos(x);
    
end

[t,y] = ode45(@(t,y) shift_fun(t,y,delta,theta) ,0:0.1:10,init');

end
function yt = shift_fun(t,y,delta,theta)
l = theta/(2*pi);
yt = delta*[y(floor(l*end)+1:end) ;y(1:floor(l*end))];

end
