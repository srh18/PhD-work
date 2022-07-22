function [y,t] = small_wall(L,R,Bo)
if nargin<3
    Bo = 1;
    if nargin<2
        R = 1;
        if nargin<1
            L = 4*pi;
        end
    end
end
x = 0:pi/128:2*pi-pi/128;
theta = atan((3*Bo*L^3*R^2)/(2*pi*(4*pi^2-L^2)));
if theta<0
    theta = theta + pi;
end
[y,t] = shift_ode(-cos(theta),theta,cos(x));
end