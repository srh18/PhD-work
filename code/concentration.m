function w = concentration()
L = 2*pi;
nz = 256;
nr = 100;

z =  0:L/nz:L-L/nz;   
r = (0:1/nr:1)';
h = 1-0.1*cos(2*pi*z/L);
w = h.^2.*(r - r.^2/2);
end