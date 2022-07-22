function h = thickness(q,r)
h = (3*q*1e-6./(2*pi*9.81*r*1e-3)).^(1/3);
end