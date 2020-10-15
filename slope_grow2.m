function slope_grow2
z =0:pi/100:2*pi;
etai = exp(-z);
detai = -exp(-z);
dddetai = -exp(z) ;
for t = 0:0.01:1
eta = t*(1/6-z-3*detai.^2 +2etai + 2*dddetai+etai*t);
clf, hold on
plot(eta,z) 
plot(etai,z)
ax = gca;
ax.YDir = 'reverse';

ylim([0 2*pi])
pause(0.01)
end
end