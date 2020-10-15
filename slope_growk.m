function slope_growk(c,h,B,A,K,t)

i = 1;

diff =zeros(1,501);
maxloc = zeros(1,501);
clf, hold on 
nz = 0:1/400:1;
plot(sin(nz*4*pi),nz,'DisplayName','Wall shape')
for k = K
    z = 0:pi/k/100:4*pi/k;
etai = A*sin(k*z);
detai = A*k*cos(k*z);
dddetai = -A*k^3*cos(k*z);
eta = (-c*z +1/3*detai.^2+1/3*etai-2/3*B*dddetai+h)*t;
diff(i)=max(eta)-min(eta);
[~,in] = max(eta);
maxloc(i) = in;
i= i+1;

plot(eta,nz,'DisplayName',sprintf('$k= %g$',k))

ax = gca;
ax.YDir ='reverse';
end



legend
end
