function slope_growb(c,h,b,A,k,t)

i = 1;

diff =zeros(1,501);
maxloc = zeros(1,501);
clf, hold on 
nz = 0:1/400:1;
plot(sin(nz*4*pi),nz,'DisplayName','Wall shape')
for B = b
    z = 0:pi/k/100:4*pi/k;
etai = A*sin(k*z);
detai = A*k*cos(k*z);
dddetai = -A*k^3*cos(k*z);
eta = (-c*z +1/3*detai.^2+1/3*etai  -2/3*B*dddetai+h)*t;
diff(i)=max(eta)-min(eta);
[~,in] = max(eta);
maxloc(i) = in;
i= i+1;

plot(eta,nz,'DisplayName',sprintf('$B= %g$',B))

ax = gca;
ax.YDir ='reverse';
end
legend



end
