function slope_grow3(c,h,B,A,k)
z = 0:pi/100:4*pi/k;
etai = A*sin(k*z);
detai = A*k*cos(k*z);
dddetai = -A*k^3*cos(k*z);
i = 1;
figure(1)
diff =zeros(1,501);
maxloc = zeros(1,501);
for t = 0:0.01:5
eta = etai +t+0.1*(-c*z +1/3*detai.^2+1/3*etai -1/3*t -2/3*B*dddetai+h)*t;
diff(i)=max(eta)-min(eta);
[~,in] = max(eta);
maxloc(i) = in;
i= i+1;
clf, hold on
plot(eta,z)
plot(etai,z)
ax = gca;
ax.YDir ='reverse';
ylim([0,4*pi/k])
xlim([-A 4*A])
pause(0.01)
end
figure(2),clf, hold on 
plot(0:0.01:5,diff)


end
