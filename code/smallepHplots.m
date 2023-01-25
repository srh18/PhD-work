function smallepHplots(Bo)
fig = gcf;
clf
fig.WindowStyle = 'normal';
fig.Position =[587 1129 800 600];

hold on 

% 
% fig.WindowStyle = 'normal';
% fig.Position = [100 100 600 500];
xlabel('$\frac{L}{\pi}$')
ylabel('$A$')

title('Amplitude of fluid disturbance for small $\epsilon$')
legend
ax = gca;
ax.YScale = 'log';
fig2 = figure;
fig2.WindowStyle = 'normal';
fig2.Position =[587 1129 800 600];
hold on 
% fig2.WindowStyle = 'normal';
% fig2.Position = [100 100 600 500];
xlabel('$\frac{L}{\pi}$')
ylabel('$\frac{\theta}{\pi}$')
title('Phase shift of fluid disturbance from wall disturbance')
legend('Location','northeast')
L = pi/8:pi/100:8*pi;

%L(176) = 2pi

ax = gca;
for R = [0.5 1 2]
  figure(fig);  
[A,theta] = smallepH(L,R,Bo);





plot(L/pi,A,'DisplayName',sprintf('$R =%g$',R) )
 


figure(fig2);


plot(L/pi,theta/pi,'DisplayName',sprintf('$R =%g$',R))

end
saveas(fig,'../plots/smallepHA','epsc')
saveas(fig2,'../plots/smallepHsmtheta','epsc')
end