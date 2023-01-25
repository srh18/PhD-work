function fill_nice_contours
l = make_nice_contours;

patch([l{1}(1,:) 0.5],[l{1}(2,:) 5],[0.2422    0.1504    0.6603])
patch(l{2}(1,:) ,l{2}(2,:),[0.2422    0.1504    0.6603])
patch([l{1}(1,:) 2 2 0],[l{1}(2,:) 5 0 0],[0.2784    0.3529    0.9763])
patch([l{2}(1,:) flip(l{3}(1,:))] ,[l{2}(2,:) flip(l{3}(2,:))],[0.2784    0.3529    0.9763])

patch([l{3}(1,:) flip(l{4}(1,:)) 2 2],[l{3}(2,:) flip(l{4}(2,:)) 0 5],[0.1540    0.5902    0.9218])
patch([l{4}(1,:) flip(l{5}(1,:)) ],[l{4}(2,:) flip(l{5}(2,:)) ],[0.0770    0.7468    0.7224])
patch([l{5}(1,:) flip(l{6}(1,:)) ],[l{5}(2,:) flip(l{6}(2,:)) ],[0.5187    0.7982    0.3374])
patch([l{6}(1,:) 10.1 flip(l{7}(1,:))  ],[l{6}(2,:) 5 flip(l{7}(2,:)) ],[0.9946    0.7407    0.2394])
patch([l{7}(1,:) 10.1 ],[l{7}(2,:) 0 ],[0.9769    0.9839    0.0805])
 text(3.25,4.25,'No steady state','FontSize',20,'Color','white')
 text(3,2.9,'Stable','FontSize',20,'Color','white')
 text(10.015,0.15,'$5$','FontSize',15)
 text(3,1.5,'$1$','FontSize',20)
text(5,1.5,'$2$','FontSize',20)
text(7,1.5,'$3$','FontSize',20)
text(9,1.5,'$4$','FontSize',20)
title('Number of unstable modes for walls $\eta = \delta\cos\left(\frac{2\pi z}{L}\right)$')
xlabel('$\frac{L}{\pi}$')
ylabel('$\delta$')
f = gcf;
f.WindowStyle = 'normal';
f.Position = [100 100 1200 750];
saveas(f,'../plots/modeshighressmooth','epsc')
end