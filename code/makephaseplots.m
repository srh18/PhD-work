function makephaseplots
folder = dir('outputs/n500/*.mat')
for i =	1:length(folder);
txt = strcat('outputs/n500/',folder(i).name);
load(txt)
clf
value.phaseplot(0,1)
savtxt = strcat('n500plots/',erase(folder(i).name,'.mat'));
saveas(gcf,savtxt,'epsc')
end
end
