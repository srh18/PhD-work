function plotsnap(vfd,vps,t)
clf,hold on 
plot(vfd.z,vfd.h(t*20,:))
plot(vps.z,vps.h(t*20,:))
end