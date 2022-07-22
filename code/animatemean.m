meanh = 0;
n= 0;
for t = 1001:4001
    clf, hold on 
    plot(value.z,value.h0,'k');
    plot(value.z,value.h(t,:));
    
    meanh = (meanh*n+value.h(t,:))/(n+1);
    n = n+1;
    plot(value.z,meanh,'--')
    pause(0.01)
end