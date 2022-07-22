function animate_mean(t0,tint,T,varargin)
filename = 0;
%filename = '123peaksmean';
record = ischar(filename);
i0 = [];
int = [];
ymin = [];
meanh = zeros(3,256);
n = 0;
if record
    v = VideoWriter(strcat('../video/',filename));
    open(v);
end
for i = 1:length(t0)
    varargin{i} = varargin{i}.get_h;
    value = varargin{i};
    
i0(i) = floor((t0(i)-value.t(1))/value.delt)+1;
int(i) = floor(tint/value.delt);
ivec = i0(i):int(i):i0(i)+T/int(i);
ymin(i) = min(min(value.h(ivec,:),[],2));
ymax(i) =max(max(value.h(ivec,:),[],2));
if int(i)<1
    int(i) =1;
end

end


clf

ax = gca;
ax.NextPlot = 'replaceChildren';
xlabel('$L$')
ylabel('$h$')





ax.YLim = [min(ymin)-0.1, max(ymax)+0.2];
ax.XLim = [0,value.L];
loops = floor(T/tint);
M(loops) = struct('cdata',[],'colormap',[]);

for i = 1:loops
    
    plot(value.z,value.h0,'k','DisplayName','Steady State')
    ax.NextPlot = 'add';
   
for j = 1:length(t0)
    value = varargin{j};
    t = i0(j) + int(j)*(i-1);
    ax.ColorOrderIndex = j;
    
    plot(value.z,value.h(t,:),'DisplayName',sprintf('Case %g',j));

    %meanh(j,:) = (meanh(j,:)*n+value.h(t,:))/(n+1);
    meanh(j,:) = mean(value.h(i0(j):t,:));
    ax.ColorOrderIndex = j;
  
    
    plot(value.z,meanh(j,:),'--','DisplayName','Time averaged thickness')
    
  
    
    

end
ax.Children =ax.Children([end 1:end-1]);
legend('NumColumns',2,'Location','north')
    drawnow
    
    
    M(i) = getframe(gcf);
    if record
        writeVideo(v,M(i));
    end
    ax.NextPlot = 'replaceChildren';
    n = n+1;
end

if record
    close(v);
end

end

