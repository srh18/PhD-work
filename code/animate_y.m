function animate_y(y,t)
loops = length(t);
            M(loops) = struct('cdata',[],'colormap',[]);
            
            for i = 1:loops

                plot(y(i,:))
                            
                drawnow
                
   
                M(i) = getframe(gcf);

                hold off
            end
end