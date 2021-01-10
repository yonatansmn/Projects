function [h1,h2,h3,h4]=plotData(frame,data,i,thermal,r,maxd)
h3=0; h4=0;
for k=1:length(i)
    f=find (data(i(k),:)==0,1);
    if f>0
        p=1:(f-1);
        p(1)=data(i(k),1); p(2)=data(i(k),2);
        for j=3:2:f-1
            p(j)=data(i(k),j);
            p(j+1)=data(i(k),j+1);
        end
    else
        p=1:(length(data(i(k),:)));
        p(1)=data(i(k),1); p(2)=data(i(k),2);
        for j=3:2:size(data,2)-1
            p(j)=data(i(k),j);
            p(j+1)=data(i(k),j+1);
        end
    end

h1=plot(frame,p(1),p(2),'xr'); hold on
h2=plot(frame,p(3:2:end-1),p(4:2:end),'xb'); hold on
if thermal
h3=line (frame,[p(1),r+maxd],[p(2),r+maxd*sind(30)]); hold on
h4=line (frame,[p(1),r+maxd],[p(2),-r-maxd*sind(30)]); hold on
end
end

% hold off
end

% 
% if f>0
%     for j=3:2:f-1
%         p=plot(frame,data(i,j),data(i,j+1),'xk'); hold on
%     end
% else
%     for j=3:2:size(data,2)-1
%         hh=plot(frame,data(i,j),data(i,j+1),'xk'); hold on
%     end
% end
% end

    