%TOPO is the 512 ppd 15x15 km of the site topography. 
%TOPOels is the 128 ppd 200x200 km 
profile on
nonius='/Users/oa/Dropbox/spaceIL/Landing_sites/site_Nonius_Lat-36.6_Lon4.03_size0.25x0.31/';
wohler='/Users/oa/Dropbox/spaceIL/Landing_sites/site_Wohler_Lat-38_Lon29.5_size0.25x0.31/';
berzelius='/Users/oa/Dropbox/spaceIL/Landing_sites/site_Berzelius-opt1_Lat37.8_Lon54.7_size0.25x0.31/';

addpath 'C:\Users\Oagroup\Dropbox (Weizmann Institute)\spaceIL\updated\Landing_sites\site_Wohler_Lat-38_Lon29.5_size5x6.3'
addpath 'C:\Users\oagroup.WISMAIN\Dropbox (Weizmann Institute)\spaceIL\updated\Landing_sites\site_Wohler_Lat-38_Lon29.5_size5x6.3';
% pth=nonius; titl='Nonius';
% pth=berzelius; titl='Berzelius';
pth=wohler; titl='Wohler';

if ~exist('TOPO','var'),
    load(['maps.mat']);
end

% clf; vimp; %waduc land

z=TOPO; ppd=512;
%z=TOPOels; ppd=128;

R=1737;
lat0=Latitude;
kmpd=2*pi*R/360;
res=kmpd/ppd;
r=15;
figure
x=(1:size(z,2))*(res*cosd(lat0)); x=x-mean(x);
y=(1:size(z,1))*(res); y=y-mean(y);

% subplot(221) 
% imagesc(x,y,z); axis tight ij equal
% set(get(colorbar,'ylabel'),'string','Topography (m)')
% xlabel('km'); ylabel('km');
th=linspace(0,2*pi); hold on; 
% plot(sin(th)*r,cos(th)*r,'--k')
% title(titl);

% subplot(222) %plot elevation angle
[yy,xx]=ndgrid(y,x);
az=atan2d(yy,xx); %azimuth
d=sqrt(xx.^2+yy.^2)*1e3; %distance in meters
h=z-z(end/2,end/2); %height relative to center [0,0]
h=h-d.^2/(2.*R*1e3); %horizon reduction
sl=atand(h./d); %slope angle
% imagesc(x,y,sl); axis tight ij equal
maxaz=30;
hold on;
% plot(sin(th)*r,cos(th)*r,'--k')
% plot([0,cosd(maxaz)]*max(x)./cosd(maxaz),[0,sind(maxaz)]*max(x)./cosd(maxaz),'--k')
% plot([0,cosd(-maxaz)]*max(x)./cosd(maxaz),[0,sind(-maxaz)]*max(x)./cosd(maxaz),'--k')
% set(get(colorbar,'ylabel'),'string','Elevation angle (deg)')
% xlabel('km'); ylabel('km');
% title([num2str(ppd),' ppd']);

% subplot(235) 
azbw=1;
azbins=-180:azbw:180;
f=find(d>3*res*1e3);
scatter(az(:),sl(:),1,d(:)/1e3,'.'); box on;
hold on;
axis tight
plot([1,1]*maxaz,ylim,'--k');
plot(-[1,1]*maxaz,ylim,'--k');
xlabel('Azimuth (deg N of E)'); ylabel('Elevation angle (deg)');
set(gca,'xtick',-180:60:180);

subs=ceil((az(:)+azbw/2+180)/azbw);
msl=accumarray(subs(f),sl(f),[length(azbins),1],@max); %max slope for every azimuth. 
hold on; plot(azbins,msl,'k')
hold on; plot(azbins,ones(length(azbins),1)*3,'--k')
set(get(colorbar,'ylabel'),'string','Distance (km)')
set(gca,'pos',get(gca,'pos')+[0,0,.3,0]);
drawnow;
%stop

%same calcultaion for every point within the ellipse

f=find(d<r*1e3);
calcf=1;
if calcf,
maxh=zeros(size(f)); meanh=zeros(size(f));
for i=1:length(f)
    %for i=[],
    % progress(i,length(f),100); 
    xxr=xx-xx(f(i)); yyr=yy-yy(f(i));
    azr=atan2d(yyr,xxr);  
    dr=sqrt(xxr.^2+yyr.^2)*1e3;
    hr=z-z(f(i));
    hr=hr-dr.^2/(2.*R*1e3);
    slr=atand(hr./dr);
    %faz=find(abs(az)<30);
    fd=find(dr>3*res*1e3);
    subs=ceil((azr(:)+azbw/2+180)/azbw);
    msl=accumarray(subs(fd),slr(fd),[length(azbins),1],@max);
    % fazbins=find(abs(azbins)<=30);
    % maxh(i)=max(msl(fazbins)); 
    % meanh(i)=mean(msl(fazbins));
    maxh(i)=max(msl);
    meanh(i)=mean(msl);
end
end

subplot(234)
s=sort(meanh);
plot(s,(1:length(s))/length(s))
hold on
s=sort(maxh);
plot(s,(1:length(s))/length(s))
legend('mean','max','Location','NorthWest'); legend('boxoff')
xlim([0,5]); set(gca,'xtick',0:5);
grid on
ylabel('CDF')
xlabel('Horizon Elevation (deg)')

% sz=size(d);
% i=sub2ind(sz,round(sz(1)/2),round(sz(2)/2));
% xxr=xx-xx(i); yyr=yy-yy(i);
% azr=atan2d(yyr,xxr);
% dr=sqrt(xxr.^2+yyr.^2)*1e3;
% hr=z-z(i);
% hr=hr-dr.^2/(2.*R*1e3);
% slr=atand(hr./dr);
% %faz=find(abs(az)<30);
% 
% fd=find(dr>3*res*1e3);
% subs=ceil((azr(:)+azbw/2+180)/azbw);
% msl=accumarray(subs(fd),slr(fd),[length(azbins),1],@max);
% fazbins=find(abs(azbins)<=30);
