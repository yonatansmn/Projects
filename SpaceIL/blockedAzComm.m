% According to section 2.3.1 in "Landing Site Requirements" document, Due to the spacecraft's minimal antenna pointing, 
% angle local slope, lunar liberation and landing mechanics, the line-of-sight to Earth shall be more than 25deg above 
% local horizon [Req 3LAN-3.1]. 
% 
% this code calculate the precentage of blocked azimuths at every point within the landing ellipse.

site='Posidonius';
% site='Wohler';

switch site
    case 'Posidonius'
        addpath posidonius 9.4x11
        addpath ('Landing_sites\site_posidoniusBig_Lat32_Lon19_size9.4x11')
    case 'Wohler'
        addpath Wohler big 7x8.9
        addpath ('Landing_sites\site_WohlerBig_Lat-38_Lon29.5_size7x8.9');
end

if ~exist('TOPO','var')
    load('maps.mat','TOPO','Latitude');
end
 
z=TOPO; ppd=512;
%z=TOPOels; ppd=128;
 
R=1737;
lat0=Latitude;
kmpd=2*pi*R/360;
res=kmpd/ppd;
r=15;
th=linspace(0,2*pi);

x=(1:size(z,2))*(res*cosd(lat0));x=x-mean(x);
y=(1:size(z,1))*(res);y=y-mean(y);

% setting map to ellipse size & find delta topography
xl=find(x<-r,1,'last'); xf=find(x>r,1); 
yl=find(y<-r,1,'last'); yf=find(y>r,1); 
zellipse=z(yl:yf,xl:xf);
zout=z; zout(yl:yf,xl:xf)=nan;
delT=(max(zout(:))-min(zellipse(:)))/1e3;
clear zout zellipse

%find max distance
ea=25; %elevation angle lower bound
maxd=delT/tand(ea);
if maxd<0, maxd=0; end %no obstacles outside the ellipse

%setting map to relevant limits
clear xl xf yl yf
xl=find(x<-r-maxd,1,'last'); xf=find(x>r+maxd,1); xp=x(xl:xf);
yl=find(y<-r-maxd,1,'last'); yf=find(y>r+maxd,1); yp=y(yl:yf);
[yy,xx]=ndgrid(yp,xp);
zz=z(yl:yf,xl:xf);
d=sqrt(xx.^2+yy.^2)*1e3;

azbw=1;
azbins=-180:azbw:180;
f=find(d<r*1e3);
zp=zeros(size(zz));
% c=1;
% maxi=1;
% m=zeros(1,1); k=zeros(length(f),0);

profile on
for i=1:length(f)
    xxr=xx-xx(f(i)); yyr=yy-yy(f(i));
    azr=atan2d(yyr,xxr);  
    dr=sqrt(xxr.^2+yyr.^2)*1e3;
    hr=zz-zz(f(i))-dr.^2/(2.*R*1e3);
    slr=atand(hr./dr);
    fd=find(dr>3*res*1e3);
    idx=find(slr(fd)>ea);
    if idx>0
%         k(c)=i;      
%         if (length(idx)>maxi)
%              m(1:c-1,maxi+1:length(idx))=0;
%              maxi=length(idx);
%         end
%         m(c,:)=padarray(idx,maxi-length(idx),'post');
%         c=c+1;
        obsaz=unique(round(azr(idx))); % in what azimuths relative to the point are obstacles
        zp(f(i))=length(obsaz)/length(azbins);
    end
end
% profile viewer

% k(k==0)=[];

% setting map to ellipse size
clear xl xf yl yf
xl=find(xp<-r,1,'last'); xf=find(xp>r,1); xr=xp(xl:xf);
yl=find(yp<-r,1,'last'); yf=find(yp>r,1); yr=yp(yl:yf);
zr=zp(yl:yf,xl:xf);

[vv,uu]=ndgrid(yr,xr);
dd=sqrt(uu.^2+vv.^2); 
ff=find(dd>r);
zr(ff)=NaN;

figure;
% suptitle ({'','','',siteName,'Blocked Azimuths due to Communication Constraint'})
% subplot 121
b=imagesc(xr,yr,zr); axis tight ij equal
set(b,'AlphaData',~isnan(zr))
hold on;
colorbar;
title 'Rate';
xlabel 'km'; ylabel 'km';
set(gca,'xtick',-r:3:r);
set(gca,'ytick',-r:3:r);

% sortzr=sort(zr(:));
% fin=find(isnan(sortzr),1);
% sortzr=sortzr(1:(fin-1));
% pos2 = [0.6 0.2 0.3 0.58];
% subplot('Position',pos2)
% plot(sortzr,(1:length(sortzr))/length(sortzr)); 
% grid on
% xlabel 'Blocked Az/Total Az'; ylabel 'CDF';
% title 'Cumulative Distribution Function';
% axis tight

%convert saved data to azimuth, distance,center point coordinates & obstacle coordinates matrices
%example: for center point [xc(1),yc(1)] there are obstacles in points
%[xobs(1,:),yobs(1,:)]. those obstacles are in azimuths azobs(1,:) and
%distance (in metres) dobs(1,:).

% azobs=zeros(size(m)); dobs=zeros(size(m));
% xc=zeros(size(m,1))'; yc=zeros(size(m,1))';
% xobs=zeros(size(m)); yobs=zeros(size(m));
% col=length(m(1,:)); row=length(m(:,1));
% 
% if k>0
% for i=1:row
%     xc(i)=xx(f(k(i))); yc(i)=yy(f(k(i)));
%     for j=1:col
%         if m(i,j)>0
%             xobs(i,j)=xx(m(i,j));
%             yobs(i,j)=yy(m(i,j));
%             delx=xobs(i,j)-xc(i);
%             dely=yobs(i,j)-yc(i);
%             dobs(i,j)=sqrt(delx^2+dely^2)*1e3;
%             azobs(i,j)=atan2d(dely,delx);
%         else
%             break;
%         end
%     end
% end
% end
% 
% azobs=azobsc; dobs=dobsc; xc=xcc; yc=ycc; xobs=xobsc; yobs=yobsc; zr=zrc; xr=xrc; yr=yrc;
% fn='dataComm2';
% 
% switch site
%     case 'Posidonius'
% save([pwd,'\Line_of_sight\posidonius 9.4x11\'...
%     ,fn],'azobs','dobs','xc','yc','xobs','yobs','zr','yr','xr');       
%     case 'Wohler'
% save([pwd,'\Line_of_sight\wohler big 7x8.9\'...
%     ,fn],'azobs','dobs','xc','yc','xobs','yobs','zr','yr','xr');       
% end
