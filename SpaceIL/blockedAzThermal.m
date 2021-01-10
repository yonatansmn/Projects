% According to section 2.1.3 in "Landing Site Requirements" document, each landing site in the landing zone shall allow a 
% non-obstructed line of sight of 3° elevation with angular sector of ±30° toward East (azimuth 0°) to prevent shading of solar panels.
% 
% this code calculate the precentage of blocked azimuths at every point within the landing ellipse.

clear
% siteName='Posidonius';
siteName='Wohler';

switch siteName
    case 'Posidonius'
         load('Landing_sites\site_posidoniusBig_Lat32_Lon19_size9.4x11\maps.mat','TOPO','Latitude');
    case 'Wohler'
        load('Landing_sites\site_WohlerBig_Lat-38_Lon29.5_size7x8.9\maps.mat','TOPO','Latitude');
end

z=TOPO;
ppd=512; %pixels per degree
lat0=Latitude;
%z=TOPOels; ppd=128;

%General params
Rmoon=1737;
kmpd=2*pi*Rmoon/360; %km per degree
res=kmpd/ppd; %km/pixels ratio
r=15; %radius [km]
th=linspace(0,2*pi);
ealb=[5]; %elevation angle lower bound (can do more than 1)
as=30; %angular sector
% X=(29.5-9148/2/ppd)+(1:size(TOPO,2))/ppd; Y=(-38-7208/2/ppd)+(1:size(TOPO,1))/ppd; figure; surf (X,Y,TOPO); shading interp; xlabel 'EAST'; ylabel ' NORTH';

% setting map to ellipse size & find (max) delta topography

% conver to km
x=(1:size(z,2))*(res*cosd(lat0)); %cosine is needed due to the surf curveness
y=(1:size(z,1))*(res);

% center around 0
x=x-mean(x); 
y=y-mean(y);

% shrink to ellipse size
xf=find(x>r,1,'first'); xl=find(x<-r,1,'last');
yf=find(y>r,1,'first'); yl=find(y<-r,1,'last');

z=z(:,xl:end); %west to the landing site is irrelevant (see explanation at the top)
zellipse=z(yl:yf,xl:xf);
zout=z; 
zout(yl:yf,xl:xf)=nan;
max_dz=(max(zout(:))-min(zellipse(:)))/1e3;
clear zout zellipse


for j=1:length(ealb)
    ea=ealb(j);

    %find max distance
    maxd=max_dz/tand(ea);
    if maxd<0, maxd=0; end %no obstacles outside the ellipse

    % setting map to relevant limits
    x_east=x(xl:end); 
    clear xl xf yf yl
    xf=find(x_east>r+maxd,1,'first'); 
    xp=x_east(1:xf);
    yf=find(y>r+maxd*sind(as),1,'first'); 
    yl=find(y<-r-maxd*sind(as),1,'last'); 
    yp=y(yl:yf);
    zz=z(yl:yf,1:xf);
    [yy,xx]=ndgrid(yp,xp);
    d=sqrt(yy.^2+xx.^2)*1e3;

    azbw=1;
    azbins=-180:azbw:180;
    f=find(d<r*1e3);
    zp=zeros(size(zz));

    % c=1;
    % maxi=1;
    % m=zeros(length(f),1);
    % k=zeros(length(f),1);
    % fn='dataThermal.txt';
    % fileID=fopen(fn,'w');
    
    %convert to meter
    Dmoon = 2*Rmoon*1e3;
    res_m = res * 1e3;
    xx = xx * 1e3;
    yy = yy * 1e3;
    
    % profile on
    [fi,fj]=ind2sub(size(d),f);
    length_f = length(f);
    for i=1:length_f
        progressbar(i,length_f);
        if (mod(fi(i),3)>0 || mod(fj(i),3)>0), continue; end
        xxr=xx-xx(f(i)); yyr=yy-yy(f(i));
        azr=atan2d(yyr,xxr);  
        dr=sqrt(xxr.^2+yyr.^2);
        hr=zz-zz(f(i))-dr.^2/Dmoon;
        slr=atand(hr./dr);
        fd=find(dr>3*res_m);
        idx=find(slr(fd)>ea);
        slraz=azr(idx); % in what azimuths relative to the point are obstacles
        idx=idx(abs(slraz)<=as);
        if idx>0
    %         k(c)=i;
    %             if (length(idx)>maxi)
    %                  m(1:c-1,maxi+1:length(idx))=0;
    %                  maxi=length(idx);
    %             end
    %         m(c,:)=padarray(idx,maxi-length(idx),'post');
    %         c=c+1;
            slraz=unique(round(slraz));
            obsaz=find(abs(slraz)<=as);
            zp(f(i))=length(obsaz)/(ceil(length(azbins)/6));
        end
    end
    % profile viewer

    % k(k==0)=[];
    % m(m(:,1)==0,:)=[];

    % setting map to ellipse size
    clear xf yf yl
    xl=find(xp<-r,1,'last'); xf=find(xp>r,1); xr=xp(xl:xf);
    yl=find(yp<-r,1,'last'); yf=find(yp>r,1); yr=yp(yl:yf);
    zr=zp(yl:yf,xl:xf);

    [vv,uu]=ndgrid(yr,xr);
    dd=sqrt(uu.^2+vv.^2); 
    ff=find(dd>r);
    zr(ff)=NaN;

    figure(j);
    subplot 121
    suptitle ({'','','',siteName,'Blocked Azimuths due to Thermal Constraint'})
    b=imagesc(xr,yr,zr); axis tight ij equal
    set(b,'AlphaData',~isnan(zr))
    hold on;
    colorbar;
    title 'Rate';
    xlabel 'km'; ylabel 'km';
    set(gca,'xtick',-r:3:r);
    set(gca,'ytick',-r:3:r);

    sortzr=sort(zr(:));
    fin=find(isnan(sortzr),1);
    sortzr=sortzr(1:(fin-1));
    pos2 = [0.6 0.2 0.3 0.58];
    subplot('Position',pos2)
    plot(sortzr,(1:length(sortzr))/length(sortzr)); 
    grid on
    xlabel 'Blocked Az/Total Az'; ylabel 'CDF';
    title 'Cumulative Distribution Function';
    axis tight

end
%convert saved data to azimuth, distance,center point coordinates & obstacle coordinates matrices
%example: for center point [xc(1),yc(1)] there are obstacles in points
%[xobs(1,:),yobs(1,:)]. those obstacles are in azimuths azobs(1,:) and distance
%dobs(1,:).
% 
% if k>0
% azobs=zeros(size(m)); dobs=zeros(size(m));
% xc=zeros(size(m,1))'; yc=zeros(size(m,1))';
% xobs=zeros(size(m)); yobs=zeros(size(m));
% col=length(m(1,:)); row=length(m(:,1));
% 
% for i=1:row
%       xc(i)=xx(f(k(i))); yc(i)=yy(f(k(i)));
%         for j=1:col
%             if m(i,j)>0
%                 xobs(i,j)=xx(m(i,j));
%                 yobs(i,j)=yy(m(i,j));
%                 delx=xobs(i,j)-xc(i);
%                 dely=yobs(i,j)-yc(i);
%                 dobs(i,j)=sqrt(delx^2+dely^2)*1e3;
%                 azobs(i,j)=atan2d(dely,delx);
%             else
%                 break;
%             end
%         end
%     end
% end

% azobs=azobst; dobs=dobst; xc=xct; yc=yct; xobs=xobst; yobs=yobst; zr=zrt; xr=xrt; yr=yrt;
% azobs=azobsc; dobs=dobsc; xc=xcc; yc=ycc; xobs=xobsc; yobs=yobsc; zr=zrc; xr=xrc; yr=yrc;
% 
% fn='dataThermal4';
% switch site
%     case 'Posidonius'
% save([pwd,'\Line_of_sight\posidonius 9.4x11\'...
%     ,fn],'azobs','dobs','xc','yc','xobs','yobs','zr','yr','xr');       
%     case 'Wohler'
% save([pwd,'\Line_of_sight\wohler big 7x8.9\'...
%     ,fn],'azobs','dobs','xc','yc','xobs','yobs','zr','yr','xr');
% end
