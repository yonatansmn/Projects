clear
site='Posidonius';
% site='Wohler';

switch site
    case 'Posidonius'
%          load('Line_of_sight\posidonius 9.4x11\dataComm'); thermal=0;
         load('Line_of_sight\posidonius 9.4x11\dataThermal'); thermal=1;
         load('Line_of_sight\posidonius 9.4x11\maps.mat','TOPO','Latitude');
    case 'Wohler'
%          load('Line_of_sight\wohler big 7x8.9\dataComm'); thermal=0;
%          load('Line_of_sight\wohler big 7x8.9\dataThermal'); thermal=1;
         load('Line_of_sight\wohler big 7x8.9\maps.mat','TOPO','Latitude');
end

z=TOPO; ppd=512;
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
ea=25; 
as=90;

if thermal
ea=3;
as=30;
end

maxd=delT/tand(ea);
if maxd<0, maxd=0; end %no obstacles outside the ellipse

%setting map to relevant limits
clear xl xf yl yf
xl=find(x<-r-maxd,1,'last'); xf=find(x>r+maxd,1); xp=x(xl:xf);
if thermal, xl=find(x<-r,1,'last'); xf=find(x>r+maxd,1); xp=x(xl:xf); end
yl=find(y<-r-maxd*sind(as),1,'last'); yf=find(y>r+maxd*sind(as),1); yp=y(yl:yf);
zz=z(yl:yf,xl:xf);

% display
if thermal, suptitle ({site,'Blocked Azimuths due to Thermal Constraint'});
else, suptitle ({site,'Blocked Azimuths due to Communication Constraint'}); end

f1=subplot (221);
imagesc(f1,xp,yp,zz); axis tight ij equal
title({'SLDEM Topography (512 ppd)'});
set(get(colorbar,'ylabel'),'string','Topography (m)')
xlabel 'km'; ylabel 'km';
hold on
plot (r.*cos(th),r.*sin(th),'k');

b=subplot (223);
bb=imagesc(b,xr,yr,zr); axis tight ij equal
title 'Rate of Blocked Azimuths';
set(bb,'AlphaData',~isnan(zr))
colorbar;
xlabel 'km'; ylabel 'km';
set(gca,'xtick',-r:5:r);
set(gca,'ytick',-r:5:r);
% linkaxes ([b,f1]);

titles={'Show','X','Y','#Blocked Az'};
X=xc(:,1); Y=yc(:,1);
numobs=zeros(1,size(xobs,1));
for i=1:length(numobs)
    numobs(i)=find(xobs(i,:)~=0,1,'last');
end

co=zeros(size(xobs,1),size(xobs,2)*2);
for j=1:size(xobs,2)
    co(:,2*j-1)=xobs(:,j);
    co(:,2*j)=yobs(:,j);
end

[numobs, or] = sort(numobs,'descend');
X=X(or,:);
Y=Y(or,:);
co=co(or,:);

data=[X Y co];

tdata=[X Y numobs'];
tdata(tdata==0)=nan;

for i1=1:size(tdata,1)
        for i2=1:size(tdata,2)
            D{i1,i2}=num2str(tdata(i1,i2),3);
        end
end
c=cell(size(tdata,1),1);c(:)={false};
uidata=[c,D];

t=uitable(figure(1),'data',uidata,'columnname',titles,'fontsize',15,'ColumnWidth',{100});
subplot(2,2,[2,4]);
ax = gca;
ax.Visible = 'off';
% pos = get(subplot(2,2,[2,4]),'position');
set(t,'units','normalized');
set(t,'position',[0.5703 0.1100 0.3447 0.7800]);
t.ColumnEditable(1)=1;

figure(3);p1=plot(1:3);close(3)
setappdata(f1.Parent,'handle',p1);
t.CellEditCallback={@(src,eventdata)callback1(src,eventdata,f1,data,thermal,r,maxd)};

%%
% [xg,yg]=ginput(1);
% xg=num2str(round(xg,6));
% yg=num2str(round(yg,6));
% fn=fopen('dataThermal','r');
% s1=strcat ('x=',xg,' y=',yg);
% tf=0; nol1=0;
% 
% while tf==0
%     s2=fgetl(fn);
%     tf=strcmp(s1,s2);
%     nol1=nol1+1;
% end
% 
% tf=0;
% nol2=nol1; 
% 
% while tf==0
%     nol2=nol2+1;
%     s3=fgetl(fn);
%     tf=strcmp('***',s3);
% end
% 
% M=dlmread(fn,'',[nol1 0 nol2 0]);
