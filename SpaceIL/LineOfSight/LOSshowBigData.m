function LOSshowBigData(point)
% unplot(1);
% [xg,yg]=ginput(1);
xg=round(point(1),6);yg=round(point(2),6);
hold on; plot (xg,yg,'xr');
xg=num2str(round(xg,6));
yg=num2str(round(yg,6));
fn='C:\Users\Oagroup\Dropbox (Weizmann Institute)\spaceIL\updated\Line_of_sight\wohler big 7x8.9\dataThermal.txt';
f=fopen(fn,'r');
s1=strcat ('x=',num2str(xg,[2,6]),' y=',num2str(yg,[2,6]));
tf=0; nol1=0;

while ~tf
    s2=fgetl(f);
    if ~ischar(s2)
    error('reached end of file before pattern was found'); end
    tf=strcmp(s1,s2);
    nol1=nol1+1;
end

tf=0;
nol2=nol1-1; 

while ~tf
    s3=fgetl(f);
    tf=strcmp('***',s3);
    nol2=nol2+1;
end

M=dlmread(fn,'',[nol1 0 nol2-1 1]);
hold on; plot(M(:,1),M(:,2),'xb');
fclose(fn);
end
