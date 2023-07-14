clear dire

dire(1).d='c:/cygwin/home/user/lesexe/les2.3/tracer06-03-04/debugtwomey/dan';
%dire(2).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_06_18AED.txt';
%dire(3).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_03_12AED.txt';
%dire(4).d='c:\documents and settings\user\my documents\HIBISCUS\baurufield\soundings\ASCII\1\04_02_02_12AED.txt';


%for i=1:size(dire,2)

n=149

fid=fopen(dire(1).d,'rb');
fseek(fid,4,'bof');

Pdat=fread(fid,n,'float64=>double'); %pressure (Pa)
fseek(fid,8,'cof');

Tdat=fread(fid,n,'float64=>double');
fseek(fid,8,'cof');

Wdat=fread(fid,n,'float64=>double'); %vapour mixing ratios
fseek(fid,8,'cof');

ALTdat=fread(fid,n,'float64=>double');
fseek(fid,8,'cof');

%capeles=fread(fid,1,'float64=>double');
%fseek(fid,8,'cof');

%cinles=fread(fid,1,'float64=>double');

fclose(fid);


%Tdat=THdat./(1e5./Pdat).^0.286;

%end

%rhstart=dat(6,1); %relative humidity at ground (%)
