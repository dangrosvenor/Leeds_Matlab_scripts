function [time,wind,vapor,cloud,rain]=read_ssmi_day_v7(data_file)
% [time,wind,vapor,cloud,rain]=read_ssmi_day_v7(data_file);
%
% this subroutine reads compressed or uncompressed RSS SSM/I and SSMIS daily byte maps
% (version-7 data released June 2012 for SSM/I and August 2011 for SSMIS)
%
% input arguments:
% data_file = the full path and name of the uncompressed data file
%
% the function returns these products:
% time = time of observation in fractional hours GMT
% wind = wind speed in meters/second
% vapor = atmospheric water vapor in millimeters
% cloud = liquid cloud water in millimeters
% rain  = rain rate in millimeters/hour
%
%  longitude is 0.25*xdim- 0.125
%  latitude  is 0.25*ydim-90.125
%
% For detailed data description, see 
% http://www.remss.com/ssmi/ssmi_description.html
%
% Remote Sensing Systems
% support@remss.com

scale=[.1,.2,.3,.01,.1];
offset=[0.,0.,0.,-0.05,0.];

xdim=1440;ydim=720;tdim=2;numvar=5;
mapsiz=xdim*ydim*tdim;

if ~exist(data_file,'file')
   disp(['file not found: ' data_file]);
   time=[];wind=[];vapor=[];cloud=[];rain=[];
   return;
end;

if ~isempty(regexp(data_file,'.gz', 'once'))
    data_file=char(gunzip(data_file));
end

fid=fopen(data_file,'rb');
data=fread(fid,mapsiz*numvar,'uint8');
fclose(fid);
disp(data_file);
map=reshape(data,[xdim ydim numvar tdim]);

for i=1:numvar
    tmp=map(:,:,i,:);
    ia=find(tmp<=250);tmp(ia)=tmp(ia)*scale(i)+offset(i);
    map(:,:,i,:)=tmp;
end;

time  = squeeze(map(:,:,1,:));
wind  = squeeze(map(:,:,2,:));
vapor = squeeze(map(:,:,3,:));
cloud = squeeze(map(:,:,4,:));
rain  = squeeze(map(:,:,5,:));

return;
end
