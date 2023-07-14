function [sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_averaged_v7(data_file)
% [sst,wspdLF,wspdMF,vapor,cloud,rain]=read_amsr_averaged_v7(data_file);
%
% this subroutine reads uncompressed AMSR-E data from Remote Sensing Systems'
% binary format  (version-7) 
%
% File name format is  	amsre_yyyymmddv7_d3d for 3-day   (average of 3 days ending on file date)
%			      amsre_yyyymmddv7	   for weekly  (start sunday, end saturday, named by saturday date)
%			      amsre_yyyymmv7	   for monthly 

% The averaged files include: 3-Day, weekly, and monthly time composites.
% These averaged time composite files all share the same data format.
%	
% input arguments:
% data_file = the full path and name of the uncompressed data file
%
% the function returns these products:
%   [sst,wspdLF,wspdMF,vapor,cloud,rain,wspdAW,wdir]
%   sst  is surface water temperature at depth of about 1 mm in deg C
%   wspdLF is 10 meter surface wind in m/s made using 10.7 GHz channel and above, low frequency channels
%   wspdMF is 10 meter surface wind in m/s made using 18.7 GHz channel and above, medium frequency channels
%   vapor is atmospheric water vapor in millimeters
%   cloud is liquid cloud water in millimeters
%   rain  is rain rate in millimeters/hour
%
%  The center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude.
%  The center of the second cell is 0.375 E longitude, -89.875 latitude.
% 		XLAT=0.25*ILAT-90.125
%		XLON=0.25*ILON-0.125
%
% For detailed data description, see 
% http://www.remss.com/amsr/amsr_data_description.html
%
% Remote Sensing Systems
% support@remss.com

xscale=[0.15,0.2,0.2,0.3, 0.01,0.1];
offset=[-3.0 ,0.0,0.0,0.0,-0.05,0.0];
xdim=1440;ydim=720;numvar=6;
mapsiz=xdim*ydim;

if ~exist(data_file,'file')
   disp(['file not found: ' data_file]);
   sst=[];wspdLF=[];wspdMF=[];vapor=[];cloud=[];rain=[];
   return;
end;

fid=fopen(data_file,'rb');
data=fread(fid,mapsiz*numvar,'uint8');
fclose(fid);
disp(data_file);
map=reshape(data,[xdim ydim numvar]);

for i=1:numvar
    tmp=map(:,:,i);
    ia=find(tmp<=250);tmp(ia)=tmp(ia)*xscale(i)+offset(i);
    map(:,:,i)=tmp;
end;

sst=squeeze(map(:,:,1));
wspdLF=squeeze(map(:,:,2));
wspdMF=squeeze(map(:,:,3));
vapor=squeeze(map(:,:,4));
cloud=squeeze(map(:,:,5));
rain=squeeze(map(:,:,6));

return;
