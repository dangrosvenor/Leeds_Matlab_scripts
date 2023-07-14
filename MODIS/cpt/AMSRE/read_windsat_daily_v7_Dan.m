function [time,sst,windLF,windMF,vapor,cloud,rain,windAW,wdir]=read_windsat_daily_v7(data_file)
% [time,sst,windLF,windMF,vapor,cloud,rain,windAW,wdir]=read_windsat_daily_v7(data_file);
% 
% this subroutine will read the compressed or uncompressed RSS WindSat daily bytemaps.
% reads version-7.0.1 files released August 2013
%
% input argument is the full path file name:
%   read_windsat_daily_v7(filename)
%
% output arguments:
%   [mingmt,sst,windLF,windMF,vapor,cloud,rain,windAW,wdir]
%   mingmt is gmt time in minutes of day
%   sst  is surface water temperature at depth of about 1 mm in deg C
%   windLF is 10 meter surface wind in m/s made using 10.7 GHz channel and above, low frequency channels
%   windMF is 10 meter surface wind in m/s made using 18.7 GHz channel and above, medium frequency channels
%   vapor is atmospheric water vapor in millimeters
%   cloud is liquid cloud water in millimeters
%   rain  is rain rate in millimeters/hour
%   windAW is 10 meter surface wind for all weather conditions made using 3 algorithms
%   wdir is wind direction oceanographic convention, blowing North = 0 in degrees
%
%  The center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude.
%  The center of the second cell is 0.375 E longitude, -89.875 latitude.
% 		XLAT=0.25*ILAT-90.125
%		XLON=0.25*ILON-0.125
%
% Details of the binary data file format are located at 
% http://www.remss.com/windsat/windsat_data_description.html#binary_data_files
%
% To contact RSS support:
% http://www.remss.com/support
%
%

xscale=[6., 0.15,0.2,0.2,0.3, 0.01,0.1,0.2,1.5];
offset=[0.,-3.0 ,0.0,0.0,0.0,-0.05,0.0,0.0,0.0];
%note the negative offset for cloud (lwp, = -0.05) - is negative because they take
%away the water vapor from the total water column. They used ot set these
%to zero, but is generally agreed now that this biases the LWP high. (from
%Rob).



xscale=[6., 0.15,0.2,0.2,0.3, 0.01,0.1,0.2,1.5];
offset=[0.,-3.0 ,0.0,0.0,0.0,-0.05,0.0,0.0,0.0];


xdim=1440;ydim=720;tdim=2;numvar=9;

if ~exist(data_file,'file'),
    disp(['file not found: ', data_file]);
    time=[];sst=[];windLF=[];windMF=[];vapor=[];cloud=[];rain=[];windAW=[];wdir=[];
    return;
end;

if ~isempty(regexp(data_file,'.gz', 'once'))
    data_file=char(gunzip(data_file));
end

fid=fopen(data_file,'rb');
map=reshape(fread(fid,xdim*ydim*numvar*tdim,'uint8'),xdim,ydim,numvar,tdim);
fclose(fid);
disp(data_file);

for i=1:numvar
    tmp=map(:,:,i,:);
    %modified by Dan to actually remove missing data (>250)
%    ia=find(tmp<=250);tmp(ia)=tmp(ia)*xscale(i)+offset(i);
%set missing data to NaN
     ib=find(tmp>250); tmp(ib)=NaN;
     %now scale all data (NaNs will remain NaN)    
     tmp=tmp*xscale(i)+offset(i);     
%    ia=find(tmp<=250);tmp(ia)=tmp(ia)*xscale(i)+offset(i);
    map(:,:,i,:)=tmp;
end;

time=squeeze(map(:,:,1,:)) / 60; %Dividing by 60 to convert from mins to hours (Windsat uses minutes rather than hours as used by
    % AMSRE, etc.)
sst=squeeze(map(:,:,2,:));
windLF=squeeze(map(:,:,3,:));
windMF=squeeze(map(:,:,4,:));
vapor=squeeze(map(:,:,5,:));
cloud=squeeze(map(:,:,6,:));
rain=squeeze(map(:,:,7,:));
windAW=squeeze(map(:,:,8,:));
wdir=squeeze(map(:,:,9,:));

return;
end

