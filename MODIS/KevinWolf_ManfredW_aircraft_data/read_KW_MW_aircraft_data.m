
txtfile = '/home/disk/eos10/d.grosvenor/KevinWolf_Manfred_aircraft_data/cloud_retrieval_20160819a_UTC_19.23-19.36_uncert03_reff15.dat'

ncols=14;

fid=fopen(txtfile,'rt');
textstr=repmat('%s',[1 ncols]);
headings = textscan(fid,textstr,1);
textstr=repmat('%f',[1 ncols]);
dat=textscan(fid,textstr);
fclose(fid);

fprintf(1,'\n Done file read.\n');

base_date = datenum('19-Aug-2016');
mat_time = datenum(base_date + dat{1}(:)/3600/24);
% >> datestr(mat_time(1))
% ans =
% 19-Aug-2016 19:24:00
% >> datestr(mat_time(end))
% ans =
% 19-Aug-2016 19:36:00

%So, only 12 mins of data here...

%Aircraft was heading NE on a straight line for this segment.
%Starts at (lat,lon) :- (16.6058, -52.3753) ends at (15.9593,-53.7544)
%So travels 
% distlatlon(dat{3}(1),dat{2}(1),dat{3}(end),dat{2}(end))
% 
% ans =
% 
%   164.0020 km



   



