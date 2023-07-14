function [dat_flt99,time_flt99,qv_flt99_fp,qv_flt99_humi,dist_flt99,X_flt99,Y_flt99]=read_BAS_Feb2010_flights_func(filename,ncol,X_latlon,Y_latlon,LATS,LONS)

%filename='/home/mbexddg5/work/manchester_flt99.txt';
%filename='Y:\BAS_flights\flight99_data\flight99_data.txt';

fid=fopen(filename,'rt');

dat_flt99 = fscanf(fid,'%f',[ncol inf]);

dat_flt99=dat_flt99'; %transpose to be as for other flights

% idat=0;
% go=1;
% while go==1
%     idat=idat+1;
%     dat_temp = fscanf(fid,'%f',15);
%     if length(dat_temp)==0
%         go=0;
%     else
%         dat_flt99(idat,1:length(dat_temp)) = dat_temp;
%     end
% 
% end


fclose(fid);



dat_flt99=sortrows(dat_flt99,1); %in case the times in the file are not in order - this puts them in order
inds=1:size(dat_flt99,1);

time_flt99=dat_flt99(:,1)/3.6e6; %convert from ms to hours

td_flt99=dat_flt99(inds,7)+273.15;  %frost point instrument ***wasn't working properly in Jan 2006***
td2_flt99=dat_flt99(inds,8)+273.15; %humicap data possibly more reliable?

f=1e6*28.97/18; 

qv_flt99_fp = SatVapPress(td_flt99,'goff','ice',100*dat_flt99(inds,6),1) / f;
qv_flt99_humi = SatVapPress(td2_flt99,'goff','ice',100*dat_flt99(inds,6),1) / f;

%calculate the distance along the flight track (total distance travelled)
%based on every 100 points throughout the data
inds=[1:100:length(dat_flt99(:,1))];
dist_flt99 = cumsum(distlatlon(dat_flt99(inds(1:end-1),2),dat_flt99(inds(1:end-1),3),dat_flt99(inds(2:end),2),dat_flt99(inds(2:end),3)),1);
%note that distlatlon function is a little out - need correcting for non-sphericity of the globe (12th Feb, 2010)

%lat_ref = -72.3142; %lat2d.var(1,1) for 3rd domain
%lon_ref = -71.9906; %lat2d.var(1,1) for 3rd domain

%find x and y positions by interpolating from the regular LAT LON grid tables for X and Y
%for domain 3 that have been previously computed

%so these are the positions from the origin of domain 3
X_flt99=interp2(LATS,LONS,X_latlon,dat_flt99(:,2),dat_flt99(:,3));
Y_flt99=interp2(LATS,LONS,Y_latlon,dat_flt99(:,2),dat_flt99(:,3));

disp('Done reading BAS Feb 2010 aircraft data');



%profiles from aircraft data
%descent starts at ~20.8 UTC and finishes at ~21 UTC

%[ibeg iend]=findheight(dat(:,1)/3.6e6,20.8,21.0);
%inds_prof=ibeg:iend;
%plot(qv2(inds_prof),dat(inds_prof,11));




