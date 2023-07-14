%indices for the first 3000m leg
i3000 = [62921:238242];

detail_level = 'Whole track, every 100 indices';
detail_level = 'Whole track, every 10000 indices';
%detail_level = 'Every 0.2 degrees on eastern part of track away from mountain';
detail_level = 'Every 0.2 degrees on eastern part of track away from mountain - only bit covered by outbound leg too';
%detail_level = 'Single point';

%make some coarse indices for faster calcs
switch detail_level
    case 'Single point'
        LONS_choose = [-65.4];
        iLONS = findheight_nearest(dat_flt19(i3000,3),LONS_choose);
        iLONS(iLONS>1e19)=length(i3000);
        i3000_2 = i3000(iLONS);
        recalc_ilat=1;
        ilat=0; ilon=0;
    case 'Whole track, every 100 indices'
        i3000_2=[i3000(1):100:i3000(end)];
        %Load in previously calculated ilat and ilon indices for the WRF
        %array that correspond to the flight track (LATS and LONS above)
        load_file_Ant='~/mat_files_various/Antarctica/ilat_ilon_3000m.mat';
        load(load_file_Ant,'ilat','ilon');
        recalc_ilat=0;
    case 'Every 0.2 degrees on eastern part of track away from mountain'
        %east of 65.4 deg lon
        %"case B"
        LONS_choose = [-65.4:0.2:-62.2];
        iLONS = findheight_nearest(dat_flt19(i3000,3),LONS_choose);
        iLONS(iLONS>1e19)=length(i3000);
        i3000_2 = i3000(iLONS);
%        recalc_ilat=1;
%        ilat=0; ilon=0;
        load_file_Ant='~/mat_files_various/Antarctica/Pressure_3000m_WRF_caseB.mat';
        load(load_file_Ant,'ilat','ilon');
        recalc_ilat=0;
    case 'Every 0.2 degrees on eastern part of track away from mountain - only bit covered by outbound leg too'
        %east of 65.4 deg lon
        %"case C"
        LONS_choose = [-63.7:0.2:-62.2];
        iLONS = findheight_nearest(dat_flt19(i3000,3),LONS_choose);
        iLONS(iLONS>1e19)=length(i3000);
        i3000_2 = i3000(iLONS);
%        recalc_ilat=1;
%        ilat=0; ilon=0;    
        load_file_Ant='~/mat_files_various/Antarctica/Pressure_3000m_WRF_caseC.mat';
        load(load_file_Ant,'ilat','ilon');
        recalc_ilat=0;
        
      case 'Whole track, every 10000 indices'
        i3000_2=[i3000(1):1e4:i3000(end)];
        %Load in previously calculated ilat and ilon indices for the WRF
        %array that correspond to the flight track (LATS and LONS above)
        recalc_ilat=1;
        ilat=0; ilon=0;    

        %load_file_Ant='~/mat_files_various/Antarctica/ilat_ilon_3000m.mat';
        %load(load_file_Ant,'ilat','ilon');
        %recalc_ilat=0;   
        
        
end



LATS = dat_flt19(i3000_2,2);
LONS = dat_flt19(i3000_2,3);



field_str = 'Pressure';
%field_str = 'Wind direction';


%zfind = 3000; %(m)
zfind = 2954.2;

%Time for WRF
        %6th Jan:- 11=6UTC, 13=12UTC, 15=18UTC, 16=21 UTC, 17=00UTC
           %Aircraft took off at 19:20 and started descent at 20:15 (flew
           %at 3000m). So, midpoint of the upper level run was around
           %19:45. Closest to 21 UTC model output.
           
%times=8:24;
times=1:25;

p3000_time = NaN*ones([length(LATS),length(times)]);           
itime=0;
for time=times  
    itime=itime+1
    [p3000_time(:,itime),HGT,iz3000,ilat,ilon] = get_field_along_flight_path(nc,lat2d,lon2d,time,LATS,LONS,zfind,field_str,recalc_ilat,ilat,ilon,DX,DY);
end

times_all = [0:3:72]; %hours after 0UTC on 6th Jan
time_utc_p3000 = times_all(times);

% average over the part away from the mountains to the east
% iav=find(LONS>-65.4);
% p3000_timser=meanNoNan(p3000_time(iav,:),1);
% figure
% plot(time_utc_p3000,p3000_timser,'bo-');


