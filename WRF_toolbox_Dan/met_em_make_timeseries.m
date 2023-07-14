%%% DON'T FORGET TO COMMENT OUT ALL THE FILE DIRECTORY SPECIFICATIONS IN LOAD_WRF_VARS.M %%%

%%%%%%%%% setting up of the location of the points for timeseries %%%%%%%%%
ih_wrf=0; %set to zero if just want to load surface values (just psfc for now)


rundir='ant_jan06_ecmwf_Nov08'; 
rundir='ant_jan06_sfRUC_v3'; 
%rundir='ant_jan06_ecmwf_ml_0.225';
%rundir='ecmwf_ml_01_18';
rundir='ncepR2_3dom_nudging';
rundir='ecmwf_ml_0.5';

cd(['Y:/WRF/' rundir]);
files = dir('met_em*d01*');

ifile=1;
fileWRF(3).file=files(ifile).name;
load_WRF_vars;

%runs through all the AWSs stored in LAT_AWS and finds the timeseries of temp, pressure, wind speed and direction at each location

 
%----   points to plot on the flight track based on the time along the flight track
 times_flight_loc = [19.84565 20.056 20.12 20.2755]; %20.12 is the time where the max was seen for aircraft data - others are just points
 %along the aircraft track
 times_flight_loc = [20.3755]; %20.3755 is the time where the max was seen for aircraft data
 clear it_flt;
 for iflt=1:length(times_flight_loc)
     it_flt(iflt) = findheight(time_flt19,times_flight_loc(iflt));
 end
% LAT_extra = dat(it_flt,2);
% LON_extra = dat(it_flt,3);
 
%LAT_extra =[];
%LON_extra =[];


%----   extra points to plot - give as x and y km
% extra_x = [575]; %for 12UTC, 6th Jan, ncep polar
% extra_y = [351];
 
% extra_x = [556]; %for 03UTC, 7th Jan, ncep polar
% extra_y = [355];
 
% extra_x = [570]; %for 12UTC, 6th Jan, ecmwf
% extra_y = [355];

%  extra_x = [600]; %random 
%  extra_y = [400];
  
%  extra_x = [275]; %lefthand side of the equiv cross sections
%  extra_y = [380];

extra_x = []; %lefthand side of the equiv cross sections
extra_y = [];
 
 
 nlat = size(lat2d.var,1);
 nlon = size(lat2d.var,2);
 dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
 dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));
 
 i_grid = dx_grid * [1:nlon];
 j_grid = dy_grid * [1:nlat];
 
 for iflt=1:length(extra_x)
     i_extra = findheight(i_grid,extra_x(iflt));
     j_extra = findheight(j_grid,extra_y(iflt));
     LAT_extra = [LAT_extra ;lat2d.var(j_extra,i_extra)];
     LON_extra = [LON_extra ;lon2d.var(j_extra,i_extra)];
 end
 
 % ascent_str='034'; %'4' is the L-shaped segments - variation of L-shaped segs is minimal
% ascent_str='03'; %0 and 3 are the first and last ascent
% ascent_str='1'; %upwind side
 ascent_str='AWS'; %AWS locations
 
 if strfind(ascent_str,'0')  %descent after going over the peninsula
     LAT=[LAT_extra' ];% lat2d(1).var(140,240)];  %places during descent plus other locations
     LON=[LON_extra'  ];% lon2d(1).var(140,240)];
     as_ds_str='descent';
 elseif strfind(ascent_str,'1')
%     LAT=[-67.55 -67.62 -67.55]; %first ascent from Rothera
%     LON=[-68.1 -67.8 -67.5];
     LAT=[-67.55 -67.62 -67.55 -66.8 -67.2]; %first ascent from Rothera
     LON=[-68.1 -67.8 -67.5 -73.9 -70.9];
     
     as_ds_str='ascent';
 elseif strfind(ascent_str,'3')  %final ascent on the way back
     LAT=[LAT_extra  ];% lat2d(1).var(140,240)];  %places during descent plus other locations
     LON=[LON_extra  ];% lon2d(1).var(140,240)];
     as_ds_str='ascent2';
 elseif strfind(ascent_str,'4')  %cumulative ascent during L-shaped flight legs
     LAT=[LAT_extra  ];% lat2d(1).var(140,240)];  %places during descent plus other locations
     LON=[LON_extra  ];% lon2d(1).var(140,240)];
     as_ds_str='ascent2';
 elseif strfind(ascent_str,'AWS')  %cumulative ascent during L-shaped flight legs
%      if wis_aws==0
%          LAT=[-68.34]; %
%          LON=[-69.01];
%      else
         for itemp=1:length(LAT_AWS)
             LAT(itemp) = LAT_AWS(itemp).dat;
             LON(itemp) = LON_AWS(itemp).dat;
         end
%         LAT(6) = -72.4053;
%         LON(6) = -61.0326;
%      end
 end


 



    
    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);
    
%%%%% end of setting of location %%%%%   



clear pressure temperature wind vapour wind
for ifile=1:length(files)
    fileWRF(3).file=files(ifile).name;
    load_WRF_vars;
    pressure(1).times(ifile) = str2num(Times(9:10)) + str2num(Times(12:13))/24; %save the times

    for iloc=1:length(ilat)
        if isnan(ilat(iloc))==1 | isnan(ilon(iloc))==1
            pressure(iloc).y(ifile)=NaN;
            temperature(iloc).y(ifile)=NaN;
            wind(iloc).y(ifile)=NaN;
            winddir(iloc).y(ifile)=NaN;
        else
            if ih_wrf==0
                %            psfc = nc{'PSFC'}(:);
                %            psfc=psfc(ilat(iloc),ilon(iloc));

                pres=nc{'PRES'}(1,:,ilat(iloc),ilon(iloc));
                temp=nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                psfc=pres(1); %surface pressure
                
                pmsl=nc{'PMSL'}(:);  %using this method we don't need SOILHGT data as use z=0 from sea level - but need sea-level temperature
                    %soilhgt method only useful when the height we want the pressure for is above SOILHGT - usually not the case
                pmsl=pmsl(ilat(iloc),ilon(iloc));
                    

                if prod(size(nc{'GHT'})) > 0      %if is NCEP analysis

                    %NOTE, might want to use this for finding pressure corrections for NCEP analysis files
                    soilhgt = nc{'GHT'}(1,:,ilat(iloc),ilon(iloc));
                    [Y,I]=sort(pres,'descend'); %sort the data in order to take into account pres(2)(=surface) being lower than first level of 1000 hPa
                    pres2=pres(I);
                    temp2=temp(I);
                    soilhgt = soilhgt(I(1));  %height of the lowest pressure level
                else                    
                    pres2 = [pmsl pres];
                    skinT = temp(1)*(pmsl/pres(1))^0.286; %temperature assuming a constant potential temperature between msl and the surface point
                    %skinT = nc{'SKINTEMP'}(:); %using skin temperature as sea level temp - could be inaccurate?
                    %skinT = skinT(ilat(iloc),ilon(iloc));
                    %    skinT=skinT+5; %to test temperature sensitivity - does not seem that senstitive - 5 degree increase led to only 0.14 hPa increase in
                    %in estimated pressure of a height below soilhgt (tried 91 m wheras soilhgt was 111m - difference was even less for lower altitudes).
                    temp2 = [skinT temp];
                    %                soilhgt=nc{'SOILHGT'}(:); soilhgt=soilhgt(ilat(iloc),ilon(iloc));  %height of PSFC level
                    soilhgt=0;
                end



                PSPAN=[pres2(1) pres2(end)]; %pressure range for integration
                [P,h] = ODE45(@hydrostatic2,PSPAN,soilhgt,[],pres2,temp2); %enter the initial height for the given the first value of PSPAN
                %note that some of the pressure levels in the PRES array will be below the

                pressure(iloc).y(ifile) = interp1(h,P,ALT_AWS(iloc).dat); %interpolate pressure curve for the height of the station

                temperature(iloc).y(ifile) = interp1(pres2,temp2,pressure(iloc).y(ifile));  %issue is that here are interpolating based on the constant potemp
                %assumption below the terrain down to sea level

                u=nc{'UU'}(1,:,ilat(iloc),ilon(iloc));
                v=nc{'VV'}(1,:,ilat(iloc),ilon(iloc));
                sp = sqrt(u.^2 + v.^2);

                if psfc>pmsl %if analysis surface is below sea level (i.e. negative GHT)
                    PSPAN2=[pres(1) pres(end)]; %pressure range for integration - just use the orginal pressre and temp
                    [P2,h2] = ODE45(@hydrostatic2,PSPAN2,0,[],pres,temp); %set initial height as zero (psfc) as just want this +5m
                    H2=interp1(P2,h2,pres); %find heights above surface for the pressures
                    HSFC=0; %setting to zero - actual height doesn't matter as are just interested in 5 m above surface
                    wind(iloc).y(ifile) = interp1(H2,sp,HSFC+5); %interpolating for 5 m above the surface (as recognised by the analysis)
                    u5=interp1(H2,u,HSFC+5);
                    v5=interp1(H2,v,HSFC+5);
                    winddir(iloc).y(ifile) = wind_dir(u5,v5,lat2d,lon2d,ilat,ilon,iloc);
                    %or could take an average over lowest 10 m perhaps?
                else
                    H2=interp1(P,h,pres);
                    HSFC = interp1(P,h,psfc); %height of the surface (analysis surface) according to the hydrostatic array
                    wind(iloc).y(ifile) = interp1(H2,sp,HSFC+5); %interpolating for 5 m above the surface (as recognised by the analysis)                    
                    u5=interp1(H2,u,HSFC+5);                      %or could take an average over lowest 10 m perhaps?
                    v5=interp1(H2,v,HSFC+5);
                    winddir(iloc).y(ifile) = wind_dir(u5,v5,lat2d,lon2d,ilat,ilon,iloc);
                end                              
                




                %            pressure(iloc).y(ifile) = psfc(ilat(iloc),ilon(iloc));
            else
                pressure(iloc).y(ifile) = get_WRF_point(nc,ih_wrf,ilat,ilon,'Pressure',iloc);
                temperature(iloc).y(ifile) = get_WRF_point(nc,ih_wrf,ilat,ilon,'Temp',iloc);
                wind(iloc).y(ifile) = get_WRF_point(nc,ih_wrf,ilat,ilon,'Wind',iloc);
                vapour(iloc).y(ifile) = get_WRF_point(nc,ih_wrf,ilat,ilon,'Vapour',iloc);
            end
        end
    end
end

for iloc=1:length(ilat)
        pressure(iloc).ilat = ilat(iloc); %diagonal of the matrix is the pairs lat2d.var(ilat(1),ilon(1)), lat2d.var(ilat(2),ilon(2)), etc.
        pressure(iloc).ilon = ilon(iloc); %diagonal of the matrix is the pairs lat2d.var(ilat(1),ilon(1)), lat2d.var(ilat(2),ilon(2)), etc.         
end


%cd(['Y:/WRF/' rundir '/AWS_positions']);
cd(['Y:/WRF/' rundir '/AWS_positions_d01']);
savename=['pressure_alt_corr_timser_level_' num2str(ih_wrf) '.mat'];
save(savename,'pressure');
savename=['temperature_timser_level_' num2str(ih_wrf) '.mat'];
save(savename,'temperature');
savename=['wind_timser_level_' num2str(ih_wrf) '.mat'];
save(savename,'wind');
savename=['winddir_timser_level_' num2str(ih_wrf) '.mat'];
save(savename,'winddir');

if ih_wrf~=0
    savename=['vapour_timser_level_' num2str(ih_wrf) '.mat'];
    save(savename,'vapour');
end
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
