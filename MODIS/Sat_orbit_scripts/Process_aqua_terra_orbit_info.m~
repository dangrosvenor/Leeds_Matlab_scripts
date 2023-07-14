%Process a bunch of TLE files to get the overpasses for each 1x1 degree
%gridbox on Earth, separated into year
%filedir='C:\Users\Dan\Documents\logbook\UW\MODIS\Orbit prediction\TLE archive\';
filedir='/home/disk/eos8/d.grosvenor/MODIS_TLE_archive/';
save_filename = 'aqua_terra_TLEs_MATLAB.mat';
%save_filename = 'aqua_TLE_TEST.mat';

iload=0;
if iload==1
    load([filedir save_filename],'tle01_aqua','tle02_aqua','tle01_terra','tle02_terra','satrec_aqua','satrec_terra');
    %    load([filedir save_filename],'tle01_aqua','tle02_aqua','tle01_terra','tle02_terra');
    %    load([filedir save_filename]);
end

%for some reason using the satrec files as freshly loaded doesn't work and
%twoline2rvMOD needs to be run first (even if the output doesn't seem to be
%used....??? Doesn't seem to matter which TLE file we run it on, as long as
%it is run...???
%clear satrec_aqua satrec_terra
for itle=1:1 %length(tle01_aqua)
    %converts the TLE into the format used by the orbit prediction
    %function
    satrec = twoline2rvMOD(tle01_aqua{itle},tle02_aqua{itle});
end

% for itle=1:1 %length(tle01_terra)
%     %converts the TLE into the format used by the orbit prediction
%     %function
%     satrec_terra(itle) = twoline2rvMOD(tle01_terra{itle},tle02_terra{itle});
% end


N_orb=16; %the number of orbits to allow for in each TLE (for array declration)

%% select the satellite and years required for processing
sat_str='terra'; years = [1999:2011];
sat_str='aqua'; years = [2002:2011];
%sat_str='aqua2'; years = [2002];
sat_str='aqua'; years = [2002];

%% copy the selected processed TLE record to satrec_process ready for
%% processing
satrec_process = eval(['satrec_' sat_str]);

%% main code
MLAT=[89.5:-1:-89.5];
MLON=[-179.5:1:179.5];

%MLAT=[50:60];
%MLON=[-4:1:4];
[LON,LAT]=meshgrid(MLON,MLAT);

clear llh
llh(1,:)=LAT(:);
llh(2,:)=LON(:);
llh(3,:)=0*ones([1 length(LAT(:))]); %zero altitude -how much difference does this make?




dtr=pi/180;
min_per_day=60*24;
%timestep for orbit trajectory calculation
dt=10/60;  %10 sec (dt is therfore in minutes)

%epoch_days and years contain all of the days and years for the whole
%record (e.g. 2002-2012 for aqua)
clear epoch_days epoch_years
for i=1:length(satrec_process)
    epoch_days(i)=satrec_process(i).epochdays;
    epoch_years(i)=satrec_process(i).epochyr;
end

epoch_years = epoch_years+2000;
epoch_years(epoch_years>=2057)=epoch_years(epoch_years>=2057)-100;
%[years,iuni] = unique(epoch_years);
datenum_year_starts = datenum(epoch_years,1,1);
[Y,epoch_month,day_of_month,H,MI,S]  = datevec(datenum(datenum_year_starts+epoch_days-1));



%% year loop
for iyear=1:length(years)
    inds_eyear2 = find(epoch_years==years(iyear));
    months_iyear = epoch_month(inds_eyear2);
    months_unique = unique(months_iyear);
    
%% month loop
    %split into months to save memory
    for imonth=3:3  %length(months_unique)
        month_now = months_unique(imonth);
        inds_emonth = find(months_iyear==month_now);
        inds_process = inds_eyear2(inds_emonth);
        
        file_savename_orbit = [num2str(years(iyear)) '_' num2str(month_now,'%.2d') '_saved_orbits_' sat_str '_' datestr(now,30)];
        
        
    
    %no. days in the year in question - only set up if are at the end of
    %the year, though (i.e. December)
    if month_now==12
        ndays_in_year = datenum(years(iyear),12,31) - datenum(years(iyear),1,1) +1;        
    elseif month_now==1
        ndays_last_year = datenum(years(iyear)-1,12,31) - datenum(years(iyear)-1,1,1) +1;
    else
        ndays_in_year = 0;
        ndays_last_year = 0;
    end


    %find the difference between each epoch and the next - but add the
    %first one for the next month to the end for differencing
    %N.B. this requires one more monthly TLE than are processing for
    %Note, days start at 1 for 1st Jan. So 1.19 of the next year is
    %equivlaent to 366.19 in a year with 365 days
    
    %need to change this to start one TLE back too - and then to only write
    %to the array if are in the correct month (as well as the 20 min thing)
    %Also process one epoch back to cover the start of the month/year
    epoch_days_thisyear = [epoch_days(inds_process(1)-1)-ndays_last_year epoch_days(inds_process) ndays_in_year+epoch_days( inds_process(end)+1 ) ];
    diff_epoch_days = diff(epoch_days_thisyear);

%    N_orb_tle = length(inds_process)-1;
    N_orb_tle = 31; %just the number of days in a month for memory allocation

    %pre-allocate the big arrays
    sat_times = NaN * ones([N_orb_tle N_orb length(MLAT)*length(MLON)]);
    sat_sza = NaN * ones([N_orb_tle N_orb length(MLAT)*length(MLON)]);
    sat_sensor = NaN * ones([N_orb_tle N_orb length(MLAT)*length(MLON)]);
    
    %separated into individual days
%    sat_times_daily = NaN*ones([31 16 length(MLAT)*length(MLON)]);
%    sat_sza_daily = NaN*ones([31 16 length(MLAT)*length(MLON)]);
%    sat_sensor_daily = NaN*ones([31 16 length(MLAT)*length(MLON)]);
    


%% TLE loop - loop through the starts of all epochs and run the sat
%% trajectory calculation forwards to the start of the next epoch. There is
%% one extra epoch in epoch_days_thisyear so that we know where to run the
%% calc forwards to
    for itle = 1:length(inds_process)
        fprintf(1,'\nProcessing tle=%d of %d',itle,length(inds_process)-1);
        %********************Run SPG4 for Multiple Orbits****************
        %NOTE:  dt, npts, and tsince offset tailored to object 32765 parameters
        %specify the time after the epoch (satrec.epochdays) that we want orbit
        %data for

        %the matlab time for this day and time (in days)
        day_mat = datenum(years(iyear),1,1) + epoch_days_thisyear(itle) -1;
%        day_mat_next = datenum(years(iyear),1,1) + epoch_days_thisyear(itle+1) -1;
        %mins since 01-Jan-0000
        mins_mat = (day_mat-1)*min_per_day;
        
%        [Ynext,Mnext,Dnext]=datevec(day_mat_next); %could prevent it from
%        doing unnecssary work at the beginning and end of each month

        %length of time that we want to span for the satellite orbit - do
        %until the next available epoch time
        tspan = diff_epoch_days(itle)*min_per_day; %convert days to mins

        %dt is in minutes
        %nhrs = 24;
        %npts=ceil(nhrs/(dt/60)); %no. of steps of dt that we want orbit info for
        tsince_offset=20; %tsince_offset is the offset before the start of
        %epoch time and after the start of the next epoch time to allow for
        %incomplete overpasses (missing the maximum elevation)

        N = ceil( (tspan + 2*tsince_offset) / dt);
        %Nref = No. steps from the reference time (01-Jan-00) for the epoch of the
        %tle minus the offset
        Nref = floor(  (mins_mat-tsince_offset) / dt  );

        %time relative to the epoch time minus the offset (in mins)
        tsince = [Nref*dt : dt : Nref*dt + N*dt] - mins_mat;
        %        tsince=-tsince_offset:dt:-tsince_offset+N*dt; %minutes
        npts = length(tsince);

        %[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,0);
        %UTsec=Ehr*3600+Emin*60+Esec;

        %calculate the satellite position and velocity for tsince array
        xsat_ecf=zeros(3,npts);
        vsat_ecf=zeros(3,npts);
        for n=1:npts
            [satrec_process(inds_process(itle)), xsat_ecf(:,n), vsat_ecf(:,n)]=spg4_ecf(satrec_process(inds_process(itle)),tsince(n));
        end
        %Scale state vectors to mks units
        xsat_ecf=xsat_ecf*1000;  %m
        vsat_ecf=vsat_ecf*1000;  %mps

        %ecf is Earth Centred Earth Fixed
        %The Earth-centered earth-fixed (ECEF or ECF) or conventional terrestrial
        %coordinate system rotates with the Earth and has its origin at the centre of the Earth.
        %The X axis passes through the equator at the prime meridian.
        %The Z axis passes through the north pole but it does not exactly coincide
        %with the instantaneous Earth rotational axis.
        %The Y axis can be determined by the right-hand rule to be passing through
        %the equator at 90� longitude.

  %convert from ECF to lat, lon, height
        %llh means  [LAT,LON,HEIGHT]. ecf2llhT produces a [3xtime] array of this for
        %the satellite
        sat_llh=ecf2llhT(xsat_ecf);            %ECF to geodetic (llh)
        % -------------------------------------------------------------------------
        % *** if just want [LAT,LON,HEIGHT] of the satellite then can stop here ***
        % *** LLH is independent of the ground location (origin) ***
        % -------------------------------------------------------------------------
        %function llh=ecf2llhT(ecf)
        %
        % Description:
        %   This function returns the [lat lon hgt] LLH coordinates
        %   of the given earth center fixed (ECF) coordinates.
        %
        % Usage:
        %   llh=ecf2llhT(ecf)
        %
        % Inputs
        %   ecf   - 3xN Matrix of Vectors in ECF coordinates (meters)
        %           where
        %              x(1,:) is Greenwich meridan; (0 lon, 0 lat)
        %              y(2,:) is Easterly
        %              z(3,:) is North Pole
        %
        % Outputs:
        %   llh   - 3xN Matrix of Vectors
        %           where
        %              llh(1,:) is latitude positive (radians)
        %              llh(2,:) is longitude positive East (radians)
        %              llh(3,:) is height (meters)
        %



        %now find the elevations and orbits for all lat lons using sat_llh
        %i.e. using
        [sat_times,sat_sza,sat_sensor] = sat_overpasses_for_latlon(llh,sat_llh,tsince,epoch_days_thisyear(itle),epoch_days_thisyear(itle+1),sat_times,sat_sza,sat_sensor,datenum(years(iyear),1,1),xsat_ecf,vsat_ecf,month_now);
        %since tsince is the time between TLEs it might span more than one day - so
        %may need more than 16 orbits in the array - but there won't be more than
        %16 per day - I think that Matlab will just add to the output array anyway

       
%        sat_times(itle,:,:) = sat_times(itle,:,:) + epoch_days_thisyear(itle);
    end




%% Reshape the arrays into lat lon grid (from a linear grid)

sat_times = reshape(sat_times,[size(sat_times,1) size(sat_times,2) length(MLAT) length(MLON)]);
sat_sza = reshape(sat_sza,[size(sat_sza,1) size(sat_sza,2) length(MLAT) length(MLON)]);
sat_sensor = reshape(sat_sensor,[size(sat_sensor,1) size(sat_sensor,2) length(MLAT) length(MLON)]);


save([filedir file_savename_orbit],'sat_times','sat_sza','sat_sensor','-V7.3');
%now have to save the results for each month separately
    %need to save the data to disk here for each month

    end %imonth
    


end %iyear


return

sat_sza2=sat_sza;
sat_sensor2=sat_sensor;
sat_times2=sat_times;

sens_thresh=65.2;
sens_thresh=67; %putting this a little high so that we include more overpasses
%sens values seems a little high anyway... why?
iremove = find(sat_sza>81.4 | sat_sensor>sens_thresh);
sat_sza2(iremove)=NaN; %remove points that aren't daylight or have too high sensorZA
sat_sensor2(iremove)=NaN;
sat_times2(iremove)=NaN;



