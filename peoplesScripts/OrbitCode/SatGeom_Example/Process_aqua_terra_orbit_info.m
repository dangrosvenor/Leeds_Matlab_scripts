%Process a bunch of TLE files to get the overpasses for each 1x1 degree
%gridbox on Earth, separated into year
filedir='C:\Users\Dan\Documents\logbook\UW\MODIS\Orbit prediction\TLE archive\';
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



sat_str='terra'; years = [1999:2011];
sat_str='aqua'; years = [2002:2011];
%sat_str='aqua2'; years = [2002];
sat_str='aqua'; years = [2002];

satrec_process = eval(['satrec_' sat_str]);

MLAT=[-89.5:1:89.5];
MLON=[-179.5:1:179.5];

%MLAT=[50:60];
%MLON=[-4:1:4];
[LAT,LON]=meshgrid(MLAT,MLON);

clear llh
llh(1,:)=LAT(:);
llh(2,:)=LON(:);
llh(3,:)=0*ones([1 length(LAT(:))]); %zero altitude -how much difference does this make?


dtr=pi/180;
min_per_day=60*24;
%timestep for orbit trajectory calculation
dt=10/60;  %10 sec (dt is therfore in minutes)

clear epoch_days epoch_years
for i=1:length(satrec_process)
    epoch_days(i)=satrec_process(i).epochdays;
    epoch_years(i)=satrec_process(i).epochyr;
end

epoch_years = epoch_years+2000;
epoch_years(epoch_years>=2057)=epoch_years(epoch_years>=2057)-100;
%[years,iuni] = unique(epoch_years);




for iyear=1:length(years)
    inds_eyear = find(epoch_years==years(iyear));
    %no. days in the year in question
    ndays_in_year = datenum(years(iyear),12,31) - datenum(years(iyear),1,1) +1;

    %find the difference between each epoch and the next - but add the
    %first one for the next year to the end for differencing
    %N.B. this requires one more TLE than are processing for
    %Note, days start at 1 for 1st Jan. So 1.19 of the next year is
    %equivlaent to 366.19 in a year with 365 days
    epoch_days_thisyear = [epoch_days(inds_eyear) ndays_in_year+epoch_days(inds_eyear(end)+1)];
    diff_epoch_days = diff(epoch_days_thisyear);

    %separated into individual days
    sat_times_daily = NaN*ones([366 16 length(MLAT)*length(MLON)]);
    sat_sza_daily = NaN*ones([366 16 length(MLAT)*length(MLON)]);
    sat_sensor_daily = NaN*ones([366 16 length(MLAT)*length(MLON)]);
    
    %separated by individual TLE periods
    sat_times = NaN * ones([length(inds_eyear)-1 64 length(MLAT)*length(MLON)]);
    sat_sza = NaN * ones([length(inds_eyear)-1 64 length(MLAT)*length(MLON)]);
    sat_sensor = NaN * ones([length(inds_eyear)-1 64 length(MLAT)*length(MLON)]);


    for itle = 1:length(inds_eyear)-1
        itle
        %********************Run SPG4 for Multiple Orbits****************
        %NOTE:  dt, npts, and tsince offset tailored to object 32765 parameters
        %specify the time after the epoch (satrec.epochdays) that we want orbit
        %data for

        %the matlab time for this day and time (in days)
        day_mat = datenum(years(iyear),1,1) + epoch_days_thisyear(itle) -1;
        %mins since
        mins_mat = (day_mat-1)*min_per_day;

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

        %time relative to the epoch time in mins
        tsince = [Nref*dt:dt:Nref*dt + N*dt] - mins_mat;
        %        tsince=-tsince_offset:dt:-tsince_offset+N*dt; %minutes
        npts = length(tsince);

        %[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,0);
        %UTsec=Ehr*3600+Emin*60+Esec;


        xsat_ecf=zeros(3,npts);
        vsat_ecf=zeros(3,npts);
        for n=1:npts
            [satrec_process(inds_eyear(itle)), xsat_ecf(:,n), vsat_ecf(:,n)]=spg4_ecf(satrec_process(inds_eyear(itle)),tsince(n));
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
        %the equator at 90° longitude.


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
        [sat_times(itle,:,:),sat_sza(itle,:,:),sat_sensor(itle,:,:)] = sat_overpasses_for_latlon(llh,sat_llh,tsince,epoch_days_thisyear(itle),epoch_days_thisyear(itle+1));
        %since tsince is the time between TLEs it might span more than one day - so
        %may need more than 16 orbits in the array - but there won't be more than
        %16 per day - I think that Matlab will just add to the output array anyway

        %convert to days since the start of the year
        sat_times(itle,:,:) = sat_times(itle,:,:) + epoch_days_thisyear(itle);
    end



for itle=1:size(sat_times,1)
    for ilat=1:size(sat_times,3)
        %sat_times are in days starting from tspan(1)
        %        sat_times = sat_times + epoch_days_thisyear(itle) + datenum(years(iyear),1,1)-1;
        %        sat_times = sat_times + epoch_days_thisyear(itle);
        D = floor(squeeze(sat_times(itle,:,ilat)));
        %make into days from the start of the current year
        %N.B. epoch_days are decimal days (i.e they include the time too)
        %        [Y,MO,D,H,MI,S] = datevec(sat_times);
        %days since the start of the current year
        %        sat_times = sat_times - datenum(years(iyear),1,1) - 1;
        D(isnan(D))='';
        day_uni = unique(D);

        for idaysat=1:length(day_uni)
            dsat = day_uni(idaysat);
            idaysati=find(D==dsat);
            sat_times_daily(dsat,1:length(idaysati),ilat) = sat_times(itle,idaysati,ilat);
            sat_sza_daily(dsat,1:length(idaysati),ilat) = sat_sza(itle,idaysati,ilat);
            sat_sensor_daily(dsat,1:length(idaysati),ilat) = sat_sensor(itle,idaysati,ilat);                   
        end

    end
end

sat_times_daily = reshape(sat_times_daily,[size(sat_times_daily,1) size(sat_times_daily,2) length(MLAT) length(MLON)]);
sat_sza_daily = reshape(sat_times_daily,[size(sat_times_daily,1) size(sat_times_daily,2) length(MLAT) length(MLON)]);
sat_sensor_daily = reshape(sat_times_daily,[size(sat_times_daily,1) size(sat_times_daily,2) length(MLAT) length(MLON)]);

%now have to save the results for each year separately

%also need to test whether running 

end



