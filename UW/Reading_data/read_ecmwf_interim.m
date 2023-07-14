%Read in ECMWF LTS, etc from Dargan's data

ecmwf_data_type = 'Dargan';
ecmwf_data_type = 'Manuel';
switch ecmwf_data_type
    case 'Dargan'

        %Monthly averaged data
        %Should have:-   3D fields: q, T, z, u, v, w, vorticity, divergence
        %                2D fields: Psurf, total column water vapour, u10, v10,
        %                T2m, Tdew
        direc_monthly = '/home/disk/eos10/dargan/ERAinterim_monthly/';

        %4xdaily data
        %Should have:-   3D fields: q, T, z, u, v
        %                2D fields: surface temp (T2??), Psurf, surface dewpoint
        %                temperature
        direc_daily = '/home/disk/eos10/dargan/ERAinterim_daily0Z/'; %and _daily6Z,  _daily12Z, _daily18Z


        nc_ecTm = netcdf([direc_monthly 't.nc'],'nowrite');

        ecTm = nc_ecTm{'t'}(:,:,:,:);
        %sclf =
        %dtemp = f{'2d'}(:,:,:);
        %dtemp = offs + sclf*dtemp

        %Here have data from 01-Jan-1989 to 01-Dec-2007
        ecTime_read = nc_ecTm{'time'}(:); %hours since 1900-01-01 00:00:0.0
        ecTime = datenum('01-Jan-1900') + ecTime_read/24;


    case 'Manuel'

        %OR Manuel has some data. I converted it from grib to netCDF using ncl
        %e.g. ncl_convert2nc t_press_01_12_2008.grib

        direc_Man = '/home/disk/eos5/d.grosvenor/ERA_Interim/ManuelZ/';
        nc_ecT_Man = netcdf([direc_Man 't_press_01_12_2008.nc'],'nowrite');
        ecTime_Man_read = nc_ecT_Man{'initial_time0_hours'}(:); %hours since 1800-01-01 00:00
        ecTime_Man = datenum('01-Jan-1800') + ecTime_Man_read/24;
%        ecDateVec_Man = datevec(ecTime_Man);    %[Y,MO,D,H,MI,S]

        ecLon_Man = nc_ecT_Man{'g0_lon_3'}(:); %runs from -180 to +178.5 (steps of 1.5 degree)
        ecLat_Man = nc_ecT_Man{'g0_lat_2'}(:); %runs from 90 to -90 (steps of 1.5 degree       


        %only has data for 6 levels, though:-
        %      200 500 700 850 925 1000
        % But for LTS only need surface and 700, so can use 1000 and 700 hPa
        % levels.
        nc_levs_Man =  netcdf('/home/disk/cresta2/mzuluaga/Data/ERA_interim/pressure/levels_for_ERApressure_v0.nc','nowrite');
        levs_Man_L6 = nc_ecT_Man{'lv_ISBL1'}(:); %pressure levels in hPa (=[200 500 700 850 925 1000])

        ecT_Man = nc_ecT_Man{'T_GDS0_ISBL'}(:,:,:,:);  %K
        
        %calculate the daily averages (average of 0,6,12 and8 hours)
        %And monthly averages for each hour - size(ecT_mon_Man) = [4    12     6   121   240]
        %    ---> [hour month level lat lon]
        [ecT_daily_Man,ecT_mon_Man,time_daily] = ecmwf_averages(ecT_Man,ecTime_Man);

        ecLTS_mon_Man = squeeze(potemp(ecT_mon_Man(:,:,3,:,:),700e2) - potemp(ecT_mon_Man(:,:,1,:,:),1000e2));
        
        

end

