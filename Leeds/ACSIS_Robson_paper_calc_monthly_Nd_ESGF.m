%ACSIS_Robson_paper_calc_monthly_Nd_ESGF.m

%Processes the daily Nd data into monthly. After this need to run the 

dat_type='CMIP';
dat_type='AMIP';

switch dat_type
    case 'CMIP'        
        esgf_dir='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_historical/';
        files=dir([esgf_dir 'r*']);
    case 'AMIP'
        esgf_dir='/home/disk/eos15/d.grosvenor/UM/UKESM/CMIP6_AMIP_data/';
        files(3).name='';
end


nd_var = 'scldncl';

%ref file to copy lat lon, etc. data from
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_all_rsut.mat';


iens_inds = [3:length(files)];


for iens=iens_inds
    iens
    
    dirUM = [esgf_dir files(iens).name '/'];            
    nc=netcdf([dirUM nd_var '/merged.nc']);
    
    var_nc=nc{nd_var};
    dims=dim(var_nc);
    nmonths = length(dims{1}) / 30; %number of days divided by 30 since uses 30 day months
    
    iday=1;
    for imon=1:nmonths
        %fprintf(1,'%d ',imon);
        
        if imon==1 & iens==iens_inds(1)
            clear opts
            opts.time_var='time';
            opts.lat_var = 'lat';
            opts.lon_var = 'lon';
            opts.time_ref = datenum('01-Jan-1850');
            opts.time_fconv = 1; %conversion multiplier to get to days
            var_UM = [nd_var '_time_only'];
            load_type = 'merged netCDF';
            pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
            dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,[],[],var_UM,opts); %data is ordered [time lat lon]. 180 times (monthly over 15 years)            
            dat_global.dat = NaN*ones([nmonths size(dat_global.dat,2) size(dat_global.dat,3)]);
        end
        
        %UKESM1 uses 30-day month calendar
        inds=iday:iday+29;
        iday=iday+30;

        nd_dat = nc{nd_var}(inds,:,:);
        inan=find(nd_dat>1e19);
        nd_dat(inan)=NaN;
        nd_mon = meanNoNan(nd_dat,1);
        
        if imon==1
            siz=size(nd_dat);
            dat_global.dat = NaN*ones([nmonths siz(2) siz(3)]);
        end
        
        dat_global.dat(imon,:,:) = nd_mon;
        
    end
    
    nc=close(nc);
    savename = [dirUM nd_var '/monthly_Nd.mat'];
    save(savename,'dat_global','-V7.3');
       
    
end