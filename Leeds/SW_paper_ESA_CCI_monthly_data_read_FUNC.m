function SW_paper_ESA_CCI_monthly_data_read_FUNC()
%Read data from the ESA CCI cloud dataset and save

dat_dir = '/home/disk/eos8/d.grosvenor/ESA_Cloud_CCI/AVHRR_PMv3_L3C_Monthly/';

var_list={'cer_liq','cot_liq','ctt'};
save_file = [dat_dir 'ESA_Cloud_CCI_Monthly_vars_for_Nd.mat'];


%var_list={'cfc','cfc_day','cfc_low','cfc_mid','cfc_high'};
%save_file = [dat_dir 'ESA_Cloud_CCI_Monthly_Cloud_Fraction.mat'];

% cf_tot = nc{'cfc'}(:); %total CF
% cf_tot_day = nc{'cfc_day'}(:); %total CF daytime
% cf_tot_night = nc{'cfc_night'}(:); %total CF night
% cf_low = nc{'cfc_low'}(:); %day and night combined?
% cf_mid = nc{'cfc_mid'}(:);
% cf_high = nc{'cfc_high'}(:);




files = dir([dat_dir '*ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA*netcdf3.nc']);

%Loop over all the files and get the times
for ifile=1:length(files)
    file = [dat_dir files(ifile).name];
    nc=netcdf(file);
    time = nc{'time'}(:); %days since 1970-1-1;
    time_matlab_in(ifile) = datenum('01-Jan-1970') + time;
    if ifile==1
        lat= nc{'lat'}(:); dlat=abs(mean(diff(lat)));
        lon = nc{'lon'}(:); dlon=abs(mean(diff(lon(1:2)))); %-180 to +180 convention, which is what we want
        %lon(lon>180)=lon(lon>180)-360;
        
        [gcm_Plon2D_CCI,gcm_Plat2D_CCI] = meshgrid(lon,lat);                
    end
    nc=close(nc);
end

%Make a dataset that has all of the months with the missing months padded
%with NaNs - for all years even if only have partial data for 1st year
[Y,M,D] = datevec(time_matlab_in);
nt = (max(Y)-min(Y)+1)*12;

%Make the blank arrays
for ivar=1:length(var_list)
    eval_str = [var_list{ivar} '=NaN*ones([nt length(lat) length(lon)]);'];
    eval(eval_str);
end

%Loop through the times covered by the files (and possibly before if starts
%midway through a year)
icount=0;
for year=min(Y):max(Y)
   for month=1:12
       icount=icount+1;
       Y_out(icount) = year;
       M_out(icount) = month;
       time_matlab(icount) = datenum(year,month,1); %output matlab time for 1st day of month
              
       im=find(M==month & Y==year);
       if length(im)==0           
          continue %No file for this month, continue and leave as NaNs
       else %open the file and copy the data to the arrays
           file = [dat_dir files(im).name];
           nc=netcdf(file);
           for ivar=1:length(var_list)                              
               eval_str = ['dat=nc{''' var_list{ivar} '''}(:);'];
               eval(eval_str);
               dat(dat<-998) = NaN; %Fill value of -999
               eval_str = [var_list{ivar} '(icount,:,:)=dat;'];
               eval(eval_str);
           end
           nc=close(nc);
       end
   end
    
end






save(save_file,'-V7.3');
%,'sw_up_toa','sw_net_coarse','M_coarse_grain','N_coarse_grain','yr_start_UM','yr_end_UM',...
%    'sw_in_TOA_UM','sw_deepc','time_matlab','lat2d_deepc','lon2d_deepc');




