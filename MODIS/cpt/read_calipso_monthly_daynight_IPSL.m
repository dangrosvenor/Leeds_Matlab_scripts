%read monthly CALIPSO data - DAYTIME only
% Use read_calipso_monthly_night_IPSL for nighttime AND average of day and
% night. But can also load day from it now too.
%data fields will be of size [nyears*12 nlat nlon]
%
%directory for data
nc_dir='/home/disk/eos8/d.grosvenor/CPT/CALIPSO_day_night_CF/day/'; calipso_daynight_label = 'DAYTIME';

time_series_type = 'CALIPSO';

gcm_time_of_day_select=0;

clear field_names

ifield=1;
% low, mid & high cloud fractions
field_names{ifield}='cllcalipso'; ifield=ifield+1; % 1000 hPa to 680 hPa
field_names{ifield}='clmcalipso'; ifield=ifield+1; % 680 hPa to 440 hPa
field_names{ifield}='clhcalipso'; ifield=ifield+1; % 440 hPa to 0 hpa


years_requested = [2007:2010]; %will leave out 2006 as it is only a partial year
%years_requested = [2007]; %will leave out 2006 as it is only a partial
%year (starts in June since the satellite only went up shortly before that)

years_calipso_str='';
if max(diff(years_requested)==1)
    years_calipso_str=[num2str(years_requested(1)) ' to ' num2str(years_requested(end))];
else
    for i=1:length(years_requested)
        years_calipso_str=[years_calipso_str ' ' num2str(years_requested(i))];
    end
end
gcm_years_loaded_str = years_calipso_str;


%% Search through all of the files in the directory to find the ones we
%% want to process
for ifield=1:length(field_names)

    files = dir([nc_dir field_names{ifield} '*']);

    clear filenames_csat years_csat months_csat iyear_csat
    j=0;
    iyear=1;
    for i=1:length(files)

        if length(strfind(files(i).name,'.nc'))>0
            year_of_file = str2num(files(i).name(60:63));

            if length(find(year_of_file==years_requested))>0
                j=j+1;
                filenames_csat{j} = files(i).name;
                months_csat(j) = str2num(files(i).name(64:65));
                years_csat(j) = year_of_file;
                %on first pass set to current year
                if j==1
                    year_old = years_csat(j);
                end

                if years_csat(j) ~= year_old
                    iyear=iyear+1;
                end
                iyear_csat(j)=iyear;
                year_old = years_csat(j);
            end

        end

    end

%% Now process the files we identified

    nc = netcdf([nc_dir filenames_csat{1}]);
    dat = nc{field_names{ifield}}(:,:,:);

    gcm_lat = nc{'lat'}(:);

    gcm_lon = -179:2:179;

    gcm_slat = -90:2:90;
    gcm_slon = -180:2:180;


    nlat=size(dat,1);
    nlon=size(dat,2);

    years_csat_unique = unique(years_csat);
    
    nyears=length(years_requested);

    %make an array of NaNs of the required size
    eval_str = [lower(field_names{ifield}) '_monthly=NaN*ones([nyears*12 nlat nlon]);'];
    eval(eval_str);


    clear month_calipso_cf year_calipso_cf

    for j=1:length(filenames_csat)

        field_name = field_names{ifield};
        ind = months_csat(j)+(years_csat(j)-years_csat(1))*12;

        eval(['dat = read_calipso_cmor(nc_dir,filenames_csat{j},''' field_name ''',0,0,''(:,:,:)'');']);
        eval([lower(field_name) '_monthly(ind,:,:) = 100*dat(:,:);']); %convert to % to be consistent with the other CALIPSO files


        month_calipso_cf(ind) = months_csat(j);
        year_calipso_cf(ind) = years_csat(j);


    end


end

daynum_timeseries3_CALIPSO_monthly = 1:length(month_calipso_cf);
gcm_time_UTC_CALIPSO_monthly = zeros([1 length(month_calipso_cf)]);

gcm_time_matlab_CALIPSO_monthly = datenum(year_calipso_cf,month_calipso_cf);


Plat=gcm_slat;
Plon=gcm_slon;


i180=find(Plon>180);
Plon2=Plon;
Plon2(i180)=Plon2(i180)-360;

Plon=Plon2;



[gcm_Plon2D_edges_CALIPSO_monthly,gcm_Plat2D_edges_CALIPSO_monthly]=meshgrid(Plon,Plat);


Plat=gcm_lat;
Plon=gcm_lon;

i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

[gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly]=meshgrid(Plon,Plat);


gcm_str = 'CALIPSO_monthly';
gcm_str_select = 'CALIPSO_monthly';
am3_dataset = [calipso_daynight_label];


fprintf(1,'\n Done read calipso monthly\n');






