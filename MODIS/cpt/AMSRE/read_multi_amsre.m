%read multiple daily amsre files
nc_dir='/home/disk/eos5/d.grosvenor/AMSRE/'; amsre_label = 'daily';

gcm_time_of_day_select=0;


years_requested = [2002:2011]; %will leave out 2006 as it is only a partial year
%years_requested = [2007]; %will leave out 2006 as it is only a partial year
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


%files = dir([nc_dir field_names{ifield} '*']);
files = dir([nc_dir 'y*']);

clear filenames_csat years_csat months_csat iyear_csat
j=0;
iyear=1;
for i=1:length(files)
    
    if length(strfind(files(i).name,'.nc'))>0       
        year_of_file = str2num(files(i).name(19:22));
        
        if length(find(year_of_file==years_requested))>0
            j=j+1;
            filenames_csat{j} = files(i).name;
            months_csat(j) = str2num(files(i).name(23:24));
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
dat = nc{field_names{1}}(:,:,:);
%lat_csat = nc{'LAT'}(:);
%lon_csat = nc{'LON'}(:);
%nasc_csat = nc{'ASCDES'}(:);

gcm_lat = nc{'latitude'}(:);
%gcm_lon = nc_inst{'lon'}(:);

gcm_lon = -179:2:179;

gcm_slat = -90:2:90;
gcm_slon = -180:2:180;


nlat=size(dat,1);
nlon=size(dat,2);

years_csat_unique = unique(years_csat);

%nyears=length(years_csat_unique);
nyears=length(years_requested);

for ifield=1:length(field_names)
     eval_str = [lower(field_names{ifield}) '_monthly_' amsre_label '=NaN*ones([nyears*12 nlat nlon]);'];
     eval(eval_str);     
end


clear month_calipso_cf year_calipso_cf

for j=1:length(filenames_csat)
        


        ind = months_csat(j)+(years_csat(j)-years_csat(1))*12;
%        low_CF_calipso = read_calipso_cmor(nc_dir,nc_inst_file,'cllcalipso',ilat,ilon,'(:,ilat,ilon)'); 

for ifield=1:length(field_names)
    field_name = field_names{ifield};
    eval(['dat = read_calipso_ipsl_night(nc_dir,filenames_csat{j},''' field_name ''',0,0,''(:,:,:)'',amsre_label);']);
    eval([lower(field_name) '_monthly_' amsre_label '(ind,:,:) = 100*dat(:,:);']); %convert to % to be consistent with the other CALIPSO files
end

%        eval([lower(field_name) '(ind+1,:,:) = dat(2,:,:);']); %descending        
%        eval([lower(field_name) '_daily(j,:,:) = 0.5*(dat(1,:,:) + dat(2,:,:));']); %daily average
%        asc_calipso_matt(ind) = 1;
%        asc_calipso_matt(ind+1) = 2;
        
        month_calipso_cf(ind) = months_csat(j);
        year_calipso_cf(ind) = years_csat(j);    
        
%        month_calipso_matt_daily(j) = months_csat(j);
%        year_calipso_matt_daily(j) = years_csat(j);    
        
            
%    end
                              
end




gcm_str = 'CALIPSO_monthly';
gcm_str_select = 'CALIPSO_monthly';
daynum_timeseries3_CALIPSO_monthly = 1:length(month_calipso_cf);
gcm_time_UTC_CALIPSO_monthly = zeros([1 length(month_calipso_cf)]);
am3_dataset = [amsre_label];


Plat=gcm_slat;
Plon=gcm_slon;


i180=find(Plon>180);
Plon2=Plon;
Plon2(i180)=Plon2(i180)-360;

%Plon=[Plon2(1:end); Plon2(1)];
Plon=Plon2;



[gcm_Plon2D_edges_CALIPSO_monthly,gcm_Plat2D_edges_CALIPSO_monthly]=meshgrid(Plon,Plat);


%cell centres - contained in gcm_lat - seems that the fields are stored as point values at gcm_lat, which runs to lat=+/-90
% - so gcm_lat are effectively the cell face positions. And we ignore the
% first and last lat value (since cells can't extend beyond 90 degree lat)
Plat=gcm_lat;
Plon=gcm_lon;

i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

[gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly]=meshgrid(Plon,Plat);

%dlon = 1.87, dlat = 2.5. Both constant throughout the grid for CAM5

% dlat = diff(gcm_lat);
% dlon = diff(gcm_lon);
% Plat=gcm_lat+[dlat; dlat(end)]/2;
% Plon=gcm_lon+[dlon; dlon(end)]/2;
% 
% i180=find(Plon>180);
% Plon(i180)=Plon(i180)-360;
% 
% [gcm_Plon2D_edges,gcm_Plat2D_edges]=meshgrid(Plon,Plat);

fprintf(1,'\n Done read amsre multi\n');






