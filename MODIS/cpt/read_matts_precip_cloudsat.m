%read Matt L's CloudSat precip data (2007-2010, monthly) in mm/hr
precip_dir='/home/disk/eos8/d.grosvenor/CPT/precip/';


ifield=1;
field_names{ifield}='RAIN_WARM_OCEAN'; ifield=ifield+1; %surface rain rate from warm clouds (mm/hr)
field_names{ifield}='RAIN_OCEAN'; ifield=ifield+1; 

field_names{ifield}='NUM_WARM_OCEAN'; ifield=ifield+1; %total no. valid warm CloudSat pixels per month
field_names{ifield}='NUM_TOTAL_OCEAN'; ifield=ifield+1; %total no. valid CloudSat pixels per month
field_names{ifield}='NUM_PRECIP_OCEAN'; ifield=ifield+1; %? - no. precipitating pixels?

field_names{ifield}='NUM_SATURATED_OCEAN'; ifield=ifield+1; %?
field_names{ifield}='NUM_SATURATED_WARM_OCEAN'; ifield=ifield+1; %?

files = dir(precip_dir);

clear filenames_csat years_csat months_csat iyear_csat month_calipso_matt year_calipso_matt month_calipso_matt_daily year_calipso_matt_daily
j=0;
iyear=1;
for i=3:length(files)
    
    if length(strfind(files(i).name,'.nc'))>0
        j=j+1;

        filenames_csat{j} = files(i).name;
        months_csat(j) = str2num(files(i).name(6:7));
        years_csat(j) = str2num(files(i).name(1:4));
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

nc = netcdf([precip_dir filenames_csat{1}]);
dat = nc{'RAIN_WARM_OCEAN'}(:,:,:);
%lat_csat = nc{'LAT'}(:);
%lon_csat = nc{'LON'}(:);
%nasc_csat = nc{'ASCDES'}(:);

nasc=size(dat,1);
nlat=size(dat,2);
nlon=size(dat,3);

years_csat_unique = unique(years_csat);
nyears=length(years_csat_unique);

 for ifield=1:length(field_names)
     eval_str = [lower(field_names{ifield}) '=NaN*ones([nasc*nyears*12 nlat nlon]);'];
     eval(eval_str);
     eval_str = [lower(field_names{ifield}) '_daily=NaN*ones([nyears*12 nlat nlon]);'];
     eval(eval_str);
 end



for j=1:length(filenames_csat)
        
    nc = netcdf([precip_dir filenames_csat{j}]);     
    
    for ifield=1:length(field_names)
        field_name = field_names{ifield};
        ind = (j-1)*2+1;
        eval(['dat = read_matts_precip_field(nc,field_name);']);        
        eval([lower(field_name) '(ind,:,:) = dat(1,:,:);']); %ascending
        eval([lower(field_name) '(ind+1,:,:) = dat(2,:,:);']); %descending        
        eval([lower(field_name) '_daily(j,:,:) = 0.5*(dat(1,:,:) + dat(2,:,:));']); %daily average
        asc_calipso_matt(ind) = 1;
        asc_calipso_matt(ind+1) = 2;
        
        month_calipso_matt(ind:ind+1) = months_csat(j);
        year_calipso_matt(ind:ind+1) = years_csat(j);    
        
        month_calipso_matt_daily(j) = months_csat(j);
        year_calipso_matt_daily(j) = years_csat(j);    
        
            
    end
                              
end


%do some screening by pixel number
npix_screen=10;
rain_warm_ocean(num_warm_ocean<npix_screen)=NaN;
rain_ocean(num_total_ocean<npix_screen)=NaN;

daynum_timeseries3_CLOUDSAT_PRECIP = 1:size(rain_warm_ocean,1);
modisyear_timeseries3_CLOUDSAT_PRECIP = 1:size(rain_warm_ocean,1); %Dummy data


%2degrees in LAT, 4 degrees in LON (90x90)
lat_matt_edges = [-90:2:90]; %91 edges
lon_matt_edges = [-180:4:180]; %91 edges

%i180=find(lon_matt_edges>180);
%lon_matt_edges(i180)=lon_matt_edges(i180)-360;

[Plon2D_matt_edges,Plat2D_matt_edges] = meshgrid(lon_matt_edges,lat_matt_edges);

lat_matt_centres = [-89:2:89]; %90 faces
lon_matt_centres = [-178:4:178]; %91 edges
%i180=find(lon_matt_centres>180);
%lon_matt_centres(i180)=lon_matt_centres(i180)-360;

[Plon2D_matt_centres,Plat2D_matt_centres] = meshgrid(lon_matt_centres,lat_matt_centres);

gcm_years_loaded_str='y2007_to_2010';
am3_dataset = '';

gcm_str = 'CLOUDSAT_PRECIP';

fprintf(1,'\n Done read precip\n');






