%Pre-processor to get everything into .mat files

%Runs :-
%   UM_quick_plot_commands.m

% Set some stuff

isave_plot=0;
iclose_plot=1;

var_UM = 'wind speed';
var_UM = 'accum_number_ukca';
%var_UM = 'accum_mass_H2SO4_ukca';
%var_UM = 'accum_mass_BC_ukca';
%var_UM = 'accum_mass_OC_ukca';
%var_UM = 'emissions_BC_biomass';
%var_UM = 'emissions_BC_biofuel';
%var_UM = 'emissions_BC_fossil';
var_UM='LWP';
%var_UM = 'SW_down_surf';
%var_UM = 'SW_down_clean_surf';
%var_UM = 'SW_down_clean_clear_surf';
%var_UM = 'Nd_lwc_weighted_UKCA';
%var_UM = 'high_cloud_amount';
%var_UM = 'mid_cloud_amount';
%var_UM = 'low_cloud_amount';


% Grid to coarsen to (if icoarse is set to one below)
coarse_grid = 'global grid'; %To match to the global model grid specified in UM_quick_plot_commands.m
coarse_grid = '1 deg'; %1x1deg regular (MODIS) grid.


itimes_DRIVER=1; %1=15UTC, 2=18UTC
itimes_DRIVER='all'; %can set to show all times or a certain one

iplot_wind_arrows=0;
iUTC=1;

%some defaults that are likely to get overwritten
file_type='nc';
%lat_var = 'grid_latitude'; lon_var = 'grid_longitude';
lat_var = 'Latitude'; lon_var = 'Longitude';
ivar_dir=0;

%um_case='u-ar365'; pole_lat=70; pole_lon=284;
%um_case='u-as103 global'; pole_lat=45.0; pole_lon=145.0;
%um_case='u-as103 nested'; pole_lat=45.0; pole_lon=145.0;
um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; %global

um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umnsaa'; icoarse=0; %nest
%um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'glm'; icoarse=0; %global
%um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'NAtlantic'; icoarse=1; %nest

um_case='u-au689'; pole_lat=45.0; pole_lon=145.0; run_type = 'umnsaa'; icoarse=0; %nest - e.g. for hi-res animation
%um_case='u-au689'; pole_lat=45.0; pole_lon=145.0; run_type = 'umnsaa'; icoarse=1; %nest

um_case='u-au652'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
%um_case='u-au536'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
um_case='u-av503'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
%um_case='u-av504'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions

if ivar_dir==1
    var_dir = ['/' var_UM '/'];
else
    var_dir = '';
end

dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case var_dir];
        
switch run_type
    case {'umglaa','umnsaa'}
        filename_multi= [dirUM '/' run_type '*' var_UM '*.nc'];  file_type='nc';
    otherwise
        filename_multi= [dirUM '/*' var_UM '*'  run_type '*.nc'];  file_type='nc';
end
        



nested_or_global='global';

%If want to plot an outline of the nest on the global map
filename_nest = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat';
nest=load(filename_nest);
irotated_pole_box=1;
              

filenames=dir(filename_multi);

if length(filenames)==0
    error('No matching files!');
end

clear date_num
for ifile=1:length(filenames)
    filename = filenames(ifile).name;
    idate = findstr(var_UM,filename);
    switch run_type
        case {'umglaa','umnsaa'}
            date_str = filename(idate-13:idate-2);
        otherwise
            date_str = filename(idate+1+length(var_UM):idate+1+length(var_UM)+11);
    end

    
    
    date_num(ifile) = str2num(date_str);    
end

[Y,I] = sort(date_num);


clear LWP_ALL time_ALL

it_ALL=0;    
for ifile=1:length(I)
    
    filename = [dirUM '/' filenames(I(ifile)).name];
    
    %Run the script
    itimes = itimes_DRIVER;
    UM_quick_plot_commands    


    for it=1:length(itimes) 
        it_ALL = it_ALL + 1;
        eval([var_UM '_ALL{it_ALL} = dat_modis_save{it};']);
        time_ALL(it_ALL) = time_save(it);
    end

end

if icoarse==1
    coarse_str='';
else
    coarse_str='native res';
end

save_file = [dirUM '/' run_type '_' var_UM '_' coarse_str '_ALL.mat']
save(save_file,[var_UM '_ALL'],'time_ALL','gcm_Plat2D_UM','gcm_Plon2D_UM','gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM');

