%Stripped down template for running plot_global_maps to produce a global
%map

% Set some stuff

var_UM = 'wind speed';
var_UM = 'accum_number_ukca';
var_UM = 'accum_mass_H2SO4_ukca';
var_UM = 'accum_mass_BC_ukca';
%var_UM = 'accum_mass_OC_ukca';
var_UM = 'emissions_BC_biomass';
var_UM = 'emissions_BC_biofuel';
%var_UM = 'emissions_BC_fossil';
var_UM='LWP';

um_case='u-ar365'; pole_lat=70; pole_lon=284;
um_case='u-as103 global'; pole_lat=45.0; pole_lon=145.0;
%um_case='u-as103 nested'; pole_lat=45.0; pole_lon=145.0;
um_case='u-at459 global'; pole_lat=45.0; pole_lon=145.0;

itimes=1; %1=15UTC, 2=18UTC
itimes='all'; %can set to show all times or a certain one

iplot_wind_arrows=0;
iUTC=1;

%soem defaults that are likely to get overwritten
file_type='nc';
%lat_var = 'grid_latitude'; lon_var = 'grid_longitude';
lat_var = 'Latitude'; lon_var = 'Longitude';
        
switch um_case
    case 'u-ar365'

        %Define lat and lon grids
        filename = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-ar365/global/u-ar365_U_200811131500___glm_saved_slice_.pp.nc';
        nc = netCDF(filename);
        %Output here is for 500m altitude
        filenameV = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-ar365/global/u-ar365_V_200811131500___glm_saved_slice_.pp.nc';
        ncV = netCDF(filenameV);
        
        iUTC=0; time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)

    case 'u-as103 global' %ukca tests
        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/glm/u-as103_' var_UM '_201608010252___glm_saved_slice_.pp.nc'];
%        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/glm/u-as103_' var_UM '_201608010252_iz=0___glm_saved_slice_.pp.nc'];        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ics/u-as103_' var_UM '_201608010252___glm_saved_slice_.pp.nc'];
        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/global_ukca_emiss_BC_biofuel_NetCDF3.nc']; file_type='nc_emiss_CMIP_global';
        nc = netCDF(filename);        
        nested_or_global='global';
        
    case 'u-as103 nested' %ukca tests        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/u-as103_' var_UM '_201608010259_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];
%        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/u-as103_' var_UM '_201608010259_iz=0_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ics/u-as103_' var_UM '_201608010000_iz=0_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ukca_emiss_BC_biomass.nc'];
        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat']; file_type='mat';
       
%        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ukca_emiss_BC_biofuel_NetCDF3.nc']; file_type='nc_emiss_CMIP';
%        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ukca_emiss_BC_fossil_NetCDF3.nc']; file_type='nc_emiss_CMIP';        
        
         nested_or_global='nested';
        
        switch file_type
            case {'nc','nc_emiss_CMIP'}
                nc = netCDF(filename);
        end
        
    case 'u-at459 global' %ukca tests
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/u-as103_' var_UM '_201608010259_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];
        %        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/u-as103_' var_UM '_201608010259_iz=0_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ics/u-as103_' var_UM '_201608010000_iz=0_NAtlantic_4p0_L70_ukv_saved_slice_.pp.nc'];
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ukca_emiss_BC_biomass.nc'];

        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat']; file_type='mat';

        %        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/ukca_emiss_BC_biofuel_NetCDF3.nc']; file_type='nc_emiss_CMIP';
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-at459/.nc']; file_type='nc_emiss_CMIP';
        
        filename= ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-at459/umglaa_p201608010252_' var_UM '_saved.nc'];  file_type='nc';      

        nested_or_global='global';


        
        filename_nest = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat';
        nest=load(filename_nest);
        irotated_pole_box=1;
              
end

%%
UM_quick_plot_commands


