% Produce maps of LWP difference between volcano on and off runs.
% Will use UM_maps_generic_time_loop_RUN_v1_SWTOA_Iceland_31Oct2016.m to
% loop through and make overall time averages.

%% Start by loading in data that was created by the indicated
%% script
% File created using UM_maps_generic_time_loop_RUN_v1_SWTOA_Iceland_31Oct2016
% (time averages of the SW_TOA field)
%loadfile = '/home/disk/eos10/d.grosvenor/UM/Iceland/saved_data_UM_maps_generic_time_loop_RUN_v1_SWTOA_Iceland_31Oct2016.mat';
loadfile = '/home/disk/eos10/d.grosvenor/UM/Iceland/saved_data_UM_maps_generic_time_loop_RUN_v1_LWP_Iceland_31Oct2016.mat';
%loadfile = '/home/disk/eos10/d.grosvenor/UM/Iceland/saved_data_UM_maps_generic_time_loop_RUN_v1_Nd_Iceland_31Oct2016.mat';
load(loadfile);

% --- settings for UM_maps_*_FUNC script ---
UM_map.isave_plot=0;
icoarsen=0; icoarsen_diff=1;  %Whether to coarse grain the data (done within this script)
iplot_maps=1;

% --- settings just for this script ---
plot_diff=1; %a difference plot?
iprc_diff=0; %percentage differnce instead of absolute?
thresh_diff=0; %min value for which to show % diff
iclear_sky_only=0;
thresh_lwp_clearsky = 5; %g/m2
plot_pdf_diff=1;


xlab_pdf_DRIVER = 'LWP';

% Indices of the sets of model runs to take the difference between
% Low aerosol cases for differences, 1st one minus second one
idiff_cases{1}(1)=3; %Low aerosol, volcano ON
idiff_cases{1}(2)=1; %Low aerosol, volcano OFF
idiff_cases{2}(1)=2; %High aerosol, volcano ON
idiff_cases{2}(2)=4; %High aerosol, volcano OFF

%Nice strings for labels for the difference performed
diff_str_DRIVER{1} = 'Low background aerosol';
diff_str_DRIVER{2} = 'High background aerosol';


%% For PDF plotting

Ybins_DRIVER = [0:20:800];
if icoarsen_diff==1
    Ybins_DRIVER = [-50:4:50];
else
    Ybins_DRIVER = [-50:1:50];    
end

logbin_norm_driver = 0; %Whether hav lognormal bins and normalise by these
i_plot_norm_driver = 1; %whether to normalise
i_div_bin_widths_driver = 1; %Whether to divide by the bin widths (also normalises)

pdf_type_driver='normal';  %normal or cumulative PDF
%pdf_type_driver='cumulative';



%--- Load the UM data

   
% Set the colour limits    
UM_map.iset_min_clim=1;
UM_map.clim_min=-20;
UM_map.iset_max_clim=1;
UM_map.clim_max=20;


UM_map.iplot_markers=0;
    
%% Plot function
% if iplot_maps==1      
%     UM_time_out = UM_maps_20141126T041216_FUNC(UM_map);
% end

gcm_Plat2D_edges_UM = UM_time_out{1}.gcm_Plat2D_edges_UM;
gcm_Plon2D_edges_UM = UM_time_out{1}.gcm_Plon2D_edges_UM;
gcm_Plat2D_UM = UM_time_out{1}.gcm_Plat2D_UM;
gcm_Plon2D_UM = UM_time_out{1}.gcm_Plon2D_UM;



UM_commands_plot_diff_PDF
