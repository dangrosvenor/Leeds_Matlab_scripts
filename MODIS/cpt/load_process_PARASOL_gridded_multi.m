%Saves the individual years of data for multiple years to a .mat file. The
% daynum_<varname>_year_str arrays are saved (one for each year).
% These are the averages over all 4 potential orbits for each location and
% day. Size is [366 180 360].

savefile_POLDER_multi = '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_Reff_2005-2012.mat';

ioverride_load_POLDER=1;

%Things that need setting
ioverall_mean=0;
years_parasol_str=''; %created below

 years_requested_multi = [2005:2012];
 
 for iyear_multi=1:length(years_requested_multi)
    
     years_requested = years_requested_multi(iyear_multi);
     
%% Run the script ------------------
     load_process_PARASOL_gridded
% ----------------------------------  

    
    if iyear_multi==1
        save(savefile_POLDER_multi,'LAT','LON','LAT_edges_POLDER','LON_edges_POLDER','gcm_Plon2D_edges_POLDER','gcm_Plon2D_POLDER','gcm_str_select','gcm_years_loaded_str','daynum_timeseries3_POLDER','modisyear_timeseries3_POLDER');
    end
    
    for ivar_multi=1:length(vars_PAR)
        new_var_name  = ['daymean_' vars_PAR{ivar_multi} '_' year_str];
        save_var_append(savefile_POLDER_multi,eval(['daymean_' vars_PAR{ivar_multi}]),new_var_name)
    end
     
 end