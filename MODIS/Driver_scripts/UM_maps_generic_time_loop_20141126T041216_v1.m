function [UM_time_out] = UM_maps_generic_time_loop_20141126T041216_v1(UM_time_in)

% Map plots/domain averages for given range of times for UM - generic
% function

%Convert all of the variable names in the input structure to actual names
%for ease of use
name_struc='UM_time_in'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end

%Just copy the variable names
UM_map = UM_time_in;


%% Set some defaults in case not supplied
%If not supplied then set a default value
if exist('time_tol')~=1
    time_tol = 2/3600/24; %set a default value
    % Can also set to 'same month' for matching by month and year
end
if exist('iclose_figs')~=1
    iclose_figs=0;
end
if exist('noplot')~=1
    noplot=0;
end
if exist('UM_time_loop_action')~=1
    UM_time_loop_action='timeseries';
end



%% For PDF plotting

Ybins_DRIVER = [0:20:800];
logbin_norm_driver = 0; %Whether hav lognormal bins and normalise by these
i_plot_norm_driver = 1; %whether to normalise
i_div_bin_widths_driver = 1; %Whether to divide by the bin widths (also normalises)

pdf_type_driver='normal';  %normal or cumulative PDF
%pdf_type_driver='cumulative';


%
% %--- Load and process the data

%varname refers to the top level variable names - this may require several
%different fields, each with a different nc_varname (as labelled in the
%NetCDf files). Need to add the option to supply several variable names,
%though. Or could deal with them in this function
% --- N.B. might want to add code in here too :- UM_maps_20141126T041216_FUNC

%% switch for variable name
switch varname
    case 'accum_num'         
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='';
        UM_map.var_units_str = 'kg^{-1}'; 
        UM_map.var_nice_name = 'Accum. aerosol number conc.'; 
        
    case 'accum_mass'         
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='';
        UM_map.var_units_str = 'kg kg^{-1}'; 
        UM_map.var_nice_name = 'Accum. aerosol mass mixing ratio';         
        
    case 'qL'         
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='';
        UM_map.var_units_str = 'kg kg^{-1}'; 
        UM_map.var_nice_name = 'qL';          
        
    case 'SW_up_TOA'
        nc_varname = 'SW_up_TOA'; %The name of the variable in the nc file
        VAR_NAME_STR='SW_up_TOA'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag=''; %set to '' to just load the variable from the .nc file
        % 'load_mat' to laod from the .mat file
        % 'calc' to calcualte e.g. LWP from qL, rho ,etc.
    case 'SW_TOA_outgoing'
        nc_varname = 'SW_TOA_outgoing'; %The name of the variable in the nc file
        VAR_NAME_STR='SW_up_TOA'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag=''; %set to '' to just load the variable from the .nc file
        % 'load_mat' to laod from the .mat file
        % 'calc' to calcualte e.g. LWP from qL, rho ,etc.        
    case 'SW_down_TOA'  %N.B. - added code in here too :- UM_maps_20141126T041216_FUNC
        nc_varname = 'SW_down_TOA'; %The name of the variable in the nc file
        VAR_NAME_STR='SW_down_TOA'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag=''; %set to '' to just load the variable from the .nc file
        % 'load_mat' to laod from the .mat file
        % 'calc' to calcualte e.g. LWP from qL, rho ,etc.
    case 'SW_down_surf_LWP_LT_0pt1' %N.B. - added code in here too :- UM_maps_20141126T041216_FUNC
        nc_varname = 'SW_down_surf_LWP_LT_0pt1'; %The name of the variable in the nc file
        VAR_NAME_STR='SW_down_surf'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag=''; %set to '' to just load the variable from the .nc file
        % 'load_mat' to laod from the .mat file
        % 'calc' to calcualte e.g. LWP from qL, rho ,etc.
    case 'Transmission_down_surf_LWP_LT_0pt1'
        %        nc_varname = 'SW_down_surf'; %The name of the variable in the nc file
        nc_varname = 'Transmission_down_surf_LWP_LT_0pt1'; %The name of the variable in the nc file
        VAR_NAME_STR='SW_down_surf'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag=''; %set to '' to just load the variable from the .nc file (or for some computations).
        % 'load_mat' to laod from the .mat file
        % 'calc' to calcualte e.g. LWP from qL, rho ,etc.
    case 'LW_up_TOA'
        nc_varname = 'LW_up_TOA'; %The name of the variable in the nc file
        VAR_NAME_STR='LW_up_TOA'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='';

    case 'Nd'
        nc_varname = 'Nd'; %The name of the variable in the nc file
        VAR_NAME_STR='Nd'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'cm^{-3}';         

    case {'LWP','LWP_incloud_20gmsq','CF_LWP_20','CF_LWP_150','LWP_incloud_150gmsq'}
        nc_varname = 'LWP'; %The name of the variable in the nc file
        VAR_NAME_STR='LWP'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'g m^{-2}';        
        UM_map.var_nice_name = 'LWP';
    case 'RWP'
        nc_varname = 'RWP'; %The name of the variable in the nc file
        VAR_NAME_STR='RWP'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'g m^{-2}';
        UM_map.var_nice_name = 'RWP';
        
    case 'CF_0pt25_LWP_4km_20'
        nc_varname = 'CF_0pt25_LWP_4km_20'; %The name of the variable in the nc file
        VAR_NAME_STR='CF_0pt25_LWP_4km_20'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        
        
    case 'accum_num_z3000'
        nc_varname = 'accum_z3000'; %The name of the variable in the nc file
        VAR_NAME_STR='accum_num_z3000'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg^{-1}';
          
    case 'accum_mass_mean_up_to_z3000'
        nc_varname = 'accum_mass_up_to_z3000'; %The name of the variable in the nc/mat file
        VAR_NAME_STR='accum_mass_mean_to_z3000'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg kg^{-1}';  
        UM_map.var_nice_name = 'Accum. mode mass MR mean to 0-3km';        
        
    case 'accum_mass_z3000'
        nc_varname = 'accum_mass_z3000'; %The name of the variable in the nc/mat file
        VAR_NAME_STR='accum_mass_z3000'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg kg^{-1}';    
        
    case 'accum_mass_ug_per_m3_mean_to_z3000'
        nc_varname = 'accum_mass_ug_per_m3_mean_to_z3000'; %The name of the variable in the nc/mat file
        VAR_NAME_STR='accum_mass_ug_per_m3_mean_to_z3000'; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = '\mug m^{-3}'; 
        UM_map.var_nice_name = 'Accum. mode mass MR mean to 0-3km';  
        
    case 'accum_mass_total_column_to_z3000'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated mass to 0-3km';   
        
    case 'accum_number_total_column_to_z3000'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated number to 0-3km';  
        
    case 'droplet_number_total_column_to_z3000'    
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Droplet number conc. column integrated number to 0-3km'; 
        
         
    case 'air_mass_total_column_to_z1500' 
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Air mass column integrated number to 0-1.5km';       
        
    case 'droplet_number_total_column_to_z1500' 
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Droplet number conc. column integrated number to 0-1.5km';
        
    case 'Total_number_aerosol_droplets_to_z1500m'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg^{-1}';
        UM_map.var_nice_name = 'Mean total aerosol+droplets number 0-1.5km';
        
     case 'droplet_number_total_column_z1500_to_top'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Droplet number conc. column integrated number 1.5km upwards';       
        
      case 'accum_number_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated number to 0-1.5km';  
        
    case 'accum_number_total_column_z1500_to_top'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated number 1.5km upwards';       
     
   case 'coarse_number_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Coarse mode column integrated number to 0-1.5km';  
        
    case 'coarse_number_total_column_z1500_to_top'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Coarse mode column integrated number 1.5km upwards';                        
        
    case 'accum_mass_column_integrated'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated mass'; 
        
    case 'accum_mass_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated mass to 0-1.5km'; 
        
    case 'coarse_mass_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Coarse mode column integrated mass to 0-1.5km';    
        
     case 'act_mass_liq_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Aerosol in droplets column integrated mass'; 
        
    case 'act_mass_rain_total_column_to_z1500'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'kg m^{-2}'; 
        UM_map.var_nice_name = 'Aerosol in rain column integrated mass';         
        
    case 'accum_number_column_integrated'
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Accum. mode column integrated number';  
        
    case 'droplet_number_column_integrated'   
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        UM_map.var_nice_name = 'Droplet number conc. column integrated number';   
        
   case {'aitken_number_column_integrated','coarse_number_column_integrated'}
       switch varname
           case 'aitken_number_column_integrated'
               UM_map.var_nice_name = 'Aitken mode column integrated number';
           otherwise
               UM_map.var_nice_name = remove_character(varname,' ');                              
       end
        nc_varname = varname; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm^{-2}'; 
        
    case 'BL height using qL'   
        nc_varname = 'zi'; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm'; 
        UM_map.var_nice_name = 'Boundary layer height using qL threshold';          
      
  case 'BL height using RH50'   
        nc_varname = 'zi'; %The name of the variable in the nc/mat file
        VAR_NAME_STR = varname; %The name part in the filename to replace VAR_NAME with
        %OR, if the VAR_NAME is not in the supplied filename then will
        %just use the supplied one - make sure it is correct!
        flag='load_mat';
        UM_map.var_units_str = 'm'; 
        UM_map.var_nice_name = 'Boundary layer height using RH<50 threshold';       
        
    otherwise
        error(['*** Need to add a case for ' varname ' in this script. ***']);
end
%% End of switch varname

if ~isfield(UM_map,'var_nice_name')
    UM_map.var_nice_name = varname;
end



%% Script to get the UM run details by providing the run set name (in
%% UM_cases)
UM_case_select_RUN

UM_case_out.VAR_NAME_STR = VAR_NAME_STR; %save this too as needed for writing the .mat file for the timeseries at the end.


%% Loop through the different UM runs
for idat_UM=1:length(fileUM)

    if iscell(dirUM)==1
        dirUM_i = dirUM{idat_UM};
    else
        dirUM_i = dirUM;
    end

    UM_map.idat_UM = 1; %needs to be set to one when just retrieving one time index at a time  %idat_UM;
    UM_map.var = nc_varname; %can set to the variable name to just read the variable
    %         UM_map.flag = flag{idat_UM};
    UM_map.flag = flag;
    UM_map.file_lwp =  remove_character([dirUM_i fileUM{idat_UM}],'VAR_NAME',VAR_NAME_STR);
    if strcmp(flag,'load_mat')==1 & length(strfind(UM_map.file_lwp(end-3:end),'.mat'))==0
        UM_map.file_lwp = [UM_map.file_lwp '.mat'];
    end
    %         UM_map.file_rho = [dirUM_i fileUM_rho{idat_UM}]; %filename_rho;
    UM_map.pole_lat = pole_lat;
    UM_map.pole_lon = pole_lon;
    
    UM_map.labs_UM = UM_case_out.labs_UM(idat_UM); %set label to be the run in question - UM_case_out comes form UM_case_select_RUN


    if exist('time_range')==1 %==1 means is a variable in the workspace

        UM_map.time_in = []; %set to [] for all times
        %         UM_map.iz = iz;



        % Get the time info only
        switch UM_map.flag
            case 'load_mat'
                load(UM_map.file_lwp,'time_matlab');
                %             case 'calc'
            otherwise
                [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(UM_map.file_lwp,[],pole_lat,pole_lon,time_tol);
        end






        its = find(time_matlab >= time_range(1) & time_matlab <= time_range(2) );
        time_specific = time_matlab(its);

        clear nc

    end %otherwise time_specifc should already be set

    %          UM_map.noplot = 0;
    %          switch UM_map.action
    %              %         case 'save all'
    %              % %            UM_time_out.all_data = NaN*ones();
    %              %             %Leave all the data in
    %              case 'timeseries domain means'
    %                  UM_map.noplot = 1; %prevent the plot from being made (speeds up processing)
    %
    %
    %          end



%% Main loop over time
    for itime=1:length(time_specific)


        %% Plot function
        %if iplot_maps==1
        UM_map.time_in = time_specific(itime);

% ------------ Get the data / do the plot ----------------------------
        UM_time_out_temp = UM_maps_20141126T041216_FUNC(UM_map);
% --------------------------------------------------------------------        

        if iclose_figs==1 & noplot==0
            close(gcf)
        end
        %end

        %     switch UM_map.action
        % %         case 'save all'
        % % %            UM_time_out.all_data = NaN*ones();
        % %             %Leave all the data in
        %         case 'timeseries domain means'
        %             %Strip out the data not needed to save memory
        %             UM_time_out_temp = rmfield(UM_time_out_temp,'P_save');
        %             UM_time_out_temp = rmfield(UM_time_out_temp,'gcm_Plat2D_edges_UM');
        %             UM_time_out_temp = rmfield(UM_time_out_temp,'gcm_Plon2D_edges_UM');
        %             UM_time_out_temp = rmfield(UM_time_out_temp,'gcm_Plat2D_UM');
        %             UM_time_out_temp = rmfield(UM_time_out_temp,'gcm_Plon2D_UM');
        %
        %
        %     end
        
        UM_time_out{idat_UM}.gcm_Plat2D_UM = UM_time_out_temp.gcm_Plat2D_UM;
        UM_time_out{idat_UM}.gcm_Plon2D_UM = UM_time_out_temp.gcm_Plon2D_UM;
        UM_time_out{idat_UM}.gcm_Plat2D_edges_UM = UM_time_out_temp.gcm_Plat2D_edges_UM;
        UM_time_out{idat_UM}.gcm_Plon2D_edges_UM = UM_time_out_temp.gcm_Plon2D_edges_UM;
         
        UM_time_out{idat_UM}.z_um = UM_time_out_temp.z_um;
         
        UM_time_out{idat_UM}.var_units_str = UM_map.var_units_str;       
        UM_time_out{idat_UM}.var_nice_name = UM_map.var_nice_name;  
        
        
        
        %% Switches for specific cases
        iclear_output=0; %flag to say whether to clear dat_UM if don't want to output full field
        switch UM_map.varname
            case 'LWP_incloud_20gmsq'
                thresh_LWP = 20; %LWP threshold to count as in cloud (g/m2)
                lwp_it = UM_time_out_temp.P_save{1};
                lwp_it(lwp_it < thresh_LWP)=NaN;
                UM_time_out{idat_UM}.timeseries(itime) = meanNoNan(lwp_it(:),1);
                UM_time_out{idat_UM}.datUM{itime} = lwp_it;  
                %            VAR_NAME_STR = 'LWP_incloud_20gmsq';
            case 'CF_LWP_20'
                thresh_LWP = 20; %LWP threshold to count as in cloud (g/m2)
                lwp_it = UM_time_out_temp.P_save{1}; %2D field
                lwp_thresh = lwp_it(lwp_it >= thresh_LWP);
                lwp_it(lwp_it < thresh_LWP)=NaN;
                UM_time_out{idat_UM}.timeseries(itime) = length(lwp_thresh(:)) ./ length(lwp_it(:));
                UM_time_out{idat_UM}.datUM{itime} = lwp_it; 
                %            VAR_NAME_STR = 'CF_LWP_20';
            case 'LWP_incloud_150gmsq'
                thresh_LWP = 150; %LWP threshold to count as in cloud (g/m2)
                lwp_it = UM_time_out_temp.P_save{1};
                lwp_thresh = lwp_it(lwp_it >= thresh_LWP);
                UM_time_out{idat_UM}.timeseries(itime) = meanNoNan(lwp_thresh(:),1);
                %            VAR_NAME_STR = 'LWP_incloud_20gmsq';
            case 'CF_LWP_150'
                thresh_LWP = 150; %LWP threshold to count as in cloud (g/m2)
                lwp_it = UM_time_out_temp.P_save{1}; %2D field
                lwp_thresh = lwp_it(lwp_it >= thresh_LWP);
                UM_time_out{idat_UM}.timeseries(itime) = length(lwp_thresh(:)) ./ length(lwp_it(:));
                %            VAR_NAME_STR = 'CF_LWP_20';
            otherwise
                iclear_output=1;
                UM_time_out{idat_UM}.datUM{itime} = UM_time_out_temp.P_save{1};
                UM_time_out{idat_UM}.timeseries(itime) = UM_time_out_temp.P_mean{1};
                switch UM_time_loop_action
                    case {'output all','output_3D'}
                        iclear_output=0;
                end
        end
        
        
        
        switch UM_time_loop_action
                    case 'time_mean_2D'
                        % Add to a running mean based on total number of
                        % times that will be used.
                        if itime==1
                           UM_time_out{idat_UM}.datUM_timemean{1} = 0; 
                           %UM_time_out{idat_UM}.datUM_timemean_Ntime{1} = zeros(size( UM_time_out_temp.P_save{1} ));
                           UM_time_out{idat_UM}.datUM_timemean_Ntime{1} = zeros(size( UM_time_out{idat_UM}.datUM{itime} ));
                        end
                        
                        %increment the locations with data for working out the mean
                        %inot_nan=find(isnan(UM_time_out_temp.P_save{1})==0);   
                        inot_nan=find(isnan(UM_time_out{idat_UM}.datUM{itime})==0);
                        
                        UM_time_out{idat_UM}.datUM_timemean_Ntime{1}(inot_nan) = UM_time_out{idat_UM}.datUM_timemean_Ntime{1}(inot_nan) + 1;
                        
                        %temp = UM_time_out_temp.P_save{1};
                        temp = UM_time_out{idat_UM}.datUM{itime};
                        inan=find(isnan(temp)==1);                        
                        temp(inan)=0;
                        
                        UM_time_out{idat_UM}.datUM_timemean{1} = UM_time_out{idat_UM}.datUM_timemean{1} + temp;
                        
                        
        end
                
        if iclear_output==1
            UM_time_out{idat_UM}.datUM{itime} = 0; %set to zero to save memory/diskspace if not wanted
        end

        end
        %% End of switches for specific variable case



        %    UM_time_out{idat_UM} = UM_time_out_temp;

    
    switch UM_time_loop_action
        case 'time_mean_2D'
            UM_time_out{idat_UM}.datUM_timemean{1} = UM_time_out{idat_UM}.datUM_timemean{1} ./ UM_time_out{idat_UM}.datUM_timemean_Ntime{1};
        otherwise
            UM_time_out{idat_UM}.datUM_timemean{1} = 0;
    end
    UM_time_out{idat_UM}.time = time_specific;
end

UM_time_out{1}.UM_case_out = UM_case_out;

