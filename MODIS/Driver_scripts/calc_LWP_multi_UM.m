% Note - this loops through time inside the get_LWP_RWP_UM routine (since
% set vars_in.flag to 'calc')


% --- Provide the name of the set of runs to process (see :-
%       UM_case_select_runs
% function, which is run via UM_case_select_RUN)
UM_cases = '12th Nov case, as of May 2016';
UM_cases = '12th Nov case, as of May 2016 adhoc';
%UM_cases = '12th Nov case, as of May 2016 adhoc eos10';
UM_cases = '12th Nov case, as of May 2016 adhoc multi-dirUM';
%UM_cases = 'Iceland_9day_runs_Nov2016'; %All the runs
%UM_cases = 'Iceland_9day_runs_Nov2016 adhoc'; %All the runs
%UM_cases = '12th Nov case, as of May 2016 processing runs multi-dirUM';
%UM_cases = '12th Nov case, as of May 2016 processing runs PLOTS multi-dirUM';
UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM';
%UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM Processing OFF';
UM_cases = '12th Nov case, as of May 2017 processing runs surface fluxes, delaying processing, etc.';

clear var_to_calc_multi
 i=1;
  var_to_calc_multi{i}='LWP'; i=i+1;
  var_to_calc_multi{i}='RWP'; i=i+1;
  var_to_calc_multi{i}='Nd'; i=i+1; %At the location of max LWC and for z<3.5km
%  
 var_to_calc_multi{i}='accum_number_total_column_to_z1500'; i=i+1; % total aerosol number in the column up to 3km (/m2).
 var_to_calc_multi{i}='droplet_number_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).
 var_to_calc_multi{i}='coarse_number_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).
% % 
 var_to_calc_multi{i}='act_mass_liq_total_column_to_z1500'; i=i+1; % total aerosol number in the column up to 3km (/m2).
 var_to_calc_multi{i}='act_mass_rain_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).
% 
var_to_calc_multi{i}='accum_mass_total_column_to_z1500'; i=i+1; % total aerosol number in the column up to 3km (/m2).
var_to_calc_multi{i}='coarse_mass_total_column_to_z1500'; i=i+1; % total droplet number in the column up to 3km (/m2).

var_to_calc_multi{i}='air_mass_total_column_to_z1500'; i=i+1; % total mass of air in the column up to xkm (kg/m2). Allows calc of mean values from the integrated ones.


var_to_calc_multi{i}='Total_number_aerosol_droplets_to_z1500m'; i=i+1;
var_to_calc_multi{i}='BL height using qL'; i=i+1;
var_to_calc_multi{i}='BL height using RH50'; i=i+1;

%  Also calculates W at these locations.
%var_to_calc_multi{i}='dBZ'; i=i+1;
%var_to_calc_multi{i}='Nd_GCM'; i=i+1; %At the location of max LWC and for z<3.5km
%var_to_calc_multi{i}='CF_0pt25_LWP_4km_20'; i=i+1;
%var_to_calc_multi{i}='LWP_incloud_20gmsq'; i=i+1;   % -- N.B. wont' need to use this one, unless wnat to save the 2D
%fields just for LWP>20 (but prob will just want to process to timeseries).
%var_to_calc_multi{i}='accum_num_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_mass_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_mass_mean_to_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_mass_ug_per_m3_mean_to_z3000'; i=i+1;
%var_to_calc_multi{i}='accum_mass_total_column_to_z3000'; i=i+1; % total aerosol mass in the column up to 3km (kg/m2). Like LWP, but for aerosol
%var_to_calc_multi{i}='accum_number_total_column_to_z3000'; i=i+1; % total aerosol number in the column up to 3km (/m2).
%var_to_calc_multi{i}='droplet_number_total_column_to_z3000'; i=i+1; % total droplet number in the column up to 3km (/m2).


%var_to_calc_multi{i}='accum_number_total_column_z1500_to_top'; i=i+1; % total aerosol number in the column up to 3km (/m2).
%var_to_calc_multi{i}='droplet_number_total_column_z1500_to_top'; i=i+1; % total droplet number in the column up to 3km (/m2).
%var_to_calc_multi{i}='coarse_number_total_column_z1500_to_top'; i=i+1; % total droplet number in the column up to 3km (/m2).



%var_to_calc_multi{i}='accum_mass_total_column_z1500_to_top'; i=i+1; % total aerosol number in the column up to 3km (/m2).
%var_to_calc_multi{i}='coarse_mass_total_column_z1500_to_top'; i=i+1; % total droplet number in the column up to 3km (/m2).



%var_to_calc_multi{i}='accum_mass_column_integrated'; i=i+1; % total aerosol mass in the column (kg/m2). Like LWP, but for aerosol

%var_to_calc_multi{i}='accum_number_column_integrated'; i=i+1; % total aerosol number in the column (/m2).
%var_to_calc_multi{i}='droplet_number_column_integrated'; i=i+1; % total droplet number in the column (/m2).
%var_to_calc_multi{i}='aitken_number_column_integrated'; i=i+1; % total aerosol number in the column (/m2).
%var_to_calc_multi{i}='coarse_number_column_integrated'; i=i+1; % total aerosol number in the column (/m2).


%var_to_calc_multi{i}='U_wind_10m';
%var_to_calc_multi{i}='V_wind_10m';
%var_to_calc_multi{i}='W_wind_10m';


for ivar_multi=1:length(var_to_calc_multi)

    var_to_calc = var_to_calc_multi{ivar_multi};


    vars_in.Nmulti_out = 1; %default
    vars_in.z_accum = 1e99; %default height up to which to integrate quantities (e.g. LWP)  

    clear fileUM xdat_import ydat_import flag


    %% Script to get the UM run details by providing the run set name
    %% Provide the case in UM_case_select_runs
    UM_case_select_RUN  %runs UM_case_select_runs


    for idat_UM=1:length(fileUM)
        if iscell(dirUM)==1
            dirUM_i = dirUM{idat_UM};
        else
            dirUM_i = dirUM;
        end

        filename = [dirUM_i fileUM{idat_UM}];
        filename_qL = remove_character(filename,'VAR_NAME','qL');
        filename_rho = remove_character(filename,'VAR_NAME','rho');
        %    filename_rho = [dirUM_i fileUM_rho{idat_UM}];

        filename_save_var = [remove_character(filename,'VAR_NAME',var_to_calc) '.mat'];

        switch var_to_calc
            case 'LWP'

                vars_in.var = 'LWP';
                vars_in.flag = 'calc';
                vars_in.file_qL = filename_qL;
                vars_in.file_lwp = filename;  %Needs thi setting for some reason??
                vars_in.file_rho = filename_rho;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                % Added this 14/07/15 - older uses of this will have saved as kg/m2
                lwp = lwp*1e3;
                %    rwp = rwp*1e3;

                %    save(filename_save,'lwp','rwp','-V7.3');
                save(filename_save_var,'lwp','-V7.3');
                %    save(filename_save_RWP,'rwp','-V7.3');

            case 'RWP'

                vars_in.var = 'RWP';
                vars_in.flag = 'calc';
                %            vars_in.file_qL = filename_qL;
                vars_in.file_qR = remove_character(filename,'VAR_NAME','qR'); %
                vars_in.file_lwp = filename;  %Needs thi setting for some reason??
                vars_in.file_rho = filename_rho;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                % Added this 14/07/15 - older uses of this will have saved as kg/m2
                rwp = rwp*1e3;
                %    rwp = rwp*1e3;

                %    save(filename_save,'lwp','rwp','-V7.3');
                save(filename_save_var,'rwp','-V7.3');
                %    save(filename_save_RWP,'rwp','-V7.3');


            case 'Nd'


                clear vars_in

                vars_in.var = 'Nd max_LWC';  %N.B. - this is also for z<3km
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.file_qL = filename_qL;
                vars_in.file_rho = filename_rho;
                vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                %            vars_in.file_W = remove_character(filename,'VAR_NAME','UVW_component_3D'); %

                vars_in.Nmulti_out = 4;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                Nd = var_out{1};
                W = var_out{2};
                wmax = var_out{3};
                %               wind_speed = var_out{4};    %Can't do this yet, since need to re-grid the
                % staggered U grid - could do this with Iris on
                % postproc and write out?

                %Outputs in #/cm3
                save(filename_save_var,'Nd','W','wmax','-V7.3');

            case 'dBZ'

                vars_in.var = 'Radar dBZ';  %N.B. - this is also for z<3km
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','qR'); %
                vars_in.file_qR = vars_in.file_lwp;
                vars_in.file_rho = filename_rho;
                vars_in.file_NR = remove_character(filename,'VAR_NAME','NR'); %
                vars_in.file_qL = remove_character(filename,'VAR_NAME','qL'); %
                vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [dBZ,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                %                Nd = var_out{1};
                %                W = var_out{2};
                %                wmax = var_out{3};
                % %               wind_speed = var_out{4};    %Can't do this yet, since need to re-grid the
                %                          % staggered U grid - could do this with Iris on
                %                          % postproc and write out?

                %Outputs
                save(filename_save_var,'dBZ','-V7.3');
                
            case 'BL height using RH50'

                vars_in.var = 'BL height using RH50';  %N.B. - this is also for z<3km
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','qv'); %
                vars_in.file_qv = vars_in.file_lwp;
                vars_in.file_rho = filename_rho;
                vars_in.file_th = remove_character(filename,'VAR_NAME','th'); %
                vars_in.file_P = remove_character(filename,'VAR_NAME','P'); %
                

                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [zi,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                %Outputs
                save(filename_save_var,'zi','-V7.3');    
                
                
            case 'BL height using qL'

                vars_in.var = 'BL height using qL';  %N.B. - this is also for z<3km
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','qL'); %
                vars_in.file_qL = vars_in.file_lwp;
                vars_in.file_rho = filename_rho;

                

                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [zi,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                %Outputs
                save(filename_save_var,'zi','-V7.3');
                
            case 'Total_number_aerosol_droplets_to_z1500m' %N.B. - will use pre-existing 2D fields of integrated aerosols and air mass column

                filename_timser = [filename '.mat'];
                
                vars_in.var = 'Total_number_aerosol_droplets';
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename_timser,'VAR_NAME','air_mass_total_column_to_z1500'); %
                vars_in.file_airmass = vars_in.file_lwp;
                vars_in.file_rho = filename_rho;
                vars_in.file_accum = remove_character(filename_timser,'VAR_NAME','accum_number_total_column_to_z1500'); %
                vars_in.file_coarse = remove_character(filename_timser,'VAR_NAME','coarse_number_total_column_to_z1500'); %
                vars_in.file_droplets = remove_character(filename_timser,'VAR_NAME','droplet_number_total_column_to_z1500'); %
                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [Total_number_aerosol_droplets_to_z1500m,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                %                Nd = var_out{1};
                %                W = var_out{2};
                %                wmax = var_out{3};
                % %               wind_speed = var_out{4};    %Can't do this yet, since need to re-grid the
                %                          % staggered U grid - could do this with Iris on
                %                          % postproc and write out?

                %Outputs
                save(filename_save_var,'Total_number_aerosol_droplets_to_z1500m','-V7.3');                


            case 'Nd_GCM'

                filename_rho = [remove_character(filename,'Nd_times_LYR','Nd_times_LYR')];
                %            filename_save_var = [remove_character(filename,'VAR_NAME','Nd_no_min_qL') '.mat'];
                filename_save_var = [remove_character(filename,'Nd_times_LYR','Nd_GCM') '.mat'];

                clear vars_in

                %% Nd_times_LYR is already in cm3 !!

                vars_in.var = 'Nd_GCM';  %N.B. - this is also for z<3.2km
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'',''); %
                vars_in.file_Nd_times_LYR = remove_character(filename,'',''); %
                vars_in.file_LYR_weight = remove_character(filename,'Nd_times_LYR','LYR_CLD_WEIGHT'); %
                %            vars_in.file_P = remove_character(filename,'Nd_times_LYR','P'); %
                %            vars_in.file_potemp = remove_character(filename,'Nd_times_LYR','potemp'); %%
                %            vars_in.file_qL = filename_qL;
                vars_in.file_rho = filename_rho;
                %            vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                %            vars_in.file_W = remove_character(filename,'VAR_NAME','UVW_component_3D'); %

                vars_in.Nmulti_out = 3;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                Nd_thresh_weight_0pt001 = var_out{1};
                Nd_zweight = var_out{2};
                Nd_thresh_weight_0pt004 = var_out{3};
                %               W = var_out{2};
                %               wmax = var_out{3};
                %               wind_speed = var_out{4};    %Can't do this yet, since need to re-grid the
                % staggered U grid - could do this with Iris on
                % postproc and write out?

                %Outputs in #/cm3
                save(filename_save_var,'Nd_thresh_weight_0pt001','Nd_zweight','Nd_thresh_weight_0pt004','-V7.3');


            case 'CF_0pt25_LWP_4km_20_OLD'
                thresh_LWP_driver=20; %g/m2
                target_dlat = 0.25;  %area over which to calculate the CF
                target_dlon = 0.25;
                % Currently corase graining LWP to GOES (4km) resolution first.

                vars_in.var = 'lwp'; %name of variable that will read from .mat file
                vars_in.flag = 'load_mat';
                vars_in.file_lwp = [remove_character(filename,'VAR_NAME','LWP') '.mat']; %
                %            vars_in.file_qR = vars_in.file_lwp;
                %            vars_in.file_rho = filename_rho;
                %            vars_in.file_NR = remove_character(filename,'VAR_NAME','NR'); %
                %            vars_in.file_qL = remove_character(filename,'VAR_NAME','qL'); %
                %            vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_ed
                %    ges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

                %--- script to load in the GOES dlat and dlon (dlat_GOES etc.)
                GOES_load_lat_lon

                nt = size(lwp,1);
                for it=1:nt

                    % Coarse grain the LWP data to GOES resolution (4km) for comparison
                    % to that
                    dlat_target = dlat_GOES;
                    dlon_target = dlon_GOES;
                    [lwp_coarse,gcm_Plat2D_UM2,gcm_Plon2D_UM2,gcm_Plat2D_edges_UM2,gcm_Plon2D_edges_UM2] = coarse_grain(squeeze(lwp(it,:,:)),gcm_Plat2D_UM,gcm_Plon2D_UM,dlat_target,dlon_target);

                    % --- Now calc the cloud fraction

                    %Work out the number of points to average over to get 0.25x0.25 degree
                    %resolution

                    d=diff(gcm_Plat2D_UM2,[],1);
                    dlat_UM = meanNoNan(meanNoNan(d,1),1);
                    N = ceil(abs(target_dlat/dlat_UM));

                    d=diff(gcm_Plon2D_UM2,[],2);
                    dlon_UM = meanNoNan(meanNoNan(d,1),1);
                    M = ceil(abs(target_dlon/dlon_UM));

                    %Cloud fraction will be defined as fraction of points within each N*M
                    %box with an LWP greater than a threshold. Total no. points =N*M
                    %Make an array of ones and make the ones that we don't want to
                    %count zero
                    Nlwp = zeros(size(lwp_coarse));
                    Nlwp(lwp_coarse>=thresh_LWP_driver)=1;
                    %Now coarse grain (avearge) the Nlwp array - this will now be our
                    %cloud fraction
                    cf_GOES_0pt25 = reduce_matrix_subsample_mean(Nlwp,N,M);


                    if it==1
                        CF_0pt25_LWP_4km_20 = NaN*ones([nt size(cf_GOES_0pt25)]);
                    end

                    CF_0pt25_LWP_4km_20(it,:,:) = cf_GOES_0pt25;
                end
                %Outputs
                save(filename_save_var,'CF_0pt25_LWP_4km_20','-V7.3');

                gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM2,N,M);
                gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM2,N,M);
                gcm_Plat2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plat2D_edges_UM2,N,M);
                gcm_Plon2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plon2D_edges_UM2,N,M);

            case 'CF_0pt25_LWP_4km_20'

                vars_in.var = 'lwp'; %name of variable that will read from .mat file
                vars_in.flag = 'load_mat';
                vars_in.file_lwp = [remove_character(filename,'VAR_NAME','LWP') '.mat']; %
                %            vars_in.file_qR = vars_in.file_lwp;
                %            vars_in.file_rho = filename_rho;
                %            vars_in.file_NR = remove_character(filename,'VAR_NAME','NR'); %
                %            vars_in.file_qL = remove_character(filename,'VAR_NAME','qL'); %
                %            vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_ed
                %    ges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                clear CF_in
                % --- Options for calculating the the cloud fraction

                % Coarse grain the LWP data to GOES resolution (4km) for comparison
                % to that
                %--- script to load in the GOES dlat and dlon (dlat_GOES etc.)
                GOES_load_lat_lon

                % ---Options for coarse graining
                CF_in.icoarsen=1;
                CF_in.dlat_target_coarsen = dlat_GOES; %not needed if not coarsening
                CF_in.dlon_target_coarsen = dlon_GOES;

                % -- Options for CF calculation
                CF_in.thresh_LWP_driver=20; %g/m2 - threshold to define cloud
                CF_in.target_dlat = 0.25;  %area over which to calculate the CF
                CF_in.target_dlon = 0.25;

                %--- Lat and lon grid
                CF_in.lat = gcm_Plat2D_UM;
                CF_in.lon = gcm_Plon2D_UM
                CF_in.lat_edges = gcm_Plat2D_edges_UM;
                CF_in.lon_edges = gcm_Plon2D_edges_UM;

                %--- data for CF_in set in the time loop below
                %CF_in.dat = <set below>

                nt = size(lwp,1);
                for it=1:nt

                    CF_in.dat = squeeze(lwp(it,:,:));
                    CF_out = UM_calc_coarsen_calc_CF(CF_in);


                    if it==1
                        CF_0pt25_LWP_4km_20 = NaN*ones([nt size(CF_out.cf)]);
                    end

                    CF_0pt25_LWP_4km_20(it,:,:) = CF_out.cf;
                end  %time loop
                %Outputs
                save(filename_save_var,'CF_0pt25_LWP_4km_20','-V7.3');

                gcm_Plat2D_UM = CF_out.lat;
                gcm_Plon2D_UM = CF_out.lon;
                gcm_Plat2D_edges_UM = CF_out.lat_edges;
                gcm_Plon2D_edges_UM =  CF_out.lon_edges;

                %            gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM2,N,M);
                %            gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM2,N,M);
                %            gcm_Plat2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plat2D_edges_UM2,N,M);
                %            gcm_Plon2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plon2D_edges_UM2,N,M);

            case 'LWP_incloud_20gmsq'

                % In-cloud LWP based on a cloud mask using LWP>20 gmsq
                thresh_LWP=20;

                vars_in.var = 'lwp'; %name of variable that will read from .mat file
                vars_in.flag = 'load_mat';
                vars_in.file_lwp = [remove_character(filename,'VAR_NAME','LWP') '.mat']; %
                %            vars_in.file_qR = vars_in.file_lwp;
                %            vars_in.file_rho = filename_rho;
                %            vars_in.file_NR = remove_character(filename,'VAR_NAME','NR'); %
                %            vars_in.file_qL = remove_character(filename,'VAR_NAME','qL'); %
                %            vars_in.file_Nd = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.Nmulti_out = 1;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;
                %    vars_in.time_in = time_select;


                %    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                %    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_ed
                %    ges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
                [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                %--- data for CF_in set in the time loop below
                %CF_in.dat = <set below>

                nt = size(lwp,1);
                for it=1:nt
                    fprintf(1,'\nit=%d, ',it);

                    lwp_it = squeeze(lwp(it,:,:));
                    lwp_thresh = lwp_it(lwp_it>thresh_LWP);

                    if it==1
                        LWP_incloud_20gmsq= NaN*ones([nt size(lwp_it)]);
                    end

                    LWP_incloud_20gmsq(it,:,:) = meanNoNan(lwp_thresh(:),1);
                end  %time loop
                %Outputs
                save(filename_save_var,'LWP_incloud_20gmsq','-V7.3');

                %            gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM2,N,M);
                %            gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM2,N,M);
                %            gcm_Plat2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plat2D_edges_UM2,N,M);
                %            gcm_Plon2D_edges_UM = reduce_matrix_subsample_mean(gcm_Plon2D_edges_UM2,N,M);


            case 'accum_num_z3000'

                clear vars_in

                vars_in.var = 'accum_num_at_z';  %
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_num'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

                accum_z3000 = var_out;

                save(filename_save_var,'accum_z3000','-V7.3');

            case 'accum_mass_z3000'

                clear vars_in

                vars_in.var = 'accum_mass_at_z';  %
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_mass'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

                accum_mass_z3000 = var_out;

                save(filename_save_var,'accum_mass_z3000','-V7.3');


            case 'accum_mass_mean_to_z3000'

                clear vars_in

                vars_in.var = 'accum_mass_mean_up_to_z';  %
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_mass'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                accum_mass_up_to_z3000 = var_out;

                save(filename_save_var,'accum_mass_up_to_z3000','-V7.3');

            case 'accum_mass_ug_per_m3_mean_to_z3000'

                clear vars_in

                vars_in.var = 'accum_mass_ug_per_m3_mean_up_to_z';
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_mass'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                accum_mass_ug_per_m3_mean_to_z3000  = var_out;

                save(filename_save_var,'accum_mass_ug_per_m3_mean_to_z3000','-V7.3');

            case 'accum_mass_total_column_to_z3000'

                clear vars_in

                vars_in.var = 'accum_mass_total_column_to_z3000'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_mass'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');

            case 'accum_number_total_column_to_z3000'

                clear vars_in

                vars_in.var = 'accum_number_total_column_to_z3000'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_num'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');
                
                
    
                
                 

            case 'droplet_number_total_column_to_z3000'

                clear vars_in

                vars_in.var = 'droplet_number_total_column_to_z3000'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 3000;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');

                

            case 'accum_mass_column_integrated'

                clear vars_in

                vars_in.var = 'accum_mass_column_integrated'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_mass'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 1e99;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');

            case 'accum_number_column_integrated'

                clear vars_in

                vars_in.var = 'accum_number_column_integrated'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','accum_num'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 1e99;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');

            case 'droplet_number_column_integrated'

                clear vars_in

                vars_in.var = 'droplet_number_column_integrated'
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME','Nd'); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;
                vars_in.z_accum = 1e99;
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');


                % Generic case for 2D plane at a specific height
            case {'U_wind_10m','V_wind_10m','W_wind_10m'}
                clear vars_in
                vars_in.z_accum = [-1]; %default to first available height
                FILE_NAME = VAR_NAME; %default - change if filename is different from variabile name
                
                %Just need to set VAR_NAME to the sname of the variable in the
                %filename (which should also be the same as the var name in
                %the .nc file)
                switch var_to_calc
                    case 'U_wind_10m'
                        VAR_NAME = 'UVW_components_3D';
                        vars_in.z_plane = [10];
                    case 'V_wind_10m'
                        VAR_NAME = 'y-wind';
                        FILE_NAME = 'UVW_components_3D';
                        vars_in.z_plane = [10];
                        vars_in.istaggered = 1;
                    case 'W_wind_10m'
                        VAR_NAME = 'dz_dt';
                        FILE_NAME = 'UVW_components_3D';
                        vars_in.z_plane = [10];
                    otherwise
                        error('Also need to set a switch for the case here!');
                end

                

                vars_in.var = 'generic_horiz_plane'; %generic name
                vars_in.VAR_NAME = VAR_NAME;
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME',FILE_NAME); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;                
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz,UM_info] = get_LWP_RWP_UM(vars_in);

                z_level = z_um(iz);
                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'z_level','z_um','iz','UM_info','-V7.3');

                % Generic case for all column integrated quantities
            case {'aitken_number_column_integrated','coarse_number_column_integrated','coarse_number_total_column_to_z1500','accum_number_total_column_to_z1500','droplet_number_total_column_to_z1500'...
                    'accum_number_total_column_z1500_to_top','coarse_number_total_column_z1500_to_top','droplet_number_total_column_z1500_to_top',...
                    'accum_mass_total_column_to_z1500','droplet_mass_total_column_to_z1500','coarse_mass_total_column_to_z1500','act_mass_liq_total_column_to_z1500','act_mass_rain_total_column_to_z1500',...
                    'air_mass_total_column_to_z1500'}
                clear vars_in
                vars_in.z_accum = [-1 1e99]; %default to all heights
                
                %Just need to set VAR_NAME to the sname of the variable in the
                %filename (which should also be the same as the var name in
                %the .nc file)
                switch var_to_calc
                    case 'aitken_number_column_integrated'
                        VAR_NAME = 'aitken_num';
                    case 'coarse_number_column_integrated'
                        VAR_NAME = 'coarse_num';
                    case 'air_mass_total_column_to_z1500'
                        VAR_NAME = 'air_mass'; %added a special case for this - just need to treat like with an MR field of 1 kg/kg
                        vars_in.z_accum = [-1 1500];                        
                    case 'coarse_number_total_column_to_z1500'
                        VAR_NAME = 'coarse_num';
                        vars_in.z_accum = [-1 1500];
                    case 'accum_number_total_column_to_z1500'
                        VAR_NAME = 'accum_num';
                        vars_in.z_accum = [-1 1500];                        
                    case 'droplet_number_total_column_to_z1500'
                        VAR_NAME = 'Nd';
                        vars_in.z_accum = [-1 1500];    
                    case 'accum_mass_total_column_to_z1500'
                        VAR_NAME = 'accum_mass';
                        vars_in.z_accum = [-1 1500];                        
                    case 'coarse_mass_total_column_to_z1500'
                        VAR_NAME = 'coarse_mass';
                        vars_in.z_accum = [-1 1500];     
                    case 'act_mass_liq_total_column_to_z1500'
                         VAR_NAME = 'act_mass_in_liq';
                        vars_in.z_accum = [-1 1500]; 
                    case 'act_mass_rain_total_column_to_z1500'    
                         VAR_NAME = 'act_mass_in_rain';
                        vars_in.z_accum = [-1 1500];                         
                    case 'accum_number_total_column_z1500_to_top'
                        VAR_NAME = 'accum_num';
                        vars_in.z_accum = [1500 1e9];
                    case 'droplet_number_total_column_z1500_to_top'
                        VAR_NAME = 'Nd';
                        vars_in.z_accum = [1500 1e9];
                    case 'coarse_number_total_column_z1500_to_top'
                        VAR_NAME = 'coarse_num';
                        vars_in.z_accum = [1500 1e9];
                    otherwise
                        error('Also need to set a switch for the case here!');
                end

                

                vars_in.var = 'generic_column_integrated'; %generic name
                vars_in.VAR_NAME = VAR_NAME;
                vars_in.flag = 'calc';
                vars_in.file_lwp = remove_character(filename,'VAR_NAME',VAR_NAME); %
                vars_in.file_rho = filename_rho;

                vars_in.Nmulti_out = 1;                
                vars_in.pole_lat = pole_lat;
                vars_in.pole_lon = pole_lon;

                [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);


                eval([var_to_calc '  = var_out;']);
                save(filename_save_var,var_to_calc,'-V7.3');
                
            otherwise
                
                error('Need to set a case here!')

        end %switch var_to_calc



        clear var_list
        i=1;
        var_list{i} = 'time_matlab'; i=i+1;
        var_list{i} = 'gcm_Plat2D_UM'; i=i+1;
        var_list{i} = 'gcm_Plon2D_UM'; i=i+1;
        var_list{i} = 'gcm_Plat2D_edges_UM'; i=i+1;
        var_list{i} = 'gcm_Plon2D_edges_UM'; i=i+1;
        var_list{i} = 'it'; i=i+1;
        var_list{i} = 'daynum_timeseries3_UM'; i=i+1;
        var_list{i} = 'modisyear_timeseries3_UM'; i=i+1;
        var_list{i} = 'gcm_time_UTC_UM'; i=i+1;
        var_list{i} = 'gcm_time_matlab_UM'; i=i+1;

        for i=1:length(var_list)
            eval(['save(''' filename_save_var ''',''' var_list{i} ''',''-V7.3'',''-APPEND'');']);
        end


        %    save(filename_save,'time','time_matlab','lat','lon','-V7.3','-APPEND');

    end

end






















