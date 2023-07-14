% This is for SAVING the transects for later plotting for the data (GCM or satellite) that is
%   loaded. Use *** saveload_plot_MULTI_lon_transect *** (or
%   saveload_lon_transect_data) for LOADING for many GCMs/data
%   sources.
% Make sure that iplot_sat_maps is set to zero if want longitude
%   transects (otherwise will just plot maps)
% Runs case 115 of watervapourMay to give the transects - select the
%   required variables in that below:-
% Runs  -- saveload_lon_transect_data --
% *** NEED to set am3_dataset properly to tell this whether to label as day
%   or night data! For CALIPSO data it is done in the read program since day
%   and night are read in separately. For AMSRE it is done in
%   plot_global_maps. Need to run this script separately to save data for both day and night
% For satellite data also make sure to set times_required2 to only one
%   value

gcm_case_multi={'Using saved GCM data','Using saved GCM data','LWP','Nd','LWP2'};
var_choose_multi={'Precip rate','low CF','','',''};

gcm_case_multi={'LWP2'};
var_choose_multi={''};

gcm_case_multi={'LWP_MOD35CF'}; %MODIS LWP multiplied by the MOD35 CF
var_choose_multi={''};

gcm_case_multi={'LWP_MOD06CF'}; %MODIS LWP multiplied by the MOD06 CF
var_choose_multi={''};

%gcm_case_multi={'Using saved GCM data'};
%var_choose_multi={'Precip rate'};

%gcm_case_multi={'Using saved GCM data'};
%var_choose_multi={'low CF'};

%gcm_case_multi={'Reff_COSP','Tau_COSP','Nd_COSP'};
%gcm_case_multi={'Reff_COSP_allCF'};
%var_choose_multi={'','',''};


%gcm_case_multi={'Reff_COSP_allCF','Reff_COSP','Tau_COSP','Nd_COSP','REFFL_maxlayer','REFFL_max_noCF','REFFL_maxliq'};
%var_choose_multi={'','','','','','',''};

%gcm_case_multi={'Reff_COSP_allCF','Reff_COSP','Tau_COSP','Nd_COSP','REFFL_maxlayer','REFFL_max_noCF','REFFL_maxliq','Using saved GCM data','Using saved GCM data','LWP','Nd','LWP2'};
%var_choose_multi={'','','','','','','','Precip rate','low CF','','',''};

%gcm_case_multi={'Reff_COSP_allCF','Reff_COSP','Tau_COSP','Nd_COSP','Using saved GCM data','Using saved GCM data','LWP','Nd','LWP2'};
%var_choose_multi={'','','','','Precip rate','low CF','','',''};

%gcm_case_multi={'Using saved GCM data','Using saved GCM data','LWP','Nd','LWP2'};
%var_choose_multi={'Precip rate','low CF','','',''};

%gcm_case_multi={'Using saved GCM data','Using saved GCM data','LWP','Nd','LWP2'};
%var_choose_multi={'Precip rate','low CF','','',''};


%gcm_case_multi={'TLWP'}; %grid-box mean LWP+RWP for comparison to AMSRE
%var_choose_multi={''};

%gcm_case_multi={'LWP_COSP'};
%var_choose_multi={''};



%gcm_case_multi={'Reff_COSP'};
%var_choose_multi={''};

%gcm_case_multi={'Re16'};
%var_choose_multi={''};

%gcm_case_multi={'Re21'};
%var_choose_multi={''};

%gcm_case_multi={'Re37'};
%var_choose_multi={''};

%gcm_case_multi={'Re16_minus_Re21'};
%var_choose_multi={''};

%treat the different wavelengths as different models for plotting
%gcm_case_multi={'ReffMODIS'}; gcm_str_select='MODIS_16';
%var_choose_multi={''};

%gcm_case_multi={'ReffMODIS'}; gcm_str_select='MODIS_21';
%var_choose_multi={''};

%gcm_case_multi={'ReffMODIS'}; gcm_str_select='MODIS_37';
%var_choose_multi={''};

%gcm_case_multi={'Re37_minus_Re21'};
%var_choose_multi={''};

%gcm_case_multi={'qv700'}; 
%var_choose_multi={''};

%gcm_case_multi={'LTS'}; 
%var_choose_multi={''};

%gcm_case_multi={'LTS1000'}; 
%var_choose_multi={''};

%gcm_case_multi={'LTS_daily'}; 
%var_choose_multi={''};

%gcm_case_multi={'qv700','LTS'}; 
%var_choose_multi={'',''};

iplot_sat_maps=0; %flag to say we just want to plot maps for satellite data

%gcm_case_multi={'Using saved GCM data'};
%var_choose_multi={'low CF'};



imulti_year_med_data=0; %flag to say whether we are using the muti-year MODIS median/PDF
%data - may be overwritten below if are using model data

clear gcm_strs
gcm_strs{1} = gcm_str_select;  %*** check that these are set correctly by the load_gcm_processed

%             gcm_case = 'LWP';
% %            gcm_case = 'Nd';
%            
% 
% %gcm_case = 'Using saved GCM data';
%             
%             var_choose ='Precip rate'; %chose this for CloudSat too
% %            var_choose =  'low CF';
% %            var_choose =  'high CF';  




time_mean_str2 = {'DJF','MAM','JJA','SON','ALL'};
%time_mean_str2 = {'ALL'};
%time_mean_str2 = {'ASON'};


times_required2 = {[0:24],[0:3],[12:15],[21:23 0:3],[9:15],[0:3 12:15],[9:12]};
times_required2 = {[0:24]}; %just set to [0 24] for MODIS
%times_required2 = {[9:12],[21:23 0]};
%times_required2 = {[0:24],[9:15],[21:23 0:3]};

LAT_val2={[-32.74 -28],[-22.74 -18;],[-12.74 -8;],[-32.74 -8]};
LON_val2={[-104 -71.25],[-104 -71.25],[-104 -71.25],[-104 -71.25]};

%LAT_val2={[-22.74 -18;]};
%LON_val2={[-104 -71.25]};

%LAT_val2={[-32.74 -8]};
%LON_val2={[-104 -71.25]};



if iplot_sat_maps==1
    times_required2 = {[0:24]};
    LAT_val2={[-32.74 -28]};
    LON_val2={[-104 -71.25]};
end



%                times_required = [0:24]; %need to specify all hours (in 1 hour increments)
                 %24 not needed, but 0-24 sounds better for time_UTC_str
%                times_required = [15:21];
%                times_required = [11:16]; 
                
%                times_required = [9:12]; %Terra daytime
%                times_required = [21:23 0]; %Terra nighttime                
                
%                times_required = [12:15];  %Aqua daytime              
%                times_required = [0:3]; %Aqua nighttime
      

%                 times_required = [9:15];  %Aqua/Terra daytime   
%                 times_required = [21:23 0:3];  %Aqua/Terra nighttime  

for ivar_lon = 1:length(gcm_case_multi)
    

for idays=1:length(time_mean_str2)
    for ihours=1:length(times_required2)
        for ilat_multi=1:length(LAT_val2)

            LAT_val = LAT_val2{ilat_multi};  LON_val = LON_val2{ilat_multi};
            time_mean_str=time_mean_str2{idays};
            time_mean_str_orig_multi_transect = time_mean_str;
            times_required=times_required2{ihours};

            ioverride_time_selection=1;
            ioverride_lat_vals=1;
                           
            %default
            switch time_mean_str
                case 'DJF'
                    days_required_for_mean = [336:366 1:60];
                case  'MAM'    
                    days_required_for_mean = [61:152];
                case 'JJA'
                    days_required_for_mean = [153:244]; 
                case 'SON'
                    days_required_for_mean = [245:335]; 
                case 'ASON'
                    days_required_for_mean = [214:335];                     
            end
            
            switch gcm_strs{1}
                case {'MODIS','CALIPSO_monthly','CLOUDSAT_PRECIP','AMSRE','POLDER','ERAInt'}
                    
                    switch gcm_strs{1}
                        case {'MODIS','ERAInt'}
                            am3_dataset ='';
                    end
                    
                    switch gcm_strs{1}
                        case {'CLOUDSAT_PRECIP','AMSRE'}
                            %am3_dataset is used to label as day or night
                            %data and is done elsewhere (see top of this script for
                            %where)
%                            am3_dataset ='';
                            imulti_year_med_data=1;
                    end
                    
                    
                    if imulti_year_med_data==1
                        switch time_mean_str
                            case 'ALL'
                                imonths = [1:12];
                            case 'DJF'
                                imonths = [12 1 2];
                            case  'MAM'
                                imonths = [3 4 5];
                            case 'JJA'
                                imonths = [6 7 8];
                            case 'SON'
                                imonths = [9 10 11];
                        end
                        
                        months_required_for_mean = imonths; %for cloudsat precip
                        
                        %set this so that time_inds_modisL3_timeseries3
                        %works properly (correct month is chosen in
                        %plot_global)

                        switch gcm_strs{1}
                            case 'MODIS'
%                                time_mean_str='ALL';
                        end
             

                    end
                    
                    if ilat_multi==1
                        ioverride_month_select_multi_years=1;
                        plot_global_maps
%                        time_mean_str=''; %to avoid it being concatenated with itself in time_inds_modisL3_timeseries3, as calle din watervap?
                        time_mean_str=time_mean_str2{idays}; %reset this to avoid it being concatenated with itself in time_inds_modisL3_timeseries3, as calle din watervap?
                        saveas_ps_fig_emf(gcf,savename,'',0,0);
                    end
                    
                otherwise
                    imulti_year_med_data=0; %otherwise the time_mean_str string is not set properly
            end
            
            if iplot_sat_maps==0

                noplot=1;
                man_choose_water_graph=1;
                graph=115;
                ioverride_longitude_transect_options=1;
                gcm_case = gcm_case_multi{ivar_lon}; %need to reset this each time, as it gets changed in watervap
                var_choose = var_choose_multi{ivar_lon};
                ichoose_styles=0;

                ioverride_time_selection=1; %override this again since time_inds_modisL3_timeseries3 is run from
                %watervap too

                waterVapourMay2005  %time_mean_str_orig changed in 
                %time_inds_modisL3_timeseries3
                
                if imulti_year_med_data==1
                    time_mean_str = time_mean_str_orig; %swap back to the original
                end
                ioverride_saveload_transect=1;
                saveorload='save';
        % --------------------------------------------                
                saveload_lon_transect_data
        % --------------------------------------------
            end

        end


    end
       
end

end



