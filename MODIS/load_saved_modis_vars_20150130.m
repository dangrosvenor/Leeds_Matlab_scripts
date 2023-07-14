% Loads in daily L3 (standard directly from MODIS) variables.
% Also creates Nd, etc. using filtering_data_get script.

try
%  *** N.B. - for the multiple year SZA Artic box data use MODIS_multi_day_processL3L2 ***

% modis_select_timeseries3_files  -- select what years of data to load here based
%                                    upon data_type (choose below)
% filename_choose_saved_MODIS_files -- where the filenames are listed
%    (this is run from load_saved_modis_vars_func)
% For NORMAL timeseries3 data the variables are chosen BELOW
%     ,but for mock L3 choose in  --- make_mockL3_variables  ---

irestrict_latlon = 0;
lat_restrict = [62 63]; lon_restrict = [27 28]; %Puijo_Sami : 62.91N, 27.66 E
lat_restrict = [-50 20]; lon_restrict = [-160 -60]; %VOCALS
%lat_restrict = [0 25]; lon_restrict = [-180 -150]; %Amy's Hawaii study
%lat_restrict = [-60 -40]; lon_restrict = [-180 180]; %Southern Ocean

lat_restrict = [-40 0]; lon_restrict = [-140 -68]; %VOCALS to match mock L3 data

lat_restrict = [40 90]; lon_restrict = [-40 10]; %For Iceland volcano (Jan 2015, Florent)
lat_restrict = [40 90]; lon_restrict = [-180 180]; %For Iceland volcano (Jan 2015, Florent)
lat_restrict = [0 90]; lon_restrict = [-180 180]; %For Iceland volcano (Jan 2015, Florent) NH only

%lat_restrict = [-22 -13]; lon_restrict = [-130 -121]; %South Pacific Koren region for JJA 2007 (March, 2015)

%to restrict the latitude long
%MLAT_standard = [-89.5:1:89.5];
%MLON_standard = [-179.5:1:179.5];
%MLAT_standard_edges = [-90:1:90];
%MLON_standard_edges = [-180:1:180];

MLAT_standard = [89.5:-1:-89.5];
MLON_standard = [-179.5:1:179.5];
MLAT_standard_edges = [90:-1:-90];
MLON_standard_edges = [-180:1:180];

thresh_LAT = lat_restrict;
thresh_LON = lon_restrict; %Needed in plot_global_maps.






%N.B. - the standard ordering for MODIS appears to become (at least for my
%files) that the lat dimension is ordered from 89.5 to -89.5, whereas LON
%runs from -179.5 to 179.5
%So, below use > min(lat) and <= max(lat) for lat restriction and <= & >
%for lon...
    ilat_restrict = find(MLAT_standard_edges>min(lat_restrict) & MLAT_standard_edges<=max(lat_restrict));
    ilon_restrict = find(MLON_standard_edges>=min(lon_restrict) & MLON_standard_edges<max(lon_restrict));

%make some lat and long for saving with AMSRE
 lat_AMSRE_time3 = MLAT_standard(ilat_restrict);
 lat_AMSRE_time3_edges = [MLAT_standard_edges(ilat_restrict) MLAT_standard_edges(ilat_restrict(end))-1];   
            %is subtract one for lon because it goes from north to south
 lon_AMSRE_time3 = MLON_standard(ilon_restrict);
 lon_AMSRE_time3_edges = [MLON_standard_edges(ilon_restrict) MLON_standard_edges(ilon_restrict(end))+1];   
    
 [Plon2D_AMSRE_time3, Plat2D_AMSRE_time3] = meshgrid(lon_AMSRE_time3,lat_AMSRE_time3);
 [Plon2D_AMSRE_time3_edges, Plat2D_AMSRE_time3_edges] = meshgrid(lon_AMSRE_time3_edges,lat_AMSRE_time3_edges);
 
%    

% *** select data type here ***
%override_loadsave=0;
if ~exist('override_loadsave') | override_loadsave==0
    data_type = 'L3 processed data';
    %data_type = 'Nd vs SZA profiles at constant LAT';
%    data_type = 'Monthly averages etc from timeseries3'; %i.e. from modis_make_monthly_averages_multi_year
    %this is the data from multiple years saved as PDFs
%    data_type = 'Mock L3 timeseries3'; 
    i_multi_year_av=0;
    savemem=0;
end

switch data_type
    case 'Mock L3 timeseries3'
        savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L2/';
    otherwise
        %savedir_var='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L3/';
        savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L3/'; %also works for C6.1 data.

end


%N.B. - choose the variables to load in below


comp='UWchallenger';

%loads multiple processed MODIS datasets and strings them together
%just leave one uncommented to load in only one year
clear modis_data_case ; imod=1;

if ~exist('i_multi_year_av')
    i_multi_year_av=0;
end







%flag to do special averaging as a function of latitude
high_lat_correction=0;

now_str=datestr(now,30);

switch data_type
    case 'Monthly averages etc from timeseries3'
% *** external script to select the required files :-
    modis_select_timeseries3_files
% ****************************************************  
    modis_make_monthly_avs_saveload_variables
    
    case 'Mock L3 timeseries3'
        % *** external script to select the required savefiles to load :-
             modis_select_timeseries3_files
        % ****************************************************
        %puts the variables that we want to load in modis_var cell array
        %This is the same as is used to determine what to save when making
        %mock L3 in MODIS_process_multiple_L2_files.m
        save_or_load = 'load';
        make_mockL3_variables
        


    case 'L3 processed data'
% *** external script to select the required files :-
    modis_select_timeseries3_files
% ****************************************************    


% data_required = 'ALL'; %load all the data
% data_required = 'SELCETED';
% 
% %flag to say whether to read in the core data - otherwise just read in some
% %extra fields
% load_core='no';
% load_core='yes';



%below will be 'core' data - data that will prob need each time, so just
%load once
% cf_time3=Cloud_Fraction_Liquid.timeseries3(ilat,ilon,itime);
% NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3(ilat,ilon,itime)./cf_time3;
% sensZA_time3 = Sensor_Zenith_Mean.timeseries3(ilat,ilon,itime);
% max_sensZA_time3 = Sensor_Zenith_Maximum.timeseries3(ilat,ilon,itime);
% maxSZA_time3 = Solar_Zenith_Maximum.timeseries3(ilat,ilon,itime);
% minSZA_time3 = Solar_Zenith_Minimum.timeseries3(ilat,ilon,itime);
% 
% tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilat,ilon,itime);
% reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3(ilat,ilon,itime)*1e-6; %convert to metres


clear modis_var field_string
istring=1;


switch savemem
    case -1  %Use a specified script
        eval(modis_vars_script_name);
    case 1
        %MOD06 Cloud Fractions
modis_var{istring}='Cloud_Fraction_Liquid'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
modis_var{istring}='Cloud_Fraction_Liquid_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=21)  %N.B. this is exactly the same as totN
modis_var{istring}='Cloud_Fraction_Combined'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=440)

modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)

modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=337)

%modis_var{istring}='Sensor_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 


%clear modis_var; istring=1;
modis_var{istring}='Cloud_Top_Temperature_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Cloud_Top_Temperature_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 

%modis_var{istring}='Cloud_Top_Pressure_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
        
modis_var{istring}='Solar_Zenith_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 

case 0


%MOD06 Cloud Fractions
modis_var{istring}='Cloud_Fraction_Liquid'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
modis_var{istring}='Cloud_Fraction_Liquid_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=21)  %N.B. this is exactly the same as totN
%modis_var{istring}='Cloud_Fraction_Ice';
%field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=436)
%modis_var{istring}='Cloud_Fraction_Undetermined'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=438)
modis_var{istring}='Cloud_Fraction_Combined'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=440)

%MOD35 Cloud Fractions (MODIS Cloud mask - is less conservative and does
%not remove the cloud edges)
modis_var{istring}='Cloud_Fraction_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
%modis_var{istring}='Cloud_Fraction_Night_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
%modis_var{istring}='Cloud_Fraction_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
%modis_var{istring}='Cloud_Fraction_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
%modis_var{istring}='Cloud_Fraction_Night_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
%modis_var{istring}='Cloud_Fraction_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);

modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Effective_Radius_Liquid_QA_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %
modis_var{istring}='Cloud_Effective_Radius_Liquid_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Effective_Radius_Liquid_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Effective_Radius_Liquid_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Effective_Radius_Liquid_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)

modis_var{istring}='Cloud_Effective_Radius_37_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)

modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=337)
%modis_var{istring}='Cloud_Optical_Thickness_Liquid_QA_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; %
modis_var{istring}='Cloud_Optical_Thickness_Liquid_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Optical_Thickness_Liquid_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Optical_Thickness_Liquid_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)
%modis_var{istring}='Cloud_Optical_Thickness_Liquid_Mean_Uncertainty'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=397)


%modis_var{istring}='Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=352)


   
%clear modis_var; istring=1;
 modis_var{istring}='Solar_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
 modis_var{istring}='Solar_Zenith_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
 modis_var{istring}='Solar_Zenith_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
 modis_var{istring}='Solar_Zenith_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Solar_Azimuth_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
 modis_var{istring}='Sensor_Zenith_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
% modis_var{istring}='Sensor_Zenith_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
% modis_var{istring}='Sensor_Zenith_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Sensor_Zenith_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Sensor_Azimuth_Mean'; field_string{istring}='.timeseries3'; istring=istring+1;
%modis_var{istring}='Scattering_Angle_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Scattering_Angle_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Scattering_Angle_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Scattering_Angle_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Scattering_Angle_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1; 


%clear modis_var; istring=1;
%keep this - required for condensation rate calculation
modis_var{istring}='Cloud_Top_Temperature_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Cloud_Top_Temperature_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Cloud_Top_Temperature_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Temperature_Day_Pixel_Counts'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Cloud_Top_Pressure_Day_Mean'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Minimum'; field_string{istring}='.timeseries3'; istring=istring+1; 
%modis_var{istring}='Cloud_Top_Pressure_Day_Standard_Deviation'; field_string{istring}='.timeseries3'; istring=istring+1; 
modis_var{istring}='Cloud_Top_Pressure_Day_Maximum'; field_string{istring}='.timeseries3'; istring=istring+1; 

%clear modis_var; istring=1;
%modis_var{istring}='Nd_timeseries'; field_string{istring}='.mean'; istring=istring+1; 
%modis_var{istring}='Nd_timeseries'; field_string{istring}='.std_dev'; istring=istring+1; 
%modis_var{istring}='W_timeseries'; field_string{istring}='.mean'; istring=istring+1; 
%modis_var{istring}='W_timeseries'; field_string{istring}='.std_dev'; istring=istring+1; 
%modis_var{istring}='H_timeseries'; field_string{istring}='.mean'; istring=istring+1; 
%modis_var{istring}='H_timeseries'; field_string{istring}='.std_dev'; istring=istring+1; 
%modis_var{istring}='LWC_timeseries'; field_string{istring}='.mean'; istring=istring+1; 
%modis_var{istring}='LWC_timeseries'; field_string{istring}='.std_dev'; istring=istring+1; 

end

    case 'Nd vs SZA profiles at constant LAT'

        modis_data_case{imod} = {'y2005_AQUA_TERRA_days_1-182'}; imod=imod+1; %
        modis_data_case{imod} = {'y2005_AQUA_TERRA_days_183-365'}; imod=imod+1; %
        
        clear modis_var
        istring=1;

        %
        modis_var{istring}='Nd_sza'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
        modis_var{istring}='Np_sza'; field_string{istring}='.timeseries3'; istring=istring+1; %(nvar=434);
        
    case 'Monthly averages etc from timeseries3'
        modis_make_monthly_avs_saveload_variables  %script
        
        

end


if i_multi_year_av==0 %if not doing a multi-year average pass all of the files to load_saved_modis_vars_func.m
    clear modis_data_case2
    for iy_multi=1:length(modis_data_case)
        modis_data_case2{iy_multi} = modis_data_case{iy_multi};
    end
end

Nd_PDF_multi=0;
%loop through the variables to load
for ivar=1:length(modis_var)
    fprintf(1,'\nLoading variable %d of %d',ivar,length(modis_var));
    if ivar==1
        iknow_size=0;
        dat_size=NaN;
    else
        iknow_size=1;
    end
    

    
    
    switch data_type
        case 'L3 processed data'
            eval_str = ['[' modis_var{ivar} field_string{ivar} ',dat_size,MLAT,MLON,daynum_timeseries3,modisyear_timeseries3,aqua_terra_timeseries3] = load_saved_modis_vars_func(data_type,modis_var{ivar},savedir_var,modis_data_case2,iknow_size,dat_size,field_string{ivar},high_lat_correction,now_str,irestrict_latlon,ilat_restrict,ilon_restrict);'];
            %Uses load_saved_modis_vars_func function above
        case 'Mock L3 timeseries3'
            eval_str = ['[' modis_var{ivar} field_string{ivar} ',dat_size,MLAT,MLON,daynum_timeseries3,modisyear_timeseries3,aqua_terra_timeseries3] = load_saved_modis_vars_func(data_type,modis_var{ivar},savedir_var,modis_data_case2,iknow_size,dat_size,field_string{ivar},high_lat_correction,now_str);'];
        case 'Nd vs SZA profiles at constant LAT'
            eval_str = ['[' modis_var{ivar} ',dat_size,LAT_sza,mid_sza] = load_saved_modis_vars_func(data_type,modis_var{ivar},savedir_var,modis_data_case,iknow_size,dat_size);'];
        case 'Monthly averages etc from timeseries3'
            eval_str = ['[' modis_var{ivar} ',dat_size,MLAT,MLON,daynum_timeseries3,modisyear_timeseries3,aqua_terra_timeseries3] = load_saved_modis_vars_func(data_type,modis_var{ivar},savedir_var,modis_data_case2,0,dat_size,'''',high_lat_correction,now_str);'];            
    end
    
    eval(eval_str);
    
%     switch data_type
%         case 'Monthly averages etc from timeseries3'
%             switch modis_var{ivar}
%                 case 'Nd_PDF_multi'
%                     Nd_PDF_multi2=Nd_PDF_multi + Nd_PDF_multi2;
%             end
%     end
        
        
end


switch data_type
    case {'L3 processed data','Mock L3 timeseries3'}      
        
        
        
        MLAT=MLAT.MLAT;
        MLON=MLON.MLON;
        
        if irestrict_latlon==1
           MLAT =  MLAT_standard(ilat_restrict);
           MLON =  MLON_standard(ilon_restrict);
        end
        
        imod=1;
        daynum_timeseries3_bk = daynum_timeseries3;
        modisyear_timeseries3_bk = modisyear_timeseries3;        
        aqua_terra_timeseries3_bk = aqua_terra_timeseries3;
        
        filename_choose_saved_MODIS_files
        
        daynum_timeseries3 = daynum_timeseries3_bk;
        modisyear_timeseries3 = modisyear_timeseries3_bk;
        aqua_terra_timeseries3 = aqua_terra_timeseries3_bk;
        
        switch data_type
            case 'L3 processed data'

                CF_temp = load(savevarname{1},'Cloud_Fraction_Liquid');
                LAT = CF_temp(1).Cloud_Fraction_Liquid.timeseries3_LAT;
                LON = CF_temp(1).Cloud_Fraction_Liquid.timeseries3_LON;
                
                %The following is to correct a bug - in make_timeseries_MODIS.m I was
                %using the values in timLAT to select the latitudes (all) required. This
                %was set to run from -89.5 to 89.5. But MLAT is actually ordered from 89.5
                %to -89.5. This causes no problem with the data. But, when I saved the
                %latitudes I saved MLAT instead of timLAT. So the order of the saved
                %latitude array (LAT) just needs to be reversed
                LAT=flipdim(LAT,2);
                
                daynum_timeseries3_MODIS=daynum_timeseries3;
                modisyear_timeseries3_MODIS = modisyear_timeseries3;
                
                LAT=MLAT;
                LON=MLON;
                
                LAT_MODIS = MLAT;
                LON_MODIS = MLON;
        
            case 'Mock L3 timeseries3'
%                MLAT = 0.5*(MLAT(1:end-1) + MLAT(2:end) );
%                MLON = 0.5*(MLON(1:end-1) + MLON(2:end) );  
                LAT=MLAT;
                LON=MLON;
        end
        
       
        
        if exist('Nd_timeseries')
            Nd_timeseries.mean = flipdim(Nd_timeseries.mean,1);
            Nd_timeseries.std_dev = flipdim(Nd_timeseries.std_dev,1);
        end
        
        if exist('H_timeseries')
            H_timeseries.mean = flipdim(H_timeseries.mean,1);
            H_timeseries.std_dev = flipdim(H_timeseries.std_dev,1);
        end
        
        if exist('W_timeseries')
            W_timeseries.mean = flipdim(W_timeseries.mean,1);
            W_timeseries.std_dev = flipdim(W_timeseries.std_dev,1);
        end
        
        if exist('LWC_timeseries')
            LWC_timeseries.mean = flipdim(LWC_timeseries.mean,1);
            LWC_timeseries.std_dev = flipdim(LWC_timeseries.std_dev,1);
        end
        
        mod_data_type = 'timeseries3';
        
        
        switch data_type
            case 'Mock L3 timeseries3'
                mod_data_type_bk = mod_data_type;
                mod_data_type = 'timeseries3';
                %run this here as it is slow & only need to do once for each
                %timseries3 dataset
                filtering_data_get
                
                [N_time3_16]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_16_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N'); 
                [N_time3_37]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');                 

                mod_data_type = mod_data_type_bk;
            otherwise
                %run this here as it is slow & only need to do once for each
                %timseries3 dataset
                filtering_data_get

        end
        

        
        

        
    case 'Monthly averages etc from timeseries3'
        MLAT=MLAT.MLAT;
        MLON=MLON.MLON;
        
         LAT=MLAT;
         LON=MLON;
         
         
        daynum_timeseries3=1; %fudge this to avoid problems in time_inds_modisL3_timeseries3
        daynum_timeseries3_MODIS=daynum_timeseries3;
           gcm_years_loaded_str='y'; 
           
           
           if max(diff(modisyear_timeseries3))==1
               gcm_years_loaded_str = ['y_' num2str(modisyear_timeseries3(1)) '_to_' num2str(modisyear_timeseries3(end))];
           else
               for iyearstr=1:length(modisyear_timeseries3)
                   gcm_years_loaded_str = [gcm_years_loaded_str ' ' num2str(modisyear_timeseries3(iyearstr))];
               end

           end
        

        
        
end

gcm_str = 'MODIS';

clear override_loadsave i_multi_year_av
catch load_save_ERROR
    clear i_multi_year_av override_loadsave
    rethrow(load_save_ERROR)
end

    
    

