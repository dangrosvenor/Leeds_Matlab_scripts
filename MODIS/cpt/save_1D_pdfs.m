%Saves 1D PDFs from watervap. Plot using case 140 of watervap.
try

save_var_dir = '/home/disk/eos1/d.grosvenor/saved_vars/';
%Filename will be var_str_file set below or before this routine.

clear var_str

if ~exist('ioverride_savePDF') | ioverride_savePDF==0
    clear var_str_file 
else
    var_str_file_final = var_str_file;
end

switch axis1D
    case 'x'
        var_switch = x_axis_vals;
    case 'y'
        var_switch = y_axis_vals;
    end

switch var_switch
    case {'Precip rate DAYTIME GCM','Cloudsat precip DAYTIME'}        
        var_str = 'Precip_day';
        var_app = [''];  
        if ~exist('ioverride_savePDF') | ioverride_savePDF==0
            var_str_file = 'Precip_rate_scatter';
        end
        
    case {'Precip rate NIGHTTIME GCM','Cloudsat precip NIGHTTIME'}
        var_str = 'Precip_night';
        var_app = ['']; 
         if ~exist('ioverride_savePDF') | ioverride_savePDF==0
             var_str_file = 'Precip_rate_scatter';
         end
                
    case {'LWP NIGHTTIME GCM, grid-box average'}
        var_str = 'LWP_night';
        var_app = [''];
    case {'LWP DAYTIME GCM, grid-box average'}
        var_str = 'LWP_day';
        var_app = [''];        
    case {'AMSRE TLWP NIGHTTIME','TLWP NIGHTTIME GCM, grid-box average'}
%        LWP_or_TLWP
        var_str = [LWP_or_TLWP '_night'];
        var_app = [''];
%        var_str_file = 'TLWP_PDFs_CAPT';

    case {'AMSRE TLWP DAYTIME','TLWP DAYTIME GCM, grid-box average'}
        var_str = [LWP_or_TLWP '_day'];
        var_app = [''];
%        var_str_file = 'TLWP_PDFs_CAPT';
     
    case 'Precip rate'
        var_str = 'Precip_rate';
        var_app = [''];
        var_str_file = 'Precip_rate';
        
    case 'LWP removal timescale'
        var_str = 'LWP_removal_timescale';
        var_app = [''];
        var_str_file = 'LWP_removal_timescale';
        
    case 'CALIPSO cloud fraction'
        var_str = ['CALISPO_' daynight_str];
        var_app = [''];
        
    case {'MOD35 Cloud Fraction from grid vals timeseries3'}
        var_str = 'MOD35_CF';
        var_app = [''];
        
    case 'Cloud Fraction from grid vals timeseries3'
        var_str = 'MOD06_CF';
        var_app = [''];
        
    case {'CF COSP-CALIPSO GCM'}
        var_str = 'COSP_CALIPSO_CF';
        var_app = [''];
        var_str_file = 'Cloud_Fraction_various';
        
    case {'CF COSP-MODIS GCM'}
        var_str = 'COSP_MODIS_CF';
        var_app = [''];
        var_str_file = 'Cloud_Fraction_various';
        
    case 'CF GCM'
        var_str = 'Model_CF';
        var_app = [''];
        var_str_file = 'Cloud_Fraction_various';
        
    case {'In-cloud averaged AMSRE TLWP'}
        var_str = 'LWP_in_cloud';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'In_cloud_LWP';
        
      case {'LWP GCM'}
        var_str = 'LWP_in_cloud_COSP_MODIS_CF';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'In_cloud_LWP';   
        
    case {'In-cloud averaged (using daytime MOD35) AMSRE TLWP'}
        var_str = 'LWP_in_cloud_AMSRE_MOD35';
        
%        if CF_gcm_thresh(1)==0.05
%            cf_str_lower = '10';
%        else
           cf_str_lower = num2str(100*CF_gcm_thresh(1),'%2.0f');
%        end
%        if CF_gcm_thresh(1)==0.05
%            cf_str_upper = '30';
%        else
            cf_str_upper = num2str(100*CF_gcm_thresh(2),'%2.0f')
%        end
        
        var_app = ['CF_GTE_' cf_str_lower '_AND_LT_' cf_str_upper '_percent'];
        var_str_file = 'In_cloud_LWP';

case {'In-cloud LWP normalised by GCM CF'}
        var_str = 'LWP_in_cloud_model_CF';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'In_cloud_LWP';
         
case {'In-cloud LWP normalised by COSP-CALIPSO CF'}
        var_str = 'LWP_in_cloud_calipso_CF';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
         var_str_file = 'In_cloud_LWP';
         
    case {'In-cloud TLWP GCM norm by COSP-MODIS-CF'}
        var_str = 'TLWP_in_cloud_COSP_MODIS_CF';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'In_cloud_LWP';
        
    case 'LWP GCM grid-box mean'
        var_str = 'LWP_grid_box_mean';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'Grid_box_mean_LWP';
        
    case 'LWP+RWP GCM grid-box mean'
        var_str = 'TLWP_grid_box_mean';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'Grid_box_mean_LWP';
        
    case 'Grid-box mean AMSRE TLWP'
        var_str = 'LWP_grid_box_mean';
        var_app = ['CF_GTE_' num2str(100*CF_gcm_thresh(1),'%2.0f') '_AND_LT_' num2str(100*CF_gcm_thresh(2),'%2.0f') '_percent'];
        var_str_file = 'Grid_box_mean_LWP';
        
    case {'Re COSP GCM'}

        var_str = 'Re';        
        var_app = ['CF_GT_' num2str(100*CF_gcm_thresh,'%2.0f') '_percent'];
        
    case {'CDR Polder2'}
        var_str = 'Re';        
        var_app = ['CF_GT_' num2str(100*0.8,'%2.0f') '_percent'];        
        
    case {'R_{eff 2.1 \mum} (\mum)'}

        var_str = 'Re';        
        var_app = ['CF_GT_' num2str(100*thresh_CF(1),'%2.0f') '_percent'];    
        
end



if ~exist('ioverride_savePDF') | ioverride_savePDF==0
    var_str_file_final = var_str_file;  %otherwise was set before
end


%set the filename for the file to save to
%switch x_axis_vals
%    case {'MOD35 Cloud Fraction from grid vals timeseries3','CF COSP-CALIPSO GCM','CF COSP-MODIS GCM','CF GCM','Cloud Fraction from grid vals timeseries3','CALIPSO cloud fraction'}
%        var_str_file = 'Cloud_Fraction_various';
%end

 eval(['midXbins_' var_str '_' gcm_str '= xdat(1).x;']);
 eval(['PDF_' var_str '_' gcm_str '= ydat(1).y;']);
 eval(['Y_mean_overall_' var_str '_' gcm_str '= Y_mean_overall;']);
 eval(['X_mean_overall_' var_str '_' gcm_str '= X_mean_overall;']);
 eval(['thresh_str_' var_str '_' gcm_str '= thresh_str;']);
 
switch axis1D
    case 'x'
        eval(['Xbins_' var_str '_' gcm_str '= Xbins;']);
    case 'y'
        eval(['Xbins_' var_str '_' gcm_str '= Ybins;']);
end

%timestamp = datestr(now,30);

filesavename = [save_var_dir var_str_file_final '_' var_app '1Dpdfs.mat'];

e = exist(filesavename);

%will need to change so that does not append if file does not exist
if e==2
    app_str = [',''-APPEND'')'];
else
    app_str = [')'];
end
eval(['save(''' filesavename ''',''PDF_' var_str '_' gcm_str ''',''midXbins_' var_str '_' gcm_str ''',''Y_mean_overall_' var_str '_' gcm_str ''',''X_mean_overall_' var_str '_' gcm_str ''',''Xbins_' var_str '_' gcm_str ''',''time_mean_str'',''LAT_str'',''LON_str'',''thresh_str_' var_str '_' gcm_str '''  ', app_str]);

'Done PDF save'


clear ioverride_savePDF
catch errorPDF
   clear ioverride_savePDF
   rethrow(errorPDF);
end
