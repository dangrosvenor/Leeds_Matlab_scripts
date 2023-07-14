filedir_gcm_load = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/';

switch_str = [gcm_year_case '_' gcm_model_case '_' gcm_month_case];

switch switch_str
    %N.B. for new type of data will only use the otherwise entry
    case 'y1990_AM3_01'
%        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120403T161555.mat';
%        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120403T165851.mat';
%        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120404T184111.mat';  
%        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120510T174146.mat';          
%        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120515T195348.mat';
        loadname_gcm = 'saved_gcm_processed_data_AM3_19900101.atmos_inst.tile5.nc_20120517T113533.mat';

    case 'y1990_AM3_with_ice_01'
        loadname_gcm = 'saved_gcm_processed_data_AM3_with_ice_19900101.atmos_inst.tile5.nc_20120601T185328.mat';
        
    case 'y2006_AM3_01'
        loadname_gcm = 'saved_gcm_processed_data_AM3_20060101.atmos_inst.tile5.nc_20120608T095057.mat';
        
    case 'y2006_AM3_07'
        loadname_gcm = 'saved_gcm_processed_data_AM3_20060701.atmos_inst.tile5.nc_20120608T100032.mat';
        
    case 'y2007_AM3_01'
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070101.atmos_inst.tile5.nc_20120608T124200.mat';        
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070101.atmos_inst.tile5.nc_20120613T155442.mat';        
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070101.atmos_inst.tile5.nc_20120613T180556.mat';
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070101.atmos_inst.tile5.nc_20120613T184856.mat';
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070101.atmos_inst.tile5.nc_20120614T135839.mat';

    case 'y2007_AM3_07'
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070701.atmos_inst.tile5.nc_20120608T124411.mat';        
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070701.atmos_inst.tile5.nc_20120613T140824.mat';
        loadname_gcm = 'saved_gcm_processed_data_AM3_20070701.atmos_inst.tile5.nc_20120614T134023.mat';
        
    case 'y2001_CAM5_01'
%        loadname_gcm = 'saved_gcm_processed_data_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120403T155412.mat';
%        loadname_gcm = 'saved_gcm_processed_data_CAM5_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120403T170224.mat';
%        loadname_gcm = 'saved_gcm_processed_data_CAM5_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120404T181423.mat';
%        loadname_gcm = 'saved_gcm_processed_data_CAM5_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120510T172243.mat';        
%        loadname_gcm = 'saved_gcm_processed_data_CAM5_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120515T192941.mat';
%        loadname_gcm = 'saved_gcm_processed_data_CAM5_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120517T113957.mat';        
        loadname_gcm = 'saved_gcm_processed_data_CAM5_CAM5_2deg_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120703T151457.mat';
        
    case 'y2001_CAM5_with_ice'
        loadname_gcm = 'saved_gcm_processed_data_CAM5_with_ice_cam5_1_17_dat.cam.h1.0001-01-01-00000.nc_20120603T084241.mat';
        
    otherwise
        if exist('am3_filepaths_chunk')
            names = fieldnames(am3_filepaths_chunk);
            loadname_gcm = eval(['am3_filepaths_chunk.' names{igcm_type}]);
        else
            disp('*** CANNOT find entry for data selected in gcm_filename_for_year ***');
            return
        end
        
end

filename_gcm =[filedir_gcm_load loadname_gcm];
%filename_gcm =[loadname_gcm];
