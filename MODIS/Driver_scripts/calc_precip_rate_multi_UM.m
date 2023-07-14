

%--- Load and process the data
dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
clear fileUM xdat_import ydat_import flag vars_in
for idat=1:99
   flag{idat} = 'nc';
   fileUM_rho{idat} = ''; 
end

idat=1;
%fileUM{idat} = 'xkqkf_qL_qR_.pp.nc'; labs_import(idat).l = '(xkqkf) 1000cm^{-3} No cloud-scheme'; flag{idat}='calc'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqko'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; flag{idat}='calc'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkl'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; flag{idat}='calc'; pole_lat=70; pole_lon=278; idat=idat+1;

idat=1;
vars_UM{idat}='LWP_RWP'; idat=idat+1;
%vars_UM{idat}='Nd'; idat=idat+1;
vars_UM{idat}='rho'; idat=idat+1;
vars_UM{idat}='qR'; idat=idat+1;
vars_UM{idat}='NR'; idat=idat+1;

for idat_UM=1:length(fileUM)
      

    for ivar=1:length(vars_UM)
        filename = [dirUM filename_start '_' vars_UM{ivar} '_.pp.nc'];
        switch flag{idat_UM}
            case 'load_mat'
                filename = [filename '.mat'];                
        end
        
        eval(['vars_in.file_' vars_UM{ivar} ' = filename;']);
        
    end
    
    vars_in.var = 'LWP';  %set as LWP just to get the times for the file (quicker than Nd)
    vars_in.flag = '';
%    vars_in.file_lwp = file_LWP_RWP;    
%    vars_in.file_rho = ''; %filename_rho;
    vars_in.pole_lat = pole_lat;
    vars_in.pole_lon = pole_lon;
    vars_in.time_in = [];

    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
    time_select = time_matlab; 
    

%    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('LWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
%    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM('RWP',flag{idat},filename,filename_rho,pole_lat,pole_lon);
    
    precip_rate = NaN*ones(size(lwp));
    
     for it_driver=1:length(time_select) 
         
         
          vars_in.var = 'Precip_rate';
          vars_in.time_in = time_select(it_driver);
          [dat,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);                
        
          precip_rate(it_driver,:,:) = 24*3600*dat;
          time_matlab(it_driver) = time_out;
         
         
         
     end
     
     filename = [dirUM fileUM{idat_UM} '_' vars_in.var '_.pp.nc'];
%    filename_rho = [dirUM fileUM_rho{idat}];
   filename_save = [filename '.mat'];
    
    filename_start = fileUM{idat_UM};
    
    
    save(filename_save,'precip_rate','-V7.3');

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
            eval(['save(filename_save,''' var_list{i} ''',''-V7.3'',''-APPEND'');']);
        end
        

%    save(filename_save,'time','time_matlab','lat','lon','-V7.3','-APPEND');    
         
end
                  

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
