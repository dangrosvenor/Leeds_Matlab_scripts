% LWP map plots for all times for UM
% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov


LAT_val_DRIVER = [-18]; LON_val_DRIVER =[-75]; %Location of the cross section (going across lon here, but put a long in to find a point on map)
LAT_val_DRIVER = [-20]; LON_val_DRIVER =[-75]; %Location of the cross section (going across lon here, but put a long in to find a point on map)

X_section_type = 'vs lon';
X_section_type = 'vs lat';
d_ind = 25; %no. of grid-cells either side of the chosen lat to average over

irestrict_domain_DRIVER=1;

UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM';
UM_cases = '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM Processing OFF';

 %% Script to get the UM run details by providing the run set name
    %% Provide the case in UM_case_select_runs
    UM_case_select_RUN  %runs UM_case_select_runs

% -- For option setting see inside the loops

zlevs_file = '/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/zlevs_orig_L70_40';
[N_levs,z_levs,dz_levs] = read_zlevs_UM(zlevs_file);

        

for idat_UM=3:3 %1:length(fileUM)
     if iscell(dirUM)==1
            dirUM_i = dirUM{idat_UM};
        else
            dirUM_i = dirUM;
        end
%     filename = [dirUM fileUM{idat}];
%     nc = netcdf(filename);
%     
%     time=nc{'t'}(:);
%     nt_driver=length(time);
%     t0_str=nc{'t'}.time_origin{1};
%     t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%     time_driver = datenum(t0_str2) + time;
%     
%     lon_UM = nc{'x'}(:);
%     lat_UM = nc{'y'}(:);
%     [lon2d,lat2d]=meshgrid(lon_UM,lat_UM);
%     %Convert to normal lat lon from rotated pole coords - make sure the
%     %rotated pole lat and lon are given correctly above
%     [lats2D_driver,lons2D_driver]=em2gm(lat2d,lon2d,pole_lat,pole_lon);  %From Annette M - see email for an example of how she uses it

    
    
    %------- Calculate the data to plot
         %read in the UM data for the specific time
%        time = time_select;
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]

%Load in the LWC
        flag{idat_UM} = '';
        
         vars_in.var = 'qL'; %can set to the variable name to just read the variable
         vars_in.flag = flag{idat_UM}; %can set to '' for just a load of a given variable

         vars_in.file_lwp =  remove_character([dirUM_i fileUM{idat_UM}],'VAR_NAME','qL');          
         vars_in.file_lwp =  remove_character(vars_in.file_lwp,'.mat','');              
         vars_in.file_rho = [dirUM_i fileUM{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = [datenum('13-Nov-2008 21:00')]; 
    
     [qL,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

%Load in the Nd
%         flag{idat_UM} = '';
%         
%          vars_in.var = 'Nd_3D'; %set up a special case to just load this in.
%          vars_in.flag = flag{idat_UM}; %can set to '' for just a load of a given variable
% 
%          vars_in.file_lwp =  remove_character([dirUM_i fileUM{idat_UM}],'VAR_NAME','Nd');          
%          vars_in.file_lwp =  remove_character(vars_in.file_lwp,'.mat','');              
%          vars_in.file_Nd =  remove_character(vars_in.file_lwp,'.mat','');                       
%          vars_in.file_rho = [dirUM_i fileUM{idat_UM}]; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = [datenum('13-Nov-2008 19:00')]; %set to this for all times
%     
%      [Nd,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

     
     
     
     
[ilat,ilon] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_val_DRIVER,LON_val_DRIVER,0.1);    

nz=20; %max height index to plot up to

figure
switch X_section_type
    case 'vs lon'
        dat_qL = meanNoNan(squeeze(qL(1:nz, ilat-d_ind : ilat+d_ind , :)),2);
        dpcolor(gcm_Plon2D_edges_UM(ilat,:),z_levs(1:nz+1),dat_qL); shading flat; colorbar
        % figure
        % dat_Nd = squeeze(Nd(1:nz,ilat,:));
        % dpcolor(gcm_Plon2D_edges_UM(ilat,:),z_levs(1:nz+1),dat_Nd); shading flat; colorbar
    case 'vs lat'
        dat_qL = meanNoNan(squeeze(qL(1:nz,: , ilon-d_ind : ilon+d_ind)),3);
        dpcolor(gcm_Plat2D_edges_UM(:,ilon),z_levs(1:nz+1),dat_qL); shading flat; colorbar
end
% 
% figure
% plot(dat_qL(:),dat_Nd(:),'bo');
% 
% % dat_qL2 = max(qL,[],1);
% % dat_Nd2 = max(Nd,[],1);
% % figure
% % plot(dat_qL2(:),dat_Nd2(:),'bo');

for i=1:600
    ii=find(dat_qL(:,i)>0.2e-4);
    if length(ii)>0
        zi(i)=max(ii);
    end
end

end
     