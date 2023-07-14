function [var_out,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz,UM_info] = get_LWP_RWP_UM(vars_in)
%Gets or calculates LWP or RWP field from the UM
%UM nc file is ordered [time height lat lon]
% Set time_in to [] for all times.
%Added the UM_info field as a variable to contain any additional info
%(model levels used, filename, etc. etc.)

UM_info.info=['get_LWP_RW_UM.m run on ' datestr(now)];

version_mat=version;
vfind = strfind(version_mat,'R2007b');    
if length(vfind)>0; version_mat='R2007b'; end

var = vars_in.var;
flag = vars_in.flag;
if isfield(vars_in,'time_in')
    time_in = vars_in.time_in;
else
    time_in=[];
end
if isfield(vars_in,'file_lwp')
    file_lwp = vars_in.file_lwp;
elseif isfield(vars_in,'file_LWP_RWP')
    file_lwp = vars_in.file_LWP_RWP;
end

if isfield(vars_in,'file_Nd')
    file_Nd = vars_in.file_Nd;
end
if isfield(vars_in,'file_qL')
    file_qL = vars_in.file_qL;
end
if isfield(vars_in,'file_qR')
    file_qR = vars_in.file_qR;
end
if isfield(vars_in,'file_NR')
    file_NR = vars_in.file_NR;
end
if isfield(vars_in,'file_qv')
    file_qv = vars_in.file_qv;
end
if isfield(vars_in,'file_exner')
    file_exner = vars_in.file_exner;
end
if isfield(vars_in,'file_theta')
    file_theta = vars_in.file_theta;
end
if isfield(vars_in,'file_rho')
    file_rho = vars_in.file_rho;
end
if isfield(vars_in,'file_W')
    file_W = vars_in.file_W;
    if exist(file_W)~=2
        error('file_W does not exist...');
    end
end
if isfield(vars_in,'isave_calc_vals')
    isave_calc_vals = vars_in.isave_calc_vals;
else
    isave_calc_vals=0;
end
if isfield(vars_in,'Nmulti_out')
     Nmulti_out = vars_in.Nmulti_out;
else
     Nmulti_out = 1;
end
if isfield(vars_in,'time_tol')
    time_tol = vars_in.time_tol;
else
    time_tol = 1/3600/24; %set as 1 minute by default
    % Can also set to 'same month' for matching by month and year
end
if isfield(vars_in,'irestrict_region')
    irestrict_region = vars_in.irestrict_region;
    lat_restrict = vars_in.lat_restrict;
    lon_restrict = vars_in.lon_restrict;    
else
    irestrict_region=0;
end          



% indices for z dimension
if isfield(vars_in,'iz')
    iz = vars_in.iz;
else
    iz=-1;
end
pole_lat = vars_in.pole_lat;
pole_lon = vars_in.pole_lon;



%default output
it=0; 
z_um=0;

switch flag
    case 'calc'
%        nc = netcdf(file_rho);
%         time=nc{'t'}(:);
%         t0_str=nc{'t'}.time_origin{1};
%         t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%         time_matlab =  datenum(t0_str2) + time;  
        

        [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_rho,[],pole_lat,pole_lon,time_tol);
        
        if isfield(vars_in,'istaggered') & vars_in.istaggered==1           
                s_grid = [size(gcm_Plat2D_UM,1)+1  size(gcm_Plat2D_UM,1)];
        else
                s_grid = size(gcm_Plat2D_UM);
        end
        
        if Nmulti_out==1
            var_out = NaN*ones([length(time_matlab) s_grid]);
        else
            for imulti_out=1:Nmulti_out
                var_out{imulti_out} = NaN*ones([length(time_matlab) s_grid]);
            end
        end
        
        for it2=1:length(time_matlab)
            fprintf(1,'\nCalculating for it2=%d...',it2);
            
            switch var
                case 'Total_number_aerosol_droplets'
                    %Avoid the z_um parts below as using pre-existing .mat
                    %files (2D)
                otherwise
                    
                    
                    %            z_um = nc{'hybrid_ht'}(:);
                    rho = nc{'rho'}(it2,:,:,:); %Actually = rho*R*R where R is...
                    %rho = netcdf_Dan(nc,'rho',['(' num2str(it2) ',:,:,:)']); %Actually = rho*R*R where R is...
                    %rho = netcdf_Dan2(nc,'rho',[1 1 1 it2(1)],[Inf Inf Inf length(it2)]);
                    
                    R = 6371000; %Radius of the Earth!
                    rho = rho/R/R;
                    
                    %z_um = netcdf_Dan(nc,'hybrid_ht','(:)');
%                    z_um = netcdf_Dan2(nc,'hybrid_ht');
                    z_um = nc{'hybrid_ht'}(:);                    
                    %z_um = repmat(z_um',[1 size(rho,2) size(rho,3)]);
                    
            end


            switch var  %column integrated quantites
                case {'generic_column_integrated','LWP','RWP','accum_mass_total_column_to_z3000','accum_number_total_column_to_z3000','droplet_number_total_column_to_z3000',...
                        'accum_mass_column_integrated','accum_number_column_integrated','droplet_number_column_integrated'}
                    
                    if length(vars_in.z_accum)==1
                        vars_in.z_accum = [-1 vars_in.z_accum(1)];
                    end
                    iz = find(z_um>vars_in.z_accum(1) & z_um<=vars_in.z_accum(2)); %do up to a limiting height in some cases
                     
                    switch var
                        case {'generic_column_integrated'}
                            
                            switch vars_in.VAR_NAME
                                case 'air_mass'
                                    q = ones(size(rho(iz,:,:))); %create an array of ones to multiply by the density.
                                otherwise
                                    nc2 = netcdf(file_lwp);
                                    q = nc2{vars_in.VAR_NAME}(it2,iz,:,:);
                                    %q = netcdf_Dan2(file_lwp,vars_in.VAR_NAME,[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
                            end                            
                            rho = rho(iz,:,:);
                            z_um = z_um(iz,:,:);
                            
                       
                        case 'LWP'
                            nc2 = netcdf(file_qL);
                            q = nc2{'qL'}(it2,:,:,:);                            
                            %q = netcdf_Dan2(file_qL,'qL',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
                        case {'accum_mass_total_column_to_z3000','accum_mass_column_integrated'}
                            nc2 = netcdf(file_lwp);
                            q = nc2{'accum_mass'}(it2,iz,:,:);
                            %nc2=file_lwp;
                            %q = netcdf_Dan2(nc2,'accum_mass',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
                            
                            
                            rho = rho(iz,:,:);
                            z_um = z_um(iz,:,:);
                        case {'accum_number_total_column_to_z3000','accum_number_column_integrated'}
                            nc2 = netcdf(file_lwp);
                            q = nc2{'accum_num'}(it2,iz,:,:);
                            
                            %nc2=file_lwp;
                            % N.B. ncread reads in the data in the opposite
                            % direction to the old method. netcdf_Dan2
                            % fixes this. Note, can only read in start end
                            % stride form now too...
                            %q = double(ncread(nc2,'accum_num',[1 1 iz(1) it2(1)],[Inf Inf it2(end) iz(end)]);
                            %q = permute(q,[4 3 2 1]);
                            %q = netcdf_Dan2(nc2,'accum_num',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
                            rho = rho(iz,:,:);
                            z_um = z_um(iz,:,:);
                            %z_um = z_um(iz);
                        case {'droplet_number_total_column_to_z3000','droplet_number_column_integrated'}
                            nc2 = netcdf(file_lwp);
                            q = nc2{'Nd'}(it2,iz,:,:);  
                            
                            %nc2 = file_lwp;
                            %q = netcdf_Dan2(nc2,'Nd',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
                            
                            
                            rho = rho(iz,:,:);
                            z_um = z_um(iz,:,:);
                        case 'RWP'
                            nc2 = netcdf(file_qR);
                            q = nc2{'qR'}(it2,:,:,:); 
                            
                            %nc2 = file_qR;
                            %q = netcdf_Dan2(nc2,'qR',[1 1 1 it2(1)],[Inf Inf Inf length(it2)]);
                    end
                    
                    clear nc2
                    
% ---------------------- do the calculation
                    if size(rho,1) < size(q,1)
                        rho(end+1,:,:)=rho(end,:,:);
                        z_um(end+1)=z_um(end)+z_um(end)-z_um(end-1);                        
                    end
                    var_out(it2,:,:) = calc_LWP(q,rho,z_um);
                    
                case {'generic_horiz_plane'}
                    [minval,iz] = min(abs(z_um-vars_in.z_plane(1)));
                    %iz = find(z_um>vars_in.z_accum(1) & z_um<=vars_in.z_accum(2));
                    nc2 = netcdf(file_lwp);
                    var_out(it2,:,:) = nc2{vars_in.VAR_NAME}(it2,iz,:,:);
                    
          

                case 'Nd max_LWC'
                    if ~exist('nc2'); nc2 = netcdf(file_Nd); end
%                    nc_W = netcdf(file_W);                    
                    if ~exist('nc_qL'); nc_qL = netcdf(file_qL); end
                    
                    
                    
                    iz = find(z_um<=3e3);
%Old NetCDf version
                    Nd=rho(iz,:,:).*nc2{'Nd'}(it2,iz,:,:)/1e6; %convert to #/cm3, use only the heights below xx km.
                    
%                    W=nc_W{'W'}(it2,iz,:,:);
%                    U=nc_W{'U'}(it2,iz,:,:);
%                    V=nc_W{'V'}(it2,iz,:,:);                    
                    qL=nc_qL{'qL'}(it2,iz,:,:);
% New NetCDF                    
%                    Nd = rho(iz,:,:)/1e6.*netcdf_Dan2(file_Nd,'Nd',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
%                    qL = netcdf_Dan2(file_qL,'qL',[1 1 iz(1) it2(1)],[Inf Inf length(iz) length(it2)]);
       

                    %Remove points where the lwc<0.05 g/m2
                    ilow_lwc = find(qL<0.05e-3); %only use values where LWC is above 0.05 g/kg
                    Nd(ilow_lwc)=NaN;
%                    W(ilow_lwc)=NaN;                    
%                    U(ilow_lwc)=NaN;   
%                    V(ilow_lwc)=NaN;                       
                    %Nd at max LWC over height - wrote a function to return
                    %the values of the second array at the max of the frist
                    %over the 1st dimension
                    [ilwc,Nd2]=max_column_inds(qL,Nd);
%                    [ilwc,W2]=max_column_inds(qL,W);                                        
%                    [ilwc,U2]=max_column_inds(U,W); 
%                    [ilwc,V2]=max_column_inds(V,W);                     
                    
                    var_out{1}(it2,:,:) = Nd2;
%                    var_out{2}(it2,:,:) = W2;                    
%                    var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   
                       %Can't do this yet, since need to re-grid the
                       %staggered U grid - could do this with Iris on
                       %postproc and write out?
%                    var_out{4}(it2,:,:) = sqrt(U2.^2+V2.^2); %Horizontal wind speed   

                case 'accum_num_at_z'
                    if ~exist('nc2'); nc2 = netcdf(file_lwp); end
%                    nc_W = netcdf(file_W);                    
%                    if ~exist('nc_qL'); nc_qL = netcdf(file_qL); end
                    
                    [temp,iz] = min(abs(z_um-vars_in.z_accum));
%                    iz = find(z_um<=3e3);

                    accum_num=nc2{'accum_num'}(it2,iz,:,:); %convert to #/cm3, use only the heights below xx km.
%                    W=nc_W{'W'}(it2,iz,:,:);
%                    U=nc_W{'U'}(it2,iz,:,:);
%                    V=nc_W{'V'}(it2,iz,:,:);                    
%                    qL=nc_qL{'qL'}(it2,iz,:,:);
       

                    %Remove points where the lwc<0.05 g/m2
%                    ilow_lwc = find(qL<0.05e-3); %only use values where LWC is above 0.05 g/kg
%                    Nd(ilow_lwc)=NaN;
%                    W(ilow_lwc)=NaN;                    
%                    U(ilow_lwc)=NaN;   
%                    V(ilow_lwc)=NaN;                       
                    %Nd at max LWC over height - wrote a function to return
                    %the values of the second array at the max of the frist
                    %over the 1st dimension
%                    [ilwc,Nd2]=max_column_inds(qL,Nd);
%                    [ilwc,W2]=max_column_inds(qL,W);                                        
%                    [ilwc,U2]=max_column_inds(U,W); 
%                    [ilwc,V2]=max_column_inds(V,W);                     
                    
                    var_out(it2,:,:) = accum_num;
%                    var_out{2}(it2,:,:) = W2;                    
%                    var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   
                       %Can't do this yet, since need to re-grid the
                       %staggered U grid - could do this with Iris on
                       %postproc and write out?
%                    var_out{4}(it2,:,:) = sqrt(U2.^2+V2.^2); %Horizontal wind speed   

                case 'accum_mass_mean_up_to_z'
                    if ~exist('nc2'); nc2 = netcdf(file_lwp); end
                    
                    [temp,iz] = min(abs(z_um-vars_in.z_accum));
%                    iz = find(z_um<=3e3);

                    accum_mass=nc2{'accum_mass'}(it2,1:iz,:,:); %convert to #/cm3, use only the heights below xx km.
                    
                    var_out(it2,:,:) = meanNoNan(accum_mass,1);
                    
                  case 'accum_mass_ug_per_m3_mean_up_to_z'
                    if ~exist('nc2'); nc2 = netcdf(file_lwp); end
                    
                    [temp,iz] = min(abs(z_um-vars_in.z_accum));
%                    iz = find(z_um<=3e3);

                    accum_mass=nc2{'accum_mass'}(it2,1:iz,:,:) .* rho(1:iz,:,:) *1e9; %convert to ug/m3, use only the heights below xx km.
                    
                    var_out(it2,:,:) = meanNoNan(accum_mass,1);
                    
                    
                case 'accum_mass_ug_per_m3_mean_up_to_z'
                    if ~exist('nc2'); nc2 = netcdf(file_lwp); end
                    
                    [temp,iz] = min(abs(z_um-vars_in.z_accum));
%                    iz = find(z_um<=3e3);

                    accum_mass=nc2{'accum_mass'}(it2,1:iz,:,:) .* rho(1:iz,:,:);
                    
                    var_out(it2,:,:) = meanNoNan(accum_mass,1);
                    
                    
                case 'accum_mass_at_z'
                    if ~exist('nc2'); nc2 = netcdf(file_lwp); end
                    
                    [temp,iz] = min(abs(z_um-vars_in.z_accum));
%                    iz = find(z_um<=3e3);

                    var_out(it2,:,:)=nc2{'accum_mass'}(it2,iz,:,:); 
                                 
                    
                case 'W max_LWC'
                    nc2 = netcdf(file_W);
                    nc_qL = netcdf(file_qL);
                    iz = find(z_um<=3e3);

                    W=nc2{'W'}(it2,iz,:,:);
                    qL=nc_qL{'qL'}(it2,iz,:,:);
       

                    %Remove points where the lwc<0.05 g/m2
                    ilow_lwc = find(qL<0.05e-3); %only use values where LWC is above 0.05 g/kg
                    W(ilow_lwc)=NaN;
                    %Nd at max LWC over height - wrote a function to return
                    %the values of the second array at the max of the first
                    %over the 1st dimension
                    [ilwc,W2]=max_column_inds(qL,W);
                    
                    var_out(it2,:,:) = W2;
                    
                case 'Total_number_aerosol_droplets'


                    
                    %Save time by restricting to lower 3km
                    %iz = find(z_um<=3e3);


                    
                    if it2==1
                        mat_accum = load(vars_in.file_accum);
                        mat_coarse = load(vars_in.file_coarse);
                        mat_droplets = load(vars_in.file_droplets);
                        mat_airmass = load(vars_in.file_airmass);
                        
                        var_out = NaN*ones([length(time_matlab) size(gcm_Plat2D_UM)]);
                    end
                    
                    dat = (...
                        mat_accum.accum_number_total_column_to_z1500(it2,:,:) +...
                        mat_coarse.coarse_number_total_column_to_z1500(it2,:,:) +...
                        mat_droplets.droplet_number_total_column_to_z1500(it2,:,:) ...
                        );             
                    dat = squeeze( dat ./ mat_airmass.air_mass_total_column_to_z1500(it2,:,:) );
       

                    var_out(it2,:,:) = dat;            
                    
%                     var_out{1}(it2,:,:) = Nd2;
%                     var_out{2}(it2,:,:) = W2;                    
%                     var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   


%clear nc_qR nc_NR nc_qL nc_Nd nc_rho   % not sure if this closes the file or not...?

                    
                    
                case 'Radar dBZ'
                    nc_qR = netcdf(file_qR);
                    nc_NR = netcdf(file_NR); 
                    nc_qL = netcdf(file_qL);
                    nc_Nd = netcdf(file_Nd);
                    nc_rho = netcdf(file_rho);
                    
                    %Save time by restricting to lower 3km
                    iz = find(z_um<=3e3);

                    q_dBZ{1} = nc_qR{'qR'}(it2,iz,:,:); %Rain MR in kg/kg                                     
                    n_dBZ{1} = nc_NR{'NR'}(it2,iz,:,:); %Rain num in #/kg
                    
                    q_dBZ{5} = nc_qL{'qL'}(it2,iz,:,:); %Liquid MR in kg/kg                                     
                    n_dBZ{5} = nc_Nd{'Nd'}(it2,iz,:,:); %Liquid num in #/kg                                        
                    
                    q_dBZ{2} = 0; n_dBZ{2} = 0; %Snow
                    q_dBZ{3} = 0; n_dBZ{3} = 0; %Graupel
                    q_dBZ{4} = 0; n_dBZ{4} = 0; %Ice                                       
                    
                    RHO = rho(iz,:,:);
                    
                    if it2==1
                       var_out = NaN*ones([length(time_matlab) length(iz) size(gcm_Plat2D_UM)]);                        
                    end
       

                    model = 'UM';
%                    model = 'LEM';                    
                    ztot = calc_radar_UM_LEM_etc(model,q_dBZ,n_dBZ,RHO);
                    var_out(it2,:,:,:) = 10.*log10(ztot);                  
                    
%                     var_out{1}(it2,:,:) = Nd2;
%                     var_out{2}(it2,:,:) = W2;                    
%                     var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   


clear nc_qR nc_NR nc_qL nc_Nd nc_rho   % not sure if this closes the file or not...?

                case 'BL height using qL'
                    nc_qL = netcdf(file_qL);

                    
                    %Save time by restricting to lower 3km
                    iz = find(z_um<=3e3);
                    
                    %qL = netcdf_Dan2(vars_in.file_qL,'qL',[1 1 1 17],[Inf Inf iz( 1]);

                    qL = nc_qL{'qL'}(it2,iz,:,:); %Rain MR in kg/kg                                     
                                                       
                    if it2==1
                       var_out = NaN*ones([length(time_matlab) size(gcm_Plat2D_UM)]);                        
                    end
       
                    %array of heights for all columns
                    zarr3D = repmat(z_um(iz),[1 size(gcm_Plat2D_UM)]);
                    %NaN all ones that are below the threshold
                    zarr3D(qL<0.2e-4) = NaN;
                    %Now the max for each column will be the max height at
                    %which our threshold qL exists
                    [zi] = max(zarr3D,[],1);
                    inan=find(isnan(zi)==1);
                    zi(inan)=0;                    
                    var_out(it2,:,:) = zi;
                                     
                    
%                     var_out{1}(it2,:,:) = Nd2;
%                     var_out{2}(it2,:,:) = W2;                    
%                     var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   


clear nc_qL   % not sure if this closes the file or not...?


                case 'BL height using RH50'
                    nc_qv = netcdf(vars_in.file_qv);
                    nc_th = netcdf(vars_in.file_th);
                    nc_P = netcdf(vars_in.file_P);
                    
                    %Save time by restricting to lower 3km
                    iz = find(z_um<=3e3);
                    
                    %qL = netcdf_Dan2(vars_in.file_qL,'qL',[1 1 1 17],[Inf Inf iz( 1]);

                    qv = nc_qv{'qv'}(it2,iz,:,:); %MR in kg/kg                                     
                    th = nc_th{'th'}(it2,iz,:,:); %MR in kg/kg                                     
                    P = nc_P{'P'}(it2,iz,:,:); %MR in kg/kg                                     

                    if it2==1
                       var_out = NaN*ones([length(time_matlab) size(gcm_Plat2D_UM)]);                        
                    end
                    
                    T = th2t(P,th);
                    rh=calc_RH(qv,P,T);
                    
       
                    %array of heights for all columns
                    zarr3D = repmat(z_um(iz),[1 size(gcm_Plat2D_UM)]);
                    %NaN all ones that are below the threshold
                    zarr3D(rh<0.5) = NaN;
                    %Now the max for each column will be the max height at
                    %which our threshold qL exists
                    [zi] = max(zarr3D,[],1);
                    inan=find(isnan(zi)==1);
                    zi(inan)=0;                    
                    var_out(it2,:,:) = zi;
                                     
                    
%                     var_out{1}(it2,:,:) = Nd2;
%                     var_out{2}(it2,:,:) = W2;                    
%                     var_out{3}(it2,:,:) = squeeze(max(W,[],1));  %max W in whole column   


clear nc_qv nc_P nc_th   % not sure if this closes the file or not...?



                case 'Nd_GCM'
                    if ~exist('nc_Nd');nc_Nd = netcdf(vars_in.file_Nd_times_LYR); end
                    if ~exist('nc_cf');nc_cf = netcdf(vars_in.file_LYR_weight); end %akin to cloud fraction
                    % Need to divide the Nd times LYR by this
                    
%                    nc_W = netcdf(file_W);                    
%                    nc_qL = netcdf(file_qL);
                    
                    iz = find(z_um>400 & z_um<=3.2e3); %Keep 400m above the surface too - N.b is prob above surface
                    %, not sea level - just stick to oceans!

%                    Nd=rho(iz,:,:).*nc2{'Nd'}(it2,iz,:,:)/1e6; %convert to #/cm3, use only the heights below xx km.
                    Nd_weight = nc_Nd{'Nd_times_LYR_CLD_WEIGHT'}(it2,iz,:,:);
                    weight = nc_cf{'LYR_CLD_WEIGHT'}(it2,iz,:,:);
                    weight2 = weight;
                    Nd_weight2 = Nd_weight;
                    % Ony include points above a minimum weight to avoid
                    % divide by zero
                    iw = find(weight<1e-3);                    
                    weight2(iw)=NaN;
                    Nd_weight2(iw)=NaN;
                  
%                    W=nc_W{'W'}(it2,iz,:,:);
%                    U=nc_W{'U'}(it2,iz,:,:);
%                    V=nc_W{'V'}(it2,iz,:,:);                    
%                    qL=nc_qL{'qL'}(it2,iz,:,:);
       
% Could screen by the LYR_WEIGHT - i.e. if not that much cloud then ignore.
% Or just do an average over the bottom 3 km

                    Nd2 = meanNoNan(Nd_weight ./ weight2,1); %already in #/cm3. use only the heights below xx km. 
                    
%Trying another methdo where we weight by the mean weight of each level - so that levels with more cloud count
% more towards the mean.
%Need to make sure that the Nd*weight array has the NaNs added too here, or will be a different set of data to that used for hte mean of 
%the weights. Doesn't matter for the other method since are divideing by
%teh NaNs in the weight array.
                    Nd_zweight = meanNoNan(Nd_weight2,1) ./ meanNoNan(weight2,1);
                    

                    
                    weight3 = weight;
                    % Ony include points above a minimum weight to avoid
                    % divide by zero
                    iw = find(weight<4e-3);                    
                    weight3(iw)=NaN;            
                    
                    Nd3 = meanNoNan(Nd_weight ./ weight3,1); %already in #/cm3. use only the heights below xx km. 

                    %Remove points where the lwc<0.05 g/m2
%                    ilow_lwc = find(qL<0.05e-3); %only use values where LWC is above 0.05 g/kg
%                    Nd(ilow_lwc)=NaN;
%                    W(ilow_lwc)=NaN;                    
%                    U(ilow_lwc)=NaN;   
%                    V(ilow_lwc)=NaN;                       
                    %Nd at max LWC over height - wrote a function to return
                    %the values of the second array at the max of the frist
                    %over the 1st dimension
%                    [ilwc,Nd2]=max_column_inds(qL,Nd);
%                    [ilwc,W2]=max_column_inds(qL,W);                                        
%                    [ilwc,U2]=max_column_inds(U,W); 
%                    [ilwc,V2]=max_column_inds(V,W);                     
                    
                    var_out{1}(it2,:,:) = Nd2;
                    var_out{2}(it2,:,:) = Nd_zweight;                    
                    var_out{3}(it2,:,:) = Nd3;
                       %Can't do this yet, since need to re-grid the
                       %staggered U grid - could do this with Iris on
                       %postproc and write out?
%                    var_out{4}(it2,:,:) = sqrt(U2.^2+V2.^2); %Horizontal
%                    wind speed   
%                         
              
                    
            end
            
       


        end
        
        if isave_calc_vals==1
            eval_str=[lower(var) ' = var_out;']; eval(eval_str);
            
            eval_str=['save(filename_save,''' lower(var) ''',''-V7.3'')']; eval(eval_str);            
            
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

         
        end

%%        
        
    case 'load_mat'
        
        dat = load(file_lwp);
        time_driver = dat.time_matlab;
        
        switch var
            case 'Nd_zweight'
                var_str='Nd_zweight'; 
             case 'Nd'
                var_str='Nd';    
            case 'W'
                var_str='W';
            case 'dBZ'
                var_str = 'dBZ';
            otherwise        
                var_str = var;
                var_str_lower=lower(var);  %Fix this so no longer need to use lower case!
        end
        
        try
            siz = eval(['size(dat.' var_str ');']);
        catch
            var_str = var_str_lower;
            siz = eval(['size(dat.' var_str ');']);
        end
        
        if irestrict_region==1
            %find the linear indices since will return as a linear list
            %rather than a block
           iregion = find(dat.gcm_Plat2D_UM >= lat_restrict(1) & dat.gcm_Plat2D_UM <= lat_restrict(2)...
               & dat.gcm_Plon2D_UM >= lon_restrict(1) & dat.gcm_Plon2D_UM <= lon_restrict(2) );   
           
            if length(siz)==4
                var_out = NaN*ones([length(time_in) siz(2) length(iregion)]);
            end
            
        else
            if length(siz)==4
                var_out = NaN*ones([length(time_in) siz(2) siz(3)*siz(4)]);
            else
                var_out = NaN*ones([length(time_in) siz(2) siz(3)]);
            end

        end
        
        if length(time_in)==0
            var_out = eval(['dat.' var_str '(:,:,:);']);
            time_matlab = time_driver;
            daynum_timeseries3_UM = dat.daynum_timeseries3_UM;
            modisyear_timeseries3_UM = dat.modisyear_timeseries3_UM;
            gcm_time_UTC_UM = dat.gcm_time_UTC_UM;
            gcm_time_matlab_UM = dat.gcm_time_matlab_UM;
        else
        

            clear it

        
            
            for i=1:length(time_in)
                
                if length(strfind(time_tol,'same month'))>0
                    [Y,M,D]=datevec(time_driver);
                    [Y2,M2,D2]=datevec(time_in(i));
                    it_val = find(Y==Y2 & M==M2);
                else
                    it_val = find(abs(time_driver-time_in(i)) < time_tol );                   
                end
                
                if length(it_val)>0
                    it(i) = it_val(1);  %for some reason the Iceland Nd files had duplicate times
                       %Only the first time had useful data
                else
                    error('*** DPG - cannot locate requested time - try reducing tolerance? ***');
                end
                if irestrict_region==1
                    if length(siz)==4
                        var_out(i,:,:) = eval(['dat.' var_str '(it(i),:,iregion);']);
                    else
                        var_out(i,:,:) = eval(['dat.' var_str '(it(i),:,:);']);
                    end
                else
                    var_out(i,:,:) = eval(['dat.' var_str '(it(i),:,:);']);
                end
                
                time_matlab(i) = time_driver(it(i));    
                daynum_timeseries3_UM(i) = dat.daynum_timeseries3_UM(i);
                modisyear_timeseries3_UM(i) = dat.modisyear_timeseries3_UM(i);
                gcm_time_UTC_UM(i) = dat.gcm_time_UTC_UM(i);
                gcm_time_matlab_UM(i) = dat.gcm_time_matlab_UM(i);
            end
                                   
        end
        


        i=1;
        %var_list{i} = 'time_matlab'; i=i+1;
        var_list{i} = 'gcm_Plat2D_UM'; i=i+1;
        var_list{i} = 'gcm_Plon2D_UM'; i=i+1;
        var_list{i} = 'gcm_Plat2D_edges_UM'; i=i+1;
        var_list{i} = ' gcm_Plon2D_edges_UM'; i=i+1;
%        var_list{i} = 'it'; i=i+1;
%        var_list{i} = 'daynum_timeseries3_UM'; i=i+1;
%        var_list{i} = 'modisyear_timeseries3_UM'; i=i+1;
%        var_list{i} = 'gcm_time_UTC_UM'; i=i+1;
%        var_list{i} = 'gcm_time_matlab_UM'; i=i+1;
        
        for i=1:length(var_list)
            eval([var_list{i} ' = dat.' var_list{i} ';']);
        end
        
%       time = dat.time;
%        time_matlab = dat.time_matlab;        
%        lat = dat.lat;
%        lon = dat.lon;    
        
        
%% Normal netCDF case (flag set to '')
    otherwise   
        switch var
            case 'Nd mat_LWP'
                %get from ,mat file instead
                
        end
        %Use LWP file just to get the times:-
        [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_lwp,time_in,pole_lat,pole_lon,time_tol);
%        var_out = NaN*ones([length(it) size(gcm_Plat2D_UM)]);
        switch var
            case {'LWP_10min','RWP_10min'}              
%                 time=nc{'t_1'}(:);
%                 t0_str=nc{'t_1'}.time_origin{1};
%                 t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%                 time_driver = datenum(t0_str2) + time;
%                 time_matlab = time_driver;
%                it = 1:length(time_driver);
                [temp,time_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,it] ...
                    = read_time_UM(nc,'t_1',time_in);
                

        end

                    
        
        switch var
            case 'LWP'
                var_out=nc{'LWP'}(it,:,:);
            case 'RWP'
                var_out=nc{'RWP'}(it,:,:);
            case {'Nd max_LWC'}
                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);                
                [nc_qL,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_qL,time_in,pole_lat,pole_lon,time_tol);                                
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc2{'hybrid_ht'}(:);
                    iz = find(z_um<=3e3);
                    try
                        Nd=nc2{'Nd'}(it3,iz,:,:); %only the heights below xx km here.
                        qL=nc_qL{'qL'}(it3,iz,:,:);                        
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end

                    %Remove points where the lwc<0.05 g/m2
                    ilow_lwc = find(qL<0.05e-3); %only use values where LWC is above 0.05 g/kg
                    Nd(ilow_lwc)=NaN;
                    %Nd at max LWC over height - wrote a function to return
                    %the values of the second array at the max of the frist
                    %over the 1st dimension
                    [ilwc,Nd2]=max_column_inds(qL,Nd);


                    
                    if length(it)>1
                        var_out(it2) = meanNoNan(Nd2(:),1);
                    else
                        var_out = Nd2;  %output the 2D field
                    end
                end                
            case {'Nd_max','Nd mat_LWP'}
                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);                
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc2{'hybrid_ht'}(:);
                    iz = find(z_um<=2e3);
                    try
                        Nd=nc2{var}(it3,iz,:,:);
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end


                    %Max over height
%                    Nd2 = meanNoNan(Nd,1);
                    Nd2 = max(Nd,[],1);  
                    switch var
                        case 'Nd_max'
                            lwp = nc{'LWP'}(it3,:,:);
                        case 'Nd mat_LWP'
                            dat = load(file_lwp);
                            lwp=dat.lwp(it3,:,:);
                    end
                    ilwp = find(lwp*1e3<10);
                    %Remove the columns where the LWP is too low
                    Nd2(ilwp) = NaN;
                    if length(it)>1
                        var_out(it2) = meanNoNan(Nd2(:),1);
                    else
                        var_out = Nd2;  %output the 2D field
                    end
                end
                
         case {'Nd_2km_minLWP'}
                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);                
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc2{'hybrid_ht'}(:);
                    iz = find(z_um<=2e3);
                    try
                        Nd=nc2{var}(it3,iz,:,:);
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end

                    switch var
                        case 'Nd'
                            lwp = nc{'LWP'}(it3,:,:);
                        case 'Nd mat_LWP'
                            dat = load(file_lwp);
                            lwp=dat.lwp(it3,:,:);
                    end
                    ilwp = find(lwp*1e3<10);
                    %Remove the columns where the LWP is too low
                    Nd(:,ilwp) = NaN;
                    var_out = Nd2;  %output the 2D field

                end        
                
        case 'Tau'
            %This is just to get the times I think
                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);
                [nc_rho,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_rho,time_in,pole_lat,pole_lon,time_tol);
                [nc_qL,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_qL,time_in,pole_lat,pole_lon,time_tol);
                
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc_rho{'hybrid_ht'}(:);
%                    iz = find(z_um<=2e3);
                    iz = 1:length(z_um);
                    
                    rho=nc_rho{'rho'}(it3,iz,:,:);  %Actually = rho*R*R where R is...
                    R = 6371000; %Radius of the Earth!
                    rho = rho/R/R;
                    qL=nc_qL{'qL'}(it3,iz,:,:);
                    qL(qL<0)=0;

                    try
                        if length(file_Nd)==0 %For the case where we have no Nd data (e.g. old mphys)
                            Nd = 100e6 * ones(size(rho));
                            fprintf(1,'\n*** WARNING - no Nd data associated with %s***',file_qL);
                        else
                            Nd=nc2{'Nd'}(it3,iz,:,:);
                        end
                        
                        Nd(Nd<0)=0;
                        Nd_t0 = nc2{'Nd'}(1,1,:,:);
  
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end


                    dz = diff(cat(1,0,z_um));
                    dz2 = repmat(dz,[1 size(Nd_t0,1) size(Nd_t0,2)]);
                    
                    Q=2;
                    k=0.8;
                    Tau = sum(pi*Q.*(Nd.*rho.*k).^(1/3) .* (3*rho.*qL./4/pi/1e3).^(2/3) .* dz2);
                    
                   
                    if length(it)>1
                        var_out(it2) = meanNoNan(Tau(:),1);
                    else
                        var_out = Tau;  %output the 2D field
                    end
                end
        
                case 'Precip_rate'
%                [nc_Nd,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);
                [nc_rho,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_rho,time_in,pole_lat,pole_lon,time_tol);
                [nc_qR,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_qR,time_in,pole_lat,pole_lon,time_tol);
                [nc_NR,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_NR,time_in,pole_lat,pole_lon,time_tol);                
                
                var_t0 = nc_qR{'qR'}(1,1,:,:);
                
                if length(it)>1
                    var_out(it2) = NaN*ones(size(it));
                else
                    var_out = NaN*ones([size(var_t0)]);  %output the 2D field
                end
                
                
                 
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc_rho{'hybrid_ht'}(:);
%                    iz = find(z_um<=2e3);
                    iz = 1:length(z_um);
                    

                    try
                        
                        rho=nc_rho{'rho'}(it3,iz,:,:);  %Actually = rho*R*R where R is...
                        R = 6371000; %Radius of the Earth!
                        rho = rho/R/R;
                        qR=nc_qR{'qR'}(it3,iz,:,:);
                        min_qR=1e-10;
                        qR(qR<min_qR)=min_qR;
                        

                        NR=nc_NR{'NR'}(it3,iz,:,:);
                        min_NR = 0.001e6; % per cc
                        NR(NR<min_NR)=NaN;
                        
                        [V,lambda,re] = fall_speed_mass_weighted(qR,-999,NR,rho);
                                               
                    
%                         if length(file_Nd)==0 %For the case where we have no Nd data (e.g. old mphys)
%                             Nd = 100e6 * ones(size(rho));
%                             fprintf(1,'\n*** WARNING - no Nd data associated with %s***',file_qL);
%                         else
%                             Nd=nc2{'Nd'}(it3,iz,:,:);
%                         end
                        
  
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end


                    precip_rate = rho.*V.*qR;
                    max_PR = max(precip_rate,[],1);
                    
                   
                    if length(it)>1
                        var_out(it2) = meanNoNan(max_PR(:),1);
                    else
                        var_out = max_PR;  %output the 2D field
                    end
                end
                
         case 'RH'
            %This is just to get the times I think
%                [nc_Nd,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_Nd,time_in,pole_lat,pole_lon,time_tol);
                [nc_rho,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_rho,time_in,pole_lat,pole_lon,time_tol);
                [nc_qv,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_qv,time_in,pole_lat,pole_lon,time_tol);
                [nc_exner,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_exner,time_in,pole_lat,pole_lon,time_tol);                
                [nc_theta,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_theta,time_in,pole_lat,pole_lon,time_tol);                
                
                for it2=1:length(it)
                    it3 = it(it2);
                    z_um = nc_rho{'hybrid_ht'}(:);
                    iz = find(z_um<=1e3);
%                    iz = 1:length(z_um);
                    

                    try
                        
                        rho=nc_rho{'rho'}(it3,iz,:,:);  %Actually = rho*R*R where R is...
                        R = 6371000; %Radius of the Earth!
                        rho = rho/R/R;
                        
                        qv=nc_qv{'qv'}(it3,iz,:,:);
                        th=nc_theta{'theta'}(it3,iz,:,:);
%                        min_qR=1e-10;
%                        qR(qR<min_qR)=min_qR;
                        
                        ex=nc_exner{'exner'}(it3,iz,:,:); %T is theta*exner
%                        NR(NR<0)=0;

                       

                        

                    
%                         if length(file_Nd)==0 %For the case where we have no Nd data (e.g. old mphys)
%                             Nd = 100e6 * ones(size(rho));
%                             fprintf(1,'\n*** WARNING - no Nd data associated with %s***',file_qL);
%                         else
%                             Nd=nc2{'Nd'}(it3,iz,:,:);
%                         end
                        
                        var_t0 = nc_exner{'exner'}(1,1,:,:);
  
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end

                    T = th.*ex;
                    P = ex.^(1/0.286).*1000e2;
                    qs = SatVapPress(T,'goff','liq',P,1);
                    f = 1.6094e+06;
                    qs = qs/f;
                    RH = qv./qs;
                    
                    
                    RH_mean = meanNoNan(RH,1);
%                    max_PR = max(precip_rate,[],1);
                    
                   
                    if length(it)>1
                        var_out(it2) = meanNoNan(RH_mean(:),1);
                    else
                        var_out = RH_mean;  %output the 2D field
                    end
                end
                
            case 'time'
                %exit now
                
                
            case {'Transmission_down_surf_LWP_LT_0pt1'}
                file_lwp_mat = [remove_character(vars_in.fileUM{1},'VAR_NAME','LWP') '.mat'];
%                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_temp,time_in,pole_lat,pole_lon,time_tol);                
                
                file_temp = remove_character([vars_in.dirUM vars_in.fileUM{1}],'VAR_NAME','SW_down_TOA');
                [nc3,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_temp,time_in,pole_lat,pole_lon,time_tol);                                
                file_temp = remove_character([vars_in.dirUM vars_in.fileUM{1}],'VAR_NAME','SW_down_surf');
                [nc4,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_temp,time_in,pole_lat,pole_lon,time_tol);                                                
                for it2=1:length(it)
                    it3 = it(it2);
                    try
                        SW_down_surf=nc4{'SW_down_surf'}(it3,:,:); %only the heights below xx km here.                        
                        lwp=load(file_lwp_mat,'lwp'); lwp=lwp.lwp(it3,:,:);
                        SW_TOA=nc3{'SW_down_TOA'}(it3,:,:);                        
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end

                    %Remove points where the lwp>0.1 g/m2
                    ilwp = find(lwp>0.1);
                    SW_TOA(ilwp)=NaN;
                    trans =  SW_down_surf./SW_TOA;


                    
                    if length(it)>1
                        var_out(it2) = meanNoNan(trans(:),1);
                    else
                        var_out = trans;  %output the 2D field
                    end
                end  
                
            case {'SW_down_surf_LWP_LT_0pt1'}
                file_lwp_mat = [remove_character(vars_in.fileUM{1},'VAR_NAME','LWP') '.mat'];
%                [nc2,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_temp,time_in,pole_lat,pole_lon,time_tol);                
                

                file_temp = remove_character([vars_in.dirUM vars_in.fileUM{1}],'VAR_NAME','SW_down_surf');
                [nc4,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(file_temp,time_in,pole_lat,pole_lon,time_tol);                                                
                for it2=1:length(it)
                    it3 = it(it2);
                    try
                        SW_down_surf=nc4{'SW_down_surf'}(it3,:,:); %only the heights below xx km here.                        
                        lwp=load(file_lwp_mat,'lwp'); lwp=lwp.lwp(it3,:,:);                       
                    catch
                        break  %For cases where the netCDF does not contain all of the times listed
                    end

                    %Remove points where the lwp>0.1 g/m2
                    ilwp = find(lwp>0.1);
                    SW_down_surf(ilwp)=NaN;
                    
                    if length(it)>1
                        var_out(it2) = meanNoNan(SW_down_surf(:),1);
                    else
                        var_out = SW_down_surf;  %output the 2D field
                    end
                end                

                
                
            otherwise
                
                switch var
                    case 'Nd_3D'
                        %Special case when just want to laod in 3D Nd field
                        var='Nd';
                end
                
                
%                switch version_mat
%                   
%                    case 'R2007b'
                        
                        eval_str = ['var_nc = nc{''' var '''};'];  eval(eval_str); %Store the NetCFD variable in var_nc
                        eval_str = ['var_ncz = nc{''hybrid_ht''};'];  eval(eval_str); %Store the NetCFD variable in var_nc    
                        z_um = var_ncz(:);                        
                        
                        %Can now do e.g. dim(var_nc) to get the info for each
                        %dimension in a cell array
                        %                eval_str = ['var_out = nc{''' var '''}(it,:,:);']; eval(eval_str);
                        
                        %Number of dimensions
                        ndims = length(dim(var_nc));
                        
                        if length(it)==0
                            error('***DPG Problem: length(it)==0');
                        end
                        
                        if ndims==3
                            var_out = var_nc(it,:,:);
%                            z_um = var_ncz(it,:,:);                            
                            iz=-1;
                        elseif ndims==4
                            if iz==-1
                                var_out = var_nc(it,:,:,:);
%                                z_um = var_ncz(it,:,:,:);                                
                            else
                                var_out = var_nc(it,iz,:,:);
%                                z_um = var_ncz(it,iz,:,:);                                
                            end
                        end
                        
%                     otherwise
%                         
%                         % N.B. ncread reads in the data in the opposite
%                             % direction to the old method. netcdf_Dan2
%                             % fixes this. Note, can only read in start end
%                             % stride form now too...
%                             %q = double(ncread(nc2,'accum_num',[1 1 iz(1) it2(1)],[Inf Inf it2(end) iz(end)]);
%                             %q = permute(q,[4 3 2 1]);
%                             
%                             %var_nc = netcdf_Dan2(nc,var,[1 1 1 1],[Inf Inf Inf Inf]);
%                             var_nc = netcdf_Dan2(nc,var);
%                             var_nc = squeeze(var_nc);
%                             
%                             
%                         %Number of dimensions
%                         ndims = length(size(var_nc));
%                         
%                         if length(it)==0
%                             error('***DPG Problem: length(it)==0');
%                         end
%                         
%                         if ndims==3
%                             var_out = var_nc(it,:,:);
%                             iz=-1;
%                         elseif ndims==4
%                             if iz==-1
%                                 var_out = var_nc(it,:,:,:);
%                             else
%                                 var_out = var_nc(it,iz,:,:);
%                             end
%                         end
%                         
%                         
%                 end
                
        
       
                
                
                
                %z_um = nc{'hybrid_ht'}(:);
                
                %z_um = netcdf_Dan2(nc,'hybrid_ht');

                
        end
end

var_out=squeeze(var_out);