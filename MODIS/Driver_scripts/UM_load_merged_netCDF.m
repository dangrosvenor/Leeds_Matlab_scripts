function [dat_global] = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,pole_lat,pole_lon,var_UM2,opts)

if exist('opts')
    %Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='opts'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
        eval(eval_str);
    end
    
end
       
if ~exist('isort_time')
    isort_time=0;
end


if ~exist('var_UM2')
    var_UM2 = var_UM;
end

istr=strfind(var_UM,'_time_only');
if length(istr)==0
    time_only=0;
else
    time_only=1;
    var_UM=var_UM(1:istr-1);
end

scale_fac=1;

switch load_type
    case 'mat'
        filename = [dirUM '/' var_UM '/' run_type '_' var_UM '_native res_ALL.mat'];        
        dat_global = load(filename);
        
    case {'merged netCDF','individual netCDF files','named netCDF'}
        dir_UM2 = [dirUM '/' var_UM];
        switch load_type
            case 'merged netCDF'
                filenames(1).name = ['merged.nc'];
                
            case 'individual netCDF files'
                filenames = [dir([dirUM '/' var_UM '/*.nc'])];
                
            case 'named netCDF'
                filenames(1).name = named_file;
                dir_UM2 = named_dir;
        end
        
        if length(filenames)>1
           if ~exist('cat_dim')
               error('Need to set cat_dim if supplying multiple NetCDF files')
           end
        else
            cat_dim = 1; %This should just return the data as is
        end
        
        for ifile=1:length(filenames)            
            filename = [dir_UM2 '/' filenames(ifile).name];
            
            if ~exist(filename)
                filename = [dirUM '/' var_UM '/' var_UM '_merged.nc'];
                if ~exist(filename)
                    error(['*** Filename ' filename ' or merged.nc does not exist!']);
                end
            end
            
            nc = netcdf(filename);
            if time_only==1
                dat = NaN;
            else
                dat = nc{var_UM2}(:);
            end
            
            if ifile==1
                
                if ~exist('lat_var')
                    lat_var = 'Latitude'; lon_var = 'Longitude';
                end
                
                %Slices have the lat lon converted from rotated coords for the nest
                %already.
                gcm_Plat2D_UM = nc{lat_var}(:);
                gcm_Plon2D_UM = nc{lon_var}(:);
                
                if length(gcm_Plat2D_UM)==0
                    error(['*** No lat/lon data - check var name ***']);
                end
                
                %Consider using UM_lat_lon_intialise_FUNC.m to sort out the
                %lat/lon, make the edges, area, etc.
                
                %Deal with emissions files and others that have single lat and lon
                %vectors
                if size(gcm_Plat2D_UM,2)==1
                    [gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(gcm_Plon2D_UM,gcm_Plat2D_UM);
                end
                
                % N.B. - Values from Python scripts are already rotated to proper lat and lons!
                %         switch run_type
                %             case 'nested'
                %                 [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(gcm_Plat2D_UM,gcm_Plo
                %                 n2D_UM,pole_lat,pole_lon);
                %        end
                
                
                i180 = find(gcm_Plon2D_UM>180);
                gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;
                
                if ~exist('time_var')
                    time_var = 'Time';
                end
                
            end
            
            time_glm2 = nc{time_var}(:);
            if length(time_glm2)==0
                time_glm2 = nc{'time'}(:);
            end
            
            %close the netcdf file
            nc=close(nc);
            
            if ifile==1
                if cat_dim==0 %set to =0 when there is e.g. one timestep per file and we want to make a big array
                    %with the time in the 1st dimension
                    dat_global.dat = NaN*ones([length(filenames) size(dat)]);
                    dat_global.dat(1,:) = dat(:);                    
                else                                              
                    dat_global.dat = dat;                    
                end
                time_glm = time_glm2;
                
            else                
                if cat_dim==0
                    dat_global.dat(ifile,:) = dat(:);
                else                    
                    dat_global.dat = cat(cat_dim,dat_global.dat,dat);
                end
                time_glm = cat(1,time_glm,time_glm2);
            end
            
        end
        
        
        
        
        %isort_time=1; %add this option in?
        if isort_time==1
%             sdat = size(dat_global.dat);
%             %sdat = 1980         144         192  [time,lat,lon]
%             rep_inds=1;
%             for idim=2:length(sdat)
%                 rep_inds=[rep_inds sdat(idim)];
%             end
%             %time_nd = repmat(time_glm,rep_inds);
%             %[B,I] = sort(time_nd,1);
%             %dat_global.dat = dat_global.dat(I);
%             [time_glm,I]=sort(time_glm);
%             
%             
%             %Get subscripts for the whole 3d array                       
%             [i1,i2,i3] = ind2sub(sdat,[1:prod(sdat)]);
%             %New indices for just the time dim using I values for
%             %sorting - replicate to 3D array
%             i1_new = repmat(I,rep_inds);
%             %Get indices using new values
%             inew = sub2ind(sdat,i1_new(:)',i2,i3);
%             %Use to sort dat and re-shape
%             dat_global.dat = dat_global.dat(inew);
%             dat_global.dat = reshape(dat_global.dat,sdat);  

            if length(size(dat_global.dat))>4
                [dat_global.dat,time_glm]=sort_multidim_array_time(dat_global.dat,time_glm,2); %chunk along 2nd dim
            else
                [dat_global.dat,time_glm]=sort_multidim_array_time(dat_global.dat,time_glm);
            end
        end
        
        if ~exist('time_ref')
            time_ref = datenum('01-Jan-1970');
            time_fconv = 1/24; %conversion multiplier to get to days
        end
        
        time_glm_matlab = time_ref + time_glm*time_fconv;
        dat_global.time_glm = time_glm;
        
        if exist('fill_range')
            inan=find(dat_global.dat>fill_range(1) & dat_global.dat<fill_range(2));
            dat_global.dat(inan)=NaN;
        end
        if exist('grid_data_to_UKESM')
            loadfile = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_clt.mat';
            grid_in = load(loadfile,'gcm_Plon2D_UM','gcm_Plat2D_UM');
            ntime = size(dat_global.dat,1);
            
            dat_new = NaN*ones([ntime size(grid_in.gcm_Plon2D_UM,1) size(grid_in.gcm_Plon2D_UM,2)]);
            for itime=1:ntime
                fprintf(1,'Processing %i of %i\n',itime,ntime);
                d=griddata(gcm_Plon2D_UM,gcm_Plat2D_UM,squeeze(dat_global.dat(itime,:,:)),grid_in.gcm_Plon2D_UM,grid_in.gcm_Plat2D_UM);
                dat_new(itime,:,:) = d;
            end
            
            gcm_Plon2D_UM = grid_in.gcm_Plon2D_UM;
            gcm_Plat2D_UM = grid_in.gcm_Plat2D_UM;
            dat_global.dat = dat_new;
        end
        
        
        switch run_type
            case 'nested ignore lat lon'
                dat_global.gcm_Plat2D_edges_UM=NaN;
                dat_global.gcm_Plon2D_edges_UM=NaN;
            otherwise
                %[dat_global.gcm_Plon2D_edges_UM,dat_global.gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
                [dat_global.gcm_Plat2D_edges_UM,dat_global.gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        end
        
        time_shift=0;
        
        
        time_LST = time_glm_matlab + time_shift;
        
        [Y,M,D,HH,MM,SS]=datevec(time_LST);
        if SS>59
            SS=0;
            MM=MM+1;
        end
        dat_global.time_ALL = datenum(Y,M,D,HH,MM,SS);
        
        dat_global.gcm_Plat2D_UM = gcm_Plat2D_UM;
        dat_global.gcm_Plon2D_UM = gcm_Plon2D_UM;
        
        switch run_type
            case 'nested ignore lat lon'
                
            otherwise                
                [dat_global.gcm_area_UM] = calc_area_lat_lon2d(dat_global.gcm_Plat2D_edges_UM,dat_global.gcm_Plon2D_edges_UM);
        end
        

        
%% cam6     
    case 'cam6'
        [var_UM,scale_fac]=convert_UM_varnames_cam(var_UM);
        
        filenames = dir([dirUM '/CAM6*.nc']);
        filename=[dirUM '/' filenames(1).name];
        nc = netcdf(filename);
        if time_only==1
            dat_global.dat = NaN;
        else
            dat_global.dat = nc{var_UM}(:);
        end
        
        lat_var = 'lat'; lon_var = 'lon';
        
        %Slices have the lat lon converted from rotated coords for the nest
        %already.
        lat1D = nc{lat_var}(:); %[192 1]
        lon1D = nc{lon_var}(:); %[288 1]
        
        [gcm_Plon2D_UM, gcm_Plat2D_UM] = meshgrid(lon1D,lat1D);
        
        % N.B. - Values from Python scripts are already rotated to proper lat and lons!
        %         switch run_type
        %             case 'nested'
        %                 [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(gcm_Plat2D_UM,gcm_Plo
        %                 n2D_UM,pole_lat,pole_lon);
        %        end
        
        
        i180 = find(gcm_Plon2D_UM>180);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;
        
        
        
        
        time_glm = nc{'time'}(:);
        if length(time_glm)==0
            time_glm = nc{'time'}(:);
        end
        time_glm = time_glm*24; %convert to hours
        time_glm_matlab = datenum('01-Mar-2009') + time_glm/24;
        
        [dat_global.gcm_Plon2D_edges_UM,dat_global.gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
        
        
        time_shift=0;
        
        
        time_LST = time_glm_matlab + time_shift;
        
        [Y,M,D,HH,MM,SS]=datevec(time_LST);
        if SS>59
            SS=0;
            MM=MM+1;
        end
        dat_global.time_ALL = datenum(Y,M,D,HH,MM,SS);
        
        dat_global.gcm_Plat2D_UM = gcm_Plat2D_UM;
        dat_global.gcm_Plon2D_UM = gcm_Plon2D_UM;
               
   
        
        
end

%dat_global.dat=dat_global.dat*scale_fac;


function [varname_out,scale_fac]=convert_UM_varnames_cam(var_UM)

scale_fac=1;

switch var_UM
    case 'SW_down_clean_surf'
        varname_out = 'FSDS';
        
    case 'SW_down_clean_clear_surf'
        varname_out = 'FSDSC_d1';
        
    case 'low_cloud_amount'
        varname_out = 'CLDLOW';
        
    case 'LWP_sec30'
        varname_out = 'TGCLDLWP';
        
    case {'Nd_lwc_weighted_UKCA','Nd_lwc_in_cloud_weighted_UKCA'}
        varname_out = 'NUMAVG';        
        
    otherwise
        varname_out = var_UM;
        
        
end

