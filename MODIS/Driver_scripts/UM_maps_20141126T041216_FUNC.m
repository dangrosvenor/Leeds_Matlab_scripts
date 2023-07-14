function [vars_map_out]=Boutle_case_12Nov2008_LWP_maps_20141126T041216_FUNC(vars_map)

% LWP map plots for all times for UM

%% Take all the input variables from the var_maps input structure
% names = fieldnames(vars_map);
% for i=1:length(names)
%     eval_str = [names{i} ' = vars_map.' names{i} ';'];
%     eval(eval_str);
% end

name_struc='vars_map'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end

% Also copy all to vars_in
vars_in = vars_map;

% %Get the names of the UM files/runs required
% [vars_UM_cases]=UM_case_select_runs(UM_cases);
% names = fieldnames(vars_UM_cases);
% for i=1:length(names)
%     eval_str = [names{i} ' = vars_UM_cases.' names{i} ';'];
%     eval(eval_str);
% end


% %varname refers to the top level variable names - this may require several
% %different fields, each with a different nc_varname (as labelled in the
% %NetCDf files). Need to add the option to supply several variable names,
% %though. Or could deal with them in this function
% switch varname
%     case 'SW TOA out'
%         nc_varname = 'SW_TOA_out'; %The name of the variable in the nc file
%         VAR_NAME_STR='SWLW_TOA_outgoing'; %The name part in the filename to replace VAR_NAME with
%         %OR, if the VAR_NAME is not in the supplied filename then will
%         %just use the supplied one - make sure it is correct!
%     case 'LW TOA out'
%         nc_varname = 'LW_TOA_out'; %The name of the variable in the nc file
%         VAR_NAME_STR='SWLW_TOA_outgoing'; %The name part in the filename to replace VAR_NAME with
%         %OR, if the VAR_NAME is not in the supplied filename then will
%         %just use the supplied one - make sure it is correct!
% end



if ~exist('optional_saveas')
    optional_saveas.dummy=NaN;
end

if ~exist('iz')
    iz=-1;
end
if ~exist('icoarsen')
    icoarsen=0;
end


if i_mask_low_LWP==1
    fprintf(1,'\n*** WARNING i_mask_low_LWP is set to one - check this is requried!! ***\n\n');
end

%% Do the plotting
%for idat_UM=length(fileUM)
for itemp=1:1  %Deal with this in UM_maps_generic_time_loop now.
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
    
    
    
    %          vars_in.var = nc_varname; %can set to the variable name to just read the variable
    % %         vars_in.flag = flag{idat_UM};
    %          vars_in.flag = flag2;          %replace the old flag system - flag likely to go with variable rather than UM run
    %          vars_in.file_lwp =  remove_character([dirUM fileUM{idat_UM}],'VAR_NAME',VAR_NAME_STR);
    %          if strcmp(flag2,'load_mat')==1 & length(strfind(vars_in.file_lwp(end-3:end),'.mat'))==0
    %              vars_in.file_lwp = [vars_in.file_lwp '.mat'];
    %          end
    % %         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
    %          vars_in.pole_lat = pole_lat;
    %          vars_in.pole_lon = pole_lon;
    %          vars_in.time_in = time_range; %set to [] for all times
    vars_in.iz = iz;
    
    
    [dat_UM,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz] = get_LWP_RWP_UM(vars_in);
    
    
    % Can add extra load calls here if need more than one variable for
    % something.
    nvars_out = 1; %default - changed below for cases where output more than one var
    switch varname
        %     case 'SW_down_surf_LWP_LT_0pt1'
        %         nvars_out=nvars_out+1;
        %         vars_in.var = 'LWP';
        %         [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz] = get_LWP_RWP_UM(vars_in);
        %         ilwp = find(lwp>0.1);
        %         dat_UM(ilwp)=NaN; %Make NaN any non-clear grid points
        %    case 'Transmission_down_surf_LWP_LT_0pt1'  --- decided to move to
        %                                                   get_LWP_RWP_UM
        %         nvars_out=nvars_out+1; %not sure if need this?
        %
        %         vars_in.file_lwp
        %         vars_in.var = 'LWP';
        %         [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz] = get_LWP_RWP_UM(vars_in);
        %         ilwp = find(lwp>0.1);
        %         dat_UM(ilwp)=NaN; %Make NaN any non-clear grid points - now have the clar-sky SW_down at surf
        %         %Get SW_down_TOA
        %         vars_in.var = 'SW_down_TOA';
        %         [SW_toa,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,z_um,iz] = get_LWP_RWP_UM(vars_in);
        %         dat_UM = dat_UM ./ SW_toa; %Transmission fraction of incoming TOA to surface
        %
        %         vars_in.var = 'Transmission_down_surf_LWP_LT_0pt1';
    end
    
    
    %UM global data runs from zero longitude and back round, which seems
    %to cause issues with the plotting. Think it is the fact that the lons are not in order
    % (e.g. flip from 180 to -180 at the date line) that causes problems
    %Will re-arrange the data to run from -180 to +180.
    %Should be able to safely do this for LAM runs too since just puts the
    %data in lon order
    
    %N.B. can also have in this order if wanted:- i.e. from 0 to 360 -
    %      gcm_Plon2D_UM(gcm_Plon2D_UM<0) =  gcm_Plon2D_UM( gcm_Plon2D_UM<0)+360;
    %just need to re-do the edges
    
    [Y2,I2]=sort(gcm_Plon2D_UM(1,:));
    gcm_Plon2D_UM = gcm_Plon2D_UM(:,I2);
    gcm_Plat2D_UM = gcm_Plat2D_UM(:,I2);
    if length(size(dat_UM))==2
        dat_UM = dat_UM(:,I2);
    elseif length(size(dat_UM))==3
        dat_UM = dat_UM(:,:,I2);
    end
    
    %re-do the edges, so that the edge is at -180
    [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
    
    
    if icoarsen==1
        %need to supply dlat_targer and dlon_target
        %         dlat_target = dlat_CERES;
        %         dlon_target = dlon_CERES;
        
        %average to the coarser resolution of goes
        
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_target/dlat_UM));
        
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_target/dlon_UM));
        
        
        
        
        dat_UM= reduce_matrix_subsample_mean(dat_UM,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
    end
    
    vars_map_out.z_um{idat_UM} = z_um;
    
    switch UM_time_loop_action
        case 'output_3D'
            vars_map_out.P_save{idat_UM} = dat_UM;
            vars_map_out.P_mean{idat_UM} = meanNoNan(dat_UM(:),1);
            vars_map_out.time_str = ''; %time_str;
            vars_map_out.z_str = ''; %z_str;
        otherwise
            
            %% Make the plot if needed
            
            %if iplot_maps==1
            
            %% select time index here
            for it_driver=1:length(it_driver_out) %22:22 %1:length(time_matlab)  %4:nt_driver
                %12=23:00 UTC on 12th, 16 = 07:00 on 13th;  22=19:00 on 13th;
                %--- run the file to set up the defaults
                plot_global_maps_defaults
                
                %re-do this since they will be overwritten in plot_global_maps_defaults
                %This will mean that if it is not supplied in vars_map then the
                %default will be used, so don't need to worry about supplying
                %all flags unless needed to change from default.
                names = fieldnames(vars_map);
                for i=1:length(names)
                    eval_str = [names{i} ' = vars_map.' names{i} ';'];
                    eval(eval_str);
                end
                
                if noplot==1
                    inew_figure=0;
                end
                
                irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not
                
                thresh_LAT = LAT_val_DRIVER;
                thresh_LON = LON_val_DRIVER;
                
                %--- set some options for these particular plot loops
                set_screening = {'none'};
                modis_data_plot = 'Map of 2D data from outside driver script';
                
                %         iset_min_clim=1;
                %         clim_min=0;
                %         iset_max_clim=1;
                %         clim_max=300;
                %
                %         isave_plot=0;
                %         iplot_markers=0;
                
                %Calculate the data to plot
                
                if length(size(dat_UM))==2
                    dat_modis = dat_UM; %already in g/m2
                else
                    dat_modis = squeeze(dat_UM(it_driver,:,:)); %already in g/m2
                end
                
                if i_mask_low_LWP==1
                    dat_modis(dat_modis<thresh_LWP_mask)=NaN;
                end
                
                %Set various things
                time = time_matlab(it_driver);
                %Round to the nearest minute as sometimes get 18:59:59
                time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM');
                if iz>-1
                    z_str = ['at ' num2str(z_um(iz)/1000,'%.2f') ' km '];
                else
                    z_str ='';
                end
                
                titlenam_driver = [varname ' for ' time_str ' for ' labs_UM(idat_UM).l ' ' z_str];
                units_str_plot = var_units_str;
                
                mod_data_type='AMSRE';
                gcm_str_select='UM';
                %        daynum_timeseries3_UM = [1:length(time)];
                %        gcm_time_UTC_UM = [1:length(time)];
                
                %        gcm_Plon2D_UM = lons2D_driver;
                %        gcm_Plat2D_UM = lats2D_driver;
                
                %        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(lats2D_driver,lons2D_driver);
                
                
                
                %        i_dpcolor=1;
                ifull_swath=0;
                igcm_screen=0;
                
                
                
                
                %--- Apply override flags
                ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
                % iocean_only=1;
                ioverride_time_selection=0; %Override the times to include
                ioverride_plotglobal_loc=1; %Override the location of the plot window
                ioverride_years_time_screen=0; %Override years for screening?
                iover_ride_plot_global=1; %overrides inew_figure=1; supress_colorbar=0; i_increase_font_size_map_figures_OFF = 0;
                %(all set in plot_global_maps_defaults)
                
                %---  Run plot script and save
                plot_global_maps
                %-------------------------------
                if isave_plot==1
                    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,optional_saveas);
                    %            close(gcf);
                end
                %
                %         if imake_anim==1
                %
                %             if ==1
                %                 mov = avifile('/home/disk/eos1/d.grosvenor/modis_work/plots/MOVIE_LWP_xlhgw.avi','fps',1);
                %             end
                %
                %
                %             animfr(jj)=getframe(hh(jj).h);
                %             mov = addframe(mov,animfr(jj));
                %
                %
                %         end
                
            end
            
            if iplot_markers==1
                
                lat_mark01 = -20; lon_mark01 = -75;  %The ship
                %        [ilat,ilon] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,lat_mark01,lon_mark01,0.1);
                m_plot(lon_mark01,lat_mark01,'wo','markersize',15,'markerfacecolor','k');
                
                lat_mark02 = -20.8543; lon_mark02 = -76.4911; %The location where the UM matches pretty well
                m_plot(lon_mark02,lat_mark02,'w^','markersize',15,'markerfacecolor','k');
                
                
                
            end
            
            
            vars_map_out.P_save{idat_UM} = P;
            vars_map_out.P_mean{idat_UM} = Pmean;
            vars_map_out.time_str = time_str;
            vars_map_out.z_str = z_str;
    end
    
    
    
    % else
    %     vars_map_out.P_mean{idat_UM} = meanNoNan(datUM(:),1);
    %end
    
    
    
end


vars_map_out.gcm_Plat2D_edges_UM = gcm_Plat2D_edges_UM;
vars_map_out.gcm_Plon2D_edges_UM = gcm_Plon2D_edges_UM;
vars_map_out.gcm_Plat2D_UM = gcm_Plat2D_UM;
vars_map_out.gcm_Plon2D_UM = gcm_Plon2D_UM;




