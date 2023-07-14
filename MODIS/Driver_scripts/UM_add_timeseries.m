%script to add timeseries together, but taking into account that some may
%be different lengths

ydat_import_incs(idat_driver).y = 0;
ydat_import_flux(idat_driver).y = 0;
ydat_import_flux2(idat_driver).y = 0;


%Read a first time to get the lengths of all the timeseries
minL=9e9;
 for ivar_aero=1:length(VAR_NAMES)
     file_tim = [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME',[VAR_NAMES{ivar_aero} '']) '_timeseries.mat'];
     if i_restricted_domain==1
         file_tim = remove_character(file_tim,'ukv_.pp.nc','ukv_restricted_domain_.pp.nc');
     end     
     dat_UM_timser = load_UM_timser(file_tim,append_str_timser{idat_UM});
     minL = min(minL,length(dat_UM_timser.timeseries_UM));
 end
 
 %Read again and sum together the ones required
 for ivar_aero=1:length(VAR_NAMES)
     file_tim = [dirUM_i remove_character(fileUM{idat_UM},'VAR_NAME',[VAR_NAMES{ivar_aero} '']) '_timeseries.mat'];
     if i_restricted_domain==1
         file_tim = remove_character(file_tim,'ukv_.pp.nc','ukv_restricted_domain_.pp.nc');
     end     
     [dat_UM_timser flux_UM_timser] = load_UM_timser(file_tim,append_str_timser{idat_UM}); %Also load in fluxes here
     
     if length(flux_UM_timser.timeseries_UM)>0
         dt=diff(dat_UM_timser.time_UM)*24*3600;
         dx=0.009*111.3195e3; %metres
         tot_flux_cum_mean = dt(1)*cumsum(flux_UM_timser.timeseries_UM(1:minL,1)) / (dx*dx*600*600); %Total change due to fluxes over dt averaged over the domain
         %area - cumsum these since then we will obtain the values we need to subtract as a function of time.
         tot_flux_mean = dt(1)*flux_UM_timser.timeseries_UM(1:minL,1) / (dx*dx*600*600); %Increments averaged over the domain.
     end
     
     
%     ydat_import_flux(idat_driver).y = ydat_import_flux(idat_driver).y + tot_flux_cum_mean;
%     ydat_import_flux2(idat_driver).y = ydat_import_flux2(idat_driver).y + tot_flux_mean;
     

 time_matlab = dat_UM_timser.time_UM(1:minL);
 
 %Sum the contributions from all aerosol species.
 if iplot_flux_contribution==1
     ydat_import(idat_driver).y = ydat_import(idat_driver).y + tot_flux_cum_mean;
     %ydat_import(idat_driver).y = ydat_import_flux2(idat_driver).y + tot_flux_mean;

 elseif iplot_flux_increments==1
     %Moving average over 3 blocks for central difference
     windowSize=3;         
     [time_matlab,mov_av] = window_average(dat_UM_timser.time_UM(1:minL),tot_flux_mean,windowSize,'mean');
     
     ydat_import(idat_driver).y = ydat_import(idat_driver).y + windowSize*mov_av; %multiplying by the windowSize since are doing increments over 3 times
      %whereas dt above was just for one time increment
     
     diffs = dat_UM_timser.timeseries_UM(3:minL) - dat_UM_timser.timeseries_UM(1:minL-2);
     ydat_import_incs(idat_driver).y = ydat_import_incs(idat_driver).y + diffs;
     

 else
     ydat_import(idat_driver).y = ydat_import(idat_driver).y + dat_UM_timser.timeseries_UM(1:minL);
 end

     
 end

 
 if inormalise_initial_absolute==1
     ydat_import(idat_driver).y = ydat_import(idat_driver).y - ydat_import(idat_driver).y(2);
 end
 if inormalise_initial==1
     ydat_import(idat_driver).y = ydat_import(idat_driver).y ./ ydat_import(idat_driver).y(2);
 end
 