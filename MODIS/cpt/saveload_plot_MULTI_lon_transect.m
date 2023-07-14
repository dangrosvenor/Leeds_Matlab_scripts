% LOAD and PLOT multiple lon transects (at different lats)
% *** saveload_lon_transect_data *** is the main script executed from here
%   - choose the models to load from there. This is also where the variable
%   name to be loaded can be changed to something different for particular
%   datasets (e.g. LTS1000 for ECMWF becomes LTS_daily)
% Data is stored in filename_lon_save file
% So can do whos('-FILE',filename_lon_save) to check the variables
% Or can be more specific in the variable name, e.g.
% whos('-FILE',filename_lon_save,'Precip*rate*CAM*CLUBBv2*')
% Uses case 116 of watervap
  lon_load_multi={'30S','20S','10S','10to30S'};
%    lon_load_multi={'30S','20S','10S'};
    lon_load_multi={'10to30S'};
%    lon_load_multi={'20S'};


 time_of_day_multi = {'daytime','nighttime'};
 time_of_day_multi = {'average'};
  time_of_day_multi = {'daytime'};
  
%        var_load_multi = {'low_CF','Nd','LWP','Precip_rate','LWP2'};
%        var_load_multi = {'Precip_rate'};
        var_load_multi = {'low_CF'};        
%        var_load_multi ={'LWP2'}; %use this for AMSRE - grid-box mean
        var_load_multi ={'TLWP'}; %actually better to use this for AMSRE  (if available)  
                                   %is grid-box mean LWP+RWP for comparison to AMSRE
                                   
           var_load_multi ={'LWP_MOD35CF'}; %LWP_MOD35CF is the MODIS LWP divided by MOD35 CF                       
%        var_load_multi ={'Precip_rate'};
%        var_load_multi ={'Nd'}; 
%         var_load_multi ={'Reff_COSP'};
%        var_load_multi ={'LTS'};         
%        var_load_multi ={'LTS1000'};                 
%        var_load_multi ={'qv700'};                 
%         var_load_multi={'Reff_COSP_allCF','Reff_COSP','REFFL_maxlayer','REFFL_max_noCF','REFFL_maxliq'};
%         var_load_multi={'Reff_COSP_allCF'};
%{'Reff_COSP','Tau_COSP','Nd_COSP'};
%         var_load_multi={'LWP'};
%         var_load_multi={'LWP_COSPCF'};
%         var_load_multi={'LWP_COSP'};         

season_str_multi={'ALL'}; %or could be {'ALL'} or {'DJF'}, etc.
%season_str_multi={'DJF','MAM','JJA','SON','ALL','ASON'};

icalc_diurnal=1;

if icalc_diurnal==1
  clear value_array
end



for ilon_multi=1:length(lon_load_multi)    
    for ivar_multi=1:length(var_load_multi)
        for itime_multi=1:length(time_of_day_multi)
            for iseason_multi=1:length(season_str_multi)
                
                season_str = season_str_multi{iseason_multi};


            time_of_day = time_of_day_multi{itime_multi};
            lon_load = lon_load_multi{ilon_multi};
            var_load = var_load_multi{ivar_multi};

            figname = [var_load ' ' lon_load];

            imulti_lonplot_override=1;
            
            %% ---------------   main script ---------------------------------------            
                saveload_lon_transect_data          
            %% ---------------------------------------  
            
          
            
            
            man_choose_water_graph=1;
            %% settings for watervap ---------------------------------------
            graph=116;
            time_highlight_path=[];

            iytick_relabel=0; %flag to say whether to relabel the y-axis ticks (e.g. for log plot)
            y_axis_type=''; %default
            x_axis_type='';
            i_set_dateticks=0;
            iadd_nums_above=0;
            iovr_leg_line=0; %flag which plots a '-' line in the legend even if linestyle has been set to 'none'

            xlims=0;
            fsize=22;

            idatetick=0; %flag to put times as HH:MM instead of decimal time

            noplot=0;

            iwrite_text_dat=0; %

            ichoose_styles=0; %flag to say whether we want to specifiy the specific line patterns and colours
            line_pattern(1).p=NaN;
            line_colour(1).c=NaN;
            marker_style(1).m=NaN;
            
             %% end of settings for watervap ---------------------------------------
    
             if itime_multi==1
                 scrsz=get(0,'ScreenSize');
                 posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];
                 posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.8];
                 %posit=[9 50 scrsz(3)/2.5 scrsz(4)/2.5];
                 posit=scrsz;
                 hf=figure('name',figname,'Position',posit);
             end
             subplot(2,1,itime_multi);
             ioverride_watervap_newfig=1;
             waterVapourMay2005
             
             
             if icalc_diurnal==1
                 %               if itime_multi==1
                 %                  time0_dat = ydat;
                 %              else
                 for idat=1:length(xdat)
                     value_array{ilon_multi,ivar_multi,itime_multi}(idat).x = xdat(idat).x;
                     value_array{ilon_multi,ivar_multi,itime_multi}(idat).y = ydat(idat).y;  
                     value_array{ilon_multi,ivar_multi,itime_multi}(idat).lab = labs(idat).l;                       
                 end
                 %             end

             end
             
             
         
            

            

        
        %make it full-screen and save
%        set(gcf,'Position',scrsz);
%        saveas_ps_fig_emf(gcf,savename,'',0);
        savename2 = [savedir 'CPT/' var_load '_' lon_load '_' time_of_day '_' season_str];
        saveas_ps_fig_emf(gcf,savename2,'',0,1); %0 at the end prevents the date being added to the filename
        
            end
        
            end
        



     end

 end



