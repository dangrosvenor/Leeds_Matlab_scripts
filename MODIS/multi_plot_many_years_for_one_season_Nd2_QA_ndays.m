%Load in the timeseries3 data for each year and plot one (or several) month(s) for each year and plot one (or several) month(s) for each
%year
%% N.B. The region loaded in is restricted in load_saved_modis_vars_20150130 and not
% overridden by the driver script as yet

iload_DRIVER=1;  %Setting to =0 prevents loading (for testing purposes, e.g. if have alrady loaded)
i_ocean_only_SCREEN = 0;
iplot_DRIVER=0; %whether to do the big multi-figure sub-plot or not.

%The plan here is to just be able to change the var selection below to use
%this same code for any variable.
var_str_DRIVER = 'Nd';
%var_str_DRIVER = 'reff';

save_file_DRIVER = ['/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_' var_str_DRIVER '_' datestr(now,30) '.mat'];
%/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150205T055754.mat

isave_plot=0; %This just determines whether to save the individual plots (the one with lots of subplots is automatically saved)

m_subplot=3;
n_subplot=5;

irestrict_domain = 0;
proj_type='global oval'; %select this for anything but polar stereographic
%proj_type='polar'; irestrict_domain = 0; 

stereo_str01 = ['lat_polar= 90; %Arctic'];
stereo_str02 = ['m_proj(''stereographic'',''lat'',lat_polar,''lon'',0,''rad'',40);'];
  %the 'rad' option is how mnay degrees to plot from pole - i.e. if set to
  %50 will plot down to 40N.



switch var_str_DRIVER
    case 'Nd'
%        modis_data_plots={'Number of droplets cell values time mean - specific days'};
        modis_data_plots={'Number of droplets QA pixel weighted cell values time mean - specific days'};        
%        modis_data_plots={'Latest or earliest day after screening'};
        modis_data_plots={'Number of droplets QA cell values time mean - specific days'};
    case 'reff'
%        modis_data_plots={'Cloud Effective Radius timeseries3 mean from selected days'};
        modis_data_plots={'Cloud Effective Radius QA pixel weighted timeseries3 mean from selected days'};        
end

%        zero_CF_out_of_range=0;
time_mean_str ='ALL';
%        cont_dat_choose = 'calipso mid+highCF';

%% Set the screening

set_screening = {'none','NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
set_screening = {'none','NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
set_screening = {'none'};
%set_screening = {'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW'};
set_screening = {'none',...
    'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW',...
    'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW',...
    'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW',...
    'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW',...
    'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff',...    
    'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff',...
    'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff',...
    'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff',...    
    };

%high low SZA, warm/cold

thresh_ndays=100; %threshold no. days
thresh_SZA=[0 90];
thresh_SZA=[0 70];

thresh_SZA_multi={[0 90],[0 90],[0 90],[0 70],[0 70],[0 90],[0 90],[0 70],[0 70]};

thresh_CF=[0.8 1.00001];
thresh_CF=[0.0 1.00001];

thresh_NP=50;
thresh_sensZA=[0 90];
thresh_CTT = [173 273+100];
thresh_CTT = [268 273+100];

minT = 268;
minT = 273.15;
thresh_CTT_multi = {[173 273+100],[173 273+100],[minT 273+100],[173 273+100],[minT 273+100],[173 273+100],[minT 273+100],[173 273+100],[minT 273+100]};


thresh_relAZ = [0 180];
thresh_CTH = [-0.01 3.2];
thresh_CTH = [-0.01 1e9];

thresh_stdW = [0 1e9];
thresh_sigCTT = [0 1e9];
thresh_reff=[0 30];





ifilter_ndays=0;


switch var_str_DRIVER
    case 'Nd'
inew_cticks=1;  %If setting to =1 then set iset_min and max below to =0
iset_min_clim=0;
clim_min=20;
iset_max_clim=0;
%        clim_max=225;
clim_max=150;

%         plot_num_datapoints=1;
%         inew_cticks=0;
%         iset_min_clim=1;
%         clim_min=0;
%         iset_max_clim=1;
%         clim_max=15;

    case 'reff'
          plot_num_datapoints=0;
          inew_cticks=0;
          iset_min_clim=1;
          clim_min=8;
          iset_max_clim=1;
          clim_max=20;
end

%% Choose the dates

%        set_month = {9,10}; %Choose required months
%        set_year = {'y2007_TERRA','y2008_TERRA'}; %

        set_month = {9}; %Choose required months
        set_year = {{'y2014_TERRA','y2014_AQUA'}}; %

% All years 2000 to 2014
% clear set_year
% set_month = {8,9,10}; %Choose required months
% set_year = {{'y2000_TERRA'},{'y2001_TERRA'}}; %Will add to these below
% years=[2002:2014];
% for iy=1:length(years)
%     if years(iy)==2006
%         year_cell = {['y' num2str(years(iy)) '_TERRA'],['y' num2str(years(iy)) '_AQUA_updated']};
%     else
%         year_cell = {['y' num2str(years(iy)) '_TERRA'],['y' num2str(years(iy)) '_AQUA']};        
%     end        
%     set_year{iy+2} = year_cell;
% end

%TERRA only
clear set_year
%set_month = {8,9,10}; %Choose required months
set_month = {1,2,3,4,5,6,7,8,9,10,11,12}; %Choose required months
set_year = {{'y2000_TERRA'},{'y2001_TERRA'}}; %Will add to these below
years=[2002:2014];
for iy=1:length(years)
    year_cell = {['y' num2str(years(iy)) '_TERRA']};    
    set_year{iy+2} = year_cell;
end

%AQUA only
% clear set_year
% set_month = {8,9,10}; %Choose required months
% %set_year = {{'y2000_TERRA'},{'y2001_TERRA'}}; %Will add to these below
% years=[2002:2014];
% for iy=1:length(years)
%         if years(iy)==2006
%             year_cell = {['y' num2str(years(iy)) '_AQUA_updated']};
%         else
%             year_cell = {['y' num2str(years(iy)) '_AQUA']};
%         end
%     set_year{iy} = year_cell;
% end


% set_month = {9}; %Choose required months
% set_year = {{'y2001_TERRA','y2001_AQUA'},{'y2007_TERRA','y2007_AQUA'},{'y2014_TERRA','y2014_AQUA'}}; %

%set_month = {9,10}; %Choose required months
%set_year = {{'y2002_TERRA','y2002_AQUA'}};  %,{'y2007_TERRA','y2007_AQUA'}}; %Load Aqua and Teera together


clear RMSE_DRIVER mean_bias_DRIVER Psave Nsave Nd_std_dev_save
icount_DRIVER=0;

for iloop_year=1:length(set_year)
    icount_subplot=0;
    year_str = remove_problem_chars(set_year{iloop_year}{1}(1:5));

    %% Load in the data for the required year

    %Can prevent loading to save time when testing here
    if iload_DRIVER==1

        i_multi_year_av=0; %flag for load_saved_modis_vars.m
        data_type='L3 processed data';
        override_loadsave=1;
        L3_case = 'override';
        savemem=-1;
        modis_vars_script_name='modis_var_L3_variables_usual_screening';

        clear modis_data_case
        s='Y'; modis_data_case_choose = set_year{iloop_year};

        %Run the load script
        load_saved_modis_vars_20150130  %modified version with date so that the old one is preserved.

    end

    for iloop_plot=1:length(modis_data_plots)
        modis_data_plot = modis_data_plots{iloop_plot};

        for iscreen_set=1:length(set_screening)
            screen_type = set_screening{iscreen_set};
            thresh_SZA = thresh_SZA_multi{iscreen_set};
            thresh_CTT = thresh_CTT_multi{iscreen_set};

            for iloop_MOD=1:length(set_month)
                month_modis = set_month{iloop_MOD}
                icount_DRIVER = icount_DRIVER+1;
                icount_subplot = icount_subplot+1;
                
                
                              
                %Pick out the correct day nunbers of the month in
                %question
                [days_required_for_mean,time_mean_str] = days_of_month(month_modis);
                
                

               

%% Now plot
 


                ioverride_plotglobal_thresh=1;
                %                      iocean_only=1;
                ioverride_time_selection=1;
                ioverride_plotglobal_loc=1;
                mod_data_type='timeseries3 lambert';
                data_select='specific_modis_data';
                plot_global_maps
%                Psave{icount_subplot}{iloop_year} = P_save;
%                Nsave{icount_subplot}{iloop_year} = Npoints;                
                Psave{icount_subplot} = P_save;
                Nsave{icount_subplot} = Npoints;
                Nd_std_dev_save{icount_subplot} = Nd_std_dev;
                %Make the title here as also use for saving
                mtitle{icount_subplot} = [time_mean_str ' ' thresh_str];
                
                %                      mean_bias_DRIVER(icount_DRIVER)=mean_bias_MODIS;
                %                      RMSE_DRIVER(icount_DRIVER)=RMSE_MODIS;

                if isave_vals==1
                    %                     fprintf(fid_mtab,'%f ',mean_bias_MODIS);
                    %                     fprintf(fid_Rtab,'%f ',RMSE_MODIS);
                end


                if isave_plot==1
                    saveas_ps_fig_emf(gcf,[savename],'',0,1);
                end
                
                ax_save{icount_DRIVER} = gca;
                
                if iplot_DRIVER==1
                    %Make the big subfigure window
                    if iloop_year==1
                        hfig_sub{icount_subplot} = figure;
                        set(hfig_sub{icount_subplot},'units','centimeters');
                        pos_fig = [0 0 50 30];
                        set(hfig_sub{icount_subplot},'position',pos_fig);

                        %Maximise the figure window
                        %                    pause(0.00001);
                        %                    frame_h = get(handle(gcf),'JavaFrame');
                        %                    set(frame_h,'Maximized',1);

                        %                    set(hfig_sub{icount_subplot},'units','normalized');

                    end




                    %switch to the subplot figure
                    figure(hfig_sub{icount_subplot});

                    ax_subplot = subplot(m_subplot,n_subplot,iloop_year);
                    pos=get(ax_subplot,'Position');
                    delete(ax_subplot);
                    hax2=copyobj(ax_save{icount_DRIVER},hfig_sub{icount_subplot});
                    set(hax2, 'Position', pos);

                    %                ax_subplot = subplot(m_subplot,n_subplot,icount_DRIVER);
                    %                copyobj(allchild(ax_save{icount_DRIVER}),ax_subplot);
                    %                copyobj(allchild(ax_save{icount_DRIVER}),ax_subplot);
                    %  Year label


                    title(year_str);

                    %Remove the x and y labels
                    m_ungrid
                    set(gca,'xticklabels','');
                    set(gca,'yticklabels','');



                    if iloop_year==length(set_year)
                        %Make the title for the whole of the subplot figure
                        %This was causing quite a lot of trouble - adding a
                        %subplot axis makes the title disappear, so do this at
                        %the end.
                        %                    mtit(hfig_sub{icount_subplot},mtitle{icount_subplot},'xoff',0,'yoff',0.03,'fontsize',16);
                        t=textwrap({mtitle{icount_subplot}},250);

                        %                    pos=get(hfig_sub{icount_subplot},'position');
                        pos = [pos_fig(1) pos_fig(2) pos_fig(3) pos_fig(4)-3];
                        ax=axes('units','centimeters');
                        set(ax,'position',pos);
                        set(ax,'visible','off');
                        %                    ht=title(t,'position',[0.5 0.96 1.0]);
                        %                    ht=title(t,'position',[pos_fig(3)/2 pos_fig(4)-1 1.0]);
                        ht=title(t)
                        set(ht,'visible','on');

                        %copy the colorbar
                        hc2=copyobj(hc,hfig_sub{icount_subplot});
                        pos=get(hc2,'position');
                        set(hc2,'position',[pos(1) 0.05 pos(3) 0.03]);
                    end

                    close(hf); %close the non-subplot window

                    if iloop_year==length(set_year)
                        saveas_ps_fig_emf(gcf,[savedir mtitle{icount_subplot}],'',0,1);
                    end

                else
                    close(hf); %close the non-subplot window
                end %if iplot_DRIVER==1

            end

        end

    end
    
    
    
    
   if iloop_year==1
       save_var_DRIVER_multi_year(save_file_DRIVER,1,MLAT,'MLAT');
       save_var_DRIVER_multi_year(save_file_DRIVER,0,MLON,'MLON');
       save_var_DRIVER_multi_year(save_file_DRIVER,0,mtitle,'mtitle');
       save_var_DRIVER_multi_year(save_file_DRIVER,0,modis_data_plots,'modis_data_plots');       
   end
   
   save_var_DRIVER_multi_year(save_file_DRIVER,0,Psave,['Psave_' year_str]);
   save_var_DRIVER_multi_year(save_file_DRIVER,0,Nsave,['Nsave_' year_str]);
   save_var_DRIVER_multi_year(save_file_DRIVER,0,Nd_std_dev_save,['Nd_std_dev_save_' year_str]);
   
end %for iloop_year=1:length(set_year)
%% End of loops

% save(save_file_DRIVER,'Psave','Nsave','mtitle','MLAT','MLON');

if isave_vals==1
    fprintf(fid_mtab,'\n');
    fprintf(fid_Rtab,'\n');
end



if isave_vals==1
    fclose(fid_mtab);
    fclose(fid_Rtab);
end

%        save(save_file,'modis_data_plots','set_screening','set_month','set_year','Psave','mean_bias_DRIVER','RMSE_DRIVER')

