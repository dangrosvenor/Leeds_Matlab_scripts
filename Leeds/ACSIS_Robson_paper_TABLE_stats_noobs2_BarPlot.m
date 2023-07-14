% Main script to plot bar plots 
% Is run from :-
%    ACSIS_Robson_paper_TABLE_stats_noobs3.m
% Runs :- 
%    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script

%% Plot bar plots
vars={'TS','SW','CF','Nd','LWPic','LWP'};
ylabs={'\DeltaT (K)','\DeltaF_{SW\uparroa} (W m^{-2})','\Deltaf_c','\DeltaN_{d} (cm^{-3})','\DeltaL (g m^{-2})','\DeltaL_{all-sky} (g m^{-2})'};
ilocal=[0 1 1 0 1 0];

vars={'TS','SW','CF','Nd','LWPic'};
ylabs={'\DeltaT (K)','\DeltaF_{SW\uparrow} (W m^{-2})','\Deltaf_c','\DeltaN_{d} (cm^{-3})','\DeltaL_{ic} (g m^{-2})'};
ilocal=[0 1 1 0 1];

vars={'SW'}; ilocal=[1]; ylabs={'\DeltaF_{SW\uparrow} (W m^{-2})'};
%vars={'TS'}; ilocal=[0]; ylabs={'\DeltaT (K)'};
vars={'LWPic'}; ilocal=[1]; ylabs={'\DeltaL (g m^{-2})'};
% vars={'Nd'}; ilocal=[0]; ylabs={'\DeltaN_{d} (cm^{-3})'};
 vars={'CF'}; ilocal=[1]; ylabs={'\Deltaf_c'};
%vars={'od550aer'}; ilocal=[0]; ylabs={'\Delta\tau_{a}'};
%vars={'od550tot'}; ilocal=[0]; ylabs={'\Delta\tau_{a}'};

for ivar=1:length(vars)
    
    var_str_tab=vars{ivar};
%     switch var_str_tab
%         case {'SW','TS'}            
%             stacked_bar_or_error_bars = 'error_bar'; %new style suggested by Ken.
%         otherwise
%             stacked_bar_or_error_bars = 'stacked_bar';
%     end


    stacked_bar_or_error_bars = 'error_bar'; %new style suggested by Ken.
    
    run_str='';
    %var_str_tab='CF';
    hist_str='';
    
    %Period 1
    itr=1; %specifies the second trend in the cell
%    i0=find(trends.years_ukesm_1d==yr_start_trend_box2(itr));
    period_str='PA';
    period_lab = '1850-1970';
    
    
    
    %run script to put data into bar_dat correctly
    clear bar_dat bar_dat_UN bar_dat_sig xlabs_str
    ibar=0;
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script
    %     aer_sig_PA = aer_sig;
    %     ghg_sig_PA = ghg_sig;
    %     nat_sig_PA = nat_sig;
    %     local_sig_PA = local_sig;
    %     had_sig_PA = had_sig;
    %     net_sig_PA = net_sig;
    %     had_sig_PA = had_sig;
    
    
    %Period 2
    period_str='PB';
    period_lab = '1971-2014';
    
    %run script to put data into bar_dat correctly
    %clear bar_dat xlabs_str
    %ibar=1;
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script
    %     aer_sig_PB = aer_sig;
    %     ghg_sig_PB = ghg_sig;
    %     nat_sig_PB = nat_sig;
    %     local_sig_PB = local_sig;
    %     had_sig_PB = had_sig;
    %     net_sig_PB = net_sig;
    
    switch stacked_bar_or_error_bars
        case 'stacked_bar'
            savename=[savedir_date  stacked_bar_or_error_bars '_PLOT_DAMIP_Region4_' var_str_tab ' ' land_ocean];
            %savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
            clear opts
            %        opts.iplot_png=1;
            opts.iplot_eps=1;
            
            savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
            
    end
    
    
    %Period 3    
    stacked_bar_or_error_bars = 'none';    
    
    itr=3; %specifies the second trend in the cell
    %i0=find(trends.years_ukesm_1d==yr_start_trend_box2(itr));
    period_str='PC';
    period_lab = '1985-2014';        
    
    %run script to put data into bar_dat correctly
    clear bar_dat bar_dat_UN bar_dat_sig xlabs_str
    ibar=0;
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script
    
    
end





