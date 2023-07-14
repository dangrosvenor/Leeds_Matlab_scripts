%Loads the data from the e.g. full minus CF_constant trend line calculations
%and makes the table and bar plot.
%Need to have run
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020_CALCs.m and
% ACSIS_Robson_paper_offline_SW_calcs_Sep2020.m
% first to save the data in .mat files (and they also make the timerseries plots).

stacked_bar_or_error_bars = 'stacked bar';
stacked_bar_or_error_bars = 'error bars';

period_labs = {'1850-1970','1971-2014'};
period_strs={'P1','P2'}; period_str_chars={'PA','PB'};


inew_folder=0;
if inew_folder==1
    savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
    eval(['!mkdir ' savedir_date]);
else
    savedir_date='/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
end

ylims_bars = [-8 7]; %will prob need to change for each plot

land_ocean = 'ocean_only,_no_sea-ice';
ens_str = 'ens_mean';

filename_prestr = 'Annual_mean_SW_TOA_up_calc_for_region_4_';
titlenam = filename_prestr;

var_str_tab = 'SWcalc';

% set models
MIP = 'CMIP'; DAMIP_runs = {''};
%MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-aer','DAMIP_hist-GHG','DAMIP_hist-nat'};
%MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-GHG'};
%MIP = 'DAMIP'; DAMIP_runs = {'DAMIP_hist-aer'};
%MIP = 'HADGEM3_GC31_LL'; DAMIP_runs = {'HADGEM'}; DAMIP_runs2 = DAMIP_runs;
%MIP = 'AerChemMIP'; DAMIP_runs = {'UKESM1-AerChemMIP_control'}; DAMIP_runs2 = DAMIP_runs
%MIP = 'AerChemMIP'; DAMIP_runs = {'AerChemMIP_Aerosol_Proxy'}; DAMIP_runs2 = DAMIP_runs
%MIP = 'AerChemMIP'; DAMIP_runs = {'AerChemMIP_hist-piAer'}; DAMIP_runs2 = DAMIP_runs

yr_start_trend_box = [1850 1971]; yr_end_trend_box = [1970 2014]; %

%Other stuff
dT_P1 = (yr_end_trend_box(1) - yr_start_trend_box(1));
dT_P2 = (yr_end_trend_box(2) - yr_start_trend_box(2));


ioutput_abs=1; %whether to output the absolute changes (dSW, etc.) rather than the trends (i.e., trends * delta_time)
if ioutput_abs==0
    fdT_P1 = 1;
    fdT_P2 = 1;
    fscale=1e2;
else
    fdT_P1 = dT_P1;
    fdT_P2 = dT_P2;
    fscale=1;
end


%Loop over all models

for idamip_run=1:length(DAMIP_runs)
    
    expt_str = DAMIP_runs{idamip_run};
    %     switch expt_str
    %         case '' %UKESM
    %             %model_str = 'UKESM1_';
    %             model_str = [expt_str 'UKESM1_'];
    %         otherwise
    %             model_str = [expt_str '_']; %expt_str='DAMIP_hist-aer', etc.
    %     end
    
    switch expt_str
        case '' %UKESM
            %model_str = 'UKESM1_';
            model_str = [expt_str ''];
        case 'HADGEM'
            model_str = ['HADGEM3_GC31_LL_'];
        otherwise
            model_str = [expt_str '_']; %expt_str='DAMIP_hist-aer', etc.
    end
    
    model_str_title = model_str;
    
    
    %% Gather all the trend values and put into .tex file for table in the paper
    %tr_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only_UKESM1.mat';
    %tr_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only_UKESM1_.mat';
    tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str '.mat'];
    
    %For the aerosol proxy for AerChemMIP we calculate the trends using the
    %All-emissions UKESM (3 member sub-ensemble) minus the piAer (GHG-only
    %proxy) run.
    switch model_str
        case 'AerChemMIP_Aerosol_Proxy_'
            tr_file = remove_character(tr_file,model_str,'UKESM1-AerChemMIP_control_');
    end    
    sw_trends = load(tr_file);
    
    clear table_vals
    
    table_vals.fscale = fscale;    
    table_vals = ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_full_vals(sw_trends,fdT_P1,fdT_P2,fscale);
    
    switch model_str
        case 'AerChemMIP_Aerosol_Proxy_'
            tr_file = remove_character(tr_file,'UKESM1-AerChemMIP_control_','AerChemMIP_hist-piAer_');
            sw_trends = load(tr_file);            
            
            clear table_vals_ukesm
            
            table_vals_ukesm = ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_full_vals(sw_trends,fdT_P1,fdT_P2,fscale);
            
            table_vals.FullModelPA = table_vals.FullModelPA - table_vals_ukesm.FullModelPA;
            table_vals.FullModelPB = table_vals.FullModelPB - table_vals_ukesm.FullModelPB;
            table_vals.FullModelPAun = table_vals.FullModelPAun + table_vals_ukesm.FullModelPAun;
            table_vals.FullModelPBun = table_vals.FullModelPBun + table_vals_ukesm.FullModelPBun;
            
            %Estimated SW trend values from the offline calc
            table_vals.FullCalcPA = table_vals.FullCalcPA - table_vals_ukesm.FullCalcPA;
            table_vals.FullCalcPB = table_vals.FullCalcPB - table_vals_ukesm.FullCalcPB;
            table_vals.FullCalcPAun = table_vals.FullCalcPAun + table_vals_ukesm.FullCalcPAun;
            table_vals.FullCalcPBun = table_vals.FullCalcPBun + table_vals_ukesm.FullCalcPBun;
            
    end
    
    % -- Period 1 --
    
    %Need to set str='' for UKESM for these ones (fix in CALC routine at some
    %point)
    %     switch expt_str
    %         case '' %UKESM
    %             %model_str = 'UKESM1_';
    %             model_str = [expt_str ''];
    %         case 'HADGEM'
    %             model_str = ['HADGEM3_GC31_LL_'];
    %         otherwise
    %             model_str = [expt_str '_']; %expt_str='DAMIP_hist-aer', etc.
    %     end
    
    
    
    
    
    
    for itr=1:length(period_strs)
        
        period_str=period_strs{itr}; period_str_char=period_str_chars{itr};
        
        %Holding Nd constant
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only_N_d_constant__' period_str '.mat'];
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'N_d_constant__' period_str '.mat'];
        var_str = 'Ndconst';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding CF constant
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only_CF_constant__' period_str '.mat'];
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'CF_constant__' period_str '.mat'];
        var_str = 'CFconst';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding LWP constant
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_SW_TOA_up_calc_for_region_4__ocean_only_LWP_constant__' period_str '.mat'];
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'LWP_constant__' period_str '.mat'];
        var_str = 'LWPconst';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding everything constant except Nd
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'Nd_vary__' period_str '.mat'];
        var_str = 'Ndvary';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding everything constant except LWP
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'lwp_vary__' period_str '.mat'];
        var_str = 'LWPvary';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding everything constant except cf
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'cf_vary__' period_str '.mat'];
        var_str = 'CFvary';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding everything constant except clear-sky flux
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'clear_sky_vary__' period_str '.mat'];
        var_str = 'ClearSkyvary';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        %Holding everything constant except surface albedo (only affects cloudy
        %sky).
        tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'albedo_vary__' period_str '.mat'];
        var_str = 'Albedovary';
        %Run script to load the file and set the variables depending on the period
        ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_load_and_set
        
        
        
        % %Holding CF and LWP constant - in this case we use the trend rather than
        % %the difference to the full calc trend
        % tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'cf_LWP_constant__' period_str '.mat'];
        % sw_trends = load(tr_file);
        % %Full SW trend values from the model
        % table_vals.cfLWPconstPA = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.coeffs(2);
        % table_vals.cfLWPconstPAun = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.uncer_max;
        
        % %Holding CF and Nd constant - in this case we use the trend rather than
        % %the difference to the full calc trend
        % tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'Nd_cf_constant__' period_str '.mat'];
        % sw_trends = load(tr_file);
        % %Full SW trend values from the model
        % table_vals.NdcfconstPA = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.coeffs(2);
        % table_vals.NdcfconstPAun = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.uncer_max;
        
        % %Holding LWP and Nd constant - in this case we use the trend rather than
        % %the difference to the full calc trend
        % tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' model_str 'LWP_Nd_constant__' period_str '.mat'];
        % sw_trends = load(tr_file);
        % %Full SW trend values from the model
        % table_vals.LWPNdconstPA = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.coeffs(2);
        % table_vals.LWPNdconstPAun = fdT_P1*fscale*sw_trends.trend_dat_box{1,itr}.uncer_max;
        
    end
    
    %% Now calculate things based on the data read in.
    table_vals.CFcontPA = table_vals.FullCalcPA - table_vals.CFconstPA; %Contribution from CF - difference between full calc trend and that form holding CF constant
    table_vals.NdcontPA = table_vals.FullCalcPA - table_vals.NdconstPA;
    table_vals.LWPcontPA = table_vals.FullCalcPA - table_vals.LWPconstPA;
    table_vals.sumCalcPA = table_vals.CFcontPA + table_vals.NdcontPA + table_vals.LWPcontPA;
    table_vals.RESIDUALcontPA = table_vals.FullCalcPA  - table_vals.sumCalcPA;
    
    table_vals.altNdcontPA = table_vals.NdvaryPA;
    table_vals.altCFcontPA = table_vals.CFvaryPA;
    table_vals.altLWPcontPA = table_vals.LWPvaryPA;
    table_vals.altClearSkycontPA = table_vals.ClearSkyvaryPA;
    table_vals.altAlbedocontPA = table_vals.AlbedovaryPA;
    table_vals.altsumCalcPA = table_vals.altCFcontPA + table_vals.altNdcontPA + table_vals.altLWPcontPA + table_vals.altClearSkycontPA + table_vals.altAlbedocontPA;
    table_vals.altRESIDUALcontPA = table_vals.FullCalcPA  - table_vals.altsumCalcPA;
    
    table_vals.CFcontPAun = table_vals.FullCalcPAun + table_vals.CFconstPAun;
    table_vals.NdcontPAun = table_vals.FullCalcPAun + table_vals.NdconstPAun;
    table_vals.LWPcontPAun = table_vals.FullCalcPAun + table_vals.LWPconstPAun;
    table_vals.sumCalcPAun = table_vals.CFcontPAun + table_vals.NdcontPAun + table_vals.LWPcontPAun;
    table_vals.RESIDUALcontPAun = table_vals.FullCalcPAun  + table_vals.sumCalcPAun;
    
    table_vals.altNdcontPAun = table_vals.NdvaryPAun;
    table_vals.altCFcontPAun = table_vals.CFvaryPAun;
    table_vals.altLWPcontPAun = table_vals.LWPvaryPAun;
    table_vals.altClearSkycontPAun = table_vals.ClearSkyvaryPAun;
    table_vals.altAlbedocontPAun = table_vals.AlbedovaryPAun;
    table_vals.altsumCalcPAun = table_vals.altCFcontPAun + table_vals.altNdcontPAun + table_vals.altLWPcontPAun + table_vals.altClearSkycontPAun + table_vals.altAlbedocontPAun;
    table_vals.RESIDUALcontPAun = table_vals.FullCalcPAun  + table_vals.altsumCalcPAun;
    
    %table_vals.CFcontperPA = 100 * table_vals.CFcontPA ./ table_vals.sumCalcPA;
    %table_vals.CFcontperPA = 100 * table_vals.CFcontPA ./ table_vals.sumCalcPA;
    %table_vals.CFcontperPA = 100 * table_vals.CFcontPA ./ table_vals.sumCalcPA;
    table_vals.CFcontperPA = 100 * table_vals.CFcontPA ./ table_vals.FullCalcPA;
    table_vals.NdcontperPA = 100 * table_vals.NdcontPA ./ table_vals.FullCalcPA;
    table_vals.LWPcontperPA = 100 * table_vals.LWPcontPA ./ table_vals.FullCalcPA;
    table_vals.RESIDUALcontperPA = 100 - (table_vals.CFcontperPA  + table_vals.NdcontperPA + table_vals.LWPcontperPA);
    
    table_vals.altCFcontperPA = 100 * table_vals.altCFcontPA ./ table_vals.FullCalcPA;
    table_vals.altNdcontperPA = 100 * table_vals.altNdcontPA ./ table_vals.FullCalcPA;
    table_vals.altLWPcontperPA = 100 * table_vals.altLWPcontPA ./ table_vals.FullCalcPA;
    table_vals.altClearSkycontperPA = 100 * table_vals.altClearSkycontPA ./ table_vals.FullCalcPA;
    table_vals.altAlbedocontperPA = 100 * table_vals.altAlbedocontPA ./ table_vals.FullCalcPA;
    table_vals.altRESIDUALcontperPA = 100 - (table_vals.altCFcontperPA  + table_vals.altNdcontperPA + table_vals.altLWPcontperPA + table_vals.altClearSkycontperPA + table_vals.altAlbedocontperPA);
    
    %Period B
    
    table_vals.CFcontPB = table_vals.FullCalcPB - table_vals.CFconstPB;
    table_vals.NdcontPB = table_vals.FullCalcPB - table_vals.NdconstPB;
    table_vals.LWPcontPB = table_vals.FullCalcPB - table_vals.LWPconstPB;
    table_vals.sumCalcPB = table_vals.CFcontPB + table_vals.NdcontPB + table_vals.LWPcontPB;
    table_vals.RESIDUALcontPB = table_vals.FullCalcPB  - table_vals.sumCalcPB;
    
    table_vals.altNdcontPB = table_vals.NdvaryPB;
    table_vals.altCFcontPB = table_vals.CFvaryPB;
    table_vals.altLWPcontPB = table_vals.LWPvaryPB;
    table_vals.altClearSkycontPB = table_vals.ClearSkyvaryPB;
    table_vals.altAlbedocontPB = table_vals.AlbedovaryPB;
    table_vals.altsumCalcPB = table_vals.altCFcontPB + table_vals.altNdcontPB + table_vals.altLWPcontPB + table_vals.altClearSkycontPB + table_vals.altAlbedocontPB;
    table_vals.altRESIDUALcontPB = table_vals.FullCalcPB  - table_vals.altsumCalcPB;
    
    %Uncertainties
    table_vals.CFcontPBun = table_vals.FullCalcPBun + table_vals.CFconstPBun;
    table_vals.NdcontPBun = table_vals.FullCalcPBun + table_vals.NdconstPBun;
    table_vals.LWPcontPBun = table_vals.FullCalcPBun + table_vals.LWPconstPBun;
    table_vals.sumCalcPBun = table_vals.CFcontPBun + table_vals.NdcontPBun + table_vals.LWPcontPBun;
    table_vals.RESIDUALcontPBun = table_vals.FullCalcPBun  + table_vals.sumCalcPBun;
    
    table_vals.altNdcontPBun = table_vals.NdvaryPBun;
    table_vals.altCFcontPBun = table_vals.CFvaryPBun;
    table_vals.altLWPcontPBun = table_vals.LWPvaryPBun;
    table_vals.altClearSkycontPBun = table_vals.ClearSkyvaryPBun;
    table_vals.altAlbedocontPBun = table_vals.AlbedovaryPBun;
    table_vals.altsumCalcPBun = table_vals.altCFcontPBun + table_vals.altNdcontPBun + table_vals.altLWPcontPBun + table_vals.altClearSkycontPBun + table_vals.altAlbedocontPBun;
    table_vals.altRESIDUALcontPBun = table_vals.FullCalcPBun  + table_vals.altsumCalcPBun;
    
    %Contributions expressed as percentages
    table_vals.CFcontperPB = 100 * table_vals.CFcontPB ./ table_vals.FullCalcPB;
    table_vals.NdcontperPB = 100 * table_vals.NdcontPB ./ table_vals.FullCalcPB;
    table_vals.LWPcontperPB = 100 * table_vals.LWPcontPB ./ table_vals.FullCalcPB;
    table_vals.RESIDUALcontperPB = 100 - (table_vals.CFcontperPB  + table_vals.NdcontperPB + table_vals.LWPcontperPB);
    
    table_vals.altCFcontperPB = 100 * table_vals.altCFcontPB ./ table_vals.FullCalcPB;
    table_vals.altNdcontperPB = 100 * table_vals.altNdcontPB ./ table_vals.FullCalcPB;
    table_vals.altLWPcontperPB = 100 * table_vals.altLWPcontPB ./ table_vals.FullCalcPB;
    table_vals.altClearSkycontperPB = 100 * table_vals.altClearSkycontPB ./ table_vals.FullCalcPB;
    table_vals.altAlbedocontperPB = 100 * table_vals.altAlbedocontPB ./ table_vals.FullCalcPB;
    table_vals.altRESIDUALcontperPB = 100 - (table_vals.altCFcontperPB  + table_vals.altNdcontperPB + table_vals.altLWPcontperPB...
        + table_vals.altLWPcontperPB + table_vals.altAlbedocontperPB);
    
    
    %% Make the table
    
    iappend=0;
    save_sw_table_vals = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/' filename_prestr '_' land_ocean '_' expt_str 'TABLE'];
    %var_str = '';
    var_str = expt_str; %pre-fix values with the DAMIP_hist-aer etc.
    latex_newcommand_from_structure(table_vals,var_str,save_sw_table_vals,iappend);
    
    switch stacked_bar_or_error_bars
        case 'error bars'
            for itr=1:length(period_strs)
                period_str=period_strs{itr}; period_str_char=period_str_chars{itr};
                SW_paper_dSW_from_cloud_vars_plot_horiz
            end
            
        case 'stacked bar'
            %% Bar plot
            
            %xlabs_str={'Period 1','Period 2'};
            
            %full_calc = 4.64;
            
            %Period 1
            
            clear bar_dat xlabs_str
            
            ibar=1;
            xlabs_str{ibar}='';
            % bar_dat(ibar,1) = table_vals.FullCalcPA - table_vals.NdconstPA;
            % bar_dat(ibar,2) = table_vals.FullCalcPA - table_vals.CFconstPA;
            % bar_dat(ibar,3) = table_vals.FullCalcPA - table_vals.LWPconstPA;
            bar_dat(ibar,1) = table_vals.altNdcontPA;
            bar_dat(ibar,2) = table_vals.altCFcontPA;
            bar_dat(ibar,3) = table_vals.altLWPcontPA;
            bar_dat(ibar,4) = table_vals.altClearSkycontPA;
            bar_dat(ibar,5) = 0; %Net sum
            bar_dat(ibar,6) = 0; %to represent the calculated trend
            bar_dat(ibar,7) = 0; %to represent the actual model trend
            
            %Next bar just for net of first bar
            ibar=ibar+1;
            xlabs_str{ibar}='Period 1';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = sum(bar_dat(ibar-1,:));
            bar_dat(ibar,6) = 0;
            bar_dat(ibar,7) = 0;
            
            %Next bar just for calcaulted trend
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 1';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = 0;
            bar_dat(ibar,6) = table_vals.FullCalcPA;
            bar_dat(ibar,7) = 0;
            
            %Next bar just for model (true) trend
            ibar=ibar+1;
            xlabs_str{ibar}='';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = 0;
            bar_dat(ibar,6) = 0;
            bar_dat(ibar,7) = table_vals.FullModelPA;
            
            
            % -- Period 2
            
            ibar=ibar+1;
            % bar_dat(ibar,1) = (table_vals.FullCalcPB - table_vals.NdconstPB);
            % bar_dat(ibar,2) = (table_vals.FullCalcPB - table_vals.CFconstPB);
            % bar_dat(ibar,3) = (table_vals.FullCalcPB - table_vals.LWPconstPB);
            bar_dat(ibar,1) = table_vals.altNdcontPB;
            bar_dat(ibar,2) = table_vals.altCFcontPB;
            bar_dat(ibar,3) = table_vals.altLWPcontPB;
            bar_dat(ibar,4) = table_vals.altClearSkycontPB;
            bar_dat(ibar,5) = 0;
            bar_dat(ibar,6) = 0; %to represent the calculated trend
            bar_dat(ibar,7) = 0; %to represent the actual model trend
            
            %Next bar just for net of first bar
            ibar=ibar+1;
            xlabs_str{ibar}='Period 2';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = sum(bar_dat(ibar-1,:));
            bar_dat(ibar,6) = 0;
            bar_dat(ibar,7) = 0;
            
            %Next bar just for calcalated trend
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 2';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = 0;
            bar_dat(ibar,6) = table_vals.FullCalcPB;
            bar_dat(ibar,7) = 0;
            
            
            %Next bar just for model (true) trend
            ibar=ibar+1;
            xlabs_str{ibar}='';
            bar_dat(ibar,1) = 0;
            bar_dat(ibar,2) = 0;
            bar_dat(ibar,3) = 0;
            bar_dat(ibar,4) = 0;
            bar_dat(ibar,5) = 0;
            bar_dat(ibar,6) = 0;
            bar_dat(ibar,7) = table_vals.FullModelPB;
            
            
            %% Uncertainties for bars
            clear bar_dat_UN
            
            ibar=1;
            %xlabs_str{ibar}='';
            % bar_dat_UN(ibar,1) = table_vals.NdcontPAun;
            % bar_dat_UN(ibar,2) = table_vals.CFcontPAun;
            % bar_dat_UN(ibar,3) = table_vals.LWPcontPAun;
            bar_dat_UN(ibar,1) = table_vals.altNdcontPAun;
            bar_dat_UN(ibar,2) = table_vals.altCFcontPAun;
            bar_dat_UN(ibar,3) = table_vals.altLWPcontPAun;
            bar_dat_UN(ibar,4) = table_vals.altClearSkycontPAun;
            bar_dat_UN(ibar,5) = 0; %Net sum
            bar_dat_UN(ibar,6) = 0; %to represent the calculated trend
            bar_dat_UN(ibar,7) = 0; %to represent the actual model trend
            
            %Next bar just for net of first bar
            ibar=ibar+1;
            xlabs_str{ibar}='Period 1';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = table_vals.altsumCalcPAun;
            bar_dat_UN(ibar,6) = 0;
            bar_dat_UN(ibar,7) = 0;
            
            %Next bar just for calcaulted trend
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 1';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = 0;
            bar_dat_UN(ibar,6) = table_vals.FullCalcPAun;
            bar_dat_UN(ibar,7) = 0;
            
            %Next bar just for model (true) trend
            ibar=ibar+1;
            xlabs_str{ibar}='';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = 0;
            bar_dat_UN(ibar,6) = 0;
            bar_dat_UN(ibar,7) = table_vals.FullModelPAun;
            
            
            % -- Period 2
            
            ibar=ibar+1;
            % bar_dat_UN(ibar,1) = table_vals.NdcontPBun;
            % bar_dat_UN(ibar,2) = table_vals.CFcontPBun;
            % bar_dat_UN(ibar,3) = table_vals.LWPcontPBun;
            bar_dat_UN(ibar,1) = table_vals.altNdcontPBun;
            bar_dat_UN(ibar,2) = table_vals.altCFcontPBun;
            bar_dat_UN(ibar,3) = table_vals.altLWPcontPBun;
            bar_dat_UN(ibar,4) = table_vals.altClearSkycontPBun;
            bar_dat_UN(ibar,5) = 0; %Net sum
            bar_dat_UN(ibar,6) = 0; %to represent the calculated trend
            bar_dat_UN(ibar,7) = 0; %to represent the actual model trend
            
            %Next bar just for net of first bar
            ibar=ibar+1;
            xlabs_str{ibar}='Period 2';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = table_vals.altsumCalcPBun;
            bar_dat_UN(ibar,6) = 0;
            bar_dat_UN(ibar,7) = 0;
            
            %Next bar just for calcaulted trend
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 1';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = 0;
            bar_dat_UN(ibar,6) = table_vals.FullCalcPBun;
            bar_dat_UN(ibar,7) = 0;
            
            %Next bar just for model (true) trend
            ibar=ibar+1;
            xlabs_str{ibar}='';
            bar_dat_UN(ibar,1) = 0;
            bar_dat_UN(ibar,2) = 0;
            bar_dat_UN(ibar,3) = 0;
            bar_dat_UN(ibar,4) = 0;
            bar_dat_UN(ibar,5) = 0;
            bar_dat_UN(ibar,6) = 0;
            bar_dat_UN(ibar,7) = table_vals.FullModelPBun;
            
            
            %%
            
            offset_errbars=[1 1 1 1 0 0 0];
            bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};
            
            
            %Matlab plots the stacked bars over the last one so that it covers them up
            %if the sign of the bars change. Here willl plot +ve and -ve separately and
            %then have a net bar.
            bar_dat_UN_neg = bar_dat_UN;
            bar_dat_UN(bar_dat<0)=0;
            bar_dat_UN_neg(bar_dat>0)=0;
            
            bar_dat_neg = bar_dat;
            bar_dat(bar_dat<0)=0;
            bar_dat_neg(bar_dat_neg>0)=0;
            
            
            figure
            set(gcf,'color','w');
            hbar = bar(bar_dat,'stacked'); %will give N bars each with M components for bar_dat[N,M]
            %manually set the colours to have control - otherwise not sure how to
            %get them as FaceColor is set to 'Flat'!
            for i=1:length(hbar)
                set(hbar(i),'FaceColor',bar_cols{i});
            end
            hold on
            hbar = bar(bar_dat_neg,'stacked');
            for i=1:length(hbar)
                set(hbar(i),'FaceColor',bar_cols{i});
            end
            %alpha(0.5); %Make transparent?
            increase_font_size_map_figures
            %title(['SW_{in}=' num2str(SW_in) ', A_{clear}=' num2str(A_clear)]);
            if ioutput_abs==0
                ylabel('SW_{up TOA} Trend (W m^{-2} yr^{-1}) x 10^{-2}');
            else
                ylabel('\DeltaSW_{up TOA} (W m^{-2})');
            end
            set(gca,'xticklabel',xlabs_str);
            leg_str = {'N_{d}','f_{c}','LWP_{ic}','Clear-sky','Net','Full Calculated','Actual'};
            % pad_list=[1 2];
            % for i=1:length(pad_list)
            %     iL=pad_list(i);
            %     leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
            % end
            
            pad_list=[1 2 3 4 5 6];
            pad_list=[1:size(bar_dat,2)];
            for i=1:length(pad_list)
                iL=pad_list(i);
                leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
            end
            %L=legend(leg_str,'location','SouthWest','fontsize',14);
            L=legend(leg_str,'location','BestOutside','fontsize',14);
            %set(L,'string',leg_str)
            
            title(remove_character(model_str_title,'_',' '));
            
            set(gca,'xlim',[0 ibar+1]);
            set(gca,'ylim',ylims_bars);
            
            
            %% Plot the error bars
            
            bar_dat_plot = bar_dat;
            bar_dat_UN_plot = bar_dat_UN;
            ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT
            
            bar_dat_plot = bar_dat_neg;
            bar_dat_UN_plot = bar_dat_UN_neg;
            ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT
            
            
    end
    
%     savename=[savedir_date titlenam ' ' ens_str ' ' land_ocean expt_str '_BAR_PLOT'];
%     %savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
%     clear opts
%     %        opts.iplot_png=1;
%     opts.iplot_eps=1;
%     opts.iplot_jpg=0;
%     savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
    
%For the errobar plots, they are saved individually for the different periods
    
    
    
end
