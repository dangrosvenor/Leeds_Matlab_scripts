%Run from :-    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script
%This script is run once for each time period from ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot.m

run_str = 'UKESMAerChemMIP';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'AerChemGHG';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'AerChemAerosol';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'HAD';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'HistAer';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'HistGhg';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

run_str = 'HistNat';
ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func01

%Trend from UKESM minus the trend from piAer as proxy for aerosol trend

AerChemAerosol2 = UKESMAerChemMIP - AerChemGHG;
AerChemAerosol2_un = UKESMAerChemMIP_un + AerChemGHG_un;
if (abs(AerChemAerosol2) - AerChemAerosol2_un > 0) AerChemAerosol2_sig=1; else AerChemAerosol2_sig=0; end

netAerChemMIP = AerChemGHG + AerChemAerosol;
netAerChemMIP_un = AerChemGHG_un + AerChemAerosol_un;
if (abs(netAerChemMIP) - netAerChemMIP_un > 0) netAerChemMIP_sig=1; else netAerChemMIP_sig=0; end

aer = eval(['table_vals.' var_str_tab 'HistAertrend' period_str]);
aer_un = eval(['table_vals.' var_str_tab 'HistAertrend' period_str 'un']);
if (abs(aer) - aer_un > 0) aer_sig=1; else aer_sig=0; end

aer2 = eval(['table_vals.' var_str_tab 'HistAer2trend' period_str]);
aer2_un = eval(['table_vals.' var_str_tab 'HistAer2trend' period_str 'un']);
if (abs(aer2) - aer2_un > 0) aer2_sig=1; else aer2_sig=0; end


%set some defaults
local=0;
local_indirect = 0;
local_direct = 0;

ukesm_local = 0;
ukesm_local_indirect = 0;
ukesm_local_direct = 0;

local_sig=0; %not used any more I think

if ilocal(ivar)==1
    %local here refers to the changes in the nudged run - i.e., the ACI+ARI
    %forcing (including adjustments) for constant SST.
    local = eval(['table_vals.' var_str_tab 'HistAerTrendAeroLocal' period_str]); %DAMIP HistAer
    local2 = eval(['table_vals.' var_str_tab 'HistAer2TrendAeroLocal' period_str]); %DAMIP HistAer proxy?
    ukesm_local = eval(['table_vals.' var_str_tab  'TrendAeroLocal' period_str]); %UKESM (ones with no HistAer etc. labels)
    %- do for HADGEM too? At the moment I'm just using the UKESM values for
    %the AerChemMIP aerosol-only proxy on the assumption that the change in
    %Nd for that is the same as in the UKESM (since proxy value = UKESM
    %sub-ensemble minus no-GHG run and the no-GHG chagne in Nd is zero.
    %Plus the sub-ensemble change in Nd is likely similar to the full
    %ensemble).
                                                                                    %
    switch var_str_tab
        case 'SW'
            local_indirect = eval(['table_vals.' var_str_tab  'HistAerIndirectTrendAeroLocal' period_str]);
            local_indirect2 = eval(['table_vals.' var_str_tab  'HistAer2IndirectTrendAeroLocal' period_str]);
            ukesm_local_indirect = eval(['table_vals.' var_str_tab  'IndirectTrendAeroLocal' period_str]);
            local_direct2 = local2 - local_indirect2; %assume direct makes up the total after indirect
    end
    local_direct = local - local_indirect; %assume direct makes up the total after indirect    
    ukesm_local_direct = ukesm_local - ukesm_local_indirect; %assume direct makes up the total after indirect
    %if (abs(local) - aer_un > 0) local_sig=1; else local_sig=0; end
    local_sig=1;
    
    non_local = aer - local;
    non_local2 = aer2 - local2;
    ukesm_non_local = AerChemAerosol - ukesm_local; %should really calculate a local value for the AerChemMIP
    %sub-ensemble - need to run the pre-processing script for Nd first.
    %Assuming full UKESM Nd change is similar to the sub-ensemble here.
    %Prob a good approximation since Nd doesn't vary much between
    %ensemble members.
    
    switch var_str_tab
        case {'SW','CF','LWPic'}
            
            %Feedbacks to SWTOA based on GHG only runs (or proxy for AerChemMIP)
            %DAMIP
            dSWdTS_ghg = eval(['table_vals.' var_str_tab 'feedbackHistGhg' period_str]);
            %dSWdTS_ghg = eval(['table_vals.SWfeedbackHistGhg' period_str]);
            
            dTS = eval(['table_vals.TSHistAertrend' period_str]);
            feedback_estimate_aer = dSWdTS_ghg .* dTS;
            
            %Using the proxy aeorosl (diff between full model and DAMIP-hist-ghg).
            dTS2 = eval(['table_vals.TSHistAer2trend' period_str]);
            feedback_estimate_aer2 = dSWdTS_ghg .* dTS2;
            
            dTS_all = eval(['table_vals.TSHADtrend' period_str]);
            feedback_estimate_all = dSWdTS_ghg .* dTS_all;
            
            %AerChemMIP (UKESM) - using the proxy GHG (piAer) run here.
            %dSWdTS_ghg_AerChemMIP = eval(['table_vals.SWfeedbackAerChemGHG' period_str]);
            dSWdTS_ghg_AerChemMIP = eval(['table_vals.' var_str_tab 'feedbackAerChemGHG' period_str]);
            
            
            dTS_AerChem = eval(['table_vals.TSAerChemAerosoltrend' period_str]);
            feedback_estimate_AerChemAerosol = dSWdTS_ghg_AerChemMIP .* dTS_AerChem;
            
            dTS_AerChem_all = eval(['table_vals.TSUKESMAerChemMIPtrend' period_str]);
            feedback_estimate_AerChem_all = dSWdTS_ghg_AerChemMIP .* dTS_AerChem_all;
            
        otherwise
            feedback_estimate_aer = NaN;
            feedback_estimate_aer2 = NaN;
            feedback_estimate_all = NaN;
            feedback_estimate_AerChemAerosol = NaN;
            feedback_estimate_AerChem_all = NaN;                        
            
    end
    
    dT_obs = table_vals.TSOBStrendPC * 10.^table_vals.TSOBStrendPCexp;
    dSW_obs = table_vals.SWOBStrendPC * 10.^table_vals.SWOBStrendPCexp;
    dSWdT_obs = dSW_obs/dT_obs;
    
    switch period_str
        case 'PB';
            dSW_obsdT_PC = dSWdTS_ghg * dT_obs;
    end
    
    %N.B. - for PC the vlaues are trends rather than deltas
    model = 'ukesm';
    var = 'TS'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
    var = 'SW'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02  
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func03 %calculates e.g. dSWdTS_ukesm_PB
        %and dSWdTSglobal_ukesm_PB    
    
    model = 'HAD';
    var = 'TS'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
    var = 'SW'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02   
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func03 %calculates e.g. dSWdTS_ukesm_PB
        %and dSWdTSglobal_ukesm_PB              
    
    model = 'HistGhg';
    var = 'TS'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
    var = 'SW'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02   
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func03 %calculates e.g. dSWdTS_ukesm_PB
        %and dSWdTSglobal_ukesm_PB            
           
    model = 'AerChemGHG';
    var = 'TS'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
    var = 'SW'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02  
    ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func03 %calculates e.g. dSWdTS_ukesm_PB
        %and dSWdTSglobal_ukesm_PB
    
        switch period_str
            case 'PC'
                model = 'OBS';
                var = 'TS'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
                var = 'SW'; ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func02
                ACSIS_Robson_paper_TABLE_stats_noobs2_BAR_script_func03 %calculates e.g. dSWdTS_ukesm_PB
                %and dSWdTSglobal_ukesm_PB
        end
    
else
    local = NaN; local2 = NaN;
    local_indirect = NaN; local_indirect2 = NaN; ukesm_local_indirect = NaN;
    local_direct = NaN; local_direct2 = NaN; ukesm_local_direct = NaN;
    non_local = NaN; non_local2 = NaN;
end





ghg = eval(['table_vals.' var_str_tab 'HistGhgtrend' period_str]);
ghg_un = eval(['table_vals.' var_str_tab 'HistGhgtrend' period_str 'un']);
if (abs(ghg) - ghg_un > 0) ghg_sig=1; else ghg_sig=0; end

nat = eval(['table_vals.' var_str_tab 'HistNattrend' period_str]);
nat_un = eval(['table_vals.' var_str_tab 'HistNattrend' period_str 'un']);
if (abs(nat) - nat_un > 0) nat_sig=1; else nat_sig=0; end

net = aer+ghg+nat;
net_un = aer_un + ghg_un + nat_un;
if (abs(net) - net_un > 0) net_sig=1; else net_sig=0; end

net2 = aer2+ghg+nat;
net2_un = aer2_un + ghg_un + nat_un;
if (abs(net2) - net2_un > 0) net2_sig=1; else net2_sig=0; end



had = eval(['table_vals.' var_str_tab  'HADtrend' period_str]);
had_un = eval(['table_vals.' var_str_tab  'HADtrend' period_str 'un']);
if (abs(had) - had_un > 0) had_sig=1; else had_sig=0; end




%Trend from HADGEM minus the trend from DAMIP-Hist-GHG as proxy for aerosol trend
aer3 = had - ghg;
aer3_un = had_un + ghg_un;
if (abs(aer3) - aer3_un > 0) aer3_sig=1; else aer3_sig=0; end

ukesm = eval(['table_vals.' var_str_tab  'trend' period_str]);
ukesm_un = eval(['table_vals.' var_str_tab  'trend' period_str 'un']);



exponent = eval(['table_vals.' var_str_tab 'HistAertrend' period_str 'exp']);

%% Do the plotting

switch stacked_bar_or_error_bars
    case 'error_bar'
        %SW_paper_dSW_plot2
        %SW_paper_dSW_plot_horiz
        SW_paper_dSW_plot_horiz_just_vals
        
        
        
        
        
        
    case 'stacked_bar'
        
        ibar=ibar+1;
        
        ibar_start = ibar; %save this
        
        xlabs_str{ibar}='';
        
        iblock=1;
        bar_dat(ibar,iblock) = local_indirect; iblock=iblock+1;
        bar_dat(ibar,iblock) = local_direct; iblock=iblock+1;
        bar_dat(ibar,iblock) = non_local; iblock=iblock+1;
        bar_dat(ibar,iblock) = ghg; iblock=iblock+1;%GHGs
        bar_dat(ibar,iblock) = nat; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1; %to represent the actual model trend
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for net of first bar
        ibar=ibar+1;
        xlabs_str{ibar} = ''; %period_lab;
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = net; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for calcaulted trend
        ibar=ibar+1;
        xlabs_str{ibar} = period_lab;
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = had; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        if ilocal(ivar)==1
            %UKESM bar
            ibar=ibar+1;
            xlabs_str{ibar}=''; %'Period 1';
            iblock=1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat(ibar,iblock) = ukesm_local_indirect; iblock=iblock+1;
            bar_dat(ibar,iblock) = ukesm_local_direct; iblock=iblock+1;
            bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        end
        
        %UKESM bar
        ibar=ibar+1;
        xlabs_str{ibar}=''; %'Period 1';
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = ukesm; iblock=iblock+1;
        
        
        
        
        
        ibar=ibar+1;
        xlabs_str{ibar}=''; %'Period 1';
        iblock=1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat(ibar,iblock) = 0; iblock=iblock+1;
        
        %% Now repeat for uncertrainties
        ibar=ibar_start;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = aer_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = ghg_un;  iblock=iblock+1;%GHGs
        bar_dat_UN(ibar,iblock) = nat_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;%to represent the actual model trend
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for net of first bar
        ibar=ibar+1;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = net_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        %Next bar just for calcaulted trend
        ibar=ibar+1;
        iblock=1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = had_un; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        if ilocal(ivar)==1
            iblock=1;
            ibar=ibar+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
            
        end
        
        iblock=1;
        ibar=ibar+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = ukesm_un; iblock=iblock+1;
        
        
        iblock=1;
        ibar=ibar+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0;  iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        bar_dat_UN(ibar,iblock) = 0; iblock=iblock+1;
        
        
        
        % Define the colours oursleves as otherwise can get a weird coloru
        % order.
        %bar_cols = {'b','r','g','c','k','m',};
        bar_cols = {[0.5 0.5 0.75],[0.75 0.75 1.0],[0 0 1],[0 1 0],[0 1 1],[0 0 0],[1 0 1],[1 0 0],[0.75 0.5 0.5],[1.0 0.75 0.75]};
        offset_errbars=[0 0 1 1 1 0 0 0 0 0];
        
        
        if ilocal(ivar)==1
            if strcmp(var_str_tab,'SW')
                leg_str = {'Aerosol ACI forcing','Aerosol ARI forcing','Aerosol feedback','GHG','Natural','Net','HADGEM total','UKESM aerosol ACI forcing','UKESM aerosol ARI forcing','UKESM total'};
                %leg_str = {'Local aerosol','N/A','GHG','Natural','Net','HADGEM'};
                cut_inds=[];
            else
                leg_str = {'Aerosol forcing','Aerosol feedback','GHG','Natural','Net','HADGEM total','UKESM aerosol forcing','UKESM total'};
                cut_inds = [1 8]; %cut the ACI bar (rebadge as local aerossol)
            end
        else
            %leg_str = {'N/A','All aerosol','GHG','Natural','Net','HADGEM','N/A','UKESM'};
            %leg_str = {'All aerosol','GHG','Natural','Net','HADGEM'};
            %leg_str = {'All aerosol','N/A','GHG','Natural','Net','HADGEM'};
            
            %Cut out the local aerosol stacked bars and colours
            cut_inds = [1 2 8 9];
            leg_str = {'All aerosol','GHG','Natural','Net','HADGEM','UKESM'};
            
        end
        
        bar_dat(:,cut_inds)=[];
        bar_dat_UN(:,cut_inds)=[];
        bar_cols(cut_inds)=[];
        offset_errbars(cut_inds)=[];
        
        
        pad_list = [1:size(bar_dat,2)];
        
        bar_dat_UN_neg = bar_dat_UN;
        bar_dat_UN(bar_dat<0)=0;
        bar_dat_UN_neg(bar_dat>0)=0;
        
        %Matlab plots the stacked bars over the last one so that it covers them up
        %if the sign of the bars change. Here willl plot +ve and -ve separately and
        %then have a net bar.
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
        %if ioutput_abs==0
        %    ylabel('SW_{up TOA} Trend (W m^{-2} yr^{-1}) x 10^{-2}');
        %else
        %    ylabel('\Delta fc (W m^{-2})');
        %end
        if exponent~=0
            exp_str = [' x10^{' num2str(exponent) '}'];
        else
            exp_str='';
        end
        ylabel([ylabs{ivar} exp_str]);
        
        set(gca,'xticklabel',xlabs_str);
        
        % pad_list=[1 2];
        % for i=1:length(pad_list)
        %     iL=pad_list(i);
        %     leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
        % end
        
        for i=1:length(pad_list)
            iL=pad_list(i);
            leg_str{iL}(2,:)=leg_str{iL}(1,:); leg_str{iL}(1,:)=' ';
        end
        %L=legend(leg_str,'location','SouthWest','fontsize',14);
        L=legend(leg_str,'location','BestOutside','fontsize',14);
        %set(L,'string',leg_str)
        
        
        set(gca,'xlim',[0 ibar+1]);
        
        
        %cumsum_bar = cumsum(bar_dat,2);
        %cumsum_bar_neg = cumsum(bar_dat_neg,2);
        
        % Plot the error bars
        
        bar_dat_plot = bar_dat;
        bar_dat_UN_plot = bar_dat_UN;
        ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT
        
        bar_dat_plot = bar_dat_neg;
        bar_dat_UN_plot = bar_dat_UN_neg;
        ACSIS_Robson_paper_TABLE_stats_noobs2_BarPlot_PLOT                                
        
        
end
