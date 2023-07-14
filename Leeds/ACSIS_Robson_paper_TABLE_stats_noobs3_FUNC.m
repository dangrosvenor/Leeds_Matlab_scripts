%Function/script to gather the able variables for a single period.
% Is run for each period from ACSIS_Robson_paper_TABLE_stats_noobs3.m
% Runs ACSIS_Robson_paper_TABLE_stats_noobs_FUNC.m

iregion={'4'};

% -- get the indices for the requried time period
%var_str_mat='SW_TOA_up'; fscale=1e2; 
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
%trends = load(tr_file);
var_str_mat='N_d'; 
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);
i0=find(trends.years_ukesm_1d==trends.yr_start_trend_box2(itr));

if idelta_vals==1
    deltaT = trends.yr_end_trend_box2(itr) - trends.yr_start_trend_box2(itr);
else
    deltaT = 1;
end



%% -- Nd -- *** Nd needs to be run first (because it uses the Nd change to work out the forcing ***

%UKESM
run_str='';
var_str_mat='N_d'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1e1;
end
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);
var_str_tab='Nd'; 
hist_str='';
ido_ens=1;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

if iobs_period==1
    iobs=1;
    run_str='OBS';
    %run script
    ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
    iobs=0;    
end

%HadGEM
tr_file = ['/home/disk/os1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
run_str='HAD'; 
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='Nd';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%% Surface temperature 
% Will also do surface temperature early too since might want to calculate a feedback effect
% UKESM

for iregion={'4' '0'};
    
run_str='';
var_str_mat='Surface_Temperature'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1;
end
var_str_tab='TS';


hist_str='';
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


    
if iobs_period==1
    iobs=1;
    run_str='OBS';
    %run script
    ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
    iobs=0;    
end



%HADGEM
%var_str_mat='SW_TOA_up';
run_str='HAD';
hist_str='';
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC





% DAMIP 

%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='SW';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

end
iregion={'4'};


%% -- SW TOA --
%UKESM
run_str='';
%var_str_mat='SW_TOA_up'; 
var_str_mat='F_SW_TOA_upwelling'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1e2;
end
var_str_tab='SW';

hist_str='';
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

if iobs_period==1
    iobs=1;
    run_str='OBS';
    %run script
    ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
    iobs=0;    
end


%HADGEM
%var_str_mat='SW_TOA_up';
var_str_mat='F_SW_TOA_upwelling'; 
run_str='HAD';
hist_str='';
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
var_str_mat='F_SW_TOA_upwelling'; 
run_str='UKESMAerChemMIP';
hist_str='';
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% -- DAMIP trends and uncertainties --

%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='SW';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

ido_ens=0;

%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%% -- fc --
%UKESM
run_str='';
var_str_mat='totCF';
if idelta_vals==1
    fscale=1e3;
else
    fscale=1e5;
end
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);
var_str_tab='CF'; 
hist_str='';
ido_ens=1;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

if iobs_period==1
    iobs=1;
    run_str='OBS';
    %run script
    ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
    iobs=0;    
end


%HadGEM
var_str_mat='totCF'; 
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
run_str = 'HAD';
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1_DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='CF';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC




%% -- LWPic --

%UKESM
run_str='';
var_str_mat='LWPic'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1e2;
end
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);

var_str_tab='LWPic';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%No obs for this as yet.
% if iobs_period==1
%     iobs=1;
%     run_str='OBS';
%     %run script
%     ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
%     iobs=0;    
% end


%HadGEM
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
trends = load(tr_file);
run_str='HAD'; 
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1_DAMIP_no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='LWPic';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC




% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC



%% LWP


%UKESM
run_str='';
var_str_mat='LWP'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1e2;
end
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);

var_str_tab='LWP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

if iobs_period==1
    iobs=1;
    run_str='OBS';
    %run script
    ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
    iobs=0;    
end


%HadGEM
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
run_str='HAD'; 
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1_DAMIP_no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='LWPic';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%% AOD


%UKESM
run_str='';
var_str_mat='Aerosol+Dust_Optical_Depth,550nm'; 
if idelta_vals==1
    fscale=1;
else
    fscale=1e2;
end
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1_.mat'];
trends = load(tr_file);

var_str_tab='od550tot';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

% if iobs_period==1
%     iobs=1;
%     run_str='OBS';
%     %run script
%     ACSIS_Robson_paper_TABLE_stats_noobs_FUNC
%     iobs=0;    
% end


%HadGEM
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL_.mat'];
trends = load(tr_file);
run_str='HAD'; 
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1_DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='LWPic';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-aer2 (HADEGM minus hist-GHG)
run_str=['HistAer2'];
hist_str='_hist_aer2';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% -- AerChemMIP

% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC




%% -- rsutcs (clear-sky TOA flux) --
%Still need to process for the third trend for obs period - don't have obs
%for this yet.
if itr<=2
    
%UKESM
run_str='';
var_str_mat='Clear-sky_SW_TOA_up';
if idelta_vals==1
    fscale=1e0;
else
    fscale=1e5;
end
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_UKESM1__no_obs.mat'];
trends = load(tr_file);
var_str_tab='RSUTCS'; 
hist_str='';
ido_ens=1;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%HadGEM
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4_' ens_str '_' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
trends = load(tr_file);
run_str = 'HAD';
hist_str='';
ido_ens=0;
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%DAMIP trends and uncertainties
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1_DAMIP_no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '__DAMIP_no_obs.mat'];
%Annual_mean_totCF_for_region_4__' land_ocean_str '_HADGEM-GC31-LL__no_obs
trends = load(tr_file);

ido_ens=0;
%var_str_damip_tab='CF';

%Hist-aer
run_str=['HistAer'];
hist_str='_hist_aer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

% %Hist-aer2 (HADEGM minus hist-GHG)
% run_str=['HistAer2'];
% hist_str='_hist_aer2';
% %run script
% ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-GHG
run_str=['HistGhg'];
hist_str='_hist_GHG';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-nat
run_str=['HistNat'];
hist_str='_hist_nat';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%Hist-linear (sum)
run_str=['HistLin'];
hist_str='_hist_linear';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


% AerChemMIP runs
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_4__' land_ocean_str '_UKESM1-AerChemMIP-control__no_obs.mat'];
%tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_mat '_for_region_' iregion{1} '__' land_ocean_str '_UKESM1-AerChemMIP-control_.mat'];
trends = load(tr_file);

% UKESM AerChemMIP
%var_str_mat='SW_TOA_up';
run_str='UKESMAerChemMIP';
hist_str='';
ido_ens=1;

%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

ido_ens=0;
%Aerosol proxy (UKESM minus piAer)
run_str=['AerChemAerosol'];
hist_str='_AerChemMIP_aero';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC

%GHG proxy (piAer)
run_str=['AerChemGHG'];
hist_str='_AerChemMIP_piaer';
%run script
ACSIS_Robson_paper_TABLE_stats_noobs_FUNC


end



