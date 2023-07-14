
%set obs_str='none' before running the RUN_multi script for each variable
%to save the data needed for this script.

%% Gather all the trend values and put into .tex file for table in the paper
clear table_vals

% -- Period 1
itr=1;

% -- SW TOA --
var_str_tab='SW_TOA_up'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.SWtrendPA = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.SWtrendPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.SWtrendPAexp_precision = '%d';
table_vals.SWtrendPAun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.SWtrendPAmin = minval;
table_vals.SWtrendPAminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.SWtrendPAmax = maxval;
table_vals.SWtrendPAmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;

i0=find(trends.years_ukesm_1d==yr_start_trend_box2(1));
%i1=find(trends.years_ukesm_1d==yr_end_trend_box2(1));
table_vals.SWavFirstFivePA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.SWtrendPctPA = 100*table_vals.SWtrendPA/table_vals.SWavFirstFivePA; %Pct trend - i.e. *100 and /fscale to get percentage
table_vals.SWtrendPctPA_precision='%.9f';





% -- fc --
var_str_tab='totCF'; fscale=1e5;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.CFtrendPA = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendPAexp_precision = '%d';
table_vals.CFtrendPAun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.CFtrendPAmin = minval;
table_vals.CFtrendPAminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.CFtrendPAmax = maxval;
table_vals.CFtrendPAmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;
table_vals.CFavFirstFivePA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctPA = 100*table_vals.CFtrendPA/table_vals.CFavFirstFivePA;
table_vals.CFtrendPctPA_precision='%.9f';

%DAMIP trends and uncertainties
var_str_tab='totCF'; fscale=1e4;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%Hist-aer
table_vals.CFtrendHistAerPA = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistAerPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistAerPAexp_precision = '%d';
table_vals.CFtrendHistAerPAun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistAerPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistAerPA = 100*table_vals.CFtrendHistAerPA/table_vals.CFavFirstFiveHistAerPA;
table_vals.CFtrendPctHistAerPA_precision='%.9f';

%Hist-GHG
table_vals.CFtrendHistGhgPA = fscale*trends.trend_dat_box_hist_GHG{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistGhgPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistGhgPAexp_precision = '%d';
table_vals.CFtrendHistGhgPAun = fscale*trends.trend_dat_box_hist_GHG{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistGhgPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistGhgPA = 100*table_vals.CFtrendHistGhgPA/table_vals.CFavFirstFiveHistGhgPA;
table_vals.CFtrendPctHistGhgPA_precision='%.9f';

%Hist-nat
table_vals.CFtrendHistNatPA = fscale*trends.trend_dat_box_hist_nat{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistNatPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistNatPAexp_precision = '%d';
table_vals.CFtrendHistNatPAun = fscale*trends.trend_dat_box_hist_nat{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistGhgPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistGhgPA = 100*table_vals.CFtrendHistGhgPA/table_vals.CFavFirstFiveHistGhgPA;
table_vals.CFtrendPctHistGhgPA_precision='%.9f';

% -- Nd --
var_str_tab='N_d'; fscale=1e1;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.NdtrendPA = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendPAexp_precision = '%d';
table_vals.NdtrendPAun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.NdtrendPAmin = minval;
table_vals.NdtrendPAminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.NdtrendPAmax = maxval;
table_vals.NdtrendPAmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;

table_vals.NdavFirstFivePA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctPA = 100*table_vals.NdtrendPA/table_vals.NdavFirstFivePA;
table_vals.NdtrendPctPA_precision='%.9f';


%DAMIP trends and uncertainties
var_str_tab='N_d'; fscale=1e1;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_UKESM1_DAMIP_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%Hist-aer
table_vals.NdtrendHistAerPA = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistAerPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistAerPAexp_precision = '%d';
table_vals.NdtrendHistAerPAun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistAerPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistAerPA = 100*table_vals.NdtrendHistAerPA/table_vals.NdavFirstFiveHistAerPA;
table_vals.NdtrendPctHistAerPA_precision='%.9f';

%hist-GHG
table_vals.NdtrendHistGhgPA = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistGhgPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistGhgPAexp_precision = '%d';
table_vals.NdtrendHistGhgPAun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistGhgPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistGhgPA = 100*table_vals.NdtrendHistGhgPA/table_vals.NdavFirstFiveHistGhgPA;
table_vals.NdtrendPctHistGhgPA_precision='%.9f';

%hist-nat
table_vals.NdtrendHistNatPA = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistNatPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistNatPAexp_precision = '%d';
table_vals.NdtrendHistNatPAun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistNatPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistNatPA = 100*table_vals.NdtrendHistNatPA/table_vals.NdavFirstFiveHistNatPA;
table_vals.NdtrendPctHistNatPA_precision='%.9f';

% -- LWP --
var_str_tab='LWP'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.LWPtrendPA = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPtrendPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPtrendPAexp_precision = '%d';
table_vals.LWPtrendPAun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.LWPtrendPAmin = minval;
table_vals.LWPtrendPAminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.LWPtrendPAmax = maxval;
table_vals.LWPtrendPAmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;

table_vals.LWPavFirstFivePA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPtrendPctPA = 100*table_vals.LWPtrendPA/table_vals.LWPavFirstFivePA;
table_vals.LWPtrendPctPA_precision='%.9f';

% -- LWPic --
var_str_tab='LWPic'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.LWPictrendPA = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendPAexp_precision = '%d';
table_vals.LWPictrendPAun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.LWPictrendPAmin = minval;
table_vals.LWPictrendPAminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.LWPictrendPAmax = maxval;
table_vals.LWPictrendPAmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;

table_vals.LWPicavFirstFivePA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctPA = 100*table_vals.LWPictrendPA/table_vals.LWPicavFirstFivePA;
table_vals.LWPictrendPctPA_precision='%.9f';


%DAMIP trends and uncertainties
var_str_tab='LWPic'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_UKESM1_DAMIP_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%Hist-aer
table_vals.LWPictrendHistAerPA = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistAerPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistAerPAexp_precision = '%d';
table_vals.LWPictrendHistAerPAun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistAerPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistAerPA = 100*table_vals.LWPictrendHistAerPA/table_vals.LWPicavFirstFiveHistAerPA;
table_vals.LWPictrendPctHistAerPA_precision='%.9f';

%Hist-GHG
table_vals.LWPictrendHistGhgPA = fscale*trends.trend_dat_box_hist_GHG{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistGhgPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistGhgPAexp_precision = '%d';
table_vals.LWPictrendHistGhgPAun = fscale*trends.trend_dat_box_hist_GHG{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistGhgPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistGhgPA = 100*table_vals.LWPictrendHistGhgPA/table_vals.LWPicavFirstFiveHistGhgPA;
table_vals.LWPictrendPctHistGhgPA_precision='%.9f';

%Hist-nat
table_vals.LWPictrendHistNatPA = fscale*trends.trend_dat_box_hist_nat{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistNatPAexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistNatPAexp_precision = '%d';
table_vals.LWPictrendHistNatPAun = fscale*trends.trend_dat_box_hist_nat{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistNatPA = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistNatPA = 100*table_vals.LWPictrendHistNatPA/table_vals.LWPicavFirstFiveHistNatPA;
table_vals.LWPictrendPctHistNatPA_precision='%.9f';


%% -- Period 2 --
itr=2; %specifies the second trend in the cell

% -- SW TOA --
var_str_tab='SW_TOA_up'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.SWtrendPB = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.SWtrendPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.SWtrendPBexp_precision = '%d';
table_vals.SWtrendPBun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.SWtrendPBmin = minval;
table_vals.SWtrendPBminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.SWtrendPBmax = maxval;
table_vals.SWtrendPBmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;

i0=find(trends.years_ukesm_1d==yr_start_trend_box2(2));
table_vals.SWavFirstFivePB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.SWtrendPctPB = 100*table_vals.SWtrendPB/table_vals.SWavFirstFivePB;
table_vals.SWtrendPctPB_precision='%.9f';

% -- fc --
var_str_tab='totCF'; fscale=1e5;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.CFtrendPB = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendPBexp_precision = '%d';
table_vals.CFtrendPBun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.CFtrendPBmin = minval;
table_vals.CFtrendPBminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.CFtrendPBmax = maxval;
table_vals.CFtrendPBmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;
table_vals.CFavFirstFivePB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctPB = 100*table_vals.CFtrendPB/table_vals.CFavFirstFivePB;
table_vals.CFtrendPctPB_precision='%.9f';

%DAMIP trends and uncertainties
var_str_tab='totCF'; fscale=1e4;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%Hist-aer
table_vals.CFtrendHistAerPB = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistAerPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistAerPBexp_precision = '%d';
table_vals.CFtrendHistAerPBun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistAerPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistAerPB = 100*table_vals.CFtrendHistAerPB/table_vals.CFavFirstFiveHistAerPB;
table_vals.CFtrendPctHistAerPB_precision='%.9f';

%Hist-GHG
table_vals.CFtrendHistGhgPB = fscale*trends.trend_dat_box_hist_GHG{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistGhgPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistGhgPBexp_precision = '%d';
table_vals.CFtrendHistGhgPBun = fscale*trends.trend_dat_box_hist_GHG{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistGhgPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistGhgPB = 100*table_vals.CFtrendHistGhgPB/table_vals.CFavFirstFiveHistGhgPB;
table_vals.CFtrendPctHistGhgPB_precision='%.9f';

%Hist-nat
table_vals.CFtrendHistNatPB = fscale*trends.trend_dat_box_hist_nat{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.CFtrendHistNatPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.CFtrendHistNatPBexp_precision = '%d';
table_vals.CFtrendHistNatPBun = fscale*trends.trend_dat_box_hist_nat{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.CFavFirstFiveHistNatPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.CFtrendPctHistNatPB = 100*table_vals.CFtrendHistNatPB/table_vals.CFavFirstFiveHistNatPB;
table_vals.CFtrendPctHistNatPB_precision='%.9f';

% -- Nd --
var_str_tab='N_d'; fscale=1e1;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.NdtrendPB = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendPBexp_precision = '%d';
table_vals.NdtrendPBun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.NdtrendPBmin = minval;
table_vals.NdtrendPBminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.NdtrendPBmax = maxval;
table_vals.NdtrendPBmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;
table_vals.NdavFirstFivePB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctPB = 100*table_vals.NdtrendPB/table_vals.NdavFirstFivePB;
table_vals.NdtrendPctPB_precision='%.9f';

%DAMIP trends and uncertainties
var_str_tab='N_d'; fscale=1e1;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_UKESM1_DAMIP_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%hist-aer
table_vals.NdtrendHistAerPB = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistAerPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistAerPBexp_precision = '%d';
table_vals.NdtrendHistAerPBun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistAerPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistAerPB = 100*table_vals.NdtrendHistAerPB/table_vals.NdavFirstFiveHistAerPB;
table_vals.NdtrendPctHistAerPB_precision='%.9f';

%hist-GHG
table_vals.NdtrendHistGhgPB = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistGhgPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistGhgPBexp_precision = '%d';
table_vals.NdtrendHistGhgPBun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistGhgPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistGhgPB = 100*table_vals.NdtrendHistGhgPB/table_vals.NdavFirstFiveHistGhgPB;
table_vals.NdtrendPctHistGhgPB_precision='%.9f';

%hist-nat
table_vals.NdtrendHistNatPB = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.NdtrendHistNatPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.NdtrendHistNatPBexp_precision = '%d';
table_vals.NdtrendHistNatPBun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.NdavFirstFiveHistNatPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.NdtrendPctHistNatPB = 100*table_vals.NdtrendHistNatPB/table_vals.NdavFirstFiveHistNatPB;
table_vals.NdtrendPctHistNatPB_precision='%.9f';


% -- LWP --
var_str_tab='LWP'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.LWPtrendPB = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPtrendPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPtrendPBexp_precision = '%d';
table_vals.LWPtrendPBun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.LWPtrendPBmin = minval;
table_vals.LWPtrendPBminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.LWPtrendPBmax = maxval;
table_vals.LWPtrendPBmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;
table_vals.LWPavFirstFivePB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPtrendPctPB = 100*table_vals.LWPtrendPB/table_vals.LWPavFirstFivePB;
table_vals.LWPtrendPctPB_precision='%.9f';

% -- LWPic --
var_str_tab='LWPic'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4_ensemble_mean_ocean_only_UKESM1_no_obs.mat'];
trends = load(tr_file);
table_vals.LWPictrendPB = fscale*trends.trend_dat_box{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendPBexp_precision = '%d';
table_vals.LWPictrendPBun = fscale*trends.trend_dat_box{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...
%ens min and max trends
for iens=1:Nens   
    dat(iens) = fscale * trends.trend_dat_box_ens{1,itr,iens}.coeffs(2); %trend values    
end
[minval imin]=min(dat);
table_vals.LWPictrendPBmin = minval;
table_vals.LWPictrendPBminUN = fscale * trends.trend_dat_box_ens{1,itr,imin}.uncer_max;
[maxval imax]=max(dat);
table_vals.LWPictrendPBmax = maxval;
table_vals.LWPictrendPBmaxUN = fscale * trends.trend_dat_box_ens{1,itr,imax}.uncer_max;
table_vals.LWPicavFirstFivePB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctPB = 100*table_vals.LWPictrendPB/table_vals.LWPicavFirstFivePB; %Pct trend
table_vals.LWPictrendPctPB_precision='%.9f';

%DAMIP trends and uncertainties
var_str_tab='LWPic'; fscale=1e2;
tr_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Annual_mean_' var_str_tab '_for_region_4__ocean_only_UKESM1_DAMIP_no_obs.mat'];
%Annual_mean_totCF_for_region_4__ocean_only_HADGEM-GC31-LL_no_obs
trends = load(tr_file);

%Hist-aer
table_vals.LWPictrendHistAerPB = fscale*trends.trend_dat_box_hist_aer{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistAerPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistAerPBexp_precision = '%d';
table_vals.LWPictrendHistAerPBun = fscale*trends.trend_dat_box_hist_aer{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistAerPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistAerPB = 100*table_vals.NdtrendHistAerPB/table_vals.NdavFirstFiveHistAerPB;
table_vals.LWPictrendPctHistAerPB_precision='%.9f';

%Hist-GHG
table_vals.LWPictrendHistGhgPB = fscale*trends.trend_dat_box_hist_GHG{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistGhgPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistGhgPBexp_precision = '%d';
table_vals.LWPictrendHistGhgPBun = fscale*trends.trend_dat_box_hist_GHG{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistGhgPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistGhgPB = 100*table_vals.NdtrendHistGhgPB/table_vals.NdavFirstFiveHistGhgPB;
table_vals.LWPictrendPctHistGhgPB_precision='%.9f';

%Hist-nat
table_vals.LWPictrendHistNatPB = fscale*trends.trend_dat_box_hist_nat{1,itr}.coeffs(2); %N.B. - latex newcommand names can't have numbers in...
    %So, using A and B instead.
table_vals.LWPictrendHistNatPBexp = -log10(fscale); %the exponent for the data when quoting. E.g. Y x 10^{exp}
table_vals.LWPictrendHistNatPBexp_precision = '%d';
table_vals.LWPictrendHistNatPBun = fscale*trends.trend_dat_box_hist_nat{1,itr}.uncer_max; %N.B. - latex newcommand names can't have numbers in...

table_vals.LWPicavFirstFiveHistNatPB = meanNoNan(trends.dat_annual_box_ukesm(i0+4),2);
table_vals.LWPictrendPctHistNatPB = 100*table_vals.NdtrendHistNatPB/table_vals.NdavFirstFiveHistNatPB;
table_vals.LWPictrendPctHistNatPB_precision='%.9f';

%% save to .tex file

iappend=0;
save_sw_table_vals = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/Trends_noobs_for_region_4_ocean_only_TABLE';
var_str = '';
latex_newcommand_from_structure(table_vals,var_str,save_sw_table_vals,iappend);

