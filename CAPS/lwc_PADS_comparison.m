clear PACS_counts2
%calculates the hotwire LWC directly from the voltages using my script for comparison to
%the value given in the spreadsheet by PADS. Also calculates the Manchester CAS LWC.

%run read_PACS_LWC.m and read_PACS_CAS.m first

%the hotwire LWC
T = data_LWC_PACS(9,:);  %these are all from the spreadsheet
PMB = data_LWC_PACS(10,:);
TAS = data_LWC_PACS(11,:);
V = data_LWC_PACS(2,:);
TB = data_LWC_PACS(20,:);
[LWC_calculated,DENS,VISC,VSCW,CND,CNDW,RE,PRF,PRW,TB,P,DRYP,FACT]=lwc_calc_from_DMT_PACS_comparison(T,PMB,TAS,V,TB);



%no work out the LWC from our CAS probe.

if ~exist('flt_dat')
    flt_dat=[];
end
PACS_time = data_CAS_PACS(1,:); %seconds
PACS_counts = data_CAS_PACS(103:132,:)';
PACS_TAS_time = data_LWC_PACS(1,:); %seconds - time that corresponds to the TAS
PACS_TAS_time_temp = PACS_TAS_time;

bins_PACS_temp = bins_PACS;
PACS_counts2 = PACS_counts;
TAS_PACS=TAS;

TAS_PACS=data_CAS_PACS(163,:);
PACS_TAS_time_temp = PACS_time;

%the two lines below shift the bins to the right for testing
%PACS_counts2(:,2:30) = PACS_counts(:,1:29);
%PACS_counts2(:,1) = 0;


PACS_LWC_CAS = data_CAS_PACS(165,:);
PACS_LWC = data_LWC_PACS(7,:);

if ~exist('flt_dat')
    flt_dat=[];
end

%have passed CIP_bins and CIP_counts as [] - then the function will ignore CIP

air_speed_type = 'CIP probe';
cut_off_size = 0; %cut off size for the number concentration (microns)
%TAS is set above
LWC_PACS_cutoff = [5 30]; %do LWC for all sizes

[sample_volume_PACS,sample_volume_CIP_dummy,air_speed_1D_PACS,air_speed_PACS,PACS_total_number,...
    CIP_total_number_dummy,LWC_dist_PACS,LWC_dist_cip_dummy,PACS_mode_diameter,PACS_mean_diameter...
    ,LWC_dist_PACS_cutoff]...
    =cas_sample_volume_and_stats(flt_dat,PACS_time,...
    bins_PACS_temp,PACS_counts,PACS_TAS_time_temp,[],[],air_speed_type,...
    cut_off_size,TAS_PACS,LWC_PACS_cutoff);



disp('Finished lwc_PADS_comparison.m')