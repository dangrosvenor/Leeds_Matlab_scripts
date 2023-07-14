function [output] = SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02(UM_run,file_suffix,data_dir)

% Called from SW_calcs_James_Weber_Jan2022_load_input_data.m


switch UM_run
    case {'u-ce358','u-ce419'}
        iNd_only = 1;
    otherwise
        iNd_only = 0;
end


vn_str = 'vn1105';
vn_str = 'vn1200';

%% 1) Load ensemble monthly mean data
if iNd_only==0
    
    % -- SW TOA up (to compare against the calculation) --
    % stash/dir 01217
    %UM_run = 'u-ck596'; %control run
    stash = '01217'; var=['UM_m01s01i217_' vn_str];
    SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    output.SW_up_TOA = squeeze(dat(:,end,:,:));
    
    %SW down TOA - stash/dir 01218
    %UM_run = 'u-ck596'; %control run
    stash = '01218'; var=['UM_m01s01i218_' vn_str];
    SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    output.SW_down_TOA = squeeze(dat(:,end,:,:));
    
    
    
    % --- SW TOA upwelling CLEAR-SKY - for use directly for overall SWup calculations : (1-fc) part---
    % --- Not used for the transmissivity above clouds.
    %SW up TOA clear-sky- stash/dir 01219
    %UM_run = 'u-ck596'; %control run
    stash = '01219'; var=['UM_m01s01i219_' vn_str];
    SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    output.SW_up_TOA_clear_sky = squeeze(dat(:,end,:,:));
    
    % --- SW surface upwelling CLEAR-SKY - for calculation of surface albedo (not used for transmissivity in this setup) ---
    output.SW_up_surf_clear_sky = squeeze(dat(:,1,:,:));
    
    
    % --- SW surface downwelling CLEAR-SKY - for calculation of surface albedo (not used for transmissivity in this setup) ---
    % stash/dir 01220
    %UM_run = 'u-ck596'; %control run
    stash = '01220'; var=['UM_m01s01i220_' vn_str];
    SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    output.SW_down_surf_clear_sky = squeeze(dat(:,1,:,:));
    
    
    
    
    
    % --- total CF ---
    % stash/dir 01220
    %UM_run = 'u-ck596'; %control run
    stash = '9217'; var=['UM_m01s09i217_' vn_str];
    SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    output.totCF = dat;
    
    
end

% --- LWP  ---
% stash/dir 01220
%UM_run = 'u-ck596'; %control run
stash = '1224'; var=['UM_m01s01i224_' vn_str];
SW_calcs_James_Weber_Jan2022_load_input_data_FUNC

%This is the layer LWP (presumably layer LWC multiplied by
%layer dz) multiplied by the cloud weight. The cloud weight
%deals with the time-averaging and is 0 for LWC<1e-6 and 1
%otherwise.
% Can perhaps ignore the weighting (since it will just zero the LWC values that are already
% very small) and just sum vertically to
% get the all-sky LWP - then divide by the total CF?

lwp_all_sky_levs = 1e3*dat; %convert to g/m2 from kg/m2
if iNd_only==0
    lwp_all_sky = squeeze(sum(dat,2));
    output.lwp_ic = lwp_all_sky ./ output.totCF;
    output.lwp_ic(output.totCF<0.05) = NaN;
end



% Nd
% stash/dir 01241 - CDNC * layer weight
%UM_run = 'u-ck596'; %control run

%CDNC
stash = '1241'; var=['UM_m01s01i241_' vn_str];
SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
Nd_times_wt = dat; %cm^-3

% %LWC in each layer * weight
% stash = '1242'; var=['UM_m01s01i242_' vn_str];
% SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
% LWC_times_wt = dat; %g/m3

%Weight for each layer
stash = '1223'; var=['UM_m01s01i223_' vn_str];
SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
layer_wt = dat;

%Calculate in-cloud Nd for each layer
output.Nd_ic = Nd_times_wt ./ layer_wt;
output.Nd_ic(layer_wt<1e-20) = NaN;

%Verticaly integrate in-cloud CDNC weighted by LWC
%dz = lwp_all_sky_levs ./ LWC_times_wt;
%dz(LWC_times_wt < 1e-20) = NaN;
%Actually can just weight by layer LWP that should be LWC*dz
[output.Nd,N,std_dev] = meanNoNan(output.Nd_ic,2,'',1,1,lwp_all_sky_levs);


output.lat2d = output.gcm_Plat2D_UM;
output.lat2d_edges = output.gcm_Plat2D_edges_UM;
output.lon2d = output.gcm_Plon2D_UM;
output.lon2d_edges = output.gcm_Plon2D_edges_UM;

% --- Aerosol optical depth


% --- Absorbing aerosol optical depth


% --- Dust optical depth





