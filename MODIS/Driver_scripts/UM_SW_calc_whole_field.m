%Run 
    %  UM_SW_get_whole_fields_for_calc.m        first to get data

SW_in = 1118.3;
A_clear = 0.05;
transmission_atmos = 0.75;
LWP_thresh=20; %(g/m2)

f0 = 1; %assume cf=1 for individual columns
W0 = UM_dat_multi{iLWP}{1}.datUM{1};
N0 = UM_dat_multi{iNd}{1}.datUM{1};
SW_down_surf = UM_dat_multi{iSW_down_surf_LWP_LT_0pt1}{1}.datUM{1};
SW_down_TOA = UM_dat_multi{iSW_down_TOA}{1}.datUM{1};

icloud=find(W0>LWP_thresh);

W02=NaN*ones(size(W0));
W02(icloud)=W0(icloud);
W0=W02;

W0_mean = meanNoNan(W0(icloud),1);
N0_mean = meanNoNan(N0(icloud),1);

%remove the NaNs from Nd where LWP below a threshodl for testing different
%LWP thresholds
%N0(isnan(N0)) = N0_mean; %in case Nd was based on a different threshold

%recalculate N0_mean using the values with the new ones filled in, but just
%for icloud
N0_mean = meanNoNan(N0(icloud),1);

SW0 = UM_dat_multi{1}{1}.datUM{1};
run_name_str = UM_dat_multi{1}{1}.UM_case_out.labs_UM.l;

i=1;
%[temp,temp,temp,Ac_calc,tau_calc,A_calc,SW_calc] = SW_partitioning(f0,W0,N0,1,1,1,SW_in,A_clear,transmission_atmos);
[Ac,tau,A,SW_out] = calc_SW(f0,W0,N0,SW_in,A_clear,transmission_atmos);

% A has the NaNs filled in with the clear-sky albedo, Ac has them left as
% NaN. N.B. Nd is also NaN below LWP=20 I think.


f0_mean = length(icloud) ./ length(W0(:));


[Ac_mean,tau_mean,A_mean,SW_mean_vals] = calc_SW(f0_mean,W0_mean,N0_mean,SW_in,A_clear,transmission_atmos);

SW0_mean = meanNoNan(SW0(:),1);  %The actual SW up TOA from the model
SW_out_mean = meanNoNan(SW_out(:),1); %The mean SW using the individual columns



%Estimate the surface albedo
%SW_down_surf = load('/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/xmmz-x/xmmzx_SW_down_surf_LWP_LT_0pt1_.pp.nc_timeseries.mat');
%iclear = find(W0<LWP_thresh);

%SW_down_surf is just for LWP<0.01, so is already clear sky only
A_surf_estimate = SW0./transmission_atmos ./ SW_down_surf;
A_surf_est_mean = meanNoNan(A_surf_estimate(:),1);
Tr_estimate = SW_down_surf ./ SW_down_TOA;
Tr_est_mean = meanNoNan(Tr_estimate(:),1);


SW_mean_vals./SW0_mean   %ratio of value computed using mean values to acutal SW from model
SW_out_mean./SW0_mean   %ratio of value computed using column values to acutal SW from model
SW_calc ./ SW'  %ratios from the other script

