SW_calcs_Hawaii_setup_vals %to set cf_min


f0_orig = CF_lwc_weighted_total_column_to_zdomain_top_PI_ALL;
W0_orig = LWP_PI_ALL ./ f0_orig;
inan = find(f0_orig<cf_min);
W0_orig(inan) = NaN;
N0_orig = Nd_PI_ALL/1e6; %per cc
N0_orig(inan) = NaN;


f1_orig = CF_lwc_weighted_total_column_to_zdomain_top_PD_ALL;
W1_orig = LWP_PD_ALL ./ f1_orig;
inan = find(f1_orig<cf_min);
W1_orig(inan) = NaN;
N1_orig = Nd_PD_ALL/1e6; %per cc
N1_orig(inan) = NaN;

        