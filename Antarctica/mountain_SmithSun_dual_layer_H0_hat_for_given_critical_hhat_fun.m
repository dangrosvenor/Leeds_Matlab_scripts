function y=mountain_SmithSun_dual_layer_H0_hat_for_given_critical_hhat_fun(Hb_hat,Ha_hat,hhat)

[hpeak,dA,dB]=Smith_dual_layer_find_hpeak(Ha_hat,Hb_hat,1);

y = hpeak-hhat; %needs to be solved so that =0