function [Hb_hat,hpeak]=mountain_SmithSun_dual_layer_solve_H0_hat_for_given_critical_hhat(Ha_hat,hhat)

range=[3*pi/6 10*pi/6]; %range in which solution lies

Hb_hat=fzero(@mountain_smith_H0_hat_for_given_critical_hhat,range,[],Ha_hat,hhat);
[hpeak,dA,dB]=Smith_dual_layer_find_hpeak(Ha_hat,Hb_hat,1);

