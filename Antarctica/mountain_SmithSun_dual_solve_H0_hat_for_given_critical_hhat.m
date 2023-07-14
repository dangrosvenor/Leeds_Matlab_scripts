function [Hb_hat,hpeak]=mountain_SmithSun_dual_solve_H0_hat_for_given_critical_hhat(Ha_hat,hhat)

range=[4*pi/6 10*pi/6]; %range in which solution lies

Hb_hat=fzero(@mountain_SmithSun_dual_layer_H0_hat_for_given_critical_hhat_fun,range,[],Ha_hat,hhat);
[hpeak,dA,dB]=Smith_dual_layer_find_hpeak(Ha_hat,Hb_hat,1);

