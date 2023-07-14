function y=mountain_Smith_hhat_for_blocking(hhat,L,hm,H0)


del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
H0_hat_strat=(hhat - del_hat + acos(hhat./del_hat))/L;

y=hm-hhat/L+H0_hat_strat - H0;
