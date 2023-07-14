function y=mountain_smith_H0_hat_for_given_critical_hhat(H0_hat,hhat)

del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));
y = hhat - del_hat + acos(hhat./del_hat) - H0_hat; %needs to be solved so that =0