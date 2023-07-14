function [hhat,del_hat]=mountain_smith_solve_hhat_for_given_critical_H0_hat(H0_hat,L)
%function [hhat,del_hat]=mountain_smith_solve_hhat_for_given_critical_H0_hat(H0_hat,L)

range=[0 1]; %range in which solution lies

hhat=fzero(@mountain_smith_hhat_for_given_critical_H0_hat,range,[],H0_hat);
del_hat = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4));

hhat=hhat/L;
del_hat=del_hat/L;


