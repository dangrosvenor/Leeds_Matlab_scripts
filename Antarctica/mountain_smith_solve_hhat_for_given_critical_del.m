function [hhat,H0_hat]=mountain_smith_solve_hhat_for_given_critical_del(del_hat)

range=[0 1];

hhat=fzero(@mountain_smith_hhat_for_given_critical_del,range,[],del_hat);
H0_hat=hhat - del_hat + acos(hhat./del_hat);

