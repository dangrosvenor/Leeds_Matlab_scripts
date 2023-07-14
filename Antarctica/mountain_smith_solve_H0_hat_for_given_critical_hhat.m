function [H0_hat,N]=mountain_smith_solve_H0_hat_for_given_critical_hhat(hhat)

range=[3*pi/6 9*pi/6]; %range in which solution lies

H0_hat=fzero(@mountain_smith_H0_hat_for_given_critical_hhat,range,[],hhat);
N=H0_hat*6/pi;

