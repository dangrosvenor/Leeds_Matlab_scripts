function y=mountain_smith_hhat_for_given_critical_del(hhat,del)

y = -1/sqrt(2) * sqrt(hhat.^2 + hhat.*sqrt(hhat.^2+4)) - del;