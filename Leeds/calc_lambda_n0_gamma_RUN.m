function [n0,lam,mu] = calc_lambda_n0_gamma_RUN(hm_cat,N,q,overide_mu)
% function [n0,lam,mu] = calc_lambda_n0_gamma_RUN(hm_cat,N,q)
%N.B.- the n0 here is just the number concentration - this is how it is
    %coded in the UM... 
    %Need to use the size dist in the form :-
    % ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);
    % trapz(ND)*dD should give back N (total number).
% Also, potentially useful - function [M] = calc_moment_gamma(p,N,lam, mu)    

    
switch hm_cat
    case 'rain'
        p1=3;
        p2=0;
        p3=6; %for triple moment, but not used here yet
        
        c_x = pi*1000/6; %pi*rhow/6
        d_x = 3;
        
        mu = 2.5;
        mu = 3.5;
        mu = 4.5;        
        
     case 'liq'
        p1=3;
        p2=0;
        p3=6; %for triple moment, but not used here yet
        
        c_x = pi*1000/6; %pi*rhow/6 
        d_x = 3;
        
        mu = 0;

end

if exist('overide_mu')
    mu = overide_mu;
end


% from lookup.F90 :-  m1=mass/params%c_x
m1=q./c_x;
m2 = N;

[n0,lam] = calc_lambda_n0_gamma(m1, m2, p1, p2, mu);