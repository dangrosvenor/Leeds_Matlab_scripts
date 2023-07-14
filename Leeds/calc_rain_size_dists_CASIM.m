dD=0.1e-9; D=[dD:dD:200e-6];
q=0.1e-3; % kg/kg
N=10e6; % /kg

%N.B.- the n0 here is just the number concentration - this is how it is
    %coded in the UM... 
    %Need to use the size dist in the form :-
    % ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);
    
    mu_in = 2.5; %default mu
[n0,lam,mu] = calc_lambda_n0_gamma_RUN('rain',N,q,mu_in); ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);

clear leg
iplot=1;

figure
plot(D*1e6,ND,'k'); leg{iplot}=['mu=' num2str(mu)]; iplot=iplot+1;
hold on

mu_in = 1.5;
[n0,lam,mu] = calc_lambda_n0_gamma_RUN('rain',N,q,mu_in); ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);
plot(D*1e6,ND,'b'); leg{iplot}=['mu=' num2str(mu)]; iplot=iplot+1;

mu_in = 3.5;
[n0,lam,mu] = calc_lambda_n0_gamma_RUN('rain',N,q,mu_in); ND = n0.* ( lam.^(1+mu) / gamma(1+mu) ) .*D.^mu.*exp(-lam.*D);
plot(D*1e6,ND,'r--'); leg{iplot}=['mu=' num2str(mu)]; iplot=iplot+1;

 set(gca,'xlim',[0 80]);
 legend(leg);
 title(['q=' num2str(q*1e3) ' g/kg, N=' num2str(N/1e6) ' mg^{-1}']);
 xlabel('D (\mum)');
 ylabel('N(D)');

