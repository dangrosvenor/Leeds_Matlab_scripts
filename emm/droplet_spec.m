 Dg=0.090870e-6; %mode diameter
 sig=6.223699;
 

 Nbins=10;
 N=569.64368e6; %total no. conc m^-3
 
 
Dlow=Dg/100/sig;
Dhigh=Dg*100*sig;
Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing

logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space

Dmid=(logD(2:end) + logD(1:end-1) )/2;
Dmid=exp(Dmid);



nd=lognormal(Dmid,N,sig,Dg)*Dspacing; %total no. in each bin

figure;
plot(Dmid*1e6,nd);
m=1000*nd.*4./3.*pi.*(Dmid/2).^3;


%sig=sig*2;


% dbar=[0.01:0.01:0.5]*1e-6;
% 
% for i=1:length(dbar)
%     Dg=dbar(i)
%     Dlow=Dg/100/sig;
% 	Dhigh=Dg*100*sig;
% 	Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing
% 	
% 	logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space
% 	
% 	Dmid=(logD(2:end) + logD(1:end-1) )/2;
% 	Dmid=exp(Dmid);
% 
% 
%     nd2=lognormal(Dmid,N,sig,Dg)*Dspacing; %total no. in each bin
% 	m2=1000*nd.*4./3.*pi.*(Dmid/2).^3;
%     msave(i)=sum(m2);
% end
%     
% 'done'
% 
% 
% 
% 
