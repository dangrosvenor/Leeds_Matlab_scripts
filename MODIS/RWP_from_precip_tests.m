%Test why there is a depedency on teh initial v_fall used


%Will plot vs LWP for different fall speed value and N values
N = 25; %cm^-3  - set Nd to 50 per cc for convienience
Hcb = 500; %Height of cloud base - assuming rho_air*qR is constant up to cloud base - likely not true due to rain evap
 %Fairly large uncertainty here, but can scale linearly.
LWP = [0:10:180]; %g/m2
%Estimate precip rate for various LWPs up to AMSRE threshold for precip
rv = 40e-6;
k=0.6; %k for the rain dist - not actully needed for gamma dists, etc.

P = precip_rate_Wood_2008(LWP,N); %returned in mm/hr


v_falls_init = [0.1:10:100]; %m/s

clear RWPs taus res V_falls lams out

for i=1:length(v_falls_init)
    v_fall_init = v_falls_init(i);
    %[RWPs{i},taus{i}, re{i}, V_falls{i}, lams{i}] = 
    [out{1:5}] = rwp_from_precip_rate(P,v_fall_init,Hcb,rv,k);
    RWPs(:,i) = out{1}(:);
    taus{i} = out{2};
    res{i} = out{3};    
    V_falls{i} = out{4};
    lams{i} = out{5};
end

figure
%for i=1:length(RWPs)
%    hold on
%    plot(LWP,RWPs{i}(:));
%end

plot(LWP,RWPs);