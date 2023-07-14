function [H0_Houghton,hpeak_Houghton,F0,h_H0]=Hougton_critical_values(del_peak,U,gd,H,HM)
%[H0_Houghton,hpeak_Houghton,F0,h_H0)=Hougton_critical_values(del_peak,U,gd,H,HM) - heights in metres

thi_crit = H - HM + del_peak;
u_crit = sqrt(gd*thi_crit); 
K4 = u_crit*thi_crit;
H0_Houghton = K4/U;
    %if know del at the mountain top then we know thi_crit and so u_crit from u_crit=sqrt(g*thi_crit)
    %then can calculate K4=u_crit*thi_crit, which gives hA, if we know uA, from uA*hA=K4
    %Want to try and get heff from the displacement of the streamline and also consider
    %the displacement at the mountain crest to calculated hA?

F0=U/sqrt(gd*H0_Houghton);
hpeak_Houghton = H0_Houghton .* ( 1 + 0.5*F0.^2 - 1.5*F0.^(2/3) );
h_H0 = hpeak_Houghton/H0_Houghton;