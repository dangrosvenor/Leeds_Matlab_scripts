function [hpeak,dA_peak,dB_peak]=Smith_dual_layer_find_hpeak(Ha,Hb,L)

dh=(Hb-pi/2)/500;
h=[0:dh:5];
for i=1:length(h)
    [dA,dB]=Smith_solve_dual_layer_solve_del(Ha,Hb,h(i),0);
    if isnan(dA)==1
        hpeak=h(i-1);
        [dA_peak,dB_peak]=Smith_solve_dual_layer_solve_del(Ha,Hb,h(i-1),0);
        dA_peak=dA_peak/L;
        dB_peak=dB_peak/L;        
        hpeak=hpeak/L;
        break
    end
end