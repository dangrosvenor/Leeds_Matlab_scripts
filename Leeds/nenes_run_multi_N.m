%r_single_all = [50:100:450]*1e-9; %radius in m
Naero_all = [1000:1000:10000]*1e6; %per m3

Naero_all = [10000]*1e6; %per m3
r_single = NaN; %dummy value

N_type = 'prescribed number';
%N_type = 'prescribed single aerosol mass';
        
for iab=1:length(Naero_all)
    N_aero0 = Naero_all(iab);
    nenes_run_multi
end