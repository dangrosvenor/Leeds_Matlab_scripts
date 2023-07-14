r_single_all = [50:100:450]*1e-9; %radius in m

%N_type = 'prescribed number';
N_type = 'prescribed single aerosol mass';

N_aero=1; %dummy value
        
for iab=1:length(r_single_all)
    r_single = r_single_all(iab);
    nenes_run_multi
end