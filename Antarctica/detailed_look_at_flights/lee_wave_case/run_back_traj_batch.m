
W_offsets=[-1.5:0.25:0 1.5];

for irun=1:length(W_offsets)
    W_offset=W_offsets(irun);

    backwards_profile_projection

    figure
    pcolor(X2,Zsound_pressure,equiv3');shading flat;colorbar
    set(gca,'xlim',[1.33 1.58]*1e5);
    set(gca,'ylim',[2400 4300]);
    hold on;
    plot(X2,Z2,'k','linewidth',2);
    
    W_str=num2str(W_offset);
    picname=['equiv_back_traj_Wplus_' W_str '_qvliq_cas_RH_thresh0.85_wideview'];
    
    saveas_ps_fig_emf(gcf,[filedir picname]);
        
end