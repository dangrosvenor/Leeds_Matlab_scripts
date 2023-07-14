%run Boutle_case_12Nov2008_NdLWP_joint_PDFs_20151111.m first
figure

fsize=16;
lwidth=4;

if iswap_xy_DRIVER==1
    datX = mid_Xbins;
    datY = Y_mean;
    Ndatap = sum(qh(1:end-1,1:end-1),1);
else
    datX = mid_Ybins;
    datY = X_mean;
    Ndatap = sum(qh(1:end-1,1:end-1),2);
end

datY(Ndatap<30)=NaN;
%datY(datX<20)=NaN;
plot(datX,datY,'b-','linewidth',lwidth);
grid on

xlabel('N_d (cm^{-3})');
ylabel('LWP (g m^{-2})');

fontsize_figure(gcf,gca,fsize);

%set(gca,'xscale','log');
set(gca,'xlim',[-10 300]);

savedir = '/home/disk/eos1/d.grosvenor/modis_work/lwp_vs_Nd_plots/';
savename = [savedir 'lwp_vs_nd_from_histogram'];
%saveas_ps_fig_emf(gcf,[savename],'',0,1);

