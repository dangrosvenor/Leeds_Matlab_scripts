line_cols='Rosenfeld';
nmin=100;


%% Make the CF vs Nd plots
%      cols='rbgkmrcymkgbrrrgykmc';
cols_rgb{1}=[0 0 1];
cols_rgb{2}=[1 0 1];
cols_rgb{3}=[0 0.7 0.7];
cols_rgb{4}=[1 1 0];
cols_rgb{5}=[1 0.7 0];
cols_rgb{6}=[1 0.3 0];
cols_rgb{7}=[1 0 0];
cols_rgb{8}=[0 0 1];
cols_rgb{9}=[1 0 1];
cols_rgb{10}=[0 0.7 0.7];
cols_rgb{11}=[1 1 0];
cols_rgb{12}=[1 0.7 0];
cols_rgb{13}=[1 0.3 0];
cols_rgb{14}=[1 0 0];
cols_rgb{15}=[0 0 1];
cols_rgb{16}=[1 0 1];
cols_rgb{17}=[0 0.7 0.7];
cols_rgb{18}=[1 1 0];



%patt={'-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--','-','--'};
patt={'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};

%figure('color','w');
%set(gcf,'position',[20 210 600 384])
for iLWP=1:length(LWP_bins)-1
    xvals = xvals_LWP{iLWP};
    switch line_cols
        case 'Reverse of Rosenfeld'
            ind_col = iLWP;
        case 'Rosenfeld';
            ind_col = 8-iLWP;
    end
    
    
    xvals(NY_vals_LWP{iLWP} < nmin) = NaN;
    
    %plot( xvals , yvals_LWP{iLWP} ,['x' cols(iLWP) patt{iLWP}],'linewidth',3);
    plot( xvals , yvals_LWP{iLWP} ,['o' patt{iLWP}],'linewidth',3,'color',cols_rgb{ind_col},'markerfacecolor',cols_rgb{ind_col});
    hold on
    leg_str{iLWP}=['LWP=' num2str(LWP_bins(iLWP),'%.1f') '-' num2str(LWP_bins(iLWP+1),'%.1f') ];
end

if ilegend==1
    hL = legend(leg_str,'location','northeastoutside');
end

set(gca,'xscale','log');
set(gca,'xlim',[10 300]);
set(gca,'ylim',[0 1]);
set(gca,'fontsize',16);
xlabel(xlabelstr);
ylabel(ylabelstr);
titlenam = remove_character( ['UKCA (' um_case_PD '), ' CF_type ', ' CF_type_LWP, ', restrict_height_method=' restrict_height_method ', iocean_only=' num2str(iocean_only) ', nmin=' num2str(nmin)] ,  '_',' ');
titwrapped = wrap_title_to_nlines(titlenam,50,3);
%title(titwrapped);
title(tit_str);
grid on

set(gca,'xlim',[10 500]); set(gca,'ylim',[0 1.0])

save_str = [region_choice '_2D_PDF_cf_vs_Nd'];

titlenam_driver = [um_case_PD '_Model_CFvsNd_' save_str];
savename=[savedir_date titlenam_driver];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
%saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);

isave_data=0;
if isave_data==1
    save([savename '.mat'],'-V7.3','xvals_LWP','yvals_LWP','NY_vals_LWP','titlenam');
end


%savedir='/home/disk/eos1/d.grosvenor/modis_work/ACSIS/';
%UM_save_plot(gcf,isave_mean_2D_PDF,savedir_date,[save_str] );







