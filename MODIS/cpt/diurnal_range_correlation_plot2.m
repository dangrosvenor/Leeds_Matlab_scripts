%diurnal range correlation plot
%use data stored in 
% /home/disk/eos1/d.grosvenor/saved_vars/TLWP_PDFs_CAPT_ocean_only_1Dpdfs.mat
% /home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_CAPT_ocean_only_1Dpdfs.mat
% /home/disk/eos1/d.grosvenor/saved_vars/Precip_PDFs_CAPT_ocean_only_1Dpdfs.mat

%


diurnal_case = 'Whole domain (ocean only) - LAT_val = [-40.5 10.5]; LON_val = [-140 -68]';
diurnal_case = 'Sc deck only - LAT_val = [-30.5 -5.5]; LON_val = [-105 -68]'           ;

plot_case = 'diurnal range';
%plot_case = 'mean_vals';

savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';
           
switch diurnal_case
    case 'Whole domain (ocean only) - LAT_val = [-40.5 10.5]; LON_val = [-140 -68]'
        file_TLWP = '/home/disk/eos1/d.grosvenor/saved_vars/TLWP_PDFs_CAPT_ocean_only_1Dpdfs.mat'; load(file_TLWP);
        file_LWP = '/home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_CAPT_ocean_only_1Dpdfs.mat'; load(file_LWP);
        file_precip = '/home/disk/eos1/d.grosvenor/saved_vars/Precip_rate_scatter_all_ocean_1Dpdfs.mat'; load(file_precip);
    case 'Sc deck only - LAT_val = [-30.5 -5.5]; LON_val = [-105 -68]'
        file_TLWP = '/home/disk/eos1/d.grosvenor/saved_vars/TLWP_PDFs_CAPT_ocean_only_Sc_deck_1Dpdfs.mat'; load(file_TLWP);
        file_LWP = '/home/disk/eos1/d.grosvenor/saved_vars/LWP_PDFs_CAPT_ocean_only_Sc_deck_1Dpdfs.mat'; load(file_LWP);
        file_precip = '/home/disk/eos1/d.grosvenor/saved_vars/Precip_rate_scatter_all_ocean_Sc_deck_1Dpdfs.mat'; load(file_precip);
end

figure
% *dop^
% krybmwc
syms = {'^','d','p','o','d','p','o','d','o'};
cols={'k','r','r','r','c','c','c','b','b'};
% 
% Xbins_Precip_day_AM3                                    1x2001            16008  double              
%   Xbins_Precip_day_AM3_CLUBBv1_2deg                       1x2001            16008  double              
%   Xbins_Precip_day_AM3_CLUBBv2_COSP_200km                 1x2001            16008  double              
%   Xbins_Precip_day_CAM5_CLUBB_COSP                        1x2001            16008  double              
%   Xbins_Precip_day_CAM5_prepostLWP                        1x2001            16008  double              
%   Xbins_Precip_day_CAMCLUBBv2_prepostLWP                  1x2001            16008  double              
%   Xbins_Precip_night_AM3 

gcm_models_to_load_LWP = {...
     'AMSRE_time3';...
     'CAM5_prepostLWP';...
     'CAM5_CLUBB_COSP';...
     'CAMCLUBBv2_prepostLWP';...
     'AM3';...
     'AM3_CLUBBv1_2deg';...
     'AM3_CLUBBv2_COSP_200km';...
     'CAM5_prepostLWP';...
     'CAMCLUBBv2_prepostLWP';...
    };

gcm_models_to_load_precip = gcm_models_to_load_LWP;
gcm_models_to_load_precip{1} = 'CLOUDSAT_PRECIP';

nice_labels = {...
     'OBS';...    
     'CAM5 LWP';...
     'CAM5CLUBBv1 LWP';...
     'CAM5CLUBBv2 LWP';...
     'AM3';...
     'AM3CLUBBv1 LWP';...
     'AM3CLUBBv2 LWP';...
     'CAM5 TLWP';...
     'CAMCLUBBv2 TLWP';...
    };

%lab_xoffset = 24*0.0005 * ones(length(LWP_or_TLWP));
lab_xoffset = -0.15 * ones([1 length(nice_labels)]);

lab_xoffset(1) = 0.025;
lab_xoffset(2) = lab_xoffset(2) + 0.03;
lab_xoffset(5) = lab_xoffset(5) + 0.17;
%lab_xoffset(8) = -0.15;
lab_xoffset(8) = lab_xoffset(8) + 0.025;
lab_xoffset(7) = lab_xoffset(7) + 0.07;



lab_yoffset = -0.6 * ones([1 length(nice_labels)]);

lab_yoffset(3) = lab_yoffset(3) + 2.2;
lab_yoffset(4) = lab_yoffset(4) - 1;
lab_yoffset(9) = lab_yoffset(9) - 1;
lab_yoffset(6) = lab_yoffset(6) - 1;
lab_yoffset(7) = lab_yoffset(7) - 1;


LWP_or_TLWP = {...
     'TLWP';...    
     'LWP';...
     'LWP';...
     'LWP';...
     'LWP';...
     'LWP';...
     'LWP';...
     'TLWP';...
     'TLWP';...
    };




nmodels=length(gcm_models_to_load_LWP);

for im=1:nmodels

    switch plot_case
        case 'diurnal range'
            y = eval(['Y_mean_overall_' LWP_or_TLWP{im} '_night_' gcm_models_to_load_LWP{im} ' - Y_mean_overall_' LWP_or_TLWP{im} '_day_' gcm_models_to_load_LWP{im}]);
            x = 24*eval(['Y_mean_overall_Precip_night_' gcm_models_to_load_precip{im} '- Y_mean_overall_Precip_day_' gcm_models_to_load_precip{im}]);
            xlab='Diurnal range of precipitation rate (mm day^{-1})';
            ylab='Diurnal range of LWP or TLWP (g m^{-2})';
            set(gca,'xlim',[0 0.8]);
            
        case 'mean_vals'
            ilon=3;
            y = meanNoNan(lwp_diurnal_save.ydat(im).y(ilon),2);
            x = meanNoNan(24*precip_diurnal_save.ydat(im).y(ilon),2);
            xlab='Mean Precip (mm day^{-1})';
            ylab='Mean LWP (g m^{-2})';
            set(gca,'xlim',[0 2]);
            title(['LON=' num2str(lwp_diurnal_save.xdat(im).x(ilon))]);
    end
    
    leg{im}=nice_labels{im};
%    leg{1} = 'OBS';
    h=plot(x,y,[cols{im} syms{im}]);
    set(h,'markerfacecolor',cols{im},'markeredgecolor','k','markersize',12)
    text(x+lab_xoffset(im),y+lab_yoffset(im),leg{im});
    hold on
end

xlabel(xlab);
ylabel(ylab);
fontsize_figure(gcf,gca,16);


savename=[savedir 'Diurnal_correlation_plot_CPT'];






