iload_data=0;

file_load = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/Nd_saved_2003-2015.mat';

if iload_data==1
    load(file_load);
end

rhow=1000;
LWP21 = 5/9 * rhow * Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6 .* Cloud_Optical_Thickness_Liquid_Mean.timeseries3; %kg/2
LWP37 = 5/9 * rhow * Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6 .* Cloud_Optical_Thickness_Liquid_Mean.timeseries3; %kg/m2

%% Spatial relationships (using time-averaged data)
LWP21_time_mean = meanNoNan(LWP21,3);
LWP37_time_mean = meanNoNan(LWP37,3);
Nd21_time_mean = meanNoNan(Droplet_Number_Concentration.timeseries3,3);
Nd37_time_mean = meanNoNan(Droplet_Number_Concentration_37.timeseries3,3);

%run 2D histo template script :-
X_driver = Nd37_time_mean;
Y_driver = meanNoNan(LWP37,3)*1e3; %convert to g/m2
xlabelstr='N_d (cm^{-3})';
ylabelstr='LWP (g m^{-2})';
DRIVER_2D_pdf_LWP_vs_Nd

return

plot(mid_Xbins,Y_mean,'wo');

set(gca,'ylim',[-90 10]);
xdat_save{isave_xdat}=mid_Xbins;
ydat_save{isave_xdat}=Y_mean;
xlabelstr_save{isave_xdat} = xlabelstr;
isave_xdat=isave_xdat+1;

axis1D = 'x';
waterVapourMay2005
xdat_save_pdf{isave_xdat_pdf}=xdat(1).x;
ydat_save_pdf{isave_xdat_pdf}=ydat(1).y;
xlabelstr_save_pdf{isave_xdat_pdf} = xlabelstr;
isave_xdat_pdf=isave_xdat_pdf+1;