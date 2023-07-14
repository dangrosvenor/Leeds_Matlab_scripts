%1D PDFs of optical depth for the vertical penetration depth Nd bias paper
% Run load_VOCALS_POLDER_saved_data.m first for laod file loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2005-2012_Aqua_mockL3_no_confidence_screening_with_POLDER.mat';



% PDF part based on DRIVER_template_1D_PDF.m
%Stripped down template for running plotTime* to produce a joint 1D PDF

%% Some flags. May want to set these outside the script :-
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;
%pdf_type_driver='normal';
pdf_type_driver='cumulative';
%LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9];


%% Define the bins
Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;  
Ybins_DRIVER = [1:3:1000];
Ybins_DRIVER = [0.5:0.5:150];

%% Define the data somewhere
%Y_driver = [1:1000]; %the data

% tau data based on L3 1x1 deg means (so not pixel level data).
% Have data for 2005-2012 here for whole VOCALS region (no confidence
% screening applied). Pixel level screening as for SZA paper I think,
% except for having no confidence screening applied here
%Y_driver = mockL3_no_conf.Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:);
Y_driver = tau_no_land;  %Cloud_Optical_Thickness_Liquid_Mean.timeseries3(:);
ylabelstr='Optical Depth';

%Define the time data
daynum_timeseries3=[1:size(Y_driver,3)];

%% Shouldn't need to change anything below here


%% Set the defaults to make sure all variables are set
watervap_defaults
pdf2D_defaults



i577 = 'MODIS_plot_UW';




iset_min_clim=1;
clim_min=0;
iset_max_clim=1;
clim_max=200;

logflag=0;
dlogflag=0;

isave_plot=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

Ybins=Ybins_DRIVER; ichoose_Ybins=1;

datatype = 'makeshift';
datatype = 'none';
gcm_str='UM';
LAT_UM = MLAT;
LON_UM = MLON;




x_axis_vals = 'Dummy data'; %dummy data
%y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data
y_axis_vals = 'General y-axis'; %with lat and lon screening

graph = 977; %new 1D PDF from 2D histo data - can choose either axis
%(for watervap)

axis1D = 'y';



logbin_norm = logbin_norm_driver;
i_plot_norm=i_plot_norm_driver;
i_div_bin_widths=i_div_bin_widths_driver;
pdf_type = pdf_type_driver;
LAT_val = LAT_val_DRIVER; LON_val = LON_val_DRIVER

%The code below to wrap around the boxes works, but perhaps best to do two
%boxes and combine the PDFs to avoid confusion

% %to allow wrapping around or box :-
% if LON_val(2) < LON_val(1)
%     lon0 = LON_val(1);
%     d = lon0 + 180;
%     LON_UM = MLON - d;
%     i180 = find(LON_UM<-180);
%     LON_UM(i180) = LON_UM(i180)+360;
%     
%     LON_val = LON_val - d;
%     i180 = find(LON_val<-180);
%     LON_val(i180) = LON_val(i180)+360;
%     
% end

% --------- Override flags for 2D PDF --------
ioverride_pdf=1;
%iocean_only=1;
man_choose_plotTimeHeight_graph=1;
ioverride_location_selection=1;
ioverride_pdf_varchoose = 1;

% --------- Override flags for watervap --------
man_choose_water_graph=1;    %for watervap

%---  Run plot script and save
plotTimeHeightVap3
ioverride_pdf=1;
close(gcf);
waterVapourMay2005
%close(gcf);

                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
