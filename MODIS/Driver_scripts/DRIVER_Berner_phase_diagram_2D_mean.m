%Stripped down template for running plotTime* to produce a joint 1D PDF
% N.B. - this does no lat,lon screening

%% Some flags. May want to set these outside the script :-
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;

LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9]; %%Still need to set this even though not doing lat lon screening

%% Define the bins
%Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;

% if inormalise==0
%     Xbins_DRIVER = [0:50:300]; ichoose_Ybins=1;
%     Ybins_DRIVER = [0:20:140]; ichoose_Xbins=1;
% else
% %    Xbins_DRIVER = [0:0.1:5.1]; ichoose_Ybins=1;
% %    Ybins_DRIVER = [0:0.1:5.2]; ichoose_Xbins=1;
%     
%     Xbins_DRIVER = [-10:1:10]; ichoose_Ybins=1;
%     Ybins_DRIVER = [-10:1:11]; ichoose_Xbins=1;    
% end

%% Define the data 
%Y_driver = magic(10); %the data
%X_driver = magic(10)'; %the data
%Z_driver = ones([10 10]); %[1:100];
xlabelstr = 'N_d (cm^{-3})';
ylabelstr = 'LWP (g m^{-2})';

%% Set the defaults to make sure all variables are set
%watervap_defaults
pdf2D_defaults


ndims_hist=3; %set to 3D and put the data to average on the 3rd dimension (as Z)

i577 = 'MODIS_plot_UW';

iset_min_clim=0;
clim_min=0;
iset_max_clim=0;
clim_max=200;

logflag=0;
dlogflag=0;

isave_plot=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

Ybins=Ybins_DRIVER; ichoose_Ybins=1;
Xbins_driver=Xbins_DRIVER; ichoose_Xbins=1;


datatype = 'makeshift';
gcm_str='';

%x_axis_vals = 'Dummy data'; %dummy data
x_axis_vals = 'General GCM-style x-axis simple2';
y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data
z_axis_vals = 'Z from outside script';



LAT_val = LAT_val_DRIVER; LON_val = LON_val_DRIVER; %Still need to set this, though, even though not doing lat lon screening

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
%close(gcf);
%waterVapourMay2005
%close(gcf);

                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
