%Stripped down template for running plotTime* to produce a joint 1D PDF

%% Some flags. May want to set these outside the script :-
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;
pdf_type_driver='cumsum';
LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9]; %
%LAT_val_DRIVER = [-1e9 -80]; LON_val_DRIVER = [-1e9 1e9]; %

%% Define the bins
Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;   
Ybins_DRIVER = [-92:2:98];

%% Define the lat and lon grids
LAT_UM = [-89.5:1:89.5];
LON_UM = [-179.5:1:179.5];

%% Define the data somewhere - of shape [ilat,ilon,itime]
Y_driver = repmat(LAT_UM',[1 360 10]); %the data - e.g. 10 days of data [180,360,10]
Y_driver(2,:,:)=97;
ylabelstr='LWP (g m^{-2})';

%% Define the time data (or just use indices)
daynum_timeseries3=[1:size(Y_driver,3)];



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



%datatype = 'makeshift';
datatype = 'none'; %needed to make it run the lat/lon selection
gcm_str='UM';

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
close(gcf);
waterVapourMay2005
%close(gcf);

                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
