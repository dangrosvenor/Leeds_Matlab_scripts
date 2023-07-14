%Stripped down template for running plotTime* to produce a joint 1D PDF
% N.B. - this does no lat,lon screening

%% Some flags. May want to set these outside the script :-
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;
pdf_type_driver='cumsum';
LAT_val_ALL = [-1e9 1e9]; LON_val_ALL = [-1e9 1e9]; %%Still need to set this even though not doing lat lon screening

%% Define the bins
%Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;   
Xbins_DRIVER = [-0.05:0.1:1.05]; ichoose_Xbins_DRIVER=1;
Ybins_DRIVER = [-39:2:39]; ichoose_Ybins_DRIVER=1;
Ybins_DRIVER = [-70:2:0]; ichoose_Ybins_DRIVER=0;
Ybins_DRIVER = [-870:2:850]; ichoose_Ybins_DRIVER=1;

%% Define the data somewhere
% Y_driver = magic(10); %the data
% X_driver = magic(10)'; %the data
%xlabelstr='Low cloud fraction';
ylabelstr='SW surface forcing (W m^{-2})';

%% Set the defaults to make sure all variables are set
watervap_defaults
pdf2D_defaults

i577 = 'MODIS_plot_UW';

iset_min_clim=0;
clim_min=0;
iset_max_clim=0;
clim_max=200;

logflag=0;
dlogflag=0;


isave_plot=0;
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

Ybins=Ybins_DRIVER; ichoose_Ybins = ichoose_Ybins_DRIVER;
Xbins_driver=Xbins_DRIVER; ichoose_Xbins = ichoose_Xbins_DRIVER;


datatype = 'makeshift';
gcm_str='';

%x_axis_vals = 'Dummy data'; %dummy data
x_axis_vals = 'General GCM-style x-axis simple2';
y_axis_vals = 'General y-axis no ilat simple'; %not doing any lat lon restrictions etc since point source data

graph = 977; %new 1D PDF from 2D histo data - can choose either axis
%(for watervap)

axis1D = 'y';

logbin_norm = logbin_norm_driver;
i_plot_norm=i_plot_norm_driver;
i_div_bin_widths=i_div_bin_widths_driver;
pdf_type = pdf_type_driver;
LAT_val = LAT_val_ALL; LON_val = LON_val_ALL; %Still need to set this, though, even though not doing lat lon screening

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

clims = get(gca,'clim');

if iset_min_clim==1
    clims(1) = clim_min;
end
if iset_min_clim==1
    clims(2) = clim_max;
end

caxis(clims);
        
clim_min=0;
iset_max_clim=1;
clim_max=200;

%close(gcf);
%waterVapourMay2005
%close(gcf);

                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
