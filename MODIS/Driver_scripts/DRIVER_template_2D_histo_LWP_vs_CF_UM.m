%Stripped down template for running plotTime* to produce a joint 1D PDF
% N.B. - this does no lat,lon screening

%% Some flags. May want to set these outside the script :-
% These aren't for the 2d PDFS I think - are they needed at all?
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;
pdf_type_driver='cumsum';
LAT_val_DRIVER = [-1e9 1e9]; LON_val_DRIVER = [-1e9 1e9]; %Still need to set this even though not doing lat lon screening

%% Define the bins
%Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;
%

dcf=0.1;
%dcf=0.2;
%dcf=0.01;

%N.B.- bins will have a small amount removed and added to the ends, so if
%put zero as the first bin then will include everything > -1e-9 or similar
Xbins_DRIVER = [-dcf/2:dcf:1+dcf/2]; ichoose_Xbins_DRIVER=1; %cloud fraction
Xbins_DRIVER = [0:dcf:1]; ichoose_Xbins_DRIVER=1; %cloud fraction
%Xbins_DRIVER = [0 0.05:dcf:0.95 1]; ichoose_Xbins_DRIVER=1; %cloud fraction - as for CPT VOCALS plot

switch LWP_type
    case 'Grid-box mean'
        %Ybins_DRIVER = [5:20:500 1e4]; ichoose_Ybins_DRIVER=1; %LWP in-cloud
        Ybins_DRIVER = [-20 10.^[log10(1):0.05:log10(1e4)]]; ichoose_Ybins_DRIVER=1; %LWP grid box mean
        Ybins_DRIVER = [10.^[log10(0.9):0.05:log10(1e4)]];
        Ybins_DRIVER = [10.^[log10(20):0.05:log10(1e4)]];
        %Ybins_DRIVER = [-20 10:10:200]; ichoose_Ybins_DRIVER=1; %LWP grid box mean
    case 'In-cloud'
        Ybins_DRIVER = [-20 10.^[log10(1):0.05:log10(1e4)]]; ichoose_Ybins_DRIVER=1; %LWP grid box mean
        Ybins_DRIVER = [10.^[log10(0.9):0.05:log10(1e4)]];
end


min_cf=0.001;
% Xbins_DRIVER = [-dcf/2 min_cf dcf/2:dcf:1+dcf/2]; ichoose_Xbins_DRIVER=1; %cloud fraction PI
% Ybins_DRIVER = [-dcf/2 min_cf dcf/2:dcf:1+dcf/2]; ichoose_Ybins_DRIVER=1; %cloud fraction PD



%% Define the data somewhere
% Y_driver = magic(10); %the data
% X_driver = magic(10)'; %the data
% Z_driver = ones([10 10]); %[1:100];
% ylabelstr = 'Label for y data';
% xlabelstr = 'Label for x data';

%% Set the defaults to make sure all variables are set
watervap_defaults
pdf2D_defaults


%ndims_hist=3; %set to 3D and put the data to average on the 3rd dimension (as Z)

i577 = 'MODIS_plot_UW';

iset_min_clim=0;
clim_min=1e-5;
iset_max_clim=0;
clim_max=200;

iminovr=0;
mincovOvr=1e-4;
mincovOvr=1e-3;
imaxovr=0;
maxcovOvr=1;

logflag=0; %set iminovr=1 and mincovOvr to a min value for this. Also, need to figure out how to get the brown-blue
%colour map to work.
dlogflag=0;

iuse_overall_norm=1;

icolmap=1;
cmap=lbmap(256,'brownblue'); %nice colormap for colorblind people

isave_plot=0;
%savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

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

iplot_mean_XY='y';

% --------- Override flags for watervap --------
man_choose_water_graph=1;    %for watervap

%---  Run plot script and save
plotTimeHeightVap3
%close(gcf);
%waterVapourMay2005
%close(gcf);

%Doesn't seem to work for these plots :-
%increase_font_size_map_figures  
'';

                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
