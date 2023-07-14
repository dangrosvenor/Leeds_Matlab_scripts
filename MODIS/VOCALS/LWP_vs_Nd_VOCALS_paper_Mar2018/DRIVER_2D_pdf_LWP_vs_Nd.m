%Stripped down template for running plotTime* to produce a joint 1D PDF
% N.B. - this does no lat,lon screening

X_driver_save = X_driver;
Y_driver_save = Y_driver;

%% Some flags. May want to set these outside the script :-
logbin_norm_driver=0;
i_plot_norm_driver=0;
i_div_bin_widths_driver=1;
pdf_type_driver='cumsum';
LAT_val_ALL = [-1e9 1e9]; LON_val_ALL = [-1e9 1e9]; %Still need to set this even though not doing lat lon screening
LAT_val_ALL = [-20 -10]; LON_val_ALL = [-90 -80]; %Frida's region

LAT_val_ALL = [-30 -20]; LON_val_ALL = [-80 -70]; %contains some land - screen this out, or jsut move across?
LAT_val_ALL = [-30 -20]; LON_val_ALL = [-83 -73]; %contains some land - screen this out, or jsut move across?
LAT_val_ALL = [-30 -20]; LON_val_ALL = [-93 -83]; %
LAT_val_ALL = [-30 -20]; LON_val_ALL = [-103 -93]; %
LAT_val_ALL = [-30 -20]; LON_val_ALL = [-113 -103]; %
LAT_val_ALL = [-30 -20]; LON_val_ALL = [-123 -113]; %

%% Define the bins
dx=1; Xbins_DRIVER = [-0.01 dx:dx:500]; ichoose_Xbins_DRIVER=1;

%Ybins_DRIVER = [-39:2:39]; ichoose_Ybins_DRIVER=1;
%Ybins_DRIVER = [0:2:0]; ichoose_Ybins_DRIVER=0;
%Ybins_DRIVER = [-870:2:850]; ichoose_Ybins_DRIVER=1;
Ybins_DRIVER = [-0.01 10.^[log10(10):0.02:log10(350)]]; ichoose_Ybins_DRIVER=1;   
dy=1; Ybins_DRIVER = [-0.01 dy:dy:350]; ichoose_Ybins_DRIVER=1;   

%% Define the data somewhere
% Y_driver = magic(10); %the data
% X_driver = magic(10)'; %the data
%xlabelstr='Low cloud fraction';
%ylabelstr='SW surface forcing (W m^{-2})';

%% Set the defaults to make sure all variables are set
watervap_defaults
pdf2D_defaults

i577 = 'MODIS_plot_UW';

iset_min_clim=1;
clim_min=1e-5;
iset_max_clim=0;
clim_max=200; %(in log form if using logflag=1)

logflag=0; %sets the colorbar to log values ***make sure to set a clim_min if usig logflag=1*** 
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

ilat = find(MLAT>=LAT_val(1) & MLAT<=LAT_val(2));
ilon = find(MLON>=LON_val(1) & MLON<=LON_val(2));
X_driver2 = NaN*ones(size(X_driver));
Y_driver2 = NaN*ones(size(Y_driver));
X_driver2(ilat,ilon) = X_driver(ilat,ilon);
Y_driver2(ilat,ilon) = Y_driver(ilat,ilon);

X_driver = X_driver2;
Y_driver = Y_driver2;


% --------- Override flags for 2D PDF --------
ioverride_pdf=1;
%iocean_only=1;
man_choose_plotTimeHeight_graph=1;
ioverride_location_selection=1;
ioverride_pdf_varchoose = 1;

% --------- Override flags for watervap --------
man_choose_water_graph=1;    %for watervap


if logflag==1
    iminovr=1;
    mincovOvr = clim_min;
end

%---  Run plot script and save
plotTimeHeightVap3

%corr_coeffXY is the correlation coefficient

clims = get(gca,'clim');

if logflag==0       
    if iset_min_clim==1
        clims(1) = clim_min;
    end
end

if iset_max_clim==1
    clims(2) = clim_max;
end

caxis(clims);




        
clim_min=0;
iset_max_clim=1;
clim_max=200;

%close(gcf);
%waterVapourMay2005
%close(gcf);

                            
X_driver = X_driver_save;
Y_driver = Y_driver_save;
                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
