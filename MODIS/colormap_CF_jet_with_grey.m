%colormap jet with grayscale added onto the end

%first set up a generic jet colormap with the middle (white bit) cut out
%out
Njet=1280;
JET=jet(Njet);
%first half up to the bit we want cut out
jetA=JET(1:round(Njet*43/128),:);
%second bit after the bit to be cut out (and have also cut out the dark bit
%at the end)
jetB=JET(Njet/2:end-round(Njet*2/128),:);
jetNEW = [jetA; jetB];

Njet_NEW = size(jetNEW,1);

%the size of the final colorbar - making it equal to the size of the
%colorbar ticks (minus one as need one less coluor block than tick mark)
Ncb_tot = length(x_cbar_vals)-1;

% ctick_range_choose.m does this below or similar - i.e. cbar_vals are the
% values represented on the plot and x_cbar_vals are just a linear spacing
% of these along the colorbar (between 0 and 1)
%%% cbar_vals = [0:20:300 400 500 750 1000 1e12 2e20 4e20 6e20 8e20];
%%% x_cbar_vals=linspace(0,1,length(cbar_vals));

%i_cbar_val_new_scale_start is the position of the grey part of the
%colorscale (on a 0 to 1 scale)
Njet=i_cbar_val_new_scale_start-1; %length of the jet part

%interpolate a new colormap of length Njet from jetNEW
jetNEW_2 = interp1([0:Njet_NEW-1]/(Njet_NEW-1),jetNEW,[0:Njet-1]/(Njet-1));

%now do the grey (gray) colorbar
grey0=gray(1280);
%cut out the white region
grey0=grey0(1:ceil(1280*20/30),:);
sgrey0=size(grey0,1);

Ngray = Ncb_tot - Njet;

grey = interp1([0:sgrey0-1]/(sgrey0-1),grey0,[0:Ngray-1]/(Ngray-1));
%grey=flipdim(grey,1);
%combine the colormaps with the grey one flipped
cf_colormap = [jetNEW_2; grey];

