%Tom LC paper digitised size distributions

idat=0;


%7th Dec, 1995 - 18:35
idat=idat+1;
scale = 0.8/22; % number per mm on the graph
crystals(idat).mm=[2 17 18.5 16 11 6.75 7.5 11 8.25 8 6 4 4 3 2];
crystals(idat).num=crystals(idat).mm*scale;
size_bins(idat).microns=[1:length(crystals(idat).num)]*20;