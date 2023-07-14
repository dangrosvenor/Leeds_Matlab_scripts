% CF>80, SZA<65, no CTH or tau screehing :-
load_dir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/';
%As above, but with CTH>1km & < 3.2km, meaTau>10 :-
load_dir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_1_3.2km_SZA_65_meanTau10/';

iload_monthly_data=1;

if iload_monthly_data==1
    
    %No landmask
%     loadfile = [load_dir 'Nd_saved_sea-ice-edge_Jan2017_no_landmask_20170111T040749.mat'];
%     Nd_mask0=load(loadfile);
%     Nd_mask0.tag='no landmask';
    
    %AMSRE landmask, no smoothing
    % CF>80, SZA<65, no CTH or tau screehing :-
    %loadfile = [load_dir 'Nd_saved_sea-ice-edge_Jan2017_landmask_no_smooth_20170111T040015.mat'];
    %As above, but with CTH>1km & < 3.2km, meaTau>10 :-
    loadfile = [load_dir 'Nd_saved_sea-ice-edge_Jan2017_landmask_no_smooth_20170116T040211.mat'];
    Nd_mask1=load(loadfile);
    Nd_mask1.tag='non-smoothed landmask';
    
%     %AMSRE landmask, smoothing 2x2
%     loadfile = [load_dir 'Nd_saved_sea-ice-edge_Jan2017_landmask_smooth02_20170111T035201.mat'];
%     Nd_mask2=load(loadfile);
%     Nd_mask2.tag='landmask smoothed with N=2';
%     
%     %AMSRE landmask, smoothing 4x4
%     loadfile = [load_dir 'Nd_saved_sea-ice-edge_Jan2017_landmask_smooth04_20170111T064813.mat'];
%     Nd_mask4=load(loadfile);
%     Nd_mask4.tag='landmask smoothed with N=4';
end

months={'Dec' 'Jan' 'Feb'}; %Data is for Sep thorugh Apr plotting Dec, Jan and feb

for im=1:length(months)           
    %figure
    %eval(['pcolor(MLON,MLAT,Nd_mask0.Nd_37_SH_summer_2005_to_2011_' months{im} ');'])    
    %title([months{im} ', no landmask']); caxis([0 200]); shading flat
    
%     eval(['qpcolor(Nd_mask0.Nd_37_SH_summer_2005_to_2011_' months{im} ');'])
%     title([months{im} ', no landmask']); caxis([0 200]);
    
    eval(['qpcolor(Nd_mask1.Nd_37_SH_summer_2005_to_2011_' months{im} ');'])
    title([months{im} ', non-smoothed landmask']); caxis([0 200]);
    
%     eval(['qpcolor(Nd_mask2.Nd_37_SH_summer_2005_to_2011_' months{im} ');'])
%     title([months{im} ', landmask smoothed with N=2']); caxis([0 200]);
%     
%     eval(['qpcolor(Nd_mask4.Nd_37_SH_summer_2005_to_2011_' months{im} ');'])
%     title([months{im} ', landmask smoothed with N=4']); caxis([0 200]);    
end