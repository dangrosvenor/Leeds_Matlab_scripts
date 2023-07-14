x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

[X,Y] = MESHGRID(x_grid,y_grid);

ih_wrf=4;
time=11;

u=0.5 * ( nc{'U'}(time,ih_wrf,:,1:end-1) + nc{'U'}(time,ih_wrf,:,2:end) );
v=0.5 * ( nc{'V'}(time,ih_wrf,1:end-1,:) + nc{'V'}(time,ih_wrf,2:end,:) );
wind2d = sqrt(u.^2+v.^2);

%%%%%%%%%%%% find the crest of the mountain
hgt=nc{'HGT'}(1,:);

np=80; %number of points along the line of longitude to make averages for
np=320;
%         latA=-70;
%         latB=-64.8;
latA=-69;
latB=-67.1;
LATS=latA:(latB-latA)/np:latB;



for imean=1:length(LATS)
    icons_inds = get_inds_constant_lat(LATS(imean),lat2d,lon2d); %get indices for a constant latitude slice
    ieast=find(lon2d.var(icons_inds)>-67.5); %to remove the mountains to the west of peninsula from being included
    [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain
    lon_peak(imean)=lon2d.var(icons_inds(ieast(ipeak)));
    [ilat_peak(imean),ilon_peak(imean)] = ind2sub(size(lat2d.var),icons_inds(ieast(ipeak))); %get lat and lon indices for peak
    peak_vs_lat(imean)=peak_height;

    [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain
    lee_height(imean) = peak_height*0.75;

    ilee=[0 0];
    acc=20/100;
    acc_step=5/100;
    while length(ilee)~=1
        ilee = find( hgt(icons_inds(ieast)) > lee_height(imean)*(1-acc) & hgt(icons_inds(ieast)) < lee_height(imean)*(1+acc) );
        if length(ilee)==0
            acc=acc+acc_step; %reset to the previous value
            acc_step=acc_step/2; %reduce the step
        else
            acc=acc-acc_step;
        end
    end

    lon_lee(imean)=lon2d.var(icons_inds(ieast(ilee)));
    [ilat_lee(imean),ilon_lee(imean)] = ind2sub(size(lat2d.var),icons_inds(ieast(ilee))); %get lat and lon indices for peak
    latlon_inds(imean) = icons_inds(ieast(ilee));

    [lee_max(imean) ilee_max]=max(wind2d(icons_inds(ieast)));
    lon_lee_max(imean)=lon2d.var(icons_inds(ieast(ilee_max)));
    lee_max_h(imean)=hgt(icons_inds(ieast(ilee_max)));
end

