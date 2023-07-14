% lon_region = [66 90];
% lat_region=[16 30];
% 
% iregion=find(gcm_Plat2D_SMOS>lat_region(1) & gcm_Plat2D_SMOS<lat_region(2) & gcm_Plon2D_SMOS>lon_region(1) & gcm_Plon2D_SMOS>lon_region(1));
% 
% [ilat,ilon]=ind2sub(size(gcm_Plat2D_SMOS),iregion);
% 
% 
% %figure
% timser=NaN*ones([size(Soil_Moisture.timeseries3,3) length(iregion)]);
% for ireg=1:length(iregion)
%     timser(:,ireg)=squeeze(Soil_Moisture.timeseries3(ilat(ireg),ilon(ireg),:));
% %    plot(timser);
% %    hold on
% end

mon_thresh=0.4;
Nwindow=15;


ifind_lats=0;  %set to one if have requested new lat values


lon_region = [66:1:90];
lat_region=[16:1:30];

lon_region = [75:1:88];
lat_region=[16:1:28];

if ifind_lats==1

    i=0;
    clear ilat2 ilon2
    for ilat=1:length(lat_region)
        for ilon=1:length(lon_region);
            i=i+1;
            [ilat2(i),ilon2(i)]=getind_latlon_quick(gcm_Plat2D_SMOS,gcm_Plon2D_SMOS,lat_region(ilat),lon_region(ilon),0.1);

        end
    end

end

%iregion=find(gcm_Plat2D_SMOS>lat_region(1) & gcm_Plat2D_SMOS<lat_region(2) & gcm_Plon2D_SMOS>lon_region(1) & gcm_Plon2D_SMOS>lon_region(1));

%[ilat,ilon]=ind2sub(size(gcm_Plat2D_SMOS),iregion);


ntime = size(Soil_Moisture.timeseries3,3);
%figure
timser_comp=NaN*ones([ntime*2+1 length(iregion)]);
timser=NaN*ones([ntime length(iregion)]);
icomp=0;
for ireg=1:length(ilat2)
    timser(:,ireg)=squeeze(Soil_Moisture.timeseries3(ilat2(ireg),ilon2(ireg),:));
    x = [1:ntime];
    y = timser(:,ireg);
    inan =find(isnan(y)==1);
    x(inan)=[]; y(inan)=[];
    
    if length(x)>Nwindow
        [xtim,ytim]=window_average(x,y,Nwindow,'mean');
        imon = find(ytim>mon_thresh);
    else
        imon=[];
    end
    
    if length(imon)>0
        
        icomp=icomp+1;
    
    %Index in timser where monsoon occurs
    imon2 = ceil(xtim(imon(1)));
    

    
    %Aiming to put the imon2 time at ntime in the double sized array
    %So if imon2=ntime (end of array) then want to put end of the array at
    %ntime(=imon2) and therefore the start of the array at ntime-ntime+1=1
    %If is end-1 then would put the end of the array at ntime+1 and so
    %start at ntime+1-ntime+1 = 2. So istart=ntime-mon2+1
    istart = ntime - imon2 + 1;
    
    timser_comp(istart:istart+ntime-1,icomp) = timser(:,ireg);
    
    end
    
%    plot(timser(:,ireg));
%    hold on
end


timser_meancomp = meanNoNan(timser_comp(:,1:icomp),2);
%make so that x(ntime)=0;
x = [-ntime+1:ntime+1];
figure; plot(x,timser_meancomp); grid
icomp
