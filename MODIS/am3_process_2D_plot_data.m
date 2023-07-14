%produce data for 2D colour plots from AM3 data

am3_2Dplot_case = 'P vs day';
am3_2Dplot_case = 'diurnal';
%am3_2Dplot_case = 'timeseries';

%plotting for one location only for the year
lat_am3 = -20;
lon_am3 = -75;
%find the nearest point to this in the am3 cubed-sphere grid
[ilat,ilon,dmins]=am3_map_lat_lon_points(am3_lat,am3_lon,lat_am3,lon_am3);

switch am3_2Dplot_case
    
    case 'P vs day'
        %think there is a 3-4 hour time difference for LT in Chile e.g. LT = UTC-4
        times_of_day = [0 3 6 9 12 15 18 21]; %UTC
%        times_of_day = [15 21 0];    
        
        iday_am32D = [];
        clear iday_am32D_diurnal
        for iam3=1:length(times_of_day)
            iday_am32D = [iday_am32D; find(am3_time_UTC==times_of_day(iam3))];
        end
        
        iday_am32D = sort(iday_am32D);              
        
        %shouldn't matter that the P profiles are different for the
        %different days since can do a 2D field of P and ndays
        
        Y = am3_pfull(iday_am32D,:,ilat,ilon)/100; %mb
        X = repmat( am3_decimal_days(iday_am32D) , [1 size(am3_pfull,2)] );
        
%        P=am3_drop2(iday_am32D,:,ilat,ilon);
        
%        P=am3_cf(iday_am32D,:,ilat,ilon);
        P = am3_liq2(iday_am32D,:,ilat,ilon);
        
        figure
        pcolor(X,Y,P);
        shading flat
        set(gca,'ylim',[450 1080]);
        set(gca,'ydir','reverse');
        colorbar
        
    case 'diurnal'
        times_of_day = [0 3 6 9 12 15 18 21];
%        times_of_day = [15 21 0];    
        
        iday_am32D = [];
        clear iday_am32D_diurnal
        for iam3=1:length(times_of_day)
            iday_am32D = [iday_am32D; find(am3_time_UTC==times_of_day(iam3))];
            iday_am32D_diurnal(iam3).inds = [find(am3_time_UTC==times_of_day(iam3))];
        end
        
        iday_am32D = sort(iday_am32D);
             
        %shouldn't matter that the P profiles are different for the
        %different days since can do a 2D field of P and ndays
        
%        Y = am3_pfull(iday_am32D,:,ilat,ilon)/100; %mb
%        X = repmat( am3_decimal_days(iday_am32D) , [1 size(am3_pfull,2)] );
        
clear timser
        for it=1:length(iday_am32D_diurnal)
            inds=iday_am32D_diurnal(it).inds;
%            timser(it) = meanNoNan(am3_Nd_meanliq(inds,ilat,ilon),1);
            timser(it) = meanNoNan( am3_CTH(inds,ilat,ilon) - am3_CBH(inds,ilat,ilon) , 1);            
            
        end

        
%        figure
        hold on

 
 plot(times_of_day,timser,'bo-');
 
    case 'timeseries'
        times_of_day = [0 3 6 9 12 15 18 21];
%        times_of_day = [15 21 0];    
        
        iday_am32D = [];
        clear iday_am32D_diurnal
        for iam3=1:length(times_of_day)
            iday_am32D = [iday_am32D; find(am3_time_UTC==times_of_day(iam3))];
            iday_am32D_diurnal(iam3).inds = [find(am3_time_UTC==times_of_day(iam3))];
        end
        
        iday_am32D = sort(iday_am32D);
        
      
        
        diurnal_case = 'Cloud top/base pressures';
%        diurnal_case = 'Cloud top/base heights';
        switch diurnal_case
            case 'Cloud top/base pressures'
        timser = 0.01*am3_CTP(iday_am32D,ilat,ilon);
        timser2 = 0.01*am3_CBP(iday_am32D,ilat,ilon);
            case 'Cloud top/base heights'
        
        timser = am3_CTH(iday_am32D,ilat,ilon);
        timser2 = am3_CBH(iday_am32D,ilat,ilon);
        
        end
               
        time_am3 = am3_decimal_days(iday_am32D);
        figure
hold on
        plot(time_am3,timser,'bo');
        hold on
        plot(time_am3,timser2,'rd');
        
        switch diurnal_case
            case 'Cloud top/base pressures'
        set(gca,'ylim',[450 1080]);
        set(gca,'ydir','reverse');
        end
        
        
end
        