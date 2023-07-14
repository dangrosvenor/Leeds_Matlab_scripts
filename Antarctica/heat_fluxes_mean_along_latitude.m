%plots bar graph of surface heat/melting budget

isave=0;
%iloc=1;

imethod='mean_at_lat';
%imethod='single_points';

limit_to_melting='yes'; %if yes only include points where there was melting for all of the times considered
limit_to_melting='no'; %makes little difference for the ECMWF case

scrsz=get(0,'ScreenSize');
%posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.13];
%posit=[9 50 scrsz(3)/1.9 scrsz(4)/2.1];  %used 22nd Jan, 2009
posit=[9 50 scrsz(3)/1.4 scrsz(4)/1.6];

nlat = size(lat2d.var,1);
nlon = size(lat2d.var,2);
 dx_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));     
 dy_grid = distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));
 
 i_grid = dx_grid * [1:nlon];
 j_grid = dy_grid * [1:nlat];

        extra_x = [550 600 800]; %south, and top of Larsen C and Larsen B
        extra_y = [150 300 400];

        clear LAT LON
        for iflt=1:length(extra_x)
            i_extra = findheight_nearest(i_grid,extra_x(iflt));
            j_extra = findheight_nearest(j_grid,extra_y(iflt));
            LAT(iflt) = lat2d.var(j_extra,i_extra);
            LON(iflt) = lon2d.var(j_extra,i_extra);
        end
        

        switch imethod
            case 'mean_at_lat'
                np=80; %number of points along the line of longitude to make averages for
                latA=-70;
                latB=-64.8;
%                latA=-69;
%                latB=-67.1;
                LAT=latA:(latB-latA)/np:latB;

            case 'single_points'
                np=40; %number of points along the line of longitude to make averages for
                latA=-70;
                latB=-64.8;
                LAT=latA:(latB-latA)/np:latB;
                LON=ones([1 length(LAT)])*-61;
                [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);
        end
        
        clear lonmin_dat lonmax_dat dist_dat sw_dat lw_dat sh_dat lh_dat grd_dat cond_dat sp10_dat i2_save rh_dat rh_iz_dat
        clear n_tot_dat
        

switch(imethod)
    case 'mean_at_lat'



        landmask=nc{'LANDMASK'}(1,:);
        seaice=nc{'SEAICE'}(1,:); %NOTE seaice is zero everywhere for NCEP as it wasn't put in properly on this run!
        hgt=nc{'HGT'}(1,:); %NOTE seaice is zero everywhere for NCEP as it wasn't put in properly on this run!
        %think have a run with seaice included
        
        for imean=1:length(LAT)            
            icons_inds = get_inds_constant_lat(LAT(imean),lat2d,lon2d); %get indices for a constant latitude slice
            
            ieast=find(lon2d.var(icons_inds)>-67.5); %to remove the mountains to the west of peninsula from being included
            [peak_height ipeak]=max(hgt(icons_inds(ieast))); %find the position of the peninsula mountain and keep to the east of it
            lon_peak=lon2d.var(icons_inds(ieast(ipeak)));
            peak_vs_lat(imean)=peak_height;
            
            switch limit_to_melting
                case 'no'
                    i2=find( abs(hgt(icons_inds)-0)<1e-5 & abs(landmask(icons_inds)-1)<1e-5 & abs(seaice(icons_inds)-0)<1e-5...
                        & lon2d(1).var(icons_inds)>lon_peak ); %
                case 'yes'
                    i2=find( abs(hgt(icons_inds)-0)<1e-5 & abs(landmask(icons_inds)-1)<1e-5 & abs(seaice(icons_inds)-0)<1e-5...
                        & lon2d(1).var(icons_inds)>lon_peak & n_tot2(icons_inds)==0 ); %only including points where melting is occuring at all the times requested
                    
                    i2_rejected=find( abs(hgt(icons_inds)-0)<1e-5 & abs(landmask(icons_inds)-1)<1e-5 & abs(seaice(icons_inds)-0)<1e-5...
                        & lon2d(1).var(icons_inds)>lon_peak & n_tot2(icons_inds)~=0 ); %only including points where melting is occuring at all the times requested
                    

            end
            
            if length(i2)~=0
                lonmin_dat(imean) = min(lon2d.var(icons_inds(i2)));
                lonmax_dat(imean) = max(lon2d.var(icons_inds(i2)));
                dist_dat(imean) = distlatlon(LAT(imean),lonmin_dat(imean),LAT(imean),lonmax_dat(imean));
                
%                max_melt_lat(imean) = max(max_melt(icons_inds(i2)));

%                rat=melt_tot_nighttime./melt_tot_daytime;
%                max_melt_lat(imean) = max(rat(icons_inds(i2)));
                switch limit_to_melting
                    case 'yes'
                        i2_save(1:length(i2_rejected),imean)=icons_inds(i2_rejected);    %store the indices that were used
                end
            else
                lonmin_dat(imean) = NaN;
                lonmax_dat(imean) = NaN;
                dist_dat(imean) = 0;

            end
            
            melt_dat(imean)=mean(melt_tot(icons_inds(i2)));
            sw_dat(imean)=mean(sw_tot(icons_inds(i2)));
            lw_dat(imean)=mean(lw_tot(icons_inds(i2)));
            sh_dat(imean)=mean(sh_tot(icons_inds(i2)));
            lh_dat(imean)=mean(lh_tot(icons_inds(i2)));
            grd_dat(imean)=mean(grd_tot(icons_inds(i2)));
            temp_dat(imean)=mean(grd_tot(icons_inds(i2)));
            sp10_dat(imean)=mean(grd_tot(icons_inds(i2)));
            cond_dat(imean)=mean(cond_tot(icons_inds(i2))); 
            rh_dat(imean)=mean(rh_tot(icons_inds(i2))); 
            rh_iz_dat(imean)=mean(rh_iz_tot(icons_inds(i2))); 
            n_tot_dat(imean)=mean(n_tot(icons_inds(i2))); 
        end

            
        
    case 'single_points';

        for iloc=1:length(LAT)            
            melt_dat(iloc)=melt_tot(ilat(iloc),ilon(iloc));
            sw_dat(iloc)=sw_tot(ilat(iloc),ilon(iloc));
            lw_dat(iloc)=lw_tot(ilat(iloc),ilon(iloc));
            sh_dat(iloc)=sh_tot(ilat(iloc),ilon(iloc));
            lh_dat(iloc)=lh_tot(ilat(iloc),ilon(iloc));
            grd_dat(iloc)=grd_tot(ilat(iloc),ilon(iloc));
            
            if (abs(hgt(ilat(iloc),ilon(iloc))-0)<1e-5 & abs(landmask(ilat(iloc),ilon(iloc))-1)<1e-5...
                    & abs(seaice(ilat(iloc),ilon(iloc))-0)<1e-5)
                shelf_dat(iloc)=1; %flag whether is actually on the shelf
            else
                shelf_dat(iloc)=0;
            end
        end
        


end


for iloc=1:0    %length(LAT)
    switch imethod
        case 'mean_at_lat'
            figname = ['Mean contribution to melt energy at LAT=' num2str(LAT(iloc)) ' for ' filestr];
        case 'single_points'
            figname = ['Contribution to melt energy at LAT=' num2str(LAT(iloc)) 'LON=' num2str(LON(iloc)) ' for ' filestr];
    end
    titname = figname;
    savename = figname;
    hf=figure('name',figname,'Position',posit);

    ydat=[melt_dat(iloc) sw_dat(iloc) lw_dat(iloc) sh_dat(iloc)...
        lh_dat(iloc) sh_dat(iloc)+lh_dat(iloc) grd_dat(iloc)];

    ydat=[ydat; zeros([1 length(ydat)])]; %so that it draws the bars in different colours

    h=bar(ydat);

    legend('Melt','SW','LW','SH','LH','SH+LH','GRD');

    set(gca,'fontsize',20);

    set(gca,'xlim',[0.6 1.6]);
%    set(gca,'ylim',[minALL(ydat) maxALL(ydat)]);
    set(gca,'ylim',[-4 10]);
    %set(gca,'xscale','log');
    %set(gca,'yscale','log');


    set(gca,'xticklabel',''); %remove x-axis labels

    %xlabel('Critical Supersaturation (%)');
    ylabel('Melt rate (mm day^{-1})');
    
    title(titname,'fontsize',16);

    if isave==1
        multisaveplot
    end
    
end

 shlh_dat = sh_dat + lh_dat;
 iref = findheight_nearest(LAT,-65.58); %Larsen B point (approx centre of Larsen B)
 iref2 = findheight_nearest(LAT,-69.00); %Larsen B point (approx centre of Larsen B)
 
 
 LAT_ridge=LAT;
 fprintf(1,'\nDone heat_fluxes_mean_along_latitude\n');
 
% clear sw_cont lw_cont sh_cont lh_cont shlh_cont grd_cont
% clear sw_cont2 lw_cont2 sh_cont2 lh_cont2 shlh_cont2 grd_cont2
% 
% iref=3;
% 
% ii=1:length(sw_dat);
% ii(iref)=[];
% for iind=1:length(ii)
%     iloc=ii(iind);
%     sw_cont(iind)=(sw_dat(iloc)-sw_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     lw_cont(iind)=(lw_dat(iloc)-lw_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     sh_cont(iind)=(sh_dat(iloc)-sh_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     lh_cont(iind)=(lh_dat(iloc)-lh_dat(iref))/(melt_dat(iloc)-melt_dat(iref));    
%     shlh_cont(iind)=( (sh_dat(iloc)+lh_dat(iloc)) - ( sh_dat(iref)+lh_dat(iref)) )/(melt_dat(iloc)-melt_dat(iref));
%     grd_cont(iind)=(grd_dat(iloc)-grd_dat(iref))/(melt_dat(iloc)-melt_dat(iref));    
% end
% 
% iref=2;
% 
% ii=1;
% for iind=1:length(ii)
%     iloc=ii(iind);
%     sw_cont(3)=(sw_dat(iloc)-sw_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     lw_cont(3)=(lw_dat(iloc)-lw_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     sh_cont(3)=(sh_dat(iloc)-sh_dat(iref))/(melt_dat(iloc)-melt_dat(iref));
%     lh_cont(3)=(lh_dat(iloc)-lh_dat(iref))/(melt_dat(iloc)-melt_dat(iref));    
%     shlh_cont(3)=( (sh_dat(iloc)+lh_dat(iloc)) - ( sh_dat(iref)+lh_dat(iref)) )/(melt_dat(iloc)-melt_dat(iref));
%     grd_cont(3)=(grd_dat(iloc)-grd_dat(iref))/(melt_dat(iloc)-melt_dat(iref));    
% end
% 
% iplot_rel=0;
% if iplot_rel==1
% 
%     figname = ['Relative contribution to melt energy difference between LAT=' num2str(LAT(1)) ' and LAT=' num2str(LAT(3)) ' for ' filestr];
%     titname = figname;
%     savename = figname;
%     hf=figure('name',figname,'Position',posit);
% 
%     ydat=100*[sw_cont(1) lw_cont(1) sh_cont(1)...
%         lh_cont(1) shlh_cont(1) grd_cont(1)];
% 
%     ydat=[ydat; zeros([1 length(ydat)])]; %so that it draws the bars in different colours
% 
%     h=bar(ydat);
% 
%     legend('SW','LW','SH','LH','SH+LH','GRD');
% 
%     set(gca,'fontsize',20);
% 
%     set(gca,'xlim',[0.6 1.6]);
% %    set(gca,'ylim',[minALL(ydat) maxALL(ydat)]);
% %    set(gca,'ylim',[-4 10]);
%  
%     %set(gca,'xscale','log');
%     %set(gca,'yscale','log');
% 
% 
%     set(gca,'xticklabel',''); %remove x-axis labels
% 
%     %xlabel('Critical Supersaturation (%)');
%     ylabel('Relative contribution (%)');
%     
%     title(titname,'fontsize',16);
%     
% %    multisaveplot %make sure that neither 'timh' nor 'prof' are selected for graph type in multisaveplot
% 
% end
% 
