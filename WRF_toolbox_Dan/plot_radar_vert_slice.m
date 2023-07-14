fileName=['Y:\WRF\30thNov_Min\d03'];
imagedir=['Y:\WRF\30thNov_Min\radar_plots\vert_slices'];

warning off;
%fileName='/ddisk1/mccikpc2/wrf/hpcx/run002_noradcloud_300cm3_20060206/wrfout_d02_2006-02-05_18:00:00'
%fileName='/ddisk1/zhu/hpcx/20060206/sonde_srtm_1800_cloud_zhu_myj_godrad_10CCN/wrfout_d01_2006-02-05_18:00:00'
%fileName='/scratch/a/zhu/wrfout_d01_2006-02-05_18:00:00'
%imagedir='~/temp/modelimages1'
%nc=netcdf(fileName);
%varl=WRFUserARW(nc,'XLONG',1);
%varla=WRFUserARW(nc,'XLAT',1);

varl=lon2d;
varla=lat2d;

%hf=figure;
%m_proj('mercator','long',[min(min(varl.var)) max(max(varl.var))],'lat',[min(min(varla.var)) max(max(varla.var))]);
%[X,Y]=m_ll2xy(varl.var,varla.var);


Times=nc{'Times'}(:);
[r,c]=size(Times);
%levels=[2500,5000,14000,16000,17400];

lon_slice=131.3; %position of slice in longitude
lon_slice=131.0882;
lon_slice=130.2766;
lat_slice=-67.3;

var='Z';
var='tot';
var='vap';
var='w';

%lat_slice=[]; %position at latitude - leave as [] if want longitude slice
lon_slice=[];

if length(lon_slice)==0
    ilon_slice=0;
    ilat_slice=1;
else
    ilon_slice=1;
    ilat_slice=0;
end

iread=1;

nlat=size(lat2d.var,1);
nlon=size(lat2d.var,2);

%levels=[5000];
for time=16:1:16  %23:1:r   %20 = 9:30
    
    clear z_slice Z_slice
    
  %  if iread==1
        if ilon_slice==1
            
            [ilat2,ilon2] = getind_latlon_quick(lat2d.var,lon2d.var,lat2d.var(1,1),lon_slice,0.1);
            
            for ilat2=1:nlat                
                z_slice(ilat2,:)=WRFUserARW(nc,'Z',time,ilat2,ilon2)';
                pot_slice(ilat2,:)=WRFUserARW(nc,'th',time,ilat2,ilon2)';
                
                switch var
                case 'Z'
                    Z=WRFRadarRefl_smallMEM_2(nc,time,'thompson',ilat2,ilon2); 

                    Z=real(10.*log10(Z)); 
                    Z(find(Z(:)<0))=NaN;
                case 'tot'
                    Z=nc{'QICE'}(time,:,ilat2,ilon2)+nc{'QSNOW'}(time,:,ilat2,ilon2)+nc{'QGRAUP'}(time,:,ilat2,ilon2)+nc{'QVAPOR'}(time,:,ilat2,ilon2); 
                case 'vap'
                    Z=nc{'QVAPOR'}(time,:,ilat2,ilon2);    
                case 'w'
		    Z=nc{'W'}(time,:,ilat2,ilon2);
		end
                
                Z_slice(ilat2,:)=Z';
            end
            
            
        else
            [ilat2,ilon2] = getind_latlon_quick(lat2d.var,lon2d.var,lat_slice,lon2d.var(1,1),0.1);
            
            for ilon2=1:nlon                
                z_slice(ilon2,:)=WRFUserARW(nc,'Z',time,ilat2,ilon2)';
                pot_slice(ilon2,:)=WRFUserARW(nc,'th',time,ilat2,ilon2)';
                
                switch var
                case 'Z'
                    Z=WRFRadarRefl_smallMEM_2(nc,time,'thompson',ilat2,ilon2); %are approximating here since the level of
                    %the desired height will vary across the domain                                                                         %Here are just plotting the same level everywhere
                    Z=real(10.*log10(Z)); 
                    Z(find(Z(:)<0))=NaN;
                    
                case 'tot'
                    Z=nc{'QICE'}(time,:,ilat2,ilon2)+nc{'QSNOW'}(time,:,ilat2,ilon2)+nc{'QGRAUP'}(time,:,ilat2,ilon2)+nc{'QVAPOR'}(time,:,ilat2,ilon2);
                case 'vap'
                    Z=nc{'QVAPOR'}(time,:,ilat2,ilon2);
		case 'w'
		    Z=0.5 * ( nc{'W'}(time,1:end-1,ilat2,ilon2) + nc{'W'}(time,2:end,ilat2,ilon2) );
                end
                
                Z_slice(ilat2,:)=Z';
            end
        end
        % end
        
        z_prof=mean(z_slice,1); %mean height profile (mean of all the different heights of each model level)
        
        'finished'
    
        break    

            
            %plot slice
            figure;
            m_gshhs_i('color',[0 0 0]);
            h1=gca;
            m_grid('box','fancy','tickdir','in'); 
            hold on;
            
            
            h=pcolor(X(1,:),Y(:,1)',squeeze(Z));shading flat;
            %             h=slice(repmat(X,[1 1 80]),...
            %                 repmat(Y,[1 1 80]),...
            %                 permute(varz.var,[2 3 1]),...
            %                 permute(Z,[2 3 1]),[],[],levels(j));shading flat;
            jjj=1;
            
            % catch
            
            % end
            %         kkp=find(varz.var(:,1,1)>levels(j));
            %         h=pcolor(X(1,:),...
            %             Y(:,1)',...
            %             squeeze(Z(kkp(1),:,:)));shading flat;
            title(['\it',strrep(Times(time,:),'_','-'),' reflectivity at ',num2str(levels(j)./1000.),'km']);
            
            caxis([0 60]);h2=colorbar;axes(h2);ylabel('\it reflectivity (dBZ)');
            axes(h1);
          %  m_grid('box','fancy','tickdir','in');
            pause(1);
            print(gcf,'-dpng',[imagedir,strrep(['/' num2str(time) '-refl-',Times(time,:),'-',num2str(levels(j)./1000.),'km.png'],':','_')]);
            if(jjj==1)
                delete(h);
            end
            close(hf);
            
end
%close(nc);
