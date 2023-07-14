function radarPlotsWRF_smallMEM(fileName,imagedir)
warning off;
%fileName='/ddisk1/mccikpc2/wrf/hpcx/run002_noradcloud_300cm3_20060206/wrfout_d02_2006-02-05_18:00:00'
%fileName='/ddisk1/zhu/hpcx/20060206/sonde_srtm_1800_cloud_zhu_myj_godrad_10CCN/wrfout_d01_2006-02-05_18:00:00'
%fileName='/scratch/a/zhu/wrfout_d01_2006-02-05_18:00:00'
%imagedir='~/temp/modelimages1'
nc=netcdf(fileName);
varl=WRFUserARW(nc,'XLONG',1);
varla=WRFUserARW(nc,'XLAT',1);

hf=figure;
m_proj('mercator','long',[min(min(varl.var)) max(max(varl.var))],'lat',[min(min(varla.var)) max(max(varla.var))]);
[X,Y]=m_ll2xy(varl.var,varla.var);


Times=nc{'Times'}(:);
[r,c]=size(Times);
levels=[2500,5000,14000,16000,17400];
%levels=[5000];
for time=19:1:25  %23:1:r
    varz=WRFUserARW(nc,'Z',time);
    
    for j=1:length(levels)
        jjj=99;
       % try
%             for ilat=1:size(varl.var,1)
%                 for ilon=1:size(varl.var,2)                
%                    kkp=find(varz.var(:,ilat,ilon)>levels(j)); %find correct height for each lat lon point
                    kkp=find(varz.var(:,1,1)>levels(j)); %find correct height for each lat lon point
                    
                    Z=WRFRadarRefl_smallMEM(nc,time,'thompson',kkp(1)); %are approximating here since the level of
                                                                         %the desired height will vary across the domain
                                                                         %Here are just plotting the same level everywhere
                    Z=real(10.*log10(Z));            
                                                            
                    
                    % end
                    % end
            
            Z(find(Z(:)<0))=NaN;
            
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
            print(gcf,'-dpng',[imagedir,strrep(['/refl-',Times(time,:),'-',num2str(levels(j)./1000.),'km.png'],':','_')]);
            if(jjj==1)
                delete(h);
            end
            close(hf);
    end
end
close(nc);
