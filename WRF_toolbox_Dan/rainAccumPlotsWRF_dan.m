function rainAccumPlotsWRF(fileName,imagedir,projection)
warning off;
fileName='/ddisk1/mccikpc2/wrf/hpcx/run002_noradcloud_1000cm3_20060206/wrfout_d02_2006-02-05_18:00:00'
fileName='/ddisk1/zhu/hpcx/20060206/sonde_srtm_1800_cloud_zhu_myj_godrad_1000CCN/wrfout_d01_2006-02-05_18:00:00'
%fileName='/scratch/a/zhu/wrfout_d01_2006-02-05_18:00:00'
imagedir='~/temp/modelimages4'
nc=netcdf(fileName);
varl=WRFUserARW(nc,'XLONG',1);
varla=WRFUserARW(nc,'XLAT',1);

%m_proj('mercator','long',[min(min(varl.var)) max(max(varl.var))],'lat',[min(min(varla.var)) max(max(varla.var))]);
m_proj(projection,'long',[min(min(varl.var)) max(max(varl.var))],'lat',[min(min(varla.var)) max(max(varla.var))]);
[X,Y]=m_ll2xy(varl.var,varla.var);
m_gshhs_i('color',[0 0 0]);
h1=gca;
m_grid('box','fancy','tickdir','in'); 
hold on;

Times=nc{'Times'}(:);
[r,c]=size(Times);
levels=[2500,5000,7500,10000,12500];
%levels=[5000];
for time=1:1:r
    %varz=WRFUserARW(nc,'Z',time);
    %Z=WRFRadarRefl(nc,time,'thompson');
    %Z=real(10.*log10(Z));
    %Z(find(Z(:)<0))=NaN;
    RAINNC=WRFUserARW(nc,'RAINNC',time);
    RAINC=WRFUserARW(nc,'RAINC',time);
    RAINTOT=RAINNC.var+RAINC.var;RAINTOT(find(RAINTOT(:)==0))=NaN;
    for j=1:1
        jjj=99;
        try
%            kkp=find(varz.var(:,1,1)>levels(j));
            h=pcolor(X(1,:),...
                Y(:,1)',...
                squeeze(RAINTOT));shading flat;
%             h=slice(repmat(X,[1 1 80]),...
%                 repmat(Y,[1 1 80]),...
%                 permute(varz.var,[2 3 1]),...
%                 permute(Z,[2 3 1]),[],[],levels(j));shading flat;
            jjj=1;
        catch
            
        end
%         kkp=find(varz.var(:,1,1)>levels(j));
%         h=pcolor(X(1,:),...
%             Y(:,1)',...
%             squeeze(Z(kkp(1),:,:)));shading flat;
        title(['\it',strrep(Times(time,:),'_','-'),' rain accumulation']);

        caxis([0 150]);h2=colorbar;axes(h2);ylabel('\it rain (mm)');
        axes(h1);
        m_grid('box','fancy','tickdir','in');
        pause(1);
        print(gcf,'-dpng',strrep([imagedir,'/rain-',Times(time,:),'.png'],':','_'));
        if(jjj==1)
            delete(h);
        end
    end
end
close(nc);
