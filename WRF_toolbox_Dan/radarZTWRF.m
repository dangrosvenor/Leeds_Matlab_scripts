function DATA=radarZTWRF(int)
% The function calculates the various radar products to compare with the
% analysis suggested by Peter May.
warning off;
if(int==1)
    fileName='/ddisk1/zhu/hpcx/20060206/sonde_srtm_1800_cloud_zhu_myj_godrad_300CCN/wrfout_d01_2006-02-05_18:00:00';
else
    fileName='/ddisk1/zhu/hpcx/20051116/wrfout_d01_2005-11-16_00:00:00-darwin'
end
nc=netcdf(fileName);
varl=WRFUserARW(nc,'XLONG',1);
varla=WRFUserARW(nc,'XLAT',1);
ZLEV=[10 20 30 40 50];
lat=[-12 -11];lon=[130 131.5];
DATA.ZLEV=ZLEV;

ind=find(varl.var(:)>=lon(1) & varl.var(:)<=lon(2) & varla.var(:)>=lat(1) & varla.var(:)<=lat(2));

Times=nc{'Times'}(:);
[r,c]=size(Times);
DATA.time=zeros(1,r);
[r1,c1]=size(varl.var);
flag1=-99;
h=waitbar(0,'Please wait...')'
for time=1:1:r
    if(flag1==-99)
        varz=WRFUserARW(nc,'Z',time);flag1=1;
        DATA.z=varz.var(:,1,1)./1000;
        DATA.Z1=NaN.*ones(length(DATA.z),r,length(ZLEV));
        DATA.rain=NaN.*ones(length(DATA.z),r);
        DATA.snow=NaN.*ones(length(DATA.z),r);
        DATA.graup=NaN.*ones(length(DATA.z),r);
        DATA.Zmax=NaN.*ones(length(DATA.z),r);
    end
    DATA.time(time)=datenum(str2num(Times(time,1:4)),str2num(Times(time,6:7)),str2num(Times(time,9:10)),...
        str2num(Times(time,12:13)),str2num(Times(time,15:16)),str2num(Times(time,18:19)));
    % Get the various snow field
    varr=WRFUserARW(nc,'QRAIN',time);
    vars=WRFUserARW(nc,'QSNOW',time);
    varg=WRFUserARW(nc,'QGRAUP',time);
    Z=WRFRadarRefl(nc,time,'thompson'); % calculate the radar reflectivity
    Z=real(10.*log10(Z));
    Z(find(Z(:)<0))=0;
    try
        for j=1:length(DATA.z)
            zh1=filter2(fspecial('average',[3 3]),squeeze(Z(j,:,:)));
            rain1=filter2(fspecial('average',[3 3]),squeeze(varr.var(j,:,:)));
            snow1=filter2(fspecial('average',[3 3]),squeeze(vars.var(j,:,:)));
            graup1=filter2(fspecial('average',[3 3]),squeeze(varg.var(j,:,:)));
            for k=1:length(ZLEV)
                DATA.Z1(j,time,k)=length(find(zh1(ind)>=ZLEV(k)))./length(ind);
            end
            
            DATA.Zmax(j,time)=nanmax(zh1(ind));
            DATA.rain(j,time)=length(find(rain1(ind)>=1e-6))./length(ind);
            DATA.snow(j,time)=length(find(snow1(ind)>=1e-6))./length(ind);
            DATA.graup(j,time)=length(find(graup1(ind)>=1e-6))./length(ind);
            
        end
        time
    catch
    end
    waitbar(time./r,h);
end
close(h);
close(nc);
