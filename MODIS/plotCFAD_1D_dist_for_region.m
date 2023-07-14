%LAT_val = [25 34]; LON_val = [122:129]; %China Sea

LAT_cal = gcm_Plat2D_CALIPSO_monthly(:,1);
LON_cal = gcm_Plon2D_CALIPSO_monthly(1,:);

ilat = find(LAT_cal>=LAT_val(1) & LAT_cal <LAT_val(end));
ilon = find(LON_cal>=LON_val(1) & LON_cal <LON_val(end));


iregion = find(gcm_Plon2D_CALIPSO_monthly>=LON_val(1) & gcm_Plon2D_CALIPSO_monthly<LON_val(end) & gcm_Plat2D_CALIPSO_monthly>=LAT_val(1) & gcm_Plat2D_CALIPSO_monthly<LAT_val(end));

%months start Jan 2007 and run for 4 years (48 times)
%itime=[1 2]; %Jan Feb 2007
%itime=[6 7 8]; %JJA 2007

%IT=repmat(itime,[1 length(iregion)]);
%IH=

% switch data_type_cf_cal
%     case 'CFAD'

dat = cf_CFAD_sr_CALIPSO(itime,:,ilat(1):ilat(end),ilon(1):ilon(end));
%combines the lat and lon (last 2) dimensions into one linear dimension
dat=dat(:,:,:);

%     case 'CF3D'
%         dat = cf_CFAD_sr_CALIPSO(itime,:,ilat(1):ilat(end),ilon(1):ilon(end));
%         dat=dat(:,:,:);               
% end

if length(itime)>1
    %average over time
    dat = meanNoNan(dat,1); %makes a [H LAT*LON] sized array
    dat = shiftdim(dat,-1); %turn into [1 H LAT*LON]
    
end

cf_cfad = meanNoNan(dat,3); %average over all locations

%figure
%hold on
%plot(cfad_alts_centre_CALIPSO/1e3,cf_cfad,'r--');
%set(gca,'xlim',[-2 10]);
%grid on
%title(['time=' num2str(itime)]);