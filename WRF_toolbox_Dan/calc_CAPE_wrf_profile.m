%calc CAPE from WRF output for a profile at a specified point
%uses function [CAPE,CIN,HLCL,TLCL,PLCL,CAPE_surf,CIN_surf,CAPEP,NAP]=calc_cape(P,T,R,RSAT,ALT)
% P = pressure (Pa)
% T = temperature (K)
% R = vapour mixing ratio (kg/kg)
% RSAT = saturation vapour mixing ratio over water (kg/kg)


%specify the indices of the required profile
ilat_wrf=20;
ilon_wrf=30;
time=5;

 f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
 T = WRFUserARW(nc,'tc',time,ilat_wrf,ilon_wrf) + 273.15;
 P = WRFUserARW(nc,'p',time,ilat_wrf,ilon_wrf) *100;                                                      
 R = nc{'QVAPOR'}(time,:,ilat_wrf,ilon_wrf); %cloud MR kg/kg
 RSAT = SatVapPress(T,'goff','liq',P,1) / f; %saturation vapour mixing ratio over water (kg/kg)
 ALT=WRFUserARW(nc,'Z',time,ilat_wrf,ilon_wrf);

 [CAPE,CIN,HLCL,TLCL,PLCL,CAPE_surf,CIN_surf,CAPEP,NAP]=calc_cape(P,T,R,RSAT,ALT);