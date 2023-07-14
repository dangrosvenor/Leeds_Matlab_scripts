function [dat,iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges,iLAT,iLON,iT,iLAT2,iLON2,iT2] = restrict_array3D_to_region(dat,LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM)

stime=size(dat,1);
[iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges,iLAT,iLON,iT,iLAT2,iLON2,iT2] = restrict_to_region_2D_lats(LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);

%Re-order array to put the time dimension last, as this is the case for the
%Plat3D arrays used in the above function
dat=permute(dat,[2 3 1]);
a=NaN*ones(size(dat));
a(iregion_lin)=0;
dat=dat+a;

%Re-arrange back to [time lat lon]
dat=permute(dat,[3 1 2]);