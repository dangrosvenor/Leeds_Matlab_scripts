function [dat_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = coarse_grain(dat_UM,gcm_Plat2D_UM,gcm_Plon2D_UM,dlat_target,dlon_target)
         %need to supply dlat_targer and dlon_target
%         dlat_target = dlat_CERES;
%         dlon_target = dlon_CERES;

         d=diff(gcm_Plat2D_UM,[],1);
         dlat_UM = meanNoNan(meanNoNan(d,1),1);
         N = ceil(abs(dlat_target/dlat_UM));

         d=diff(gcm_Plon2D_UM,[],2);
         dlon_UM = meanNoNan(meanNoNan(d,1),1);
         M = ceil(abs(dlon_target/dlon_UM));

         dat_out = reduce_matrix_subsample_mean(dat_UM,N,M);
         gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
         gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
         %Work out the cell edges (as halfway between the centres)
         [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);