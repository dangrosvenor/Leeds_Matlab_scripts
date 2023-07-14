%Different options for smoothing/coarse graining UM data

        switch coarsen_method
            case 'lat lon blocks'
                %coarse grain to given lat lon size
                [dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = coarse_grain(dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM,dlat_target,dlon_target);
            case 'gridpoint blocks'

                % Or do as N*M block - perhaps preferable as this will be more
                % equal area since are at high lat
%                N=64; M=N; %64 seems to work the best
                [dat_modis, N_dat_modis] = reduce_matrix_subsample_mean(dat_modis,N,M);
                gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
                gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
                %Work out the cell edges (as halfway between the centres)
                [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

            case '2d smoothing filter'
%                N=128; M=N; 
                B2d = fspecial('average',[M N]);
                dat_modis = FILTER2(B2d,dat_modis);

        end