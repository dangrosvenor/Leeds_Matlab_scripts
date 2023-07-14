function [nc_out,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time_in,pole_lat,pole_lon,time_tol)

if  exist(filename)~=2
    error([filename ' file does not exist!!']);
end

version_mat=version;
vfind = strfind(version_mat,'R2007b');    
if length(vfind)>0; version_mat='R2007b'; end
    
%switch version_mat
%    case 'R2007b'        
        nc_out = netcdf(filename);
        lon_UM = nc_out{'x'}(:);
        lat_UM = nc_out{'y'}(:);
 
        
%     otherwise        
%         nc_out = filename;     
%         lon_UM = netcdf_Dan(nc_out,'x','(:)');
%         lat_UM = netcdf_Dan(nc_out,'y','(:)');
% end
    
      [time_driver,time_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,it] ...
            = read_time_UM(nc_out,'t',time_in,time_tol); 
    


    
    
    lat_flag = 'nested';  %In my ususal nested suite runs the lat and lon are in x and y
                            %but for the global output from Florent they
                            %were in 'latitude' and 'longitude'
 

    

    if length(lat_UM)==0
       lat_flag='global';
       %lon_UM = nc_out{'longitude'}(:);
       %lat_UM = nc_out{'latitude'}(:); 
       lon_UM = netcdf_Dan(nc_out,'longitude','(:)');
       lat_UM = netcdf_Dan(nc_out,'latitude','(:)');
       
       lon_UM(lon_UM>180) = lon_UM(lon_UM>180)-360; %convert to this format for consistency with the nested grids.       
    end
    if length(lat_UM)>0
        [lon2d,lat2d]=meshgrid(lon_UM,lat_UM);
        %Convert to normal lat lon from rotated pole coords - make sure the
        %rotated pole lat and lon are given correctly above
        if strcmp(lat_flag,'nested')==1
            [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(lat2d,lon2d,pole_lat,pole_lon);  %From Annette M - see email for an example of how she uses it
        else
            gcm_Plat2D_UM = lat2d;
            gcm_Plon2D_UM = lon2d;
        end
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

    end
    
    