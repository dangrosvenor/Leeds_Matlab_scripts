function [ilat,ilon]=getind_latlon2d(lat2d,lon2d, LAT, LON )
                        

  N  = length( LAT )  ;     %number of LAT,LON pairs to find the indices for

%  lat1d  = ndtooned( lat2d )  %converts 2d array into one long 1d vector
%  lon1d  = ndtooned( lon2d )  
%  n2d    = dimsizes( lat2d )    

  for n=1:N
%     dist   = gc_latlon(LAT(n),LON(n),lat1d,lon1d, 2,2)
     for i=1:size(lat2d,1);
         for j=1:size(lat2d,2)
            dist(i,j) = distlatlon(LAT(n),LON(n),lat2d(i,j),lon2d(i,j));
        end
    end
    
    [Y,I]= minALL(dist);
    I=squeeze(I);
    
    ilat(n)=I(1);
    ilon(n)=I(2);
    
end
         
    

%   N  = dimsizes( LAT )          
%   ij = new ( (/N,2/) , "integer")
% 
%   lat1d  = ndtooned( lat2d )  
%   lon1d  = ndtooned( lon2d )  
%   n2d    = dimsizes( lat2d )    
% 
%   do n=0,N-1
%      dist   = gc_latlon(LAT(n),LON(n),lat1d,lon1d, 2,2)
%      mndist = min( dist )
%      ind1d  = ind(dist.eq.mndist)
%      if (.not.ismissing(ind1d)) then
%          ij(n,:) = ind_resolve( ind1d(0), n2d )
%      else
%          print("getind_latlon2d: lat="+ LAT(n)+"  lon="+ LON(n)+" problem")
%      end if
% 
%      delete(mndist)
%      delete(ind1d)
%   end do
%   ij@long_name = "indicies closest to specified LAT/LON coordinate pairs"
% 
%   if (.not.any(ismissing(ij))) then
%       delete(ij@_FillValue)
%   end if
%      
%   return( ij )
% 
% 
% 
