function [ilat,ilon,dist_min,isave]=getind_latlon_quick(lat2d,lon2d, LAT, LON ,dlat,dlon)
% finds the latitude longitude index for values given in LAT and LON
% function [ilat,ilon,dist_min]=getind_latlon2d_quick(lat2d,lon2d, LAT, LON,dlat,dlon)
% lat2d and lon2d are the 2D latitude and longitude values for each grid point
% dlat is the accuracy of the first guess made to reduce the number of points being considered
%  x=find(abs(lat2d-LAT(n))<dlat);  %find most likely candiates - so finds all latitudes within dlat degrees
% and then works from those so is quicker than going through all points
% default is dlat=1 if left blank - the smaller it is the longer it will take
% if it's too small the routine will fail
% distmin is an array of the distances between each LAT,LON point and the
% closest point within lat2d, lon2d

if nargin==4
    Ndiv = sqrt(prod(size(lat2d)))/10;
%    dlat = 1.1*(maxALL(lat2d)-minALL(lat2d)) /size(lat2d,1);
%    dlon = 1.1*(maxALL(lon2d)-minALL(lat2d)) /size(lon2d,1);
%    dlat = (maxALL(lat2d)-minALL(lat2d)) /100;
%    dlon = (maxALL(lon2d)-minALL(lat2d)) /100;
    dlat = 1.1*(maxALL(lat2d)-minALL(lat2d)) /Ndiv;
    dlon = 1.1*(maxALL(lon2d)-minALL(lat2d)) /Ndiv;

end

  N  = length( LAT )  ;     %number of LAT,LON pairs to find the indices for

%  lat1d  = ndtooned( lat2d )  %converts 2d array into one long 1d vector
%  lon1d  = ndtooned( lon2d )  
%  n2d    = dimsizes( lat2d )  



for n=1:N
%    n
%        x=find(abs(lat2d-LAT(n))<dlat);  %find most likely candiates
%        y=find(abs(lon2d(x)-LON(n))<dlon);


    
%    if length(y)==0
%        ilat(n)=NaN;
%        ilon(n)=NaN;
%    else
        
        %     dist   = gc_latlon(LAT(n),LON(n),lat1d,lon1d, 2,2)
%        clear distance
%         for i=1:length(y)
%             
%             distance(i) = distlatlon(LAT(n),LON(n),lat2d(x(y(i))),lon2d(x(y(i)))); 
%         end
        
        distance = distlatlon(LAT(n),LON(n),lat2d(:),lon2d(:));
        
        [dist_min(n),I]= min(distance);
        
        %[ilat(n) ilon(n)] = ind2sub(size(lat2d),x(y(I)) );
%        [ilat(n) ilon(n)] = ind2sub(size(lat2d),I);
        isave(n) = I;        
        
        if isnan(dist_min(n))==1
            ilat(n)=NaN;
            ilon(n)=NaN;
        end
        
end
    
[ilat ilon] = ind2sub(size(lat2d),isave);
    
    
    %  ilat(n)=I(1);
    %  ilon(n)=I(2);
    
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
