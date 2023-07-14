function var=WRFUserARW( nc_file, variable, time, ilat, ilon, iz)

%nc_file=netcdf('/data/pjc/WRF/WRFV2/test/em_real/wrfout_d01_2006-02-05_18:00:00')
%variable='Z'
%time=1

variable=lower(variable);

if (nargin==3)  
    profile=0;
    horiz_slice=0;
    single=0;
elseif nargin==5
    %Dan - if only want profiles at one latiude/longitude
    % 	latitude  = nc_file{'XLAT'}(time,:,:);
    % 	longitude = nc_file{'XLONG'}(time,:,:);
    %     [ilat,ilon] = getind_latlon(latitude,longitude,lat,lon); %returns index of point nearest to lat lon (i.e. closest_lat=latitude(ilat,ilon)
    %                                                              %closest_lon=longitude(ilat,ilon)
    profile=1;    
    horiz_slice=0;
    single=0;
elseif nargin==4
    horiz_slice=1;
    profile=0;  %set flag to say only want profile or not
    single=0;
    iz=ilat;
elseif nargin==6
    single=1;
    profile=0;
    horiz_slice=0;
    
end


if( strcmp(variable , 'umet') | strcmp(variable ,'vmet'))
    
    if profile==0

   vartmp = nc_file{'V'}(time,:,:,:);

   v = 0.5.*(vartmp(:,1:end-1,:)+vartmp(:,2:end,:));
   vartmp = nc_file{'U'}(time,:,:,:);
   u = 0.5.*(vartmp(:,:,1:end-1) + vartmp(:,:,2:end));

   pii = 3.14159265;
   radians_per_degree = pii/180.;
   map_projection = nc_file.MAP_PROJ(:);
   if(map_projection == 0) % no projection

     if(strcmp(variable, 'umet')) 
       var.u = u;
       var.description = 'u met velocity';
       var.units = 'm/s';
       return;
     else  % must want vmet
       var.v = v;
       var.description = 'v met velocity';
       var.units = 'm/s';
       return;
     end
   end

   cen_lat  = nc_file.CEN_LAT(:);
   if(length(nc_file.STAND_LON(:)))
       cen_long = nc_file.STAND_LON(:);
   else
       cen_long = nc_file.CEN_LON(:);
   end
   true_lat1 = nc_file.TRUELAT1(:);
   true_lat2 = nc_file.TRUELAT2(:);
   latitude  = nc_file{'XLAT'}(1,:,:);
   longitude = nc_file{'XLONG'}(1,:,:);

   cone = 1.
   if( map_projection == 1)   % Lambert Conformal mapping
     if( (abs(true_lat1 - true_lat2) > 0.1) & ...
         (abs(true_lat2 - 90. )      > 0.1)       )
         cone = 10.^(cos(true_lat1.*radians_per_degree)) ...
               -10.^(cos(true_lat2.*radians_per_degree));
         cone = cone./(10.^(tan(45. -abs(true_lat1./2.).*radians_per_degree)) - ...
                      10.^(tan(45. -abs(true_lat2./2.).*radians_per_degree)));  
     else
         cone = sin(abs(true_lat1).*radians_per_degree);
     end
   end
   if(map_projection == 2)     % polar steraographic
     cone = 1.;
   end
   if(map_projection == 3)     % Mercator
     cone = 0.;
   end

   diff = longitude - cen_long;
   dims = size(longitude)

   diff(find(diff(:)>180))=diff(find(diff(:)>180))-360.;
   diff(find(diff(:)<-180))=diff(find(diff(:)<-180))+360.;

%      alpha = diff * cone * radians_per_degree *sign(1.,latitude)
   alpha = diff;
   alpha(find(latitude(:)<0)) = -diff(find(latitude(:)<0)).*cone.*radians_per_degree;
   alpha(find(latitude(:)>=0)) = diff(find(latitude(:)>=0)).*cone.*radians_per_degree;

   dims=size(v);
   if(strcmp(variable,'umet'))
     for k=1:dims(1)
       var.v(k,:,:) = v(k,:,:).*sin(alpha) + u(k,:,:).*cos(alpha);
     end
     var.description = 'u met velocity';
     var.units = 'm/s';
   else  % must want vmet
     for k=1:dims(1)
       var.v(k,:,:) = v(k,:,:).*cos(alpha) - u(k,:,:).*sin(alpha);
     end 
     var.description = 'v met velocity';
     var.units = 'm/s';
   end
   
else  %if profile==0
    vartmp = nc_file{'V'}(time,:,ilat:ilat+1,ilon:ilon+1);

   v = 0.5.*(vartmp(:,1,1)+vartmp(:,2,1));
   vartmp = nc_file{'U'}(time,:,ilat:ilat+1,ilon:ilon+1);
   u = 0.5.*(vartmp(:,1,1) + vartmp(:,1,2));

   pii = 3.14159265;
   radians_per_degree = pii/180.;
   map_projection = nc_file.MAP_PROJ(:);
   if(map_projection == 0) % no projection

     if(strcmp(variable, 'umet')) 
       var = u;
       %var.description = 'u met velocity';
       %var.units = 'm/s';
       return;
     else  % must want vmet
       var = v;
       %var.description = 'v met velocity';
       %var.units = 'm/s';
       return;
     end
   end

   cen_lat  = nc_file.CEN_LAT(:);
   if(length(nc_file.STAND_LON(:)))
       cen_long = nc_file.STAND_LON(:);
   else
       cen_long = nc_file.CEN_LON(:);
   end
   true_lat1 = nc_file.TRUELAT1(:);
   true_lat2 = nc_file.TRUELAT2(:);
   latitude  = nc_file{'XLAT'}(time,:,:);
   longitude = nc_file{'XLONG'}(time,:,:);

   cone = 1.
   if( map_projection == 1)   % Lambert Conformal mapping
     if( (abs(true_lat1 - true_lat2) > 0.1) & ...
         (abs(true_lat2 - 90. )      > 0.1)       )
         cone = 10.^(cos(true_lat1.*radians_per_degree)) ...
               -10.^(cos(true_lat2.*radians_per_degree));
         cone = cone./(10.^(tan(45. -abs(true_lat1./2.).*radians_per_degree)) - ...
                      10.^(tan(45. -abs(true_lat2./2.).*radians_per_degree)));  
     else
         cone = sin(abs(true_lat1).*radians_per_degree);
     end
   end
   if(map_projection == 2)     % polar steraographic
     cone = 1.;
   end
   if(map_projection == 3)     % Mercator
     cone = 0.;
   end

   diff = longitude - cen_long;
   dims = size(longitude)

   diff(find(diff(:)>180))=diff(find(diff(:)>180))-360.;
   diff(find(diff(:)<-180))=diff(find(diff(:)<-180))+360.;

%      alpha = diff * cone * radians_per_degree *sign(1.,latitude)
   alpha = diff;
   alpha(find(latitude(:)<0)) = -diff(find(latitude(:)<0)).*cone.*radians_per_degree;
   alpha(find(latitude(:)>=0)) = diff(find(latitude(:)>=0)).*cone.*radians_per_degree;

   dims=size(v);
   if(strcmp(variable,'umet'))
     for k=1:dims(1)
       var(k) = v(k).*sin(alpha(ilat,ilon)) + u(k).*cos(alpha(ilat,ilon));
     end
     %var.description = 'u met velocity';
     %var.units = 'm/s';
   else  % must want vmet
     for k=1:dims(1)
       var(k) = v(k).*cos(alpha(ilat,ilon)) - u(k).*sin(alpha(ilat,ilon));
     end 
     %var.description = 'v met velocity';
     %var.units = 'm/s';
   end
   
end

 return

end




if( (strcmp(variable,'umeta')) | (strcmp(variable,'vmeta') )) 
    
    
if profile==0

   v = nc_file{'V'}(time,:,:,:);
   u = nc_file{'U'}(time,:,:,:);

   pii = 3.14159265;
   radians_per_degree = pii/180.;

   map_projection = nc_file.MAP_PROJ(:);

   if(map_projection == 0) % no projection

     if(strcmp(variable,'umeta'))
       var.u = 0.5.*(u(:,:,1:end-1) + u(:,:,2:end));
       var.description = 'u met velocity';
       var.units = 'm/s';
       return;
     else  % must want vmeta
       var.v = 0.5.*(v(:,1:end-1,:)+v(:,2:end,:));
       var.description = 'v met velocity';
       var.units = 'm/s';
       return;
     end
   end

   cen_lat  = nc_file.CEN_LAT;
   if(length(nc_file.STAND_LON(:)))
       cen_long = nc_file.STAND_LON(:);
   else
       cen_long = nc_file.CEN_LON(:);
   end
   true_lat1 = nc_file.TRUELAT1(:);
   true_lat2 = nc_file.TRUELAT2(:);
   latitude  = nc_file{'XLAT'}(1,:,:);
   longitude = nc_file{'XLONG'}(1,:,:);

   cone = 1.;
   if( map_projection == 1) % Lambert Conformal mapping
     if( (abs(true_lat1 - true_lat2) > 0.1) & ...
         (abs(true_lat2 - 90. )      > 0.1)       ) 
         cone = 10.^(cos(true_lat1.*radians_per_degree)) ...
               -10.^(cos(true_lat2.*radians_per_degree));
         cone = cone./(10.^(tan(45. -abs(true_lat1./2.).*radians_per_degree)) - ...
                      10.^(tan(45. -abs(true_lat2./2.).*radians_per_degree))  );
     else
         cone = sin(abs(true_lat1).*radians_per_degree);
     end
   end
   if(map_projection == 2)      % polar steraographic
     cone = 1.;
   end
   if(map_projection == 3)       % Mercator
     cone = 0.;
   end

   dimv=size(v);
   dimu=size(u);

   diff = longitude - cen_long;
   diff(find(diff(:))>180) = diff(find(diff(:))>180)-360.;
   diff(find(diff(:))<180) = diff(find(diff(:))<180)+360.;
    
   alpha=diff;
   alpha(find(latitude(:))<0) = -diff(find(latitude(:))<0).*cone.*radians_per_degree;
   alpha(find(latitude(:))<0) = diff(find(latitude(:))>=0).*cone.*radians_per_degree;

   diff=cos(alpha);
   alpha=sin(alpha);

   for k=1:dimv(3)
       uk=0.5.*(u(1:end-1,:,k)+u(2:end,:,k));
       vk=0.5.*(v(:,1:end-1,k)+v(:,2:end,k));
       var.uvmet(:,:,k,1)=vk.*alpha+uk.*diff;
       var.uvmet(:,:,k,2)=vk.*diff-uk.*alpha;
   end
   diff(' returned from compute ');

     var.description = ' u,v met velocity';
     var.units = 'm/s';
     
else %if profile==0
    
   v = nc_file{'V'}(time,:,ilat:ilat+1,ilon:ilon+1);
   u = nc_file{'U'}(time,:,ilat:ilat+1,ilon:ilon+1);

   pii = 3.14159265;
   radians_per_degree = pii/180.;

   map_projection = nc_file.MAP_PROJ(:);

   if(map_projection == 0) % no projection

     if(strcmp(variable,'umeta'))
      % var = 0.5.*(u(:,:,1) + u(:,:,2));
      %% var.description = 'u met velocity';
      %% var.units = 'm/s';
       return;
     else  % must want vmeta
      % var = 0.5.*(v(:,1,:)+v(:,2,:));
      %% var.description = 'v met velocity';
      %% var.units = 'm/s';
       return;
     end
   end

   cen_lat  = nc_file.CEN_LAT;
   if(length(nc_file.STAND_LON(:)))
       cen_long = nc_file.STAND_LON(:);
   else
       cen_long = nc_file.CEN_LON(:);
   end
   true_lat1 = nc_file.TRUELAT1(:);
   true_lat2 = nc_file.TRUELAT2(:);
   latitude  = nc_file{'XLAT'}(1,:,:);
   longitude = nc_file{'XLONG'}(1,:,:);

   cone = 1.;
   if( map_projection == 1) % Lambert Conformal mapping
     if( (abs(true_lat1 - true_lat2) > 0.1) & ...
         (abs(true_lat2 - 90. )      > 0.1)       ) 
         cone = 10.^(cos(true_lat1.*radians_per_degree)) ...
               -10.^(cos(true_lat2.*radians_per_degree));
         cone = cone./(10.^(tan(45. -abs(true_lat1./2.).*radians_per_degree)) - ...
                      10.^(tan(45. -abs(true_lat2./2.).*radians_per_degree))  );
     else
         cone = sin(abs(true_lat1).*radians_per_degree);
     end
   end
   if(map_projection == 2)      % polar steraographic
     cone = 1.;
   end
   if(map_projection == 3)       % Mercator
     cone = 0.;
   end

   dimv=size(v);
   dimu=size(u);

   diff = longitude - cen_long;
   diff(find(diff(:))>180) = diff(find(diff(:))>180)-360.;
   diff(find(diff(:))<-180) = diff(find(diff(:))<-180)+360.;
    
   alpha=zeros(size(diff));
   alpha(find(latitude(:)<0)) = -diff(find(latitude(:)<0)).*cone.*radians_per_degree;
   alpha(find(latitude(:)>=0)) = diff(find(latitude(:)>=0)).*cone.*radians_per_degree;

   diff=cos(alpha);
   alpha=sin(alpha);

   for k=1:dimv(1)
       uk=0.5.*(u(k,1,1)+u(k,2,1));
       vk=0.5.*(v(k,1,1)+v(k,1,2));
       tmp(k,1)=vk.*alpha(ilat,ilon)+uk.*diff(ilat,ilon);  %umeta
       tmp(k,2)=vk.*diff(ilat,ilon)-uk.*alpha(ilat,ilon);  %vmeta
   end
  % disp(' returned from compute ');

  
  if strcmp(variable,'umeta')
      var=tmp(:,1);
  else
      var=tmp(:,2);
  end
   %  var.description = ' u,v met velocity';
   %  var.units = 'm/s';
     
 end

 return;

end

if( strcmp(variable,'ua') )
   vartmp = nc_file{'U'}(time,:,:,:);
   var.var = 0.5.*(vartmp(:,:,1:end-1) + vartmp(:,:,2:end));
   var.description = 'u Velocity';
   var.units = 'm/s';

   return;

end

if( strcmp(variable,'u') ) 
    
    if profile==0
        var.var = nc_file{'U'}(time,:,:,:);
        var.description = 'u Velocity';
        var.units = 'm/s';
    else
        var=nc_file{'U'}(time,:,ilat,ilon);
    end

   return

end

if( strcmp(variable,'va' ))
   vartmp = nc_file{'V'}(time,:,:,:);
   var.var = 0.5.*(vartmp(:,1:end-1,:)+vartmp(:,2:end,:));
   var.description = 'v Velocity';
   var.units = 'm/s';

   return

end

if( strcmp(variable,'v') ) 
    
    if profile==0
        var.var = nc_file{'V'}(time,:,:,:);
        var.description = 'v Velocity';
        var.units = 'm/s';
        
    else
        var=nc_file{'V'}(time,:,ilat,ilon);
    end

   return

end

if( strcmp(variable,'wa') ) 
   vartmp = nc_file{'W'}(time,:,:,:);
   var.var = 0.5.*(vartmp(1:end-1,:,:)+vartmp(2:end,:,:));
   var.description = 'Vertical Velocity w';
   var.units = 'm/s';

   return

end

if( strcmp(variable,'w') ) 
   var.var = nc_file{'W'}(time,:,:,:);
   var.description = 'Vertical Velocity w';
   var.units = 'm/s';

   return

end

if( strcmp(variable,'eta-dot') ) 
   var.var = nc_file{'WW'}(time,:,:,:);
   var.description = 'Coordinate Vertical Velocity';
   var.units = 's-1';

   return;

end

if( strcmp(variable,'th') ) 
    if  horiz_slice==0 & profile==0 & single==0
        var.var = nc_file{'T'}(time,:,:,:);
        var.var = var.var + 300.;
        var.description = 'Potential Temperature (theta) ';
        var.units = 'K';
    elseif horiz_slice==1 & profile==0 & single==0
        var = nc_file{'T'}(time,iz,:,:);
        var = var + 300.;
    elseif horiz_slice==0 & profile==1 & single==0
        var = nc_file{'T'}(time,:,ilat,ilon);
        var = var + 300.;
    elseif horiz_slice==0 & profile==0 & single==1
        var = nc_file{'T'}(time,iz,ilat,ilon);
        var = var + 300.;
    end
    
    return
    
end

if( strcmp(variable,'p') ) 
    
    if  horiz_slice==0 & profile==0 & single==0
        var = nc_file{'P'}(time,:,:,:);
        base = nc_file{'PB'}(time,:,:,:);
        var = 0.01*(var+base); 
        
%       var.description = 'Pressure';
%       var.units = 'mb';
        
    elseif horiz_slice==1 & profile==0 & single==0
        var = nc_file{'P'}(time,iz,:,:);
        base = nc_file{'PB'}(time,iz,:,:);
        var = 0.01*(var+base);
        
    elseif horiz_slice==0 & profile==0 & single==1
        var = nc_file{'P'}(time,iz,ilat,ilon);
        base = nc_file{'PB'}(time,iz,ilat,ilon);
        var = 0.01*(var+base);
        
    elseif horiz_slice==0 & profile==1 & single==0
        var = nc_file{'P'}(time,:,ilat,ilon);
        base = nc_file{'PB'}(time,:,ilat,ilon);
        var = 0.01*(var+base);
        
    end            
    

   return;

end

if( strcmp(variable,'z') ) 
%    vartmp = nc_file{'PH'}(time,:,:,:);
%    base = nc_file{'PHB'}(time,:,:,:);
%    vartmp = (vartmp+base)./9.81;
%    var.var= 0.5.*(vartmp(1:end-1,:,:)+vartmp(2:end,:,:));
%    var.description = 'Height';
%    var.units = 'm';
   
   
%   vartmp = nc_file{'PH'}(time,:,:,:);
%   base = nc_file{'PHB'}(time,:,:,:);


   if profile==0 & horiz_slice==0 & single==0
   
       var = ( nc_file{'PH'}(time,:,:,:) + nc_file{'PHB'}(time,:,:,:) )./9.81;
       var(1:end-1,:,:) = 0.5.*(var(1:end-1,:,:)+var(2:end,:,:));
       var(end,:,:)=[];
%       var.description = 'Height';
%       var.units = 'm';
       
   elseif horiz_slice==1 & profile==0 & single==0
       var = ( nc_file{'PH'}(time,iz:iz+1,:,:) + nc_file{'PHB'}(time,iz:iz+1,:,:) )./9.81;
       var(1:end-1,:,:) = 0.5.*(var(1:end-1,:,:)+var(2:end,:,:));
       var(end,:,:)=[];
   elseif horiz_slice==0 & profile==1 & single==0       
       var = ( nc_file{'PH'}(time,:,ilat,ilon) + nc_file{'PHB'}(time,:,ilat,ilon) )./9.81;
       var(1:end-1) = 0.5.*(var(1:end-1)+var(2:end));
       var(end)=[];
%       var.description = 'Height';
%       var.units = 'm';
   elseif horiz_slice==0 & profile==0 & single==1       
       var = ( nc_file{'PH'}(time,iz:iz+1,ilat,ilon) + nc_file{'PHB'}(time,iz:iz+1,ilat,ilon) )./9.81;
       var(1:end-1) = 0.5.*(var(1:end-1)+var(2:end));
       var(end)=[];
   
   end
       

   return

end

if( strcmp(variable,'tc') ) 

% compute theta, p, and then tc

%   vartheta = nc_file{'T'}(time,:,:,:);
%   vartheta = vartheta+300.;

 %  varp = nc_file{'P'}(time,:,:,:) + nc_file{'PB'}(time,:,:,:) ;
  % base = nc_file{'PB'}(time,:,:,:);
   %varp = (varp+base);

%   varp=(varp./100000).^(287./1005);

%   varp=( ( nc_file{'P'}(time,:,:,:) + nc_file{'PB'}(time,:,:,:) ) ./100000).^(287./1005);

   
%   var.var=varp.*vartheta;

%rewritten to use as little memory as possible (fewer stored variables)

if profile==0 & horiz_slice==0
   var=   ( ( nc_file{'P'}(time,:,:,:) + nc_file{'PB'}(time,:,:,:) ) ./100000) .^(287./1005)   .*( nc_file{'T'}(time,:,:,:) + 300 );

   var = var - 273.16;
%   var.description = 'Temperature';
%   var.units = 'C';
   
elseif profile==1 & horiz_slice==0
   var=   ( ( nc_file{'P'}(time,:,ilat,ilon) + nc_file{'PB'}(time,:,ilat,ilon) ) ./100000) .^(287./1005)   .*( nc_file{'T'}(time,:,ilat,ilon) + 300 );

   var = var - 273.16;
%   var.description = 'Temperature';
%   var.units = 'C';

else %horiz slice
    var=   ( ( nc_file{'P'}(time,iz,:,:) + nc_file{'PB'}(time,iz,:,:) ) ./100000) .^(287./1005)   .*( nc_file{'T'}(time,iz,:,:) + 300 );
    var = var - 273.16;
   
end

   return

end

if( strcmp(variable,'tc2') ) 

   vartheta = nc_file{'TH2'}(time,:,:);

   varp = nc_file{'P'}(time,1,:,:);
   base = nc_file{'PB'}(time,1,:,:);
   varp = (varp+base);


   pi=(varp./100000.).^(287./1005);
   var.var=pi.*vartheta;
   var.var = var.var - 273.16;
   var.description = '2m Temperature';
   var.units = 'C';

   clear varp;
   clear vartheta;

   return

end

if( strcmp(variable,'td') ) 

% compute p, then td

   qv = nc_file{'QVAPOR'}(time,:,:,:);
   p = nc_file{'P'}(time,:,:,:);
   base = nc_file{'PB'}(time,:,:,:);
   p = 0.01*(p+base);


   qv=max(qv,0);
   tdc=qv.*p./(0.622+qv);
   tdc=max(tdc,0.001);
   var.var=(243.5.*log(tdc)-440.8)./(19.48-log(tdc));
%      qv = qv > 0.000
%      var = qv*p/(.622+qv)  ; vapor pressure 
%      var = var > 0.001            ; avoid problems near zero
%      var = (243.5/( (17.67/log(var/6.112)) - 1.0)) ; Bolton's approximation
%      var = (243.5*log(var)-440.8)/(19.48-log(var))

   var.description = 'Dewpoint Temperature';
   var.units = 'C';

   return

end

if( strcmp(variable,'td2') ) 
   qv = nc_file{'Q2'}(time,:,:);
   p = nc_file{'P'}(time,1,:,:);
   base = nc_file{'PB'}(time,1,:,:);
   p = 0.01*(p+base);

   qv=max(qv,0);
   tdc=qv.*p./(0.622+qv);
   tdc=max(tdc,0.001);
   var.var=(243.5.*log(tdc)-440.8)./(19.48-log(tdc));

   var.description = '2m Dewpoint Temperature'
   var.units = 'C';

   clear p;
   clear qv;

   return

end

if( strcmp(variable ,'iclw') ) 

% compute p, then iclw (p (or P) is in Pa)

   qc = nc_file{'QCLOUD'}(time,:,:,:);
   p = nc_file{'P'}(time,:,:,:);
   base = nc_file{'PB'}(time,:,:,:);
   p = 0.01*(p+base);

   dimv = size(qc);

   gg=1000./9.8;
   iclw=zeros([dimv(2) dimv(3)]);
   qclw=max(qc,0);
   for k=1:dimv(1)
       if(k==1)
           dp=p(:,:,k-1)-p(:,:,k);
       elseif(k==dimv(1))
           dp=p(:,:,k)-p(:,:,k+1);
       else
           dp=(p(:,:,k-1)-p(:,:,k+1))./2;
       end
       iclw=iclw+qclw.*dp.*gg;
   end
   var.var=iclw;
   var.description = 'Int Cloud Water';
   var.units = 'mm';

   return;

end

if( strcmp(variable,'slvl') ) 

% compute theta

   vartheta = nc_file{'T'}(time,:,:,:);
   vartheta = vartheta+300.;

   p = nc_file{'P'}(time,:,:,:);
   base = nc_file{'PB'}(time,:,:,:);
   p = (p+base);

   pi=(p./100000).^(287./1005);
   tk=pi.*vartheta;

   qv = nc_file{'QVAPOR'}(time,:,:,:);
   qv = max(qv,0.000);

%    surf = new( (/ dimv(1), dimv(2) /), float)
%    t_surf = new( (/ dimv(1), dimv(2) /), float)
%    tmp1 = new( (/ dimv(1), dimv(2) /), float)
%    itmp1 = new( (/ dimv(1), dimv(2) /), integer)

   zw = nc_file{'PH'}(time,:,:,:);
   zbase = nc_file{'PHB'}(time,:,:,:);
   zw = (zw + zbase)./9.81;
   z = 0.5.*(zw(1:end-1,:,:)+zw(2:end,:,:));

%        wrf_user_fortran_util_0 :: compute_seaprs( dimv(2),dimv(1),dimv(0),  \
%                                                   z, tk, p, qv,             \
%                                                   surf, t_surf,             \
%                                                   tmp1, itmp1 )


%    surf = 0.01*surf;
%    surf.description = 'Sea Level Pressure';
%    surf.units = 'mb';

   return

end

if( strcmp(variable,'rh') ) 

% compute p, tc, then rh

   qv = nc_file{'QVAPOR'}(time,:,:,:);
   qv = max(qv,0.000);
   p = nc_file{'P'}(time,:,:,:);
   base = nc_file{'PB'}(time,:,:,:);
   p = (p+base);

   vartheta = nc_file{'T'}(time,:,:,:);
   vartheta = vartheta+300.;


   pi=(p./100000).^(287./1005);
   tk=pi.*vartheta;

   es=10.*0.6112.*exp(17.67.*(tk-273.15)./(tk-29.65));
   qvs=0.622.*es./(0.01.*p-(1-0.622).*es);
   var.var=100.*(max(min(qv./qvs,1),0));

   var.description = 'Relative Humidty';
   var.units = '%';

   return

end

%  end of diagnostic variable list - we must want a variable already in
%  the file.  check variable dimensionality and pull proper time  out of file

var.var = nc_file{upper(variable)}(time,:);

return


