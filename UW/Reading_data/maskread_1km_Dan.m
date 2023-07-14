function  [mask_1km] = maskread_1km_Dan(mask1)

%==================================================================================
% PROCEDURE	:	maskread_1km
% VERSION	:	1.0
% AUTHOR	:	Robert Wood - adapted for Matlab by Daniel Grosvenor 17th May, 2011
%                                                                 
% 
% DATE		:	March 19 2001
%
% DESCRIPTION	:	Takes HDF buffers for HDF variable ID 
%			'Cloud_Mask_1km' (mask1) 
%			and extracts mask into useable format
%===================================================================================
% BITS MASK

%   bts=2^indgen(8)

% CHECK DIMENSIONS of ARRAY
%   nx=n_elements(mask1(0,*,0))
%   ny=n_elements(mask1(0,0,*))
   
   nx=length(mask1(1,:,1));
   ny=length(mask1(1,1,:));

% OUTPUTS

%   mask_1km=bytarr(9,nx,ny)
   mask_1km=zeros([9 nx ny]);

        
%===========================================================================
% CLOUD MASK FLAGS for 1 km DATA
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html  
%   
%  mask_1km(0)  :  CLOUD MASK 
%			0: undetermined% 1: determined
%  mask_1km(1)  :  CLOUD MASK QUALITY FLAG
%			0: 0-20% cloudy pixels% 1: 20-40%% 3: 40-60%% 4: 60-100%
%  mask_1km(2)  :  DAY/NIGHT FLAG
%			0: night% 1: day
%  mask_1km(3)  :  SUN GLINT FLAG
%			0: yes% 1: no
%  mask_1km(4)  :  SNOW/ICE FLAG
%			0: yes% 1: no
%  mask_1km(5)  :  LAND/WATER FLAG
%			0: water (ocean)% 1: coastal% 2: desert% 3: land
%  mask_1km(6)  :  HEAVY AEROSOL
%			0: yes% 1: no
%  mask_1km(7)  :  THIN CIRRUS DETECTED
%			0: yes% 1: no
%  mask_1km(8)  :  SHADOW FOUND
%			0: yes% 1: no
%
%---------------------------------------------------------------------------   
 
%    bins=dec2bin(mask1,8); %dec2bin converts to an N by 8 array of the bits (in string format)
%   
%    
%    bins = reshape(bins,[size(mask1) 2]);
%    Lsiz = length(size(bins));
%    bins = permute(bins,[Lsiz 1:Lsiz-1]); %put the bits as the first dimension as makes life easier
   %becuase of way Matlab rearranges the dimensions when an array has singular dimensions

   %now we can use the bits using e.g. bin2dec_array(bins(2:4,1,:)); for
   %bits 2,3 and 4 (where 1 is the first bit) of the first byte of the 10
   
   
%    mask_1km(0,*,*)=byte(mask1(0,*,*) and bts(0))
%    mask_1km(1,*,*)=byte(((mask1(0,*,*) and bts(1))/bts(1))+$
%                    ((mask1(0,*,*) and bts(2))/bts(2))*2)
%    mask_1km(2,*,*)=byte((mask1(0,*,*) and bts(3))/bts(3))
%    mask_1km(3,*,*)=byte((mask1(0,*,*) and bts(4))/bts(4))
%    mask_1km(4,*,*)=byte((mask1(0,*,*) and bts(5))/bts(5))
%    mask_1km(5,*,*)=byte(((mask1(0,*,*) and bts(6))/bts(6))+$
%                    ((mask1(0,*,*) and bts(7))/bts(7))*2)  
%    mask_1km(6,*,*)=byte((mask1(1,*,*) and bts(0)))
%    mask_1km(7,*,*)=byte((mask1(1,*,*) and bts(1))/bts(1))
%    mask_1km(8,*,*)=byte((mask1(1,*,*) and bts(2))/bts(2))


%    mask_1km(1,:,:)  = bin2dec_array( bins(8,1,:,:) );
%    mask_1km(2,:,:)  = bin2dec_array( bins(6:7,1,:,:) );
%    mask_1km(3,:,:)  = bin2dec_array( bins(5,1,:,:) );
%    mask_1km(4,:,:)  = bin2dec_array( bins(4,1,:,:) );
%    mask_1km(5,:,:)  = bin2dec_array( bins(3,1,:,:) );
%    mask_1km(6,:,:)  = bin2dec_array( bins(1:2,1,:,:) );
%    mask_1km(7,:,:)  = bin2dec_array( bins(8,2,:,:) );
%    mask_1km(8,:,:)  = bin2dec_array( bins(7,2,:,:) );
%    mask_1km(9,:,:)  = bin2dec_array( bins(6,2,:,:) );
   
   mask_1km(1,:,:)  = bitand(bitshift(mask1(1,:,:), 0),1);      
   mask_1km(2,:,:)  = bitand(bitshift(mask1(1,:,:),-1),3);   
   mask_1km(3,:,:)  = bitand(bitshift(mask1(1,:,:),-3),1);      
   mask_1km(4,:,:)  = bitand(bitshift(mask1(1,:,:),-4),1);   
   mask_1km(5,:,:)  = bitand(bitshift(mask1(1,:,:),-5),1);      
   mask_1km(6,:,:)  = bitand(bitshift(mask1(1,:,:),-6),3);   
   mask_1km(7,:,:)  = bitand(bitshift(mask1(2,:,:), 0),1);      
   mask_1km(8,:,:)  = bitand(bitshift(mask1(2,:,:),-1),1);   
   mask_1km(9,:,:)  = bitand(bitshift(mask1(2,:,:),-2),1);      



end




