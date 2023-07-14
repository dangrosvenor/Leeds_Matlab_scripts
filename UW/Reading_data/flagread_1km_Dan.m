function [qapq_1km,qapp_1km] = flagread_1km_Dan(qa1,all_flags)

%============================================================================
% PROCEDURE	:	maskread_1km
% VERSION	:	1.0
% AUTHOR	:	Robert Wood - adapted for Matlab by Daniel Grosvenor, 17th
% May, 2011.
% DATE		:	March 19 2001
%
% DESCRIPTION	:	Takes HDF buffers for HDF variable ID 
%			'Quality_Assurance_1km' (qa1) 
%			and extracts flags into useable format
%============================================================================
% BITS MASK

%   bts=2^indgen(8)

% CHECK DIMENSIONS of ARRAY
%   nx=n_elements(qa1(0,*,0))
%   ny=n_elements(qa1(0,0,*))
   
   nx=length(qa1(1,:,1));
   ny=length(qa1(1,1,:));

% OUTPUTS

%   qapq_1km=bytarr(8,nx,ny)
%   qapq_1km=bytarr(16,nx,ny)     
   
    qapq_1km=zeros([9 nx ny]);
    qapp_1km=zeros([16 nx ny]);


%===========================================================================
%  QUALITY ASSURANCE FLAGS for 1km DATA - PRODUCT QUALITY (BYTES 0 and 1)
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html
%
%  qapq_1km(0) : OPTICAL THICKNESS RETRIEVAL 
%                    0: not useful% 1 useful
%  qapq_1km(1) : OPTICAL THICKNESS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_1km(2) : OPTICAL THICKNESS OUT OF BOUNDS 
%                    0: within bounds (tau<150)
%                    1: (100<tau<150)
%                    2: (tau>150)
%                    3: surface reflectance too large
%  qapq_1km(3) : EFFECTIVE RADIUS GENERAL QA
%		     0: not useful
% 		     1: useful
%  qapq_1km(4) : EFFECTIVE RADIUS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good			
%  qapq_1km(5) : WATER PATH GENERAL QA
%		     0: not useful
% 		     1: useful
%  qapq_1km(6) : WATER PATH CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_1km(7) : CLOUD PHASE DETERMINATION
%		     0: SWIR algorithm not run
%		     1: CLEAR		
%		     2: WATER
%		     3: ICE
%		     4: MIXED PHASE or UNDETERMINED	

%  qapq_1km(8) : 1621 cloud retrieval outcome
%		     0: not attempted or unsuccessful
%		     1: successful		

%  
%---------------------------------------------------------------------------

%    bins=dec2bin(qa1,8); %dec2bin converts to an N by 8 array of the bits (in string format)
%    
%    bins = reshape(bins,[size(qa1) 8]);
%    Lsiz = length(size(bins));
%    bins = permute(bins,[Lsiz 1:Lsiz-1]); %put the bits as the first dimension as makes life easier
   %becuase of way Matlab rearranges the dimensions when an array has singular dimensions

   %now we can use the bits using e.g. bin2dec_array(bins(2:4,1,:)); for
   %bits 2,3 and 4 (where 1 is the first bit) of the first byte of the 10
       

%    qapq_1km(0,*,*)=reform(byte(qa1(0,*,*) and bts(0))) 
%    qapq_1km(1,*,*)=reform(byte(((qa1(0,*,*) and bts(1))/bts(1))+$
%                ((qa1(0,*,*) and bts(2))/bts(2))*2))
%    qapq_1km(2,*,*)=reform(byte(((qa1(0,*,*) and bts(3))/bts(3))+$
%                ((qa1(0,*,*) and bts(4))/bts(4))*2))
%    qapq_1km(3,*,*)=reform(byte(qa1(0,*,*) and bts(5))/bts(5)) 
%    qapq_1km(4,*,*)=reform(byte( ((qa1(0,*,*) and bts(6))/bts(6))+$
%                ((qa1(0,*,*) and bts(7))/bts(7))*2 ))   
%    qapq_1km(5,*,*)=reform(byte(qa1(1,*,*) and bts(0))/bts(0)) 
%    qapq_1km(6,*,*)=reform(byte( ((qa1(1,*,*) and bts(1))/bts(1))+$
%                ((qa1(1,*,*) and bts(2))/bts(2))*2 ))
%    qapq_1km(7,*,*)=reform(byte( ( (qa1(1,*,*) and bts(3)) /bts(3) )+$
%                ( (qa1(1,*,*) and bts(4)) /bts(4) )*2 + $
%                ( (qa1(1,*,*) and bts(5)) /bts(5) )*4 ))

% 
%    qapq_1km(1,:,:)  = bin2dec_array( bins(8,1,:,:) );
%    qapq_1km(2,:,:)  = bin2dec_array( bins(6:7,1,:,:) );
%    qapq_1km(3,:,:)  = bin2dec_array( bins(4:5,1,:,:) );
%    qapq_1km(4,:,:)  = bin2dec_array( bins(3,1,:,:) );
%    qapq_1km(5,:,:)  = bin2dec_array( bins(1:2,1,:,:) );
%    qapq_1km(6,:,:)  = bin2dec_array( bins(8,2,:,:) );
%    qapq_1km(7,:,:)  = bin2dec_array( bins(6:7,2,:,:) );
%    qapq_1km(8,:,:)  = bin2dec_array( bins(3:5,2,:,:) );

   
   
   qapq_1km(1,:,:)  = bitand(bitshift(qa1(1,:,:), 0),1); %COT usefulness
   qapq_1km(2,:,:)  = bitand(bitshift(qa1(1,:,:),-1),3); %COT confidence
   qapq_1km(3,:,:)  = bitand(bitshift(qa1(1,:,:),-3),3); %COT out-of-bounds
   qapq_1km(4,:,:)  = bitand(bitshift(qa1(1,:,:),-5),1); %re usefulness     
   qapq_1km(5,:,:)  = bitand(bitshift(qa1(1,:,:),-6),3); %re confidence
   qapq_1km(6,:,:)  = bitand(bitshift(qa1(2,:,:), 0),1); %CWP usefulness
   qapq_1km(7,:,:)  = bitand(bitshift(qa1(2,:,:),-1),3); %CWP confidence
   qapq_1km(8,:,:)  = bitand(bitshift(qa1(2,:,:),-3),7); %phase flag for 1621 
   qapq_1km(9,:,:)  = bitand(bitshift(qa1(2,:,:),-6),1); %1621 cloud retrieval outcome  
   
if all_flags==1

%===========================================================================
% RETRIEVAL PROCESSING FLAGS - PROCESSING PATH FLAGS for 1 km DATA (BYTE 2)
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html
%
%  qapp_1km(0)  :  CLOUD PHASE USED IN RETRIEVAL PROCESSING
%			0: cloud mask undetermine
%			1: decision tree "stop"
%			2: cloud, liquid water phase
%			3: cloud, ice phase
%			4: cloud, undetermined phase
%			5: cloud, mixed phase
%			6-7: not used
%  qapp_1km(1)  :  RETIEVAL OUTCOME
%			0: not attempted or unsuccessful
%			1: retrieval successful
%  qapp_1km(2)  :  RAYLEIGH CORRECTION
%			0: no% 1: yes
%  qapp_1km(3)  :  ATMOSPHERIC CORRECTION
% 			0: no% 1: yes
%  qapp_1km(4)  :  BAND USED FOR OPTICAL THICKNESS RETRIEVAL
%			0: retrieval not attempted
%			1: 0.645 micron (land)
%			2: 0.858 micron (water)
%			3: 1.24 micron (snow/ice)	
%

%    qapp_1km(0,*,*)=reform(byte(((qa1(2,*,*) and bts(0))/bts(0))+$
%                ((qa1(2,*,*) and bts(1))/bts(1))*2 + $
%                ((qa1(2,*,*) and bts(2))/bts(2))*4 )) 
%    qapp_1km(1,*,*)=reform(byte(qa1(2,*,*) and bts(3))) 
%    qapp_1km(2,*,*)=reform(byte(qa1(2,*,*) and bts(4)))
%    qapp_1km(3,*,*)=reform(byte(qa1(2,*,*) and bts(5)))
%    qapp_1km(4,*,*)=reform(byte(((qa1(2,*,*) and bts(6))/bts(6))+$
%                ((qa1(2,*,*) and bts(7))/bts(7))*2))
               
%    qapp_1km(1,:,:)  = bin2dec_array( bins(6:8,3,:,:) );
%    qapp_1km(2,:,:)  = bin2dec_array( bins(5,3,:,:) );
%    qapp_1km(3,:,:)  = bin2dec_array( bins(4,3,:,:) );
%    qapp_1km(4,:,:)  = bin2dec_array( bins(3,3,:,:) );
%    qapp_1km(5,:,:)  = bin2dec_array( bins(1:2,3,:,:) );

   qapp_1km(1,:,:)  = bitand(bitshift(qa1(3,:,:), 0),7); %primary phase retrieval
   qapp_1km(2,:,:)  = bitand(bitshift(qa1(3,:,:),-3),1); %outcome of above
   qapp_1km(3,:,:)  = bitand(bitshift(qa1(3,:,:),-4),1); %Rayleigh correction
   qapp_1km(4,:,:)  = bitand(bitshift(qa1(3,:,:),-5),1); %Atmospheric correction     
   qapp_1km(5,:,:)  = bitand(bitshift(qa1(3,:,:),-6),3); %Band used for optical thickness retrieval  
   
   

%===========================================================================
% RETRIEVAL PROCESSING FLAGS - INPUT DATA RESOURCE FLAGS for 1 km (BYTES 3,4)
%
   
%   qapp_1km(5)  : TOTAL PRECIPITABLE WATER
%			0: NCEP GDAS% 1: DAO% 2: MOD07-IR% 3: cld mask clear/undet
%   qapp_1km(6)  : CLOUD TOP HEIGHT
%			0: MOD06% 1: DAO% 2: other% 3: cld mask clear/undet
%   qapp_1km(7)  : TEMPERATURE PROFILE
%			0: NCEP% 1: DAO% 2: AIRS/AMSU/HSB% 3: cld mask clear/undet
%   qapp_1km(8)  : SURFACE TEMPERATURE OVER LAND
%			0: NCEP% 1: DAO% 2: MOD11% 3: cld mask clear/undet
%   qapp_1km(9)  : SURFACE TEMPERATURE OVER OCEAN
%			0: Reynolds% 1: DAO% 2:MOD28% 3: cld mask clear/undet
%   qapp_1km(10) : BRDF/Albedo
%			0: CERES/SARB% 1: DAO% 2: MOD43% 3: cld mask clear/undet
%   qapp_1km(15) : OZONE PROFILE
%			0: TOMS% 1: TOVS% 2: DAO% 3: cld mask clear/undet
%
%---------------------------------------------------------------------------
  
%    qapp_1km(9,*,*)=reform(byte( ((qa1(3,*,*) and bts(0))/bts(0)) +$
%                		        ((qa1(3,*,*) and bts(1))/bts(1))*2))
%    qapp_1km(10,*,*)=reform(byte(((qa1(3,*,*) and bts(2))/bts(2))+$
%                		((qa1(3,*,*) and bts(3))/bts(3))*2))                   
%    qapp_1km(11,*,*)=reform(byte(((qa1(3,*,*) and bts(4))/bts(4))+$
%                		((qa1(3,*,*) and bts(5))/bts(5))*2))
%    qapp_1km(12,*,*)=reform(byte(((qa1(3,*,*) and bts(6))/bts(6))+$
%                		((qa1(3,*,*) and bts(7))/bts(7))*2))
%    qapp_1km(13,*,*)=reform(byte(((qa1(4,*,*) and bts(0))/bts(0))+$
%                		((qa1(4,*,*) and bts(1))/bts(1))*2))
%    qapp_1km(14,*,*)=reform(byte(((qa1(4,*,*) and bts(2))/bts(2))+$
%                		((qa1(4,*,*) and bts(3))/bts(3))*2))                   
%    qapp_1km(15,*,*)=reform(byte(((qa1(4,*,*) and bts(4))/bts(4))+$
%                		((qa1(4,*,*) and bts(4))/bts(4))*2))               %MISTAKE?
   

%                     
%    qapp_1km(10,:,:)  = bin2dec_array( bins(7:8,4,:,:) );
%    qapp_1km(11,:,:)  = bin2dec_array( bins(5:6,4,:,:) );
%    qapp_1km(12,:,:)  = bin2dec_array( bins(3:4,4,:,:) );
%    qapp_1km(13,:,:)  = bin2dec_array( bins(1:2,4,:,:) );
%    qapp_1km(14,:,:)  = bin2dec_array( bins(7:8,5,:,:) );
%    qapp_1km(15,:,:)  = bin2dec_array( bins(5:6,5,:,:) );
%    qapp_1km(16,:,:)  = bin2dec_array( bins(3:4,5,:,:) );   %?? check ??
   
   %byte 4
   qapp_1km(10,:,:)  = bitand(bitshift(qa1(4,:,:), 0),1); %1.6/2.1 tau usefulness
   qapp_1km(11,:,:)  = bitand(bitshift(qa1(4,:,:),-1),3); %1.6/2.1 tau confidence
   qapp_1km(12,:,:)  = bitand(bitshift(qa1(4,:,:),-3),1); %1.6/2.1 re usefullness
   qapp_1km(13,:,:)  = bitand(bitshift(qa1(4,:,:),-4),3); %1.6/2.1 re confidence
   clear_sky_restoral(1,:,:) = bitand(bitshift(qa1(4,:,:),-6),3); %clear-sky restoral
   
   %byte 5
   qapp_1km(14,:,:)  = bitand(bitshift(qa1(5,:,:), 0),1);  %1.6/2.1 LWP usefulness 
   qapp_1km(15,:,:)  = bitand(bitshift(qa1(5,:,:),-1),3);  %1.6/2.1 LWP confidence    
   qapp_1km(16,:,:)  = bitand(bitshift(qa1(5,:,:),-3),7);  %Multi-layer flag 
   qapp_1km(17,:,:)  = bitand(bitshift(qa1(5,:,:),-6),1);  %Primary cloud retrieval outcome  (duplication)  
   
%===========================================================================   

   
               
           
  

end




