
function [mask_5km,qapq_5km,qapp_5km] = flagread_C6_5km_Dan(mask5_in,qa5)
%For C6 haven't changed any of the reading code - jsut listed the extra
%valeus defined for cirrus and high cloud flags

%============================================================================
% PROCEDURE	:	flagread_5km
% VERSION	:	1.0
% AUTHOR	:	Robert Wood - adapted for Matlab by Daniel Grosvenor
% DATE		:	March 19 2001
%
% DESCRIPTION	:	Takes HDF buffers for HDF variable IDs 
%			'Cloud_Mask_5km' (mask5) and 
%			'Quality_Assurance_5km' (qa5) 
%			and extracts flags into useable format
%
% NOTE that the convention for the MODIS QA is that bits run from left to
% right in a byte - i.e. the most significant bit is on the right. Matlab
% does it the other way around so have to use the bits from right to left.
%============================================================================
% BITS MASK

%   bts=2^indgen(8)  % indgen just creates an array of [0:n] (IDL uses zero
%   based indices)
%    bts=2.^[0:7]; % 0,2,4,8, etc.


% CHECK DIMENSIONS of ARRAY
%   nx=n_elements(qa5(0,*,0))
%   ny=n_elements(qa5(0,0,*))

%For C6 mask5_in is of size [2 X Y] - i.e. 2 bytes instead of one
mask5=squeeze(mask5_in(1,:,:)); %1st byte


nx = size(mask5,1);
ny = size(mask5,2);
%   nx=length(qa5(1,:,1));
 %  ny=length(qa5(1,1,:));
   

% OUTPUTS
%    mask_5km=bytarr(6,nx,ny)
%    qapq_5km=bytarr(10,nx,ny)
%    qapp_5km=bytarr(18,nx,ny)  
%these create byte arrays of zeros - I think that zeros should be the
%Matlab equivalent

    mask_5km=zeros([10 nx ny]);  %Are 10 pieces of info in C6 vs 6 for C5.1      
    qapq_5km=zeros([10 nx ny]);
    qapp_5km=zeros([18 nx ny]);

        
%===========================================================================
% CLOUD MASK FLAGS for 5 km DATA
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html  
%   
%  mask_5km(0)  :  CLOUD MASK 
%			0: undetermined% 1: determined
%  mask_5km(1)  :  CLOUD MASK QUALITY FLAG
%			0: 0-20% cloudy pixels% 1: 20-40%% 2: 40-60%% 3: 60-100%
%  mask_5km(2)  :  DAY/NIGHT FLAG
%			0: night% 1: day
%  mask_5km(3)  :  SUN GLINT FLAG
%			0: yes% 1: no
%  mask_5km(4)  :  SNOW/ICE FLAG
%			0: yes% 1: no
%  mask_5km(5)  :  LAND/WATER FLAG
%			0: water (ocean)% 1: coastal% 2: desert% 3: land
%
%---------------------------------------------------------------------------   
        
%    mask_5km(0,*,*)=byte(mask5 and bts(0))
%    mask_5km(1,*,*)=byte(((mask5 and bts(1))/bts(1))+((mask5 and bts(2))/bts(2))*2)
%    mask_5km(2,*,*)=byte((mask5 and bts(3))/bts(3))
%    mask_5km(3,*,*)=byte((mask5 and bts(4))/bts(4))
%    mask_5km(4,*,*)=byte((mask5 and bts(5))/bts(5))  
%    mask_5km(5,*,*)=byte(((mask5 and bts(6))/bts(6))+((mask5 and bts(7))/bts(7))*2)    

% tic
%    bins=dec2bin(mask5,8); %dec2bin converts to an N by 8 array of the bits (in string format)
%    %NOTE - Matlab bits are listed with the most signficant bit on the left
%    %(e.g. 100 = 4). The MODIS QA document (and Rob's script - IDL) number
%    %these bits so that 0 is the least significant and 7 the most
%    %So here I use e.g. indices 7&8 in the array of bit characters for bit
%    %nunbers 0 and 1
%    
%    mask_5km(1,:,:)=reshape( bin2dec(bins(:,8)) , size(mask5) );
%    mask_5km(2,:,:)=reshape( bin2dec(bins(:,6:7)) , size(mask5) );
%    mask_5km(3,:,:)=reshape( bin2dec(bins(:,5)) , size(mask5) );
%    mask_5km(4,:,:)=reshape( bin2dec(bins(:,4)) , size(mask5) );
%    mask_5km(5,:,:)=reshape( bin2dec(bins(:,3)) , size(mask5) );
%    mask_5km(6,:,:)=reshape( bin2dec(bins(:,1:2)) , size(mask5) ); 
%    
%    toc
   
   
   
   %tic - this is about 40x faster than the above dec2bin method!
   
   mask_5km(1,:,:) = bitand(mask5,1);  
   mask_5km(2,:,:) = bitand(bitshift(mask5,-1),3); %3 in binary = 00000011  - so when do bitand are stipping out the rightmost 
                                          %two bits. Then shift the bits two
                                          %places to the right so that we
                                          %remove those right bits
                                          %that have just been processed
                                          %ready to strip the next bits to
                                          %the left
   mask_5km(3,:,:) = bitand(bitshift(mask5,-3),1);
   mask_5km(4,:,:) = bitand(bitshift(mask5,-4),1);
   mask_5km(5,:,:) = bitand(bitshift(mask5,-5),1);
   mask_5km(6,:,:) = bitand(bitshift(mask5,-6),3);
   
   %toc
   
   %% 2nd byte
   mask5=squeeze(mask5_in(2,:,:)); %2nd byte
   
   % Bits 1-2 :- C6 sunglint flag, 0 = Fill OR CTP retreival, 1=no sunglint & CTP retrieval successs, 2=sunglint &
   % CTP retrieval successs
   % Bits 3-4 :- C6 snow/ice flag, 0=Fill OR CTP retreival fail, 1 = no
   % snow/ice & CTP retrieval success, 2 = snow/ice & CTP retrieval success
   % Bits 5-7 :- C6 surface type flag, 0=Fill OR CTP retreival fail, 1 =
   % ocean, deep lakes and rivers & CTP retrieval success, 2 = coast, shallow lakes and rivers & CTP retrieval success   
   % 3 = desert & CTP retrieval success,4= land & CTP retrieval success,
   % 5 = all other valid (non-fill) surface types& CTP retrieval success
   % Bits 8 - Day/Night flag, 0= night (or fill if status flag=0), 1=Day
   
   mask_5km(7,:,:) = bitand(mask5,3);  
   mask_5km(8,:,:) = bitand(bitshift(mask5,-2),3); 
   mask_5km(9,:,:) = bitand(bitshift(mask5,-4),7);
   mask_5km(10,:,:) = bitand(bitshift(mask5,-7),1);

   
   

% **************************************************************************
   % stop here if we don't have qa5 (as for Joint L2 5 km files)
% **************************************************************************   
   if nargin>=2
            

%qa5 contains 10 bytes of bits for the QA flags - here they are processed
%in stages
%===========================================================================
% QUALITY ASSURANCE FLAGS for 5 km DATA (BYTES 0 to 2.5)
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html
%
%  qapq_5km(0) : CLOUD TOP PRESSURE 
%                    0: not useful% 1 useful
%  qapq_5km(1) : CLOUD TOP PRESSURE CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_5km(2) : CLOUD TOP TEMPERATURE 
%                    0: not useful% 1 useful
%  qapq_5km(3) : CLOUD TOP TEMPERATURE CONFIDENCE
%		     0: bad% 1 marginal% 2 good% 3 very good 
%  qapq_5km(4) : CLOUD FRACTION
%		     0: not useful% 1 useful
%  qapq_5km(5) : CLOUD FRACTION CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good			
%  qapq_5km(6) : CLOUD EFFECTIVE EMISSIVITY
%		     0: not useful% 1 useful
%  qapq_5km(7) : CLOUD EFFECTIVE EMISSIVITY CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_5km(8) : CLOUD PHASE INFRARED
%	     	     0: not useful% 1 useful
%  qapq_5km(9): CLOUD PHASE INFRARED CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good


%    bins=dec2bin(qa5,8); %dec2bin converts to an N by 8 array of the bits (in string format)
%   
%    
%    bins = reshape(bins,[size(qa5) 8]);
%    Lsiz = length(size(bins));
%    bins = permute(bins,[Lsiz 1:Lsiz-1]); %put the bits as the first dimension as makes life easier
   %becuase of way Matlab rearranges the dimensions when an array has singular dimensions

   %now we can use the bits using e.g. bin2dec_array(bins(2:4,1,:)); for
   %bits 2,3 and 4 (where 1 is the first bit) of the first byte of the 10
   
%    qapq_5km(1,:,:)  = bin2dec_array( bins(8,1,:,:) );
%    qapq_5km(2,:,:)  = bin2dec_array( bins(5:7,1,:,:) );
%    qapq_5km(3,:,:)  = bin2dec_array( bins(4,1,:,:) );
%    qapq_5km(4,:,:)  = bin2dec_array( bins(1:3,1,:,:) );
%    qapq_5km(5,:,:)  = bin2dec_array( bins(8,2,:,:) );
%    qapq_5km(6,:,:)  = bin2dec_array( bins(5:7,2,:,:) );
%    qapq_5km(7,:,:)  = bin2dec_array( bins(4,2,:,:) );
%    qapq_5km(8,:,:)  = bin2dec_array( bins(1:3,2,:,:) );
%    qapq_5km(9,:,:)  = bin2dec_array( bins(8,3,:,:) );
%    qapq_5km(10,:,:) = bin2dec_array( bins(5:7,3,:,:) );
%    
   
   
   
   qapq_5km(1,:,:)  = bitand( qa5(1,:,:),1 );
   qapq_5km(2,:,:)  = bitand( bitshift(qa5(1,:,:),-1) , 7 ); % 7 = 0000111 in binary. i.e. stripping 3 bits
   qapq_5km(3,:,:)  = bitand( bitshift(qa5(1,:,:),-4) , 1 ); % -4 here since we are on the 4th bit here (stripped 1 and 3 bits so far)
   qapq_5km(4,:,:)  = bitand( bitshift(qa5(1,:,:),-5) , 7 );
   qapq_5km(5,:,:)  = bitand( bitshift(qa5(2,:,:), 0) , 1 );
   qapq_5km(6,:,:)  = bitand( bitshift(qa5(2,:,:),-1) , 7 );
   qapq_5km(7,:,:)  = bitand( bitshift(qa5(2,:,:),-4) , 1 );
   qapq_5km(8,:,:)  = bitand( bitshift(qa5(2,:,:),-5) , 7 );
   qapq_5km(9,:,:)  = bitand( bitshift(qa5(3,:,:), 0) , 1 );
   qapq_5km(10,:,:) = bitand( bitshift(qa5(3,:,:),-1) , 7 );

  
   
   
   
   
%    qapq_5km(0,*,*)=reform(byte(qa5(0,*,*) and bts(0))) 
%    qapq_5km(1,*,*)=reform(byte(((qa5(0,*,*) and bts(1))/bts(1))+$
%                ((qa5(0,*,*) and bts(2))/bts(2))*2 + $
%                ((qa5(0,*,*) and bts(3))/bts(3))*4 ))   
%    qapq_5km(2,*,*)=reform(byte(qa5(0,*,*) and bts(4))) 
%    qapq_5km(3,*,*)=reform(byte(((qa5(0,*,*) and bts(5))/bts(5))+$
%                ((qa5(0,*,*) and bts(6))/bts(6))*2 + $
%                ((qa5(0,*,*) and bts(7))/bts(7))*4 ))               
%    qapq_5km(4,*,*)=reform(byte(qa5(1,*,*) and bts(0))) 
%    qapq_5km(5,*,*)=reform(byte(((qa5(1,*,*) and bts(1))/bts(1))+$
%                ((qa5(1,*,*) and bts(2))/bts(2))*2 + $
%                ((qa5(1,*,*) and bts(3))/bts(3))*4 ))   
%    qapq_5km(6,*,*)=reform(byte(qa5(1,*,*) and bts(4))) 
%    qapq_5km(7,*,*)=reform(byte(((qa5(1,*,*) and bts(5))/bts(5))+$
%                ((qa5(1,*,*) and bts(6))/bts(6))*2 + $
%                ((qa5(1,*,*) and bts(7))/bts(7))*4 ))               
%    qapq_5km(8,*,*)=reform(byte(qa5(2,*,*) and bts(0))) 
%    qapq_5km(9,*,*)=reform(byte(((qa5(2,*,*) and bts(1))/bts(1))+$
%                ((qa5(2,*,*) and bts(2))/bts(2))*2 + $
%                ((qa5(2,*,*) and bts(3))/bts(3))*4 ))



%===========================================================================
% QUALITY ASSURANCE FLAGS for 5 km DATA (BYTES 2.5-6.75)               
% RETRIEVAL PROCESSING FLAGS - PROCESSING PATH FLAGS
%
%  qapp_5km(0)  :  CIRRUS LEVEL 3 FLAG  
% 			0: missing% 1: no cirrus found% 2: cirrus found 3:clear-sky
% 			(new for C6)
%  qapp_5km(1)  :  HIGH CLOUD LEVEL 3 FLAG
%			0: missing% 1: no high cld found% 2: high cld found 3:clear-sky
% 			(new for C6)
%  qapp_5km(2)  :  NUMBER of CLOUDY PIXELS WITHIN 5x5 KM BOX
%  qapp_5km(3)  :  NUMBER of CLEAR PIXELS WITHIN 5x5 KM BOX
%  qapp_5km(4)  :  NUMBER of MISSING PIXELS WITHIN 5x5 KM BOX
%  qapp_5km(5)  :  MAXIMUM LIKELYHOOD ESTIMATOR
%			0: not used% 1 used
%  qapp_5km(6)  :  CLUSTER ANALYSIS
%			0: not used% 1 used   
%  qapp_5km(7)  :  GOODNESS of FIT
%			0: 0<1% 1: 0=1     
%  qapp_5km(8)  :  CHI-SQUARED
%			0: < npts used in MLE
%			1: > npts used in MLE
%---------------------------------------------------------------------------  
   
   
%    qapp_5km(1,:,:)  = bin2dec_array( bins(3:4,3,:,:) );
%    qapp_5km(2,:,:)  = bin2dec_array( bins(1:2,3,:,:) );
%    qapp_5km(3,:,:)  = bin2dec_array( bins(1:8,4,:,:) );
%    qapp_5km(4,:,:)  = bin2dec_array( bins(1:8,5,:,:) );
%    qapp_5km(5,:,:)  = bin2dec_array( bins(1:8,6,:,:) );
%    qapp_5km(6,:,:)  = bin2dec_array( bins(8,7,:,:) );
%    qapp_5km(7,:,:)  = bin2dec_array( bins(7,7,:,:) );
%    qapp_5km(8,:,:)  = bin2dec_array( bins(6,7,:,:) );
%    qapp_5km(9,:,:)  = bin2dec_array( bins(5,7,:,:) );
   
   qapp_5km(1,:,:)  = bitand( bitshift(qa5(3,:,:),-4) , 3 );
   qapp_5km(2,:,:)  = bitand( bitshift(qa5(3,:,:),-6) , 3 );
   qapp_5km(3,:,:)  = qa5(4,:,:);  %just use the straight integers for these three
   qapp_5km(4,:,:)  = qa5(5,:,:);
   qapp_5km(5,:,:)  = qa5(6,:,:);
   qapp_5km(6,:,:)  = bitand( bitshift(qa5(7,:,:), 0) , 1 );
   qapp_5km(7,:,:)  = bitand( bitshift(qa5(7,:,:),-1) , 1 );
   qapp_5km(8,:,:)  = bitand( bitshift(qa5(7,:,:),-2) , 1 );
   qapp_5km(9,:,:)  = bitand( bitshift(qa5(7,:,:),-3) , 1 );   
   
   %Rob's IDL code
%    qapp_5km(0,*,*)=reform(byte(((qa5(2,*,*) and bts(4))/bts(4))+$
%                ((qa5(2,*,*) and bts(5))/bts(5))*2))
%    qapp_5km(1,*,*)=reform(byte(((qa5(2,*,*) and bts(6))/bts(6))+$
%                ((qa5(2,*,*) and bts(7))/bts(7))*2))                   
%    qapp_5km(2,*,*)=qa5(3,*,*)
%    qapp_5km(3,*,*)=qa5(4,*,*)
%    qapp_5km(4,*,*)=qa5(5,*,*)
%    qapp_5km(5,*,*)=reform(byte(qa5(6,*,*) and bts(0)))
%    qapp_5km(6,*,*)=reform(byte(qa5(6,*,*) and bts(1))/bts(1))
%    qapp_5km(7,*,*)=reform(byte(qa5(6,*,*) and bts(2))/bts(2))
%    qapp_5km(8,*,*)=reform(byte(qa5(6,*,*) and bts(3))/bts(3))
   
%===========================================================================   
   
              
%===========================================================================
% QUALITY ASSURANCE FLAGS for 5 km DATA (BYTES 6.75-9)               
% RETRIEVAL PROCESSING FLAGS - INPUT DATA RESOURCE FLAGS
%
   
%   qapp_5km(9)  : CLEAR RADIANCE ORIGIN
%			0: MOD35% 1: Forward% 2: other
%   qapp_5km(10) : MOISTURE PROFILE
%			0: NCEP% 1: DAO% 2: AIRS/AMSU% 3: other
%   qapp_5km(11) : TEMPERATURE PROFILE
%			0: NCEP% 1: DAO% 2: AIRS/AMSU% 3: other
%   qapp_5km(12) : SURFACE TEMPERATURE OVER LAND
%			0: NCEP% 1: DAO% 2: AIRS/AMSU% 3: other
%   qapp_5km(13) : SURFACE TEMPERATURE OVER OCEAN
%			0: Reynolds% 1: DAO% 2:other% 3: not used
%   qapp_5km(14) : SURFACE PRESSURE
%			0: NCEP% 2: DAO% 3: other% 4: not used
%   qapp_5km(15) : TOPOGRAPHY
%			0: EOS DEM% 1: Other
%   qapp_5km(16) : SURFACE EMISSIVITY
%			0: CERES% 1: MOD11	
%   qapp_5km(17) : SURFACE TYPE
%			0: LOVELAND% 1: Olson% 2: MOD12% 3: Other
%
%---------------------------------------------------------------------------

%    qapp_5km(10,:,:)  = bin2dec_array( bins(1:2,7,:,:) );
%    qapp_5km(11,:,:)  = bin2dec_array( bins(7:8,8,:,:) );
%    qapp_5km(12,:,:)  = bin2dec_array( bins(5:6,8,:,:) );
%    qapp_5km(13,:,:)  = bin2dec_array( bins(3:4,8,:,:) );
%    qapp_5km(14,:,:)  = bin2dec_array( bins(1:2,8,:,:) );
%    qapp_5km(15,:,:)  = bin2dec_array( bins(7:8,9,:,:) );
%    qapp_5km(16,:,:)  = bin2dec_array( bins(5:6,9,:,:) );
%    qapp_5km(17,:,:)  = bin2dec_array( bins(3:4,9,:,:) );
%    qapp_5km(18,:,:)  = bin2dec_array( bins(1:2,9,:,:) );
   
   
   qapp_5km(10,:,:)  = bitand( bitshift(qa5(7,:,:),-6) , 3 );  
   qapp_5km(11,:,:)  = bitand( bitshift(qa5(8,:,:), 0) , 3 );
   qapp_5km(12,:,:)  = bitand( bitshift(qa5(8,:,:),-2) , 3 );
   qapp_5km(13,:,:)  = bitand( bitshift(qa5(8,:,:),-4) , 3 );   
   qapp_5km(14,:,:)  = bitand( bitshift(qa5(8,:,:),-6) , 3 );  
   qapp_5km(15,:,:)  = bitand( bitshift(qa5(9,:,:), 0) , 3 );
   qapp_5km(16,:,:)  = bitand( bitshift(qa5(9,:,:),-2) , 3 );
   qapp_5km(17,:,:)  = bitand( bitshift(qa5(9,:,:),-4) , 3 );   
   qapp_5km(18,:,:)  = bitand( bitshift(qa5(9,:,:),-6) , 3 );  
      
   
%   Rob's IDL code   - bts(4) and bts(5) not used? Spare bits?
%   
%    qapp_5km(9,*,*)=reform(byte(((qa5(6,*,*) and bts(6))/bts(6))+$
%                ((qa5(6,*,*) and bts(7))/bts(7))*2))
%    qapp_5km(10,*,*)=reform(byte(((qa5(7,*,*) and bts(0))/bts(0))+$
%                ((qa5(7,*,*) and bts(1))/bts(1))*2))                   
%    qapp_5km(11,*,*)=reform(byte(((qa5(7,*,*) and bts(2))/bts(2))+$
%                ((qa5(7,*,*) and bts(3))/bts(3))*2))
%    qapp_5km(12,*,*)=reform(byte(((qa5(7,*,*) and bts(4))/bts(4))+$
%                ((qa5(7,*,*) and bts(5))/bts(5))*2))
%    qapp_5km(13,*,*)=reform(byte(((qa5(7,*,*) and bts(6))/bts(6))+$
%                ((qa5(7,*,*) and bts(7))/bts(7))*2))
%    qapp_5km(14,*,*)=reform(byte(((qa5(8,*,*) and bts(0))/bts(0))+$
%                ((qa5(8,*,*) and bts(1))/bts(1))*2))                   
%    qapp_5km(15,*,*)=reform(byte(((qa5(8,*,*) and bts(2))/bts(2))+$
%                ((qa5(8,*,*) and bts(3))/bts(3))*2))
%    qapp_5km(16,*,*)=reform(byte(((qa5(8,*,*) and bts(4))/bts(4))+$
%                ((qa5(8,*,*) and bts(5))/bts(5))*2))
%    qapp_5km(17,*,*)=reform(byte(((qa5(8,*,*) and bts(6))/bts(6))+$
%                ((qa5(8,*,*) and bts(7))/bts(7))*2))              
%                
%                
%===========================================================================   


   end
   
               
               



