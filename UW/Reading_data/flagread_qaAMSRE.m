function [qaAMSRE] = flagread_qaAMSRE(qa1,all_flags)

%==========================================================================
% VERSION	:	1.0
% AUTHOR	:	Daniel Grosvenor
% DATE		:	7th Sep, 2013.
%
% DESCRIPTION	:	Converts the byte format of AMSRE qualtiy assurance
% array into useful format (normal numbers)
%============================================================================
   
%    qaAMSRE=zeros([9 nx ny]);
%    qapp_1km=zeros([16 nx ny]);


%===========================================================================
%  QUALITY ASSURANCE FLAGS 1st BYTE
%  INFORMATION: 
%
%  (0) : OPTICAL THICKNESS RETRIEVAL 
%                    0: not useful% 1 useful
%  qaAMSRE(1) : OPTICAL THICKNESS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qaAMSRE(2) : OPTICAL THICKNESS OUT OF BOUNDS 
%                    0: within bounds (tau<150)
%                    1: (100<tau<150)
%                    2: (tau>150)
%                    3: surface reflectance too large
%  qaAMSRE(3) : EFFECTIVE RADIUS GENERAL QA
%		     0: not useful
% 		     1: useful
%  qaAMSRE(4) : EFFECTIVE RADIUS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good			
%  qaAMSRE(5) : WATER PATH GENERAL QA
%		     0: not useful
% 		     1: useful
%  qaAMSRE(6) : WATER PATH CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qaAMSRE(7) : CLOUD PHASE DETERMINATION
%		     0: SWIR algorithm not run
%		     1: CLEAR		
%		     2: WATER
%		     3: ICE
%		     4: MIXED PHASE or UNDETERMINED	

%  qaAMSRE(8) : 1621 cloud retrieval outcome
%		     0: not attempted or unsuccessful
%		     1: successful		

%  
%---------------------------------------------------------------------------

   qaAMSRE(1,:,:)  = bitand(bitshift(qa1(1,:,:), 0),1); %Sea-ice possibility based on climo
   qaAMSRE(2,:,:)  = bitand(bitshift(qa1(1,:,:),-1),1); %Sea-ice possibility based on Tb
   qaAMSRE(3,:,:)  = bitand(bitshift(qa1(1,:,:),-2),1); %Vlow res Tb out of bounds or missing
   qaAMSRE(4,:,:)  = bitand(bitshift(qa1(1,:,:),-3),1); %Low res Tb out of bounds or missing
   qaAMSRE(5,:,:)  = bitand(bitshift(qa1(1,:,:),-4),1); %Med res Tb out of bounds or missing
   qaAMSRE(6,:,:)  = bitand(bitshift(qa1(1,:,:),-5),1); %High res Tb out of bounds or missing
   %bits 6&7 not used
   
   %Byte 2 - data acceptability for different resolutions - 0=normal, 1=out
   %                                              of bounds, 2=no retrieval
   qaAMSRE(7,:,:)  = bitand(bitshift(qa1(2,:,:), 0),3); %Very low resolution 
   qaAMSRE(8,:,:)  = bitand(bitshift(qa1(2,:,:),-2),3); %Low resolution 
   qaAMSRE(9,:,:)  = bitand(bitshift(qa1(2,:,:),-4),3); %Medium resolution 
   qaAMSRE(10,:,:)  = bitand(bitshift(qa1(2,:,:),-6),3); %High resolution 
