%PRO read_amsr_day_v4, file_name, gmt, sst, wspd, vapor, cloud, rain

% this routine will read the AMSR-E daily bytemap files (version 4 released January 2005).
%
% arguments are:
%   file_name :  name of file to read complete with path
% 	file_names have form satname_yyyymmdd_v4.gz
%	 where satname  = name of satellite (amsre or amsr)
%	       yyyy		= year
%		   mm  		= month
%		   dd 		= day of month
%
% The routine returns:
%   gmt, sst, wspd, vapor, cloud, rain real arrays sized (1440,720,2)
%   gmt   is the mean gmt time in minutes of the observations within that grid cell
%   sst   is the sea surface temperature in degree Celcius, valid range=[-3.0,34.5]
%   wspd  is the 10 meter surface wind speed in m/s,  valid range=[0.,50.]
%   vapor is the columnar atmospheric water vapor in mm,  valid range=[0.,75.]
%   cloud is the liquid cloud water in mm, valid range = [0.,2.5]
%   rain  is the derived radiometer rain rate in mm/hr,  valid range = [0.,25.]
%
% Longitude  is 0.25*(xdim+1)-0.125		!IDL is zero based    East longitude
% Latitude   is 0.25*(ydim+1)-90.125
%
%
% Please read the data description on www.remss.com
% for infomation on the various fields
% To contact RSS support:
% http://www.remss.com/support


%filedir = '/home/disk/eos10/robwood/AMSR/bmaps_v05/y2007/m01/';
filedir = '/home/disk/eos5/d.grosvenor/AMSRE/y2002/m06/';
%file_name2 = 'amsre_200701v5.gz';
file_name2 = 'amsre_200701v7'; %these are monthly values - should be 6 maps - however they've averaged ascending and descending
file_name2 = 'amsre_20070131v5'; %thess are the daily maps (should be 7 maps - the extra one is time UTC,
%which is the first one)
%file_name2 = 'amsre_20070131v5_d3d'; %d3d are the 3 day averages
%file_name2 = 'amsre_20111001v7'; %weekly files - same as monthly (i.e. have also averaged day and night)
file_name2 = 'amsre_20020612v7';

file_name =[filedir file_name2];



if exist(file_name)~=2  %exist returns 2 when the argument is a file
    fprintf(1,'\nFILE DOES NOT EXIST!!');
    return
end

%unzip the file
%eval(['!gzip -d ' file_name]);


fid=fopen(file_name,'rb');
b=fread(fid);

  %open file, read binary data, close file
  close,2
  openr,2,file_name, error=err, /compress	%compress keyword allows reading of gzip file, remove if data already unzipped
  if (err gt 0) then begin
  	print, 'ERROR OPENING FILE: ', file_name
  endif else begin
  	readu,2,binarydata
  	close,2
  endelse

% multipliers to change binary data to real data
xscale=[6.,0.15,.2,.3,.01,.1]

% loop through asc/dsc  and all 6 variables
for iasc=0,1 do begin
   for ivar=0,5 do begin

		% extract 1 variable, scale and assign to real array
        dat=binarydata[*,*,ivar,iasc]

        case ivar of
        	0: BEGIN		%gmt time
        	     gmt[*,*,iasc] =dat*xscale[ivar]

               END

			1: BEGIN		%sea surface temperature
				ok=where(dat le 250)
				dat=float(dat)
				dat[ok]=dat[ok]*xscale[ivar]-3.0
                sst[*,*,iasc] =dat
   			   END

			2: BEGIN		%wind speed
				ok=where(dat le 250)
				dat=float(dat)
				dat[ok]=dat[ok]*xscale[ivar]
                wspd[*,*,iasc] =dat
   			   END

            3: BEGIN		%water vapor
                ok=where(dat le 250)
				dat=float(dat)
				dat[ok]=dat[ok]*xscale[ivar]
            	vapor[*,*,iasc]=dat
               END

			4: BEGIN		%cloud
			    ok=where(dat le 250)
				dat=float(dat)
				dat[ok]=dat[ok]*xscale[ivar]
				cloud[*,*,iasc]=dat
			   END

			5: BEGIN		%rain
			    ok=where(dat le 250)
				dat=float(dat)
				dat[ok]=dat[ok]*xscale[ivar]
				rain[*,*,iasc] =dat
			   END

         ENDCASE

	endfor	%ivar
endfor		%iasc

endelse


return
END
