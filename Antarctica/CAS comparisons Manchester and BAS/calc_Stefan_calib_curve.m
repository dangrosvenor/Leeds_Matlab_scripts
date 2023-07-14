function [latex0u25_17,water_stefan,water2u_105]=calc_Stefan_calib_curve

%read the latex0u25_17 file (presumably used by Stefan? - check this).
dir_calib='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\from Stefan Benz\';
file_calib = 'latex0u25_17.txt';            
latex0u25_17=readcalib([dir_calib file_calib]);

%read the WELAS_Wasser_kalib_Dan.TXT file created by Stefan for each channel
dir_calib='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons\100729_cas_intercomparison_data\welas\from Stefan Benz\';
file_calib = 'WELAS_Wasser_kalib_Dan.TXT';
water_stefan=readcalib([dir_calib file_calib]);

%read the latex2u_105.txt file used in our 29th July Welas experiments
dir_calib='Y:\BAS_flights\29thJuly2010_BAS_chamber_comparisons/';
file_calib = 'latex2u_105.txt';            
latex2u_105=readcalib([dir_calib file_calib]);

water2u_105 = interp1(latex0u25_17,water_stefan,latex2u_105)';


function dat=readcalib(filename)
fid=fopen(filename,'r');
for i=1:8
    fgetl(fid);
end
dat=fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
dat=dat(2,:);