%read in the GOES Tau-Re data from Grant's (I think?) .nc files
%at the moment this just loads in a single file

% N.B. - regarding netCDF reading - if fails to work if using a differente
% version of Matlab then might need to add e.g.
% /usr/local/MATLAB/R2013a/toolbox/mexdcf or similar to the path

% See --- read_GOES_vocals_netcdf_files_MULTI.m ---
% for reading in muliple files and savign the data in a .mat file

%Will want to make it load in all the files and process them into a 1x1
%degree product with CF included. For the CF I guess I would use the mask array.
%Looks like 1 for cloud and 0 for no cloud and -9 for NaN I guess?
%And also to calculate the SZA for each time using the sun_pos function.
%Also have the Phase field - has numbers 1-5. Ask Gallia or Grant about
%what they mean.

%data is present from 8th Nov to 2nd Dec, 2008
%goes_dir = '/home/disk/eos4/gallia/Minnis Cloud Retrievals/';
% Think the above is no longer available (been archived to tape).
%Now have the data below for all of Oct 2008

try

goes_dir = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/';

if ~exist('ioverride_goes') | ioverride_goes==0

    goes_action = 'load all and process in L3';
    goes_action = 'load a particular file';
    %goes_action = 'load one-by-one and plot';

    %Data is every 30mins at 15 ann 45 past the hour

    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261145.nc';
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810260845.nc';
    %goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810262345.nc';
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810270115.nc';
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc';
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131715.nc';    
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131745.nc';        
    %Can find on Manchester FTP - /DATAbase3/VOCALS/Data/satellite_products/cloud-products/

    ifile_inds=[1:2];

end

clear file_list
switch goes_action
    case {'load all and process in L3','load one-by-one and plot'}
        file_list = dir([goes_dir '*.nc']);

    case 'load a particular file'
        file_list(1).name = goes_file_to_load;
        ifile_inds=[1];
end

clear var_list
ivar=1;

switch goes_action 
    case 'load one-by-one and plot'
        var_list{ivar} = 'Co-Latitude'; ivar=ivar+1;
        var_list{ivar} = 'Co-Longitude'; ivar=ivar+1;
        var_list{ivar} = 'Tau'; ivar=ivar+1;
        var_list{ivar} = 'Reff-Deff'; ivar=ivar+1;
        var_list{ivar} = 'Teff'; ivar=ivar+1;
%        var_list{ivar} = 'BTCH4'; ivar=ivar+1;
        var_list{ivar} = 'mask'; ivar=ivar+1; 
        var_list{ivar} = 'Phase'; ivar=ivar+1;
        
        
    otherwise
var_list{ivar} = 'Co-Latitude'; ivar=ivar+1;
var_list{ivar} = 'Co-Longitude'; ivar=ivar+1;
var_list{ivar} = 'Ztop'; ivar=ivar+1;
var_list{ivar} = 'Vis'; ivar=ivar+1;
var_list{ivar} = 'NIR'; ivar=ivar+1;
var_list{ivar} = 'BTCH2'; ivar=ivar+1;
var_list{ivar} = 'BTCH4'; ivar=ivar+1;
var_list{ivar} = 'CTH'; ivar=ivar+1;
var_list{ivar} = 'BBA'; ivar=ivar+1;
var_list{ivar} = 'BBOLR'; ivar=ivar+1;
var_list{ivar} = 'IRemit'; ivar=ivar+1;
var_list{ivar} = 'mask'; ivar=ivar+1;
var_list{ivar} = 'Phase'; ivar=ivar+1;
var_list{ivar} = 'Tau'; ivar=ivar+1;
var_list{ivar} = 'Reff-Deff'; ivar=ivar+1;
var_list{ivar} = 'LWP'; ivar=ivar+1;
var_list{ivar} = 'Teff'; ivar=ivar+1;
var_list{ivar} = 'Ptop'; ivar=ivar+1;
var_list{ivar} = 'Peff'; ivar=ivar+1;
var_list{ivar} = 'Pbot'; ivar=ivar+1;
var_list{ivar} = 'Zeff'; ivar=ivar+1;
var_list{ivar} = 'Zbot'; ivar=ivar+1;

end


Ntimes = length(ifile_inds);


for ifile=ifile_inds  %length(file_list)

    goes_file = [goes_dir file_list(ifile).name];
    if exist(goes_file,'file')==0
        error(['*** ERROR - ' goes_file ' does not exist ***']);        
    end
    verstr=version;
    vers=str2num(verstr(1));
    if vers<8
            nc = netcdf(goes_file,'nowrite');            
%    else
%            nc = netcdf.open(goes_file,'nowrite');           
    end
    
    goes_year(ifile) = str2num(file_list(ifile).name(23:26));
    goes_month(ifile) = str2num(file_list(ifile).name(27:28));
    goes_day(ifile) = str2num(file_list(ifile).name(29:30));
    goes_hours(ifile) = str2num(file_list(ifile).name(31:32));
    goes_mins(ifile) = str2num(file_list(ifile).name(33:34));
    

    
    for ivar=1:length(var_list)
        varname = remove_character(var_list{ivar},'-','_');
%        eval(['goes_' varname ' = nc{var_list{ivar}}(:,:);']);

if vers<8
        dat = nc{var_list{ivar}}(:,:);
else
        dat = ncread(goes_file,var_list{ivar}); dat=double(dat);
end
        
%        dat(dat<-8) = NaN; %I think -9 is NaN for most data, but shoudl check

        
        switch goes_action
            case 'load one-by-one and plot'
                Ntimes = 1;

        end
        
         if ifile==1
             eval(['goes_' varname '=NaN*ones([size(dat,1) size(dat,2) Ntimes]);']);
         end
                
         eval(['goes_' varname '(:,:,ifile) = dat;']);

    end
    
%    icloud = find(goes_mask(:,:,ifile)~=1);
%    [I,J] = ind2sub(size(goes_mask(:,:,ifile)),icloud);
%    K = ifile*ones(size(I));
%    icloud2 = sub2ind(size(goes_mask),I,J,K);
%    goes_Tau(icloud2)=NaN;
end

%Older files seem to have this set to zero everywhere....
if maxALL(goes_mask)<0.01
    goes_mask(:)=1;    
end

icloud2 = find(goes_mask~=1 | goes_Tau<-8 | goes_Reff_Deff<-8);
goes_Tau(icloud2)=NaN;
goes_Reff_Deff(icloud2)=NaN;
goes_Teff(icloud2)=NaN;
goes_Phase(icloud2)=NaN;
goes_LWP(icloud2)=NaN;

goes_IRemit(goes_IRemit<-8.9)=NaN;
goes_Vis(goes_Vis<-7.9)=NaN;

goes_Reff = goes_Reff_Deff;
inotliq = find(goes_Phase~=1);
goes_Reff(inotliq) = NaN;

Plat_L2 = goes_Co_Latitude(:,:,1); 
Plon_L2 = goes_Co_Longitude(:,:,1);

%-9 seems to be a fill value for these files
Plat_L2(Plat_L2<-8)=NaN;
Plon_L2(Plon_L2<-8)=NaN;

Plat_L2 = 90 - Plat_L2;  %something is wrong with the latitudes - taking from 90 seems ot make them
%right??
Plon_L2 = Plon_L2 - 360;


%create the cell edges
Plat2_L22 = 0.5*(Plat_L2(:,1:end-1)+Plat_L2(:,2:end));
Plat2_L22 = 0.5*(Plat2_L22(1:end-1,:)+Plat2_L22(2:end,:));

Plon2_L22 = 0.5*(Plon_L2(1:end-1,:)+Plon_L2(2:end,:));
Plon2_L22 = 0.5*(Plon2_L22(:,1:end-1)+Plon2_L22(:,2:end));
%will put the very edge as NaN
Plat3_L2 = NaN*ones([size(Plat_L2,1)+1 size(Plat_L2,2)+1]);
Plat3_L2(2:end-1,2:end-1) = Plat2_L22;

Plon3_L2 = NaN*ones([size(Plon_L2,1)+1 size(Plon_L2,2)+1]);
Plon3_L2(2:end-1,2:end-1) = Plon2_L22;



gcm_Plat2D_edges_GOES = Plat3_L2;
gcm_Plon2D_edges_GOES = Plon3_L2;

gcm_Plat2D_GOES = Plat_L2;
gcm_Plon2D_GOES = Plon_L2;

%LAT_GOES = minALL(Plat2D_GOES_time3) :1: maxALL(Plat2D_GOES_time3);
%LON_GOES = minALL(Plon2D_GOES_time3) :1: maxALL(Plon2D_GOES_time3);

nT = size(goes_LWP,3);
daynum_timeseries3_GOES = goes_day(1:nT);
modisyear_timeseries3_GOES = goes_year(1:nT);
gcm_time_UTC_GOES = goes_hours(1:nT) + goes_mins(1:nT)/60;
gcm_time_matlab_GOES = datenum(goes_year(1:nT),goes_month(1:nT),goes_day(1:nT),goes_hours(1:nT),goes_mins(1:nT),0);

gcm_str = 'GOES';
gcm_str_select = 'GOES';

savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/GOES/';

disp('Done read GOES data');

%% ncdump of the output variables

% dimensions:
%         x = 1200 ;
%         y = 550 ;
% variables:
%         float Co-Latitude(y, x) ;
%                 Co-Latitude:units = "decimal degrees" ;
%                 Co-Latitude:long_name = "Mid-points of latitude bins" ;
%         float Co-Longitude(y, x) ;
%                 Co-Longitude:units = "decimal degrees" ;
%                 Co-Longitude:long_name = "Mid-points of longitude bins" ;
%         float Ztop(y, x) ;
%                 Ztop:units = "Kilometres" ;
%                 Ztop:long_name = "Cloud top height" ;
%         float Vis(y, x) ;
%                 Vis:units = "Unitless" ;
%                 Vis:long_name = "Visible reflectance" ;
%         float NIR(y, x) ;
%                 NIR:units = "Unitlesss" ;
%                 NIR:long_name = "1.6 micron reflectance" ;
%         float BTCH2(y, x) ;
%                 BTCH2:units = "Kelvin" ;
%                 BTCH2:long_name = "3.7 micron brightness temperature" ;
%         float BTCH4(y, x) ;
%                 BTCH4:units = "Kelvin" ;
%                 BTCH4:long_name = "10.8 micron brightness temperature" ;
%         float CTH(y, x) ;
%                 CTH:units = "Kelvin" ;
%                 CTH:long_name = "12 micron brightness temperature" ;
%         float BBA(y, x) ;
%                 BBA:units = "Percentage" ;
%                 BBA:long_name = "Broadband albedo" ;
%         float BBOLR(y, x) ;
%                 BBOLR:units = "Watts per square metre" ;
%                 BBOLR:long_name = "Broadband flux" ;
%         float IRemit(y, x) ;
%                 IRemit:units = "Unitless" ;
%                 IRemit:long_name = "IR emittance" ;
%         float mask(y, x) ;
%                 mask:units = "Unitless" ;
%                 mask:long_name = "Cloud Mask" ;
%         float Phase(y, x) ;
%                 Phase:units = "Unitless" ;
%                 Phase:long_name = "Cloud phase (1-5)" ;
%         float Tau(y, x) ;
%                 Tau:units = "Unitless" ;
%                 Tau:long_name = "Optical Depth" ;
%         float Reff-Deff(y, x) ;
%                 Reff-Deff:units = "Microns" ;
%                 Reff-Deff:long_name = "Cloud dropled effective radius / diameter of ice crystal" ;                                                                          
%         float LWP(y, x) ;
%                 LWP:units = "Grams per square metre" ;
%                 LWP:long_name = "Liquid (or ice) water path)" ;
%         float Teff(y, x) ;
%                 Teff:units = "Kelvins" ;
%                 Teff:long_name = "Effective temeprature" ;
%         float Ptop(y, x) ;
%                 Ptop:units = "Millibars" ;
%                 Ptop:long_name = "Cloud top pressure" ;
%         float Peff(y, x) ;
%                 Peff:units = "Millibarss" ;
%                 Peff:long_name = "Cloud effective pressure" ;
%         float Pbot(y, x) ;
%                 Pbot:units = "Millibars" ;
%                 Pbot:long_name = "Cloud bottom pressure" ;
%         float Zeff(y, x) ;
%                 Zeff:units = "Kilometres" ;
%                 Zeff:long_name = "Cloud effectiev height" ;
%         float Zbot(y, x) ;
%                 Zbot:units = "Kilometres" ;
%                 Zbot:long_name = "Cloud bottom height" ;
% 
% // global attributes:
%                 :Title = "Cloud property retrievals unpacked from Minnis data" ;
%                 :Data_contact = "Grant Allen (grant.allen@manchester.ac.uk) " ;
%                 :Campaign = "VOCALS, Arica, Chile Oct 2008 - Dec 2008" ;
%                 :Principal_Investigator = "Hugh Coe (hugh.coe@manchester.ac.uk)"
%  ;
%                 :History = "V1.0 Compiled: Mon Jun 15 20:32:36 2009 by Grant All
% en, Uni. Manchester, UK" ;
%                 :Conventions = "CF-1.0" ;
%                 :Institution = "University of Manchester, UK" ;
%                 :Comment = "This file represents cloud bulk properties as derive
% d from GOES10 radiances; contact grant.allen@manchester.ac.uk for further inform
% ation" ;
%                 :Location = "South East Pacific (SEP and South America)" ;
%                 :Mission_Purpose = "To understand the climate of the SEP and its
%  representation in models" ;
%                 :Data_Protocol = "Contact Grant Allen to discuss use of this dat
% aset" ;
%                 :Date = "20081202" ;
%                 :START_time = "200812021815 UTC" ;
%                 :END_time = "200812021815 UTC" ;
%                 :source = "GOES10 IR and Visible radiances" ;


clear ioverride_goes
catch goes_ERROR
        clear ioverride_goes
        rethrow(goes_ERROR);
end