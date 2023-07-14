%read a MODIS hdf file
%get data from http://ladsweb.nascom.nasa.gov/browse_images/l2_browser.html?form=AADS&browseType=Level+2
%note this file contains all the data products (not just e.g. cloud top temperature)
%can just double click on the file in Matlab to browse and import different fields
%is probably the easiest way

%INFO = hdfinfo('MOD06_L2.A2010039.1405.005.2010039225429.hdf');
%gets info about the file

%when have the data need to do double(data) when plotting as is in a strange format

%NOTE the data has strange offset and scale factor values
%e.g. for Cloud_Top_Temperature there is an offset of -15000 and an sf of 0.01
%so think need to be minus 15000 and multiply by 0.01
%then seems to be sensible if is in celsius
%also note that there is a fill value of -32768
CT=(double(Cloud_Top_Temperature));
CT(CT==-32768)=NaN;
CT=(CT-15000)*0.01 + 273.15;

CTP=(double(Cloud_Top_Pressure));
CTP(CTP==-32768)=NaN;
CTP=(CTP)*0.1; %zero offset 0.1 factor

CF=(double(Cloud_Fraction));
CF(CF==-32768)=NaN;
CF=(CF)*0.01; %zero offset 0.01 factor

CTP_IR=(double(Cloud_Top_Pressure_Infrared));
CTP_IR(CTP_IR==-32768)=NaN;
CTP_IR=(CTP_IR)*0.01; %zero offset 0.01 factor

CEFR=(double(Cloud_Effective_Radius));
CEFR(CEFR==-32768)=NaN;
CEFR=(CEFR)*0.1; %zero offset 0.01 factor (microns) -this one has different XY dimensions

filename='C:\Documents and Settings\dan\My Documents\MATLAB\WRF_toolbox_Dan\XY_latlon.mat';
load(filename);

X_LAT=interp2(LATS,LONS,X_latlon,Latitude(:),Longitude(:));
Y_LAT=interp2(LATS,LONS,Y_latlon,Latitude(:),Longitude(:));
XLAT=reshape(X_LAT,size(Latitude));
YLAT=reshape(Y_LAT,size(Latitude));
