%N.B. MODIS files are hdf4 not hdf5, so need to use hdfinfo, etc rather
%than hdf5info
if ~exist('comp')
    comp='uni';
    comp='laptop';
    comp='UWchallenger';
end

if exist('imodis_file_override') & imodis_file_override==1
    clear imodis_file_override  %don't set file_name_h5 and reset for next run
else
    file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000056.005.2006254074330.hdf';
    %file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000112.005.2006259161535.hdf';
    %file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000168.005.2006254095414.hdf';
    file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000280.005.2006270142446.hdf'; %213
    %file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000336.005.2006277024458.hdf';  
    
%    file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000055.005.2006254074717.hdf';  %1
%    file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000056.005.2006254074330.hdf';  %3
     file_name_h5='MODIS/D3/y2005/MOD08_D3.A2005001.005.2006203153909.hdf';
     file_name_h5='MODIS/D3/y2004/MOD08_D3.A2004208.005.2007025114245.hdf';
     file_name_h5='MODIS/D3/y2008/MOD08_D3.A2008040.005.2008042055628.hdf';
     
    

    filedir='/home/disk/eos10/robwood/';  
    
%    filedir='/home/disk/eos1/d.grosvenor/';

end

filename_h5=[filedir file_name_h5];


iday=findstr(file_name_h5,'.A');
modis_year_str=file_name_h5(iday+2:iday+5);
modis_day_str=file_name_h5(iday+6:iday+8);





switch comp
case {'uni','UWchallenger'}
    
    INFO = hdfinfo(filename_h5);
    
%    INFO.GroupHierarchy.Datasets.Name %gives the info of all the dataset names
    
    %INFO.GroupHierarchy.Datasets(1)  %this gives all the info about dataset 1
    
%    NVARS = length(INFO.GroupHierarchy.Datasets); %number of variables

%    test = hdf5read(INFO.GroupHierarchy.Datasets(NVARS-1)); %this retrieves dataset 10

    
    
end

SD_id = hdfsd('start',filename_h5,'read'); %open the file

nvar = hdfsd('nametoindex',SD_id,'XDim')+1; %XDim (LON)
MLON = double(hdfread(INFO.Vgroup(1).Vgroup(1).SDS(nvar))); 
nvar = hdfsd('nametoindex',SD_id,'YDim')+1; %YDim (LAT)
MLAT = double(hdfread(INFO.Vgroup(1).Vgroup(1).SDS(nvar))); 

disp('Done read MODIS');

