%grid the PARASOL data (after being read in from text)
%into the standard 1x1 degree lat lon grid that MODIS uses

%now the data is in a big array of [nyears 366 nMAX]
%For each day there are a number of datapoints for different locations,
%etc.

iload=1;
if iload==1
%    file_load = '/home/disk/eos5/d.grosvenor/PARASOL/CDR_L2.v01.01/saved_CDR_after_textread_20121214T224544.mat';
    file_load = '/home/disk/eos5/d.grosvenor/PARASOL/CDR_L2.v01.01/saved_CDR_after_textread_20121217T172434.mat';
    load(file_load,'Par_OrbitNum','Par_MatlabTime','Par_Year','Par_Month','Par_Day','Par_Hours','Par_Mins',...
        'Par_Secs','Par_Lat','Par_Lon','Par_CDR','Par_CDRstd','Par_QI');
end



%get the indicies within the 180*360 grid for the nearest lat, lon
%MODIS grid - cell centres
MLAT = [-89.5:1:89.5];
MLON = [-179.5:1:179.5];



sdatALL = size(Par_Lat);

for iyear=1:sdatALL(1)
    
    file_save = ['/home/disk/eos5/d.grosvenor/PARASOL/CDR_L2.v01.01/saved_CDR_processed_' num2str(minALL(Par_Year(iyear,:,:))) '_' datestr(now,30) '.mat'];    
    
    %find the closest latitude for all points
    ilat = ceil(Par_Lat(iyear,:,:))+90;
    ilat(ilat==0)=1; %ones for which lat=-90
    ilon = ceil(Par_Lon(iyear,:,:))+180;
    ilon(ilon==0)=1; %ones for which lat=-90



%make an index for each entry of the Par_XXX arrays that says where it
%should end up in the [nyear nday lat lon] array
sdat = size(Par_Lat(iyear,:,:));
%first the location of the data from the old array
iyearOUT = repmat([1:sdat(1)]',[1 sdat(2) sdat(3)]);
idayOUT = repmat([1:sdat(2)]',[1 sdat(1) sdat(3)]); idayOUT = permute(idayOUT,[2 1 3]);
%for ilatOUT ilonOUT equivalents will just use ilat and ilon

%roll out the indices into a linear array ready for ind2sub
ilat2 = ilat(:);
ilon2 = ilon(:);
iyearOUT2 = iyearOUT(:);
idayOUT2 = idayOUT(:);
%remove the NaNs
% inan=find(isnan(ilat2)==1);
% ilat2(inan)=[]; ilon2(inan)=[]; iyearOUT2(inan)=[]; idayOUT2(inan)=[];
% inan=find(isnan(ilon2)==1);
% ilat2(inan)=[]; ilon2(inan)=[]; iyearOUT2(inan)=[]; idayOUT2(inan)=[];

% %now the ones that will determine the location in the new array
% iyear = repmat([1:sdat(1)]',[1 sdat(2) 180 360]);
% iday = repmat([1:sdat(2)]',[1 sdat(1) 180 360]); iday = permute(iday,[2 1 3 4]);
% ilatALL = repmat([1:180]',[1 sdat(1) sdat(2) 360]); ilatALL = permute(ilatALL,[2 3 1 4]);
% ilonALL = repmat([1:360]',[1 sdat(1) sdat(2) 180]); ilonALL = permute(ilonALL,[2 3 4 1]);

iALL = [1:prod(size(Par_OrbitNum(iyear,:,:)))]'; %all of the indicies in the [iyear 366 2500] array
ilin = sub2ind([sdat(2) 180 360],idayOUT2,ilat2,ilon2); %the indices based on the computed lat and lon indices
%plus the day index (just 1:366)

%remove the NaNs - will be NaNs where Par_Lat is NaN
inan = find(isnan(ilin)==1);
ilin(inan)=[];
iALL(inan)=[]; %do the same for iALL
%iALL contains the indices in the orginal arrays that correspond to ilin
%indices

%sort ilin and reorder iALL to correspond to the same order. This is done
%because unique (done later) sorts the array - sorting now means that we
%can keep track of the order and apply it to iALL
[ilin,isort] = sort(ilin);
iALL = iALL(isort);
%ilin will contain NaNs because ilat and ilon do (due to Par_Lat and
%Par_Lon missing data)

[M,noverpass]=mode(ilin); %noverpass is now the max number of repitions for a year,day,lat,lon location
%I.e. the max number of datapoints per day for one location (globally).
%So, this is the max no. overpasses that we need to store
%Will therefore make the array [nyear nday nlat nlon noverpass]

Par2_OrbitNum = NaN*ones([366 180 360 noverpass]);
Par2_MatlabTime = NaN*ones([366 180 360 noverpass]);
%Par2_Year = NaN*ones([366 180 360 noverpass]);
%Par2_Month = NaN*ones([366 180 360 noverpass]);
%Par2_Day = NaN*ones([366 180 360 noverpass]);
%Par2_Hours = NaN*ones([366 180 360 noverpass]);
%Par2_Mins = NaN*ones([366 180 360 noverpass]);
%Par2_Secs = NaN*ones([366 180 360 noverpass]);
Par2_Lat= NaN*ones([366 180 360 noverpass]);
Par2_Lon = NaN*ones([366 180 360 noverpass]);
Par2_CDR = NaN*ones([366 180 360 noverpass]);
Par2_CDRstd = NaN*ones([366 180 360 noverpass]);
Par2_QI = NaN*ones([366 180 360 noverpass]);

'doing unique bit...'
%NOW done more quickly using ndHistc
%now go through all of the unique values of the indices and make new
%indices for a bigger array with an noverpass dimension if necessary (if
%there is more than one datapoint for a location per day)
[B,I,J] = unique(ilin); %B=ilin(I) and ilin=B(J)
%WARNING - unique sorts the values of ilin! So need to put i1,i2,i3 below
%in the same order as for the B values. So need to use ilin_ex (below) for
%the i1,i2,i3 values NOT this :- [i1,i2,i3] = ind2sub([sdat(2) 180
%360],ilin);

%Will make use of the B=ilin(I) feature later for iALL. Note that the I
%values are the indices for the LAST index in ilin (in the sorted order)
%that gives each B value. E.g. for [B,I,J]=unique([10 10 20 20 20]) will
%give I = [2 5].


%This next part was included since it is much faster than the for loop
%commented out. The key part is the histogram, which uses C and so is much
%faster.
%make bins that enclose the whole numbers
bins=[B-0.5; B(end)+0.5];
%get histrogram of occurrences of all the unique ilin values
qh = ndHistc_run([ilin], bins);

%this expands a 1D PDF so that all of the entries are contained in the new
%array in the correct frequencies
ilin_ex = expand_PDF(B,qh);
[i1,i2,i3] = ind2sub([sdat(2) 180 360],ilin_ex);
%This does the same as above except that it contains an array 1:N where N
%is the no. of occurences. These will be the overpass numbers for each
%time,lat,lon
ioverpass = expand_PDF_special(B,qh);

% iALL(I) correspond to locations in iALL for each unique ilin value in B
%Becuase we sorted ilin, and because the I from unique gives the LAST index for a given B value
%iALL(I-n+1:I) will correspond to the same ilin value
%for ilin values where there are n occurences for a given B value, i.e.
%ilin(I-n+1:I)=ilin(I)
%So we expand the I values to the same size as ilin_ex
I_ex = expand_PDF(I,qh);
%We want to reference iALL(I-n+1:I). ioverpass contains e.g. [1 1 1 2 1 2 3]
%If we take one off this then I_ex-(ioverpass-1) will give us this.
%since I_ex = e.g. [1092 294 565 565 1001 1001 1001], i.e. repitions of the
%I vector I_ex-(ioverpass-1) = [1092 294 565 564 1001 1000 999]
%E.g. 1001, 1000 and 999 correspond to the same values in ilin because of
%the way that unique sets I. So we reference these values from e.g.
%Par_Lat(iyear,:)
iALL_ex2 = iALL(I_ex - (ioverpass-1) );


% ilin2=NaN*ones(size(ilin));
% ind=1;
% for iuni=1:length(B)
%     ival = find(ilin==B(iuni));
%     N = length(ival);
%     ilin_tmp = sub2ind([sdat(1) sdat(2) 180 360,noverpass],i1(ival),i2(ival),i3(ival),i4(ival),[1:N]');    
%     ilin2(ind:ind+N-1) = ilin_tmp;
%     ind=ind+N;
% end

ilin2 = sub2ind([sdat(2) 180 360 noverpass],i1,i2,i3,ioverpass);    

%now put in the new array
Par2_OrbitNum(ilin2) = Par_OrbitNum(iyear,iALL_ex2);
Par2_MatlabTime(ilin2) = Par_MatlabTime(iyear,iALL_ex2);
%Par2_Year(ilin2) = Par_Year(iyear,iALL_ex2);
%Par2_Month(ilin2) = Par_Month(iyear,iALL_ex2);
%Par2_Day(ilin2) = Par_Day(iyear,iALL_ex2);
%Par2_Hours(ilin2) = Par_Hours(iyear,iALL_ex2);
%Par2_Mins(ilin2) = Par_Mins(iyear,iALL_ex2);
%Par2_Secs(ilin2) = Par_Secs(iyear,iALL_ex2);
Par2_Lat(ilin2) = Par_Lat(iyear,iALL_ex2);
Par2_Lon(ilin2) = Par_Lon(iyear,iALL_ex2);
Par2_CDR(ilin2) = Par_CDR(iyear,iALL_ex2);
Par2_CDRstd(ilin2) = Par_CDRstd(iyear,iALL_ex2);
Par2_QI(ilin2) = Par_QI(iyear,iALL_ex2);

fprintf(1,'\nSaving\n');

%save(file_save,'Par2_OrbitNum','Par2_Year','Par2_Month','Par2_Day','Par2_Hours','Par2_Mins',...
%        'Par2_Secs','Par2_Lat','Par2_Lon','Par2_CDR','Par2_CDRstd','Par2_QI');

save(file_save,'Par2_OrbitNum','Par2_MatlabTime','Par2_Lat','Par2_Lon','Par2_CDR','Par2_CDRstd','Par2_QI');

end


fprintf(1,'\n Done process PARASOL text data\n');



