%match the  POLDER retrievals to MODIS ones
%size(Par2_CDR_ALL)   [366   180   360     4]

%days_resized_MODIS;

savefile_MODIS = '/home/disk/eos5/d.grosvenor/PARASOL/modis_DateTime_for_co_location.mat';
%list of vars saved:-
%save(savefile_MODIS,'daynum_resized_MODIS','Date_Time_Swath_MODIS','modisyear_timeseries3_MODIS','daynum_timeseries3_MODIS','ilat_loaded_MODIS','ilon_loaded_MODIS','-V7.3');
load(savefile_MODIS);

Date_Time_Swath.timeseries3 = Date_Time_Swath_MODIS.timeseries3;

siz = size(Date_Time_Swath.timeseries3);

sMod = size(Date_Time_Swath_MODIS.timeseries3);

%declare the arrays to start with
for ivar=1:length(vars_PAR)
    eval_str = [vars_PAR{ivar} '_MODIS = NaN*ones([sMod(1) sMod(2) 4 sMod(4)]);']; eval(eval_str);
end

                                
%loop through all of the days
for it=1:sMod(4)
    it2 = (it-1)*20+1; %assuming here that we have blocks of 20 in daynum
    
    i = find(modisyear_timeseries3_POLDER==modisyear_timeseries3_MODIS(it2) & daynum_timeseries3_POLDER==daynum_timeseries3_MODIS(it2));
    if length(i)>0
        if length(i)>1
            error('Is more than one match per day - only include one satellite for MODIS data');
        end
        
       for ivar=1:length(vars_PAR)
          eval_str=[vars_PAR{ivar} '_MODIS(:,:,:,it) = ' vars_PAR{ivar} '_ALL(i,ilat_loaded_MODIS,ilon_loaded_MODIS,:);']; eval(eval_str);
       end
    end
    
end


%loop through all the time indices and make a POLDER equivalent array




%reshape 
Date_resize = NaN*ones([180 360 20 366]);

if siz(3)==7260
    Date_Time_Swath.timeseries3 = reshape(Date_Time_Swath.timeseries3,[siz(1) siz(2) 20 siz(3)/20]);
end

days = unique(daynum_timeseries3);    
Date_resize(:,:,:,days) = Date_Time_Swath.timeseries3(:,:,:,:);
Date_resize = permute(Date_resize,[4 1 2 3]);

siz_new = size(Date_resize);

max_dt = 0.5/24; %max allowed difference in time (days)



siz2 = size(Par2_MatlabTime_ALL);
index_match = NaN*ones(siz2);

%    imin = NaN*ones(siz_new(1) siz_new(2) size_new(3));
%    imin = NaN*ones(siz_new);
for ioverpass=1:1 %siz(4)
    minval = NaN*ones(siz_new);
    for j=1:20
        minval(:,:,:,j) = abs(Par2_MatlabTime_ALL(:,:,:,ioverpass) - Date_resize(:,:,:,j));
    end
    
    [minval2 imin2] = min(minval,[],4);
    imin2(minval2>max_dt)=NaN;
    index_match(:,:,:,ioverpass) = imin2;
end

