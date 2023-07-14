%Matches the POLDER data to the VOCALS region in the file below and saves
%(doing all the flipping, etc.)

%First load the POLDER data using load_process_PARASOL_gridded.m
%Select ioverall_mean=0;

%VOCALS region for 2007 with the POLDER data co-located onto MODIS data
%(mock L3)
savefile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only.mat';
%loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only_bkup.mat';
loadfile = savefile;
load(loadfile,'LAT_MODIS','LON_MODIS');

clear ilat ilon
for i=1:length(LAT_MODIS)
    ilat(i) = find(MLAT_POLDER==LAT_MODIS(i));
end
for i=1:length(LON_MODIS)
    ilon(i) = find(MLON_POLDER==LON_MODIS(i));
end

if ~exist('daymean_Par2_CDR_orig')
    daymean_Par2_CDR_orig = daymean_Par2_CDR;
end

clear daymean_Par2_CDR

%declare the arrays to start with

vars_PAR2 = {'daymean_Par2_CDR_orig', ...
            };

for ivar=1:length(vars_PAR2)
%    eval_str = [vars_PAR2{ivar} '_MODIS = NaN*ones([length(LAT_MODIS) length(LON_MODIS) size(daynum_timeseries3_POLDER,1)]);']; eval(eval_str);
    eval_str = [vars_PAR2{ivar} '_MODIS = permute(' vars_PAR2{ivar} '(:,ilat,ilon),[2 3 1]);']; eval(eval_str);
% No need to flipdim since have searched for the MODIS values.
end

 
%Now match to the time in the big (5 year) POLDER array in the VOCALS .mat
%file
vocals = load(loadfile,'daymean_Par2_CDR','modisyear_timeseries3_MODIS','daynum_timeseries3_MODIS');
daymean_Par2_CDR = vocals.daymean_Par2_CDR;
   
for i=1:length(daynum_timeseries3_POLDER)
   it = find(vocals.daynum_timeseries3_MODIS==daynum_timeseries3_POLDER(i) & vocals.modisyear_timeseries3_MODIS==modisyear_timeseries3_POLDER(i));
   if length(it)>0
       daymean_Par2_CDR(:,:,it) = daymean_Par2_CDR_orig_MODIS(:,:,i);
   end
end

%Then run

%  save(savefile,'daymean_Par2_CDR','-V7.3','-APPEND');



% %loop through all of the days
% for it=1:length(modisyear_timeseries3_MODIS)
%     
%     i = find(modisyear_timeseries3_POLDER==modisyear_timeseries3_block(it) & daynum_timeseries3_POLDER==daynum_timeseries3_MODIS(it));
%     if length(i)>0
%        sst_amsre_time3(:,:,it) = sst_amsre_smooth(ilat,ilon,i);
%        %lwp_amsre is larger in the 3rd dimension than sst_asmre_smooth, but
%        %the end of it will just be padded with zeros as is created as
%        %nmonths*31, whereas smooth one is only the size of the actual data
%        %since it is based on year_amsre, which grows each loop
%        if exist('lwp_amsre')
%            lwp_amsre_time3(:,:,it,:) = lwp_amsre(ilat,ilon,i,:);
%        end
%     end
%     
% end
