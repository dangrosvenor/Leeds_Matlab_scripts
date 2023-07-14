%Changed from original file on 28th Sep, 2016 (same name, but without date)
%Matches the POLDER data to the VOCALS region in the file below and saves
%(doing all the flipping, etc.)

%First load the POLDER data using load_process_PARASOL_gridded.m
%Select ioverall_mean=0;

%data_str = 'MODIS';
data_str = 'POLDER';

%Choose which vars to save
switch data_str
    case 'MODIS'
        modis_var_mockL3 = modis_var;  %e.g. if have just run MODIS... to load and
        %concatenate the mockL3 files
        
        modis_var_mockL3{end+1}='CTH';
        modis_var_mockL3{end+1}='CTH_all';
        modis_var_mockL3{end+1}='CTH_max';        
        modis_var_mockL3{end+1}='sst_amsre_time3';
        modis_var_mockL3{end+1}='lwp_amsre_time3';
        modis_var_mockL3{end+1}='N_time3';        
        modis_var_mockL3{end+1}='N_time3_16';        
        modis_var_mockL3{end+1}='N_time3_37';        
        modis_var_mockL3{end+1}='W_time3';   
        modis_var_mockL3{end+1}='N_un_time3';        
        modis_var_mockL3{end+1}='homog_time3';        
        modis_var_mockL3{end+1}='homog_time3_W';        
        modis_var_mockL3{end+1}='homog_time3_meanW';                
        modis_var_mockL3{end+1}='stdW_time3';   
        
    case 'POLDER'

        %Or select some other variables (e.g. POLDER)
        modis_var_mockL3 = {'daymean_Par2_CDR', ...
            };

end

%VOCALS region for 2007 with the POLDER data co-located onto MODIS data
%(mock L3)
%savefile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only.mat';
%savefile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/VOCALS_MODIS_Aqua_2007_no_confidence_screening.mat';
savefile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2005-2012_Aqua_mockL3_no_confidence_screening_with_POLDER.mat';
%The following was set as of 28th Sep 2016 - changed to the one above - is
%there an issue with no. pixels between 1.6, 2.1 and 3.7 for no confidence
%screening?
%savefile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2005-2012_Aqua_mockL3_confidence_screen_20um_limit_re16_37_consistent_pixels.mat';

loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2005-2012_Aqua_mockL3_no_confidence_screening_with_POLDER.mat';
%loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only_bkup.mat';
%loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only.mat';
%loadfile = savefile;

%loadfile is where the big array to be matched to is stored - the new data
%should be in the main memory
Lfile = load(loadfile,'MLAT','MLON','modisyear_timeseries3_MODIS','daynum_timeseries3_MODIS');


LAT_dat = MLAT;
LON_dat = MLON;

clear ilat ilon
for i=1:length(Lfile.MLAT)
    ilat(i) = find(LAT_dat==Lfile.MLAT(i));
end
for i=1:length(Lfile.MLON)
    ilon(i) = find(LON_dat==Lfile.MLON(i));
end

 
%Now match to the time in the big (5 year) POLDER array in the VOCALS .mat
%file

clear it_mat it_pol   
i2=0;
for i=1:length(daynum_timeseries3)
   it = find(Lfile.daynum_timeseries3_MODIS==daynum_timeseries3(i) & Lfile.modisyear_timeseries3_MODIS==modisyear_timeseries3(i));
   if length(it)>0
       %daymean_Par2_CDR(:,:,it) = daymean_Par2_CDR_orig_MODIS(:,:,i);
       i2=i2+1;
       it_mat(i2) = it; %reference in the small (from .mat file) array
       it_pol(i2) = i; %reference in the POLDER (potentially all lats and lons) array (which becomes dat in the function)
   end
end

for ivar=1:length(modis_var_mockL3)
%    eval_str = [vars_PAR2{ivar} '_MODIS = NaN*ones([length(LAT_MODIS) length(LON_MODIS) size(daynum_timeseries3_POLDER,1)]);']; eval(eval_str);
%    eval_str = [vars_PAR2{ivar} '_MODIS = permute(' vars_PAR2{ivar} '(:,ilat,ilon),[2 3 1]);']; eval(eval_str);
% No need to flipdim since have searched for the MODIS values.
    
    if length(eval(modis_var_mockL3{ivar}))>1
        varname_ext = '';
    else
        varname_ext = '.timeseries3';
    end
    eval_str = ['add_extra_mockL3_to_VOCALS_mat_file_FUNC_28Sep2016(' modis_var_mockL3{ivar} ',''' savefile ''',''' modis_var_mockL3{ivar} ''',varname_ext,ilat,ilon,it_mat,it_pol);'];
    eval(eval_str);
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
