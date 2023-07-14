%save Par2_CDR_MODIS, etc. from parasol_match_modis_datapoints.m


%Par2_CDR_MODIS is of size e.g. [3 51 20 72] = [lat lon orbit day*year]. Can reshape 
%the last 3 indices into MODIS sized array (E.g. 1440)

savefile_POLDER = '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_CDR_colocated_Arctic.mat';
savefile_POLDER = '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_CDR_colocated_VOCALS_2006-2007.mat'; %Actually only contains 2007.
savefile_POLDER = '/home/disk/eos5/d.grosvenor/PARASOL/POLDER_CDR_colocated_VOCALS_2008.mat'; %

 for ivar=1:length(vars_PAR)
          eval_str=[vars_PAR{ivar} '_coloc = reshape(' vars_PAR{ivar} '_MODIS,[sMod(1) sMod(2) sMod(3)*sMod(4)]);']; eval(eval_str);
          if ivar>1
              app_str = ',''-APPEND'');';
          else
               app_str = ');';
          end
          eval_str=['save(savefile_POLDER,''' vars_PAR{ivar} '_coloc'',''-V7.3''' app_str]; eval(eval_str);
 end
       
disp('Done save');

%load(savefile_POLDER)