isave=1;

thresh_CTH=[0 3.2];
thresh_SZA=[0 65];

modis_day_of_year = day_of_year_from_date_func('02-Aug-2016');

modis_day_of_year_str = num2str(modis_day_of_year);
%include the / at the end of filedir
%filedir = '/home/disk/eos15/d.grosvenor/eos8/MOD_L2/NAtlantic/1deg/'; %Aqua

isave_dir = 1; %which of the filedir_mutli dirs to use for saving
%Combine Aqua and Terra for all swaths
filedir_multi{1} = '/home/disk/eos10/d.grosvenor/MOD_L2/NAtlantic/01Aug2016/1deg/'; %Terra
filedir_multi{2} = '/home/disk/eos15/d.grosvenor/eos8/MOD_L2/NAtlantic/1deg/'; %Aqua

fnames=[];
ifile=0;
for idir=1:length(filedir_multi)
    filedir = filedir_multi{idir};
    fnames_new = dir([filedir '*' modis_day_of_year_str '*.mat']);
    
    for ifiles=1:length(fnames_new)
        ifile=ifile+1;
        filedir_multi_ALL{ifile} = filedir;
    end
    
    fnames = cat(1,fnames, fnames_new); %e.g. to pick out just one day
    %fnames = dir([filedir '*.mat']); %all days

end

start_file = 1; %in case want to skip some of the files, etc.
end_file = length(fnames);
files = [start_file:end_file];
%files=[1];

nfiles=length(files);
Nd37_all_t = NaN*ones([180 360 nfiles]);

imat=0;
clear Nd37_cell
for ifname=files
    imat = imat + 1;
    
    file_name_h5 = fnames(ifname).name;
    
    %load the .mat data
    load([filedir_multi_ALL{ifname} file_name_h5]);
    
%    Nd37_cell{imat} = Nd37_1deg;    
    Nd37_all_t(:,:,imat) = Nd37_1deg;
    
%    savename = [filedir '1deg/' file_name_h5 '_1deg.mat'];
%    save(savename,'Nd37_1deg','MLAT','MLON','-V7.3');
        
end


%	combine_overlapping_data_with_Nans.m - can use this to combine the data
%	from several days into one map (averaging when have more than one
%	data point).
%icombine = [1:length(Nd37_cell)];
%Nd37_combined = combine_overlapping_data_with_Nans(Nd37_cell,icombine);
[Nd37_combined, Nvals] = meanNoNan(Nd37_all_t,3);

figure
pcolor(MLON,MLAT,Nd37_combined); shading flat; colorbar;
caxis([0 500]);
set(gca,'xlim',[-80 0]);
set(gca,'ylim',[22 80]);

figure
pcolor(MLON,MLAT,Nvals); shading flat; colorbar;
%caxis([0 500]);
set(gca,'xlim',[-80 0]);
set(gca,'ylim',[22 80]);

if isave==1
    savename = [filedir_multi{isave_dir} modis_day_of_year_str remove_character(filedir_multi{isave_dir},'/','_') '.mat'];
    save(savename, 'Nd37_combined','Nvals','filedir_multi_ALL','MLAT','MLON','fnames','modis_day_of_year_str','-V7.3');    
end
    


