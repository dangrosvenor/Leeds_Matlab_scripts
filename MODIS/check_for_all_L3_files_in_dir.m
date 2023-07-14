%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/Collection_6/aqua/2008/du.out';
%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/aqua/2008/du.out';
%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/Collection_6/terra/2008/du.out';
%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/terra/2006/du.out';

tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2011/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2011_a/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2012/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2012_a/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2013/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2013_a/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2014/';
%tfile = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/y2014_a/';

%fid=fopen(tfile,'rt');
%dat=textscan(fid,'%s%s');
%fclose(fid);

dat=dir(tfile);

full=1:366;

clear day
for i=3:length(dat)
    if strcmp(dat(i).name(end-3:end),'.hdf')==1
        str01 = str2num( dat(i).name(15:17) );

        if length(str01)>0   %
            day(i-2) = str01;
        end

    end
end




missing=[];
ii=0;
for i=1:length(full)
    ifind = find(day==full(i));
    if length(ifind)==0
        ii=ii+1;
        missing(ii)=full(i);
    end
end

missing

fprintf(1,'\n Done file check.\n');

