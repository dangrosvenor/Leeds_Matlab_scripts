tfile = '/home/disk/eos5/d.grosvenor/joint_L2/Collection_6/aqua/2008/du.out';
tfile = '/home/disk/eos5/d.grosvenor/joint_L2/aqua/2011/du.out';
%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/Collection_6/terra/2008/du.out';
%tfile = '/home/disk/eos5/d.grosvenor/joint_L2/terra/2006/du.out';
tfile = '/home/disk/eos5/d.grosvenor/joint_L2/aqua/2015/du.out';

fid=fopen(tfile,'rt');
dat=textscan(fid,'%s%s');
fclose(fid);

full=1:366;

clear day
for i=1:length(dat{1})-1  %the number of folders (i.e. days) in the list
    str01 = str2num( dat{2}{i}(3:end) ); %the day of the year
    
    if length(str01)>0   %if e.g. has 032_old then its length will be zero after str2num
        day(i) = str01;
        siz(i) = str2num(dat{1}{i});
    end
end

isiz = find(siz<100);  %find the empty folders with no swaths in at all
missing2 = day(isiz);  %Days with no data in at all
%day = sort(day);

%This checks which folders are missing
missing=[];
ii=0;
for i=1:length(full)
    ifind = find(day==full(i));
    if length(ifind)==0
        ii=ii+1;
        missing(ii)=full(i);
    end
end



fprintf(1,'\n Done file check.\n');

