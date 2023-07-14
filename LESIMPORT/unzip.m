function unzip(dire,fname)

zip=0;
lis=dir(dire);
for izip=1:length(lis)
    if strcmp(lis(izip).name,fname)==1 %if unzipped file is present then use that instead
        zip=0;
        break
    end
    if strcmp(lis(izip).name,[fname '.gz'])==1 %is zipped file present prepare to use that
        zip=1;
    end
end

if zip==1
    'unzipping file ...'
    eval(['!c:/cygwin/bin/gzip -d ' [dire fname] '.gz']);
    'finished unzipping'
end