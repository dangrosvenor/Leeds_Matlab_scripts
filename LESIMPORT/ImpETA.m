%imports meso-eta data as modified by Jorge - format as in etabru10.2004020912.ctl

clear xdan thd temp;
pat='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\bin\05.02.04\etabru10.2004020512+00.bin';
fid=fopen(pat,'rb');

for j=1:34
    fseek(fid,4,'cof');
    if j>1
        fseek(fid,4,'cof');
    end
	for i=1:78
        xdan(j).x(:,i)=fread(fid,112,'float=>double');
        %fseek(fid,4,'cof');
	end
    xdan(j).x=xdan(j).x';
end 

for j=1:6
    for k=1:39
        fseek(fid,8,'cof');
	    for i=1:78
            temp(:,i)=fread(fid,112,'float=>double');
	    end
        thd(j).x(:,:,k)=temp';
    end
end 