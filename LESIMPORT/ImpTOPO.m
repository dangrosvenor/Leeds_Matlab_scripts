%imports meso-eta data as modified by Jorge - format as in etabru10.2004020912.ctl

clear xdan thd temp x;
pat='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\topo\30S50W.bin';
fid=fopen(pat,'rb');

for j=1:1200
    %fseek(fid,4,'cof');
%     if j>1
%         fseek(fid,4,'cof');
%     end
	%for i=1:78
        x(:,j)=fread(fid,1200,'float=>double');
        %x(:,j)=fread(fid,1200);
        %fseek(fid,4,'cof');
        %end
    %xdan(j).x=xdan(j).x';
end 

x=flipdim(x',1);


break;
for j=1:6
    for k=1:39
        fseek(fid,8,'cof');
	    for i=1:78
            temp(:,i)=fread(fid,112,'float=>double');
	    end
        thd(j).x(:,:,k)=temp';
    end
end 