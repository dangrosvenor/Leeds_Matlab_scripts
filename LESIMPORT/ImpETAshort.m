%imports meso-eta data as modified by Jorge - format as in etabru10.2004020912.ctl
grid.x=[-54.262:0.084:-44.938]; %112 points
grid.y=[-25.483:0.084:-19.015]; %78 points

bx=63; 
by=38; %indices for nearest to Bauru co-ords

%xdan(i).x = [78 112] , i.e. [x y]

latx=49.03; 
lony=22.36;  %long and lat for Bauru


clear xdan thd temp;
%pat='c:\documents and settings\user\my documents\HIBISCUS\baurufield\mesoeta\bin\05.02.04\etabru10.2004020512+04.bin';
fid=fopen(pat,'rb');

for j=1:14
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

% 
% for j=1:6
%     for k=1:39
%         fseek(fid,8,'cof');
% 	    for i=1:78
%             temp(:,i)=fread(fid,112,'float=>double');
% 	    end
%         thd(j).x(:,:,k)=temp';
%     end
% end 