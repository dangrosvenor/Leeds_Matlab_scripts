% to import TRMM precipitaiion data
clear x y p

if ~exist('oridefile')
    oridefile=0
end
if oridefile==0
    direc(1).d='c:/documents and settings/login/my documents/hibiscus/troccinoxdb/trmmrainfall/3B42RT_2004022500.dat';
end
% direc(1).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-TopL';
% direc(2).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-BotL';
% direc(3).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-TopR';
% direc(4).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-BotR';


for ii=1:size(direc,2)

pr(ii).p=[];    
fid=fopen(direc(ii).d);




% fseek(fid,4,'cof');
% for i=1:n
% x(1:39,i)=fread(fid,39,'float=>double');
% fseek(fid,24,'cof');
% end

% a=1;
% while a>
% 	c=fscanf(fid,'%s',1)
% 	[a b]=size(c);
%     %c=fscanf(fid,'%c',[1 1]);
% end

b=0;
while (b<10)
nowt=fscanf(fid,'%s',1);
pr(ii).p=fscanf(fid,'%g',[3 inf]);
[a b]=size(pr(ii).p);
end

p(ii).p=reshape(pr(ii).p(3,:),181,161);
x(ii).x=pr(ii).p(2,1:181);
y(ii).y=pr(ii).p(1,1:181:end);
        
    

end

oridefile=0;

