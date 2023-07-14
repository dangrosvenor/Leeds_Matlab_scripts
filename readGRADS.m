skip=0;
rev=0;
sound=0; %flag to make sounding
oned=0;
n1=78;
n2=85;


clear x psfc topo ii direc;






if rev==1
    n=n1;
    n1=n2;
    n2=n;
end

%sens07.02-15vslice contains a slice at lon=-49.03 lat=-23.86 to -20.86 (~334km) at
%time of 15:00LT (38 elements) + lat07.02... too set t 43 in GRADS
if ~exist('oridefile');oridefile=0;end
if oridefile==0
    direc(1).d='c:\program files\pcgrads\win32e\matlab\senGRID';
    direc(2).d='c:\program files\pcgrads\win32e\matlab\latGRID';
    direc(3).d='c:\program files\pcgrads\win32e\matlab\Z0GRID';
    direc(4).d='c:\program files\pcgrads\win32e\matlab\cbntGRID';
    
end
% direc(1).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-TopL';
% direc(2).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-BotL';
% direc(3).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-TopR';
% direc(4).d='c:\program files\pcgrads\win32e\30.01-00\profs\prof-16-2-BotR';


for ii=1:size(direc,2)

pr(ii).p=[];    
fid=fopen(direc(ii).d,'rb');




% fseek(fid,4,'cof');
% for i=1:n
% x(1:39,i)=fread(fid,39,'float=>double');
% fseek(fid,24,'cof');
% end

if skip==1
    fseek(fid,4,'cof');
end

if oned==1
    x=fread(fid,n1,'float=>double');
    pr(ii).p=x;
else
    if sound==1
        psfc(ii)=fread(fid,1,'float=>double')/100; %surface pressure
        topo(ii)=fread(fid,1,'float=>double'); %surface height
    end
    for i=1:n1
        %fseek(fid,4,'cof');
        xt=fread(fid,n2,'float=>double');
        pr(ii).p(1:n2,i)=xt;
%         if rev==1
%             fseek(fid,4,'cof');
%         end
    end
    
   
end

fclose(fid);

pr(ii).p=pr(ii).p';
temp=pr(ii).p(:,2:end);
pr2(ii).p=temp;

end

if sound==1
	pr(ii).p(:,3)=pr(ii).p(:,3)-273.13; %convert to C
	pr(ii).p(:,7)=log10(pr(ii).p(:,3).*(1000./pr(ii).p(:,2)).^0.286);%log theta
	pr(ii).p(:,8)=sqrt ( pr(ii).p(:,5).^2 + pr(ii).p(:,6).^2 ); %wind mag
	pr(ii).p(:,9)=bear( pr(ii).p(:,5),pr(ii).p(:,6) ); %wind bearing

    etasoundGRAD;
end

oridefile=0;


lat=[-25.483:0.084:-19.015];
lon=[-52.582:0.084:-45.526];

lon2=[-52.582+0.084:0.084:-45.526];









