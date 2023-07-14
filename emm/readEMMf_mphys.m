function [Temm,Zemm,mp]=readEMMf_mphys(fdir,fname);

max_size=5e5;
%    mp=dlmread([fdir fname],' ');
try
    mp=dlmread([fdir fname],' ',[0 0 max_size 38]);
catch
    mp=dlmread([fdir fname],' '); % in case the file is smaller than  max_size
end
%fid=fopen([fdir fname],'rt');
%mp=fscanf(fid,'%f',[38 Inf]);
    ncols=size(mp,2)-1;
    %fid=fopen([emmdir2 'mphys_diags'],'rt');
    %mp=fscanf(fid,'%e',[58 Inf]);
    
    
%    mp=mp';
    mp(end,:)=[];
    
    Temm=unique(mp(:,1));
    
	a=find(mp(:,1)==Temm(1));
	s1=length(a);
	b=size(mp,1);
	s3=size(mp,2);
	s2=floor(b/s1);
	
%	mp=mp(1:s1*s2,:);
    mp(s1*s2+1:end,:)=[];
	mp=reshape(mp,[s1 s2 s3]); %re-shape into 3-d array
    
    Zemm=squeeze(mp(:,1,2));
    
    a=find(Zemm~=0);
    Zemm(a(end):end)=[];
    mp(a(end):end,:,:)=[];






% iwczt=dlmread([fdir fname],' ');
% dt=iwczt(2:end,3)-iwczt(1:end-1,3);
% a=find(dt>0);
% da=a(2:end)-a(1:end-1);
% nt=size(iwczt,1);
% aim=(length(a)+1)*da(1);
% if (nt<aim)
%     iwczt(nt+1:aim,:)=NaN;
% end
% 
% arr=reshape(iwczt,[da(1) length(a)+1 size(iwczt,2)]);
% Temm=arr(:,1,2);
% temm=arr(1,:,3);