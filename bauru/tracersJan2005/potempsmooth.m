function [xx,zz,ilin]=potempsmooth(TwoD,Grid,pot,fac,acc);

tref=repmat(Grid(1).THREF,[1 size(TwoD(1).TH1,2)]);
dat=squeeze(sum(TwoD(1).TH1,3))+tref; %potemp

[a b]=find(dat>(pot-acc) & dat<(pot+acc));

z=Grid(1).Z(a)/1000+0.62;
x=Grid(1).Y1(b)/1000;

% dzz=diff(z,2);
% dxx=diff(x(2:end));
% 
% g=dzz./dxx;
% inan=find(isnan(g));
% g(inan)=[];
% x(inan)=[];
% z(inan)=[];
% a(inan)=[];
% b(inan)=[];
% 
% i=find(g~=Inf);
% 
% me=mean(g(i));
% 
% ii=find(abs(g)<fac*abs(me)); %keep points with 2nd derivative less than fac*mean
% 
% xx=x(ii);
% zz=z(ii);
% 
% aa=a(ii);
% bb=b(ii);

xx=x;
zz=z;
aa=a;
bb=b;

ilin=sub2ind(size(TwoD.TH1),aa,bb); %linear indices of required points for 2-d grid
