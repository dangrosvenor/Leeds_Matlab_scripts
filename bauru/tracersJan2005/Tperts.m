run='3d';
run='2d';
run='2d2';

switch run
case '2d2'        

idirs=1:3;
for idat=1:length(idirs)
    idir=idirs(idat);
	backgroundTH=median(TwoDDan(idir).TH2,2);
    backgroundT=backgroundTH ./ (1000e2./GridDan(idir).PREFN).^0.286;
    T=TwoDDan(idir).TH2 ./ repmat( (1000e2./GridDan(idir).PREFN).^0.286 , [1 size(TwoDDan(idir).TH2,2)]);
	max_pert(idat).m=max( T-repmat(backgroundT,[1 size(TwoDDan(idir).TH2,2)]) ,[],2 );
    
end


% mean7km=mean(tpertTimH_full(1).t(:,[end-6:end 1:7],:),2) - backgroundT;
% 
% ih=findheight(GridDan(1).Z,2.5e3);
% 
% hmean=mean(mean7km(2:ih,:,:),1);
% 
% backgroundQ=median(vapfull(1).v(:,:,:),2);
% vmean7km=mean(vapfull(1).v(:,[end-6:end 1:7],:),2) - backgroundQ;
% vhmean=mean(vmean7km(2:ih,:,:),1);


case '2d'        
backgroundT=median(tpertTimH_full(1).t(:,:,:),2);

mean7km=mean(tpertTimH_full(1).t(:,[end-6:end 1:7],:),2) - backgroundT;

ih=findheight(GridDan(1).Z,2.5e3);

hmean=mean(mean7km(2:ih,:,:),1);

backgroundQ=median(vapfull(1).v(:,:,:),2);
vmean7km=mean(vapfull(1).v(:,[end-6:end 1:7],:),2) - backgroundQ;
vhmean=mean(vmean7km(2:ih,:,:),1);

case '3d'
%load in theta.mat in 3d2km_30_20mins

%[ia ib ic id]=size(theta);    
%th_mean=repmat(meTH(1:ic),[ia ib 1 id]); %ic is the height index of the max height in theta (2.5 km)
%th_mean=permute(th_mean,[2 3 1 4]);   %to put in same order as theta


%th=theta-th_mean;

%mean_pert=mean(mean(mean(theta(3:9,4:10,:,:))));


clear diff
dy=diff(GridDan(1).Y1(1:2));
D=14e3;
%D=5e3;
ny=D/dy;

[iy iy2]=findheight( Grid.Y1,Grid.Y1(1)+D/2,Grid.Y1(end)-D/2 );
[ix ix2]=findheight( Grid.X1,Grid.X1(1)+D/2,Grid.X1(end)-D/2 );
[iz2]=findheight( Grid.Z , 2.52e3 ); %should be index 29

mean2 = mean(mean(   TH1( [end-8:end-1 2:9] , [end-7:end 1:8] , 1:iz2 )   ));

mean2 = mean(mean(   TH1( [end-ny-1:end-1 2:2+ny] , [end-ny:end 1:ny+1] , 1:iz2 )   ));

mean2 = mean(mean(   TH1( [2] , [end-ny:end 1:ny+1] , 1:iz2 )   ))


end


disp('done');