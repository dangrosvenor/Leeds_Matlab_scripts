%horizontal sum of depleted vapour over height
clear diff zsum

for idir=1:2;

	vap=TwoDDan(idir).Q(2:end,:,1);
	
	dy=diff( GridDan(idir).Y1(1:2) );
	air=diff( GridDan(idir).Z ) .* GridDan(idir).RHON(2:end) ;
	for i=2:size(vap,2)
		i5=find(vap(:,i)<5/f);   %height indices of air w/ MR < 5 ppmv
        if length(i5)>=1
            zsum(i)=sum( (5/f - vap(i5,i)) .* air(i5) ) * dy;
        else
            zsum(i)=0;
        end
	end 
	
	
	ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+300e3);  %for radar plots
	ix=1;       
	xinds=[length(GridDan(idir).Y1)-ix2:length(GridDan(idir).Y1) ix:ix2 ]; %
	
	zsum(1)=0;
	
	zs(idir).z = zsum(xinds);
	xs(idir).x = GridDan(idir).Y1(1:length(xinds))'/1000;
    
    cum(idir).c=cumsum(zs(idir).z) * dy;

end

c=cum(1).c(1:2:end);
c2=cum(2).c(1:end-1);

figure
plot(xs(2).x(1:end-1),c-c2);

%%%%%%%%%%%%%%%%%%   new bit for working out means across selected domain section %%%%%

for idir=1:2;

	vap=TwoDDan(idir).Q(2:end,:,1);
	tot=sum(TwoDDan(idir).Q(2:end,:,[1:6]),3);
	
	dy=diff( GridDan(idir).Y1(1:2) );
	air=diff( GridDan(idir).Z ) .* GridDan(idir).RHON(2:end) ;

	
	
	ix2=findheight(GridDan(idir).Y1,GridDan(idir).Y1(1)+600e3);  
    % areas plotted in 2-d field that values taken from
    
	ix=1;       
	xinds=[length(GridDan(idir).Y1)-ix2:length(GridDan(idir).Y1) ix:ix2 ]; %
	
	zsum(1)=0;
	
	zs(idir).z = zsum(xinds);
	xs(idir).x = GridDan(idir).Y1(1:length(xinds))'/1000;
    

    x1=-720;
    x2=100;
    
    x1=-950;
    x2=200;
    
    
    i1=findheight(xs(idir).x,x1);
    i2=findheight(xs(idir).x,x2);
    
    xinds2=xinds(i1:i2);
    
    meanvap(idir).m=mean( vap(:,xinds2) ,2);
    meantot(idir).t=mean( tot(:,xinds2) ,2);
    
end

figure
plot(f*meanvap(1).m,GridDan(1).Z(2:end)+620);
hold on;
plot(f*meanvap(2).m,GridDan(1).Z(2:end)+620,'r');
set(gca,'ylim',[14.8 18]*1e3);
grid



    
    
    
    
		'done'


