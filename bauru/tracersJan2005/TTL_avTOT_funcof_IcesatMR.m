t1=23.46;
it1=findheight(GridDan(idir).t+3,t1);
   
ih2=findheight(GridDan(1).Z+620,17e3);
ih1=findheight(GridDan(1).Z+620,15e3);

ppmvs=[0:0.02:20]; %find likely global variation in TTL ice sat MR (radiosonde papers for e.g.)

for idir=1:2
    
    for ippmv=1:length(ppmvs)
        
        tot_init = f*sum(icediagsALL(idir).i(:,1,[37:42]),3)/npess2(idir);   %total water for 1st casei0=find(xdat(idat).x>6.5);            
		tot = f*sum(icediagsALL(idir).i(:,it1,[37:42]),3)/npess2(idir);   %total water for 1st casei0=find(xdat(idat).x>6.5);            
	%            i0=find(xdat(idat).x>8.5);            
        
        i0=find( tot > ppmvs(ippmv) );
        tot(i0)=tot_init(i0); %make them the same as inital so diff=0 for tot > X ppmv        
        
        rho=GridDan(idir).RHON(:); %convert to kg/km3 as xdat in g/kg km 
        
        dz=diff(GridDan(idir).Z(ih1-1:ih2))/1000;        
        air=sum( diff( GridDan(idir).Y1([1 end]) ).*dz.*rho(ih1:ih2) );
                
        total(idir).t(ippmv)=sum( ( tot(ih1:ih2) - tot_init(ih1:ih2) ).*dz.*rho(ih1:ih2) ) / air * diff( GridDan(1).Y1([1 end]) ) ;
        total(idir).t = total(idir).t 
        
	end
    
    
end

figure
plot(ppmvs,total(1).t);
hold on
plot(ppmvs,total(2).t,'r');