
times=15:25;
%times=[15 20];

clear max_tot max_ind max_potemp

 for itime2=1:length(times)   
     itime=times(itime2);
     disp(itime);
     
     potemp=WRFUserARW(nc,'th',itime);
	 i380=find(potemp.var>380);
    
     qtotal=nc{'QICE'}(itime,:,:,:)+nc{'QSNOW'}(itime,:,:,:)+nc{'QGRAUP'}(itime,:,:,:)+nc{'QVAPOR'}(itime,:,:,:); %total water mixing ratio
     
     [a,b]=maxALL(qtotal(i380));
     
     max_tot(itime2)=a;
     
     [I,J,Z]=ind2sub(size(potemp.var),i380(b(1)));
     
     max_ind(itime2).inds=[I,J,Z];
     max_potemp(itime2)=potemp.var(I,J,Z);
     
 end