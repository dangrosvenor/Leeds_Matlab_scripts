

dx=1e3;
dy=1e3;

clear diff

%for itime=1:size(Times,1)
times=218:234;
%times=[15 20];


iread_1=1;

for itime2=1:1   %length(times)   
     itime=times(itime2);
     disp(itime);                        

    
    if iread_1==1
        qice=nc{'QICE'}(itime,:,:,:);
        qsnow=nc{'QSNOW'}(itime,:,:,:);
        qgraup=nc{'QGRAUP'}(itime,:,:,:);
        qvap=nc{'QVAPOR'}(itime,:,:,:);
        qrain=nc{'QRAIN'}(itime,:,:,:);

        potemp=WRFUserARW(nc,'th',itime);
        zwrf=WRFUserARW(nc,'Z',itime); %height profile at position of max total water
        pressure=WRFUserARW(nc,'p',itime); %pressure (mb)
        pressure.var=pressure.var*100; %convert to Pa
        temp=WRFUserARW(nc,'tc',itime); %temperature degC 
    end        
    

fileName=['/home/mbexddg5/work/wrfvars_' num2str(itime) '_' strrep(Times(itime,12:16),':','-') '.mat'];
save(fileName,'-V6','qice','qsnow','qgraup','qvap','qrain','potemp','zwrf','pressure','temp');
   

end

disp('Finished');
    
    
         
