%constant rates of change for U and V and evolve towards sounding for wind.

clear force

time(1)=(20.25-12)/24; %1st sounding at 12UTC, 2nd at 20:15 (units = day)
time(2)=(24-20.25)/24; 
stats=[2 3]; %stations to be used for 1st and 3rd soundings

for i=1:2
    diff=pro(2).p(2:end,1,i)-pro(2).p(1:end-1,1,i);
    izero=find(abs(diff-0)<0.01); %point where data ends in sounding
    if length(izero)==0; 
        izero=length(pro(2).p(:,1,i)) - 1;  
    else 
        izero=izero(1);
    end
    %[a closest]=min(  abs( pro(2).p(i,1,1)-pro(5).p(:,1,1) )  );
    force(i).q(1:izero)=0;
    force(i).T(1:izero)=0;
    
    izero2=find(abs(pro(stats(i)).p(:,4,i)-100)<0.01); %for 100% rh values in sounding
    iz=min([izero izero2(1)-1]);
    
    qv=interp1(pro(5).p(:,1,1),pro(5).p(:,10,1),pro(stats(i)).p(1:iz,1,i),'linear','extrap');
    th=theta(pro(5).p(:,3,1),pro(5).p(:,2,1));
    th2=theta(pro(stats(i)).p(:,3,1),pro(stats(i)).p(:,2,1));
    
    TH=interp1(pro(5).p(:,1,1),th,pro(stats(i)).p(1:iz,1,i),'linear','extrap');
   
    
        force(i).q(1:iz)=(3-2*i)*(qv-pro(stats(i)).p(1:iz,10,i))/time(i); %so is minus for i=2
        force(i).T(1:iz)=(3-2*i)*(TH-th2(1:iz))/time(i);
        
        force(i).V(1:iz
        
        ktm(i)=izero;
    
    
end

force(4)=force(2);
force(3)=force(2);


force(5).q(1:ktm(1))=0;
force(5).T(1:ktm(1))=0;

force(2)=force(1);

outfile='g:casestudy/force.dat';
ntm=5;
ktm=max(ktm);

ftimes=[0 time(1) time(1)+(10/3600/24) time(1)+time(2) time(1)+time(2)+(10/3600/24)].*3600*24; %convert to seconds

fid=fopen(outfile,'w');

fprintf(fid,'%g\n',3600);
fprintf(fid,'%g\n%g\n%s\n',ktm,ntm,'KTMFORntmfor');

fprintf(fid,'%g\n',pro(2).p(1:ktm,1,1)-pro(2).p(1,1,1));
fprintf(fid,'%s\n','Forcing_levels');


fprintf(fid,'%g\n',ftimes);
fprintf(fid,'%s\n','Forcing_times');

fprintf(fid,'%g\n%s\n',0,'Flag_for_U_forcing');
fprintf(fid,'%g\n%s\n',0,'Flag_for_V_forcing');

fprintf(fid,'%g\n',1);
for i=1:3
	fprintf(fid,'%g\n',force(i).T);
end
fprintf(fid,'%s\n','flag_and_TH_forcing');

fprintf(fid,'%g\n',1);
for i=1:3
	fprintf(fid,'%g\n',force(i).q);
end
fprintf(fid,'%s\n','flag_and_Q01_forcing');

for i=1:10
    fprintf(fid,'%g\n%s\n',0,'comment');
end
fclose(fid);
 