clear force

time(1)=(18-12)/24; %1st sounding at 12UTC, 2nd at 20:15 (units = day)
%time(2)=(24-20.25)/24; 
stats=[5]; %stations to be used for 1st and 3rd soundings

origstat=5;

itim=[2]; %index i for time values in pro.p(:,:,i) for each force(i)

for i=1:length(itim)
    diff=pro(stats(i)).p(2:end,1,itim(i))-pro(stats(i)).p(1:end-1,1,itim(i));
    izero=find(abs(diff-0)<0.0001); %point where data ends in sounding
    if length(izero)==0; 
        izero=length(pro(origstat).p(:,1,itim(i))) - 1;  
    else 
        izero=izero(1);
    end
    %[a closest]=min(  abs( pro(2).p(i,1,1)-pro(5).p(:,1,1) )  );
    force(i).q(1:izero)=0;
    force(i).T(1:izero)=0;
    
    
        
    izero2=find(abs(pro(stats(i)).p(:,4,itim(i))-100)<0.01); %for 100% rh values in sounding
    if length(izero2==0);izero2=9e9;end
    iz=min([izero-1 izero2(1)-1]);
    
    if iz>150
        newi=round([1:iz/150:iz]);
        newheights=pro(stats(i)).p(newi,1,itim(i));
    else
        newheights=pro(stats(i)).p(1:iz,1,itim(i));
    end
    
    qv=interp1(pro(origstat).p(:,1,1),pro(origstat).p(:,10,1),newheights,'linear','extrap'); %vap for orig sounding at new heights
    th=theta(pro(origstat).p(:,3,1),pro(origstat).p(:,2,1));
    TH2=theta(pro(stats(i)).p(:,3,itim(i)),pro(stats(i)).p(:,2,itim(i)));
    th2=interp1(pro(stats(i)).p(1:iz,1,itim(i)),TH2(1:iz),newheights,'linear','extrap');
    TH=interp1(pro(origstat).p(:,1,1),th,newheights,'linear','extrap'); %Theta for orig at new heights
   
    qv2=interp1(pro(stats(i)).p(1:iz,1,itim(i)),pro(stats(i)).p(1:iz,10,itim(i)),newheights,'linear','extrap');
    
        force(i).q= (qv2 - qv) / time(i); 
        force(i).T= (th2 - TH) / time(i);
        
        
        ktm(i)=length(newheights);
    
    
end

zlevs=newheights-pro(origstat).p(1,1,1);

force(3).q(1:ktm(1))=0;
force(3).T(1:ktm(1))=0;

force(2).q=force(1).q;
force(2).T=force(1).T;

outfile='g:casestudy/force.dat_13_02';
ntm=3;
ktm=ktm(1);

ftimes=[0 time(1) time(1)+(10/3600/24)].*3600*24; %convert to seconds

fid=fopen(outfile,'w');

fprintf(fid,'%g\n',3600);
fprintf(fid,'%g\n%g\n%s\n',ktm,ntm,'KTMFORntmfor');

fprintf(fid,'%g\n',zlevs);
fprintf(fid,'%s\n','Forcing_levels');


fprintf(fid,'%g\n',ftimes);

fprintf(fid,'%s\n%g\n','Flag_for_U_forcing',0);
fprintf(fid,'%s\n%g\n','Flag_for_V_forcing',0);

fprintf(fid,'%s\n','flag_and_TH_forcing');
fprintf(fid,'%g\n',1);
for i=1:3
	fprintf(fid,'%g\n',force(i).T);
end

fprintf(fid,'%s\n','flag_and_Q01_forcing');
fprintf(fid,'%g\n',1);
for i=1:3
	fprintf(fid,'%g\n',force(i).q);
end


for i=1:10
    fprintf(fid,'%s\n%g\n','comment',0);
end
fclose(fid);
 