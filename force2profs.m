clear force

isave=1;
kkp=150;
time(1)=(20.25-12)/24; %1st sounding at 12UTC, 2nd at 20:15 (units = day), 3rd at 24:00
time(2)=(24-12)/24; 
stats=[2 5 3]; %stations to be used for 1st and 3rd soundings
timestat=[1 1 2];

for i=1:length(stats)

    [hMAX(i) izero]=max(pro(stats(i)).p(:,1,timestat(i)));
    hMIN(i)=pro(stats(i)).p(1,1,timestat(i));
    %hMAX(i)=imax(1);
    %izero=find( abs(pro(stats(i)).p(:,1,timestat(i))-hMAX)<0.01);
    %diff=pro(stats(i)).p(2:end,1,1)-pro(stats(i)).p(1:end-1,1,1);
    %izero=find(abs(diff-0)<0.01); %point where data ends in sounding
    %if length(izero)==0; 
%         izero=length( pro(stats(i)).p(:,1,timestat(i))) - 1;
%         izerST=izero;
%     else 
%         izero=izero(1);
%       %  hMAX(i)=pro(stats(i)).p(izero,1,1); %max height point of sounding
%     end
    
    


%first reconstruct first sounding filling in gaps using Bauru sounding
	orig(i).p=pro(stats(i)).p(1:izero,:,timestat(i));    
  
  if stats(i)~=5;
	sat=find(abs(orig(i).p(:,4)-100)<0.01); %make sure that any points that are saturated are not to be kept
	for j=1:length(sat)
        ii=sat(j);
        orig(i).p(ii,10)=interp1(pro(5).p(:,1,1),pro(5).p(:,10,1),pro(stats(i)).p(ii,1,1));
        orig(i).p(ii,4)=interp1(pro(5).p(:,1,1),pro(5).p(:,4,1),pro(stats(i)).p(ii,1,1));
	end
    
    zrows=find(abs(orig(i).p(:,8))<0.01);
    for j=1:length(zrows)
        ii=zrows(j);
        orig(i).p(ii,8)=interp1(pro(5).p(:,1,1),pro(5).p(:,8,1),orig(i).p(ii,1),'linear','extrap'); %mag
        orig(i).p(ii,9)=interp1(pro(5).p(:,1,1),pro(5).p(:,9,1),orig(i).p(ii,1),'linear','extrap'); %dir
	end
    
 
  end

end %for i



h0=max(hMIN+0.01); %first height in sounding
hm=min(hMAX-0.01);
hts=[h0:(hm-h0)/(kkp-1):hm]; %heights to aim for from first sounding

for i=1:length(stats)
        force(i).q=interp1(orig(i).p(:,1),orig(i).p(:,10),hts);
        
        dircomp=135; %direction at which to take wind component (degrees)
        V=orig(i).p(:,8).*-cos((orig(i).p(:,9)-dircomp)*pi/180); %take component of speed at dircomp degrees
        force(i).V=interp1(orig(i).p(:,1),V,hts);
        
        th2=theta(orig(i).p(:,3),orig(i).p(:,2));
        force(i).TH=interp1(orig(i).p(:,1),th2,hts);
end
    
    


outfile='g:casestudy/force.dat';
ntm=3;
ktm=kkp;

ftimes=[0 time(1) time(2)].*3600*24; %convert to seconds

fid=fopen(outfile,'w');

fprintf(fid,'%g\n',1800); %relaxation time (s)
fprintf(fid,'%g\n%g\n%s\n',ktm,ntm,'KTMFORntmfor');

fprintf(fid,'%g\n',hts-h0);
fprintf(fid,'%s\n','Forcing_levels');


fprintf(fid,'%g\n',ftimes);
%fprintf(fid,'%s\n');


fprintf(fid,'%s\n%g\n','Flag_for_U_forcing',0);
fprintf(fid,'%s\n%g\n','Flag_for_V_forcing',2);
for i=1:3
	fprintf(fid,'%g\n',force(i).V);
end

fprintf(fid,'%s\n','flag_and_TH_forcing');
fprintf(fid,'%g\n',2);
for i=1:3
	fprintf(fid,'%g\n',force(i).TH);
end

fprintf(fid,'%s\n','flag_and_Q01_forcing');
fprintf(fid,'%g\n',2);
for i=1:3
	fprintf(fid,'%g\n',force(i).q);
end

for i=1:10
    fprintf(fid,'%s\n%g\n','comment',0);
end
fclose(fid);


savedir='c:\documents and settings\Login\my documents\HIBISCUS\troccibras\radiosondesnetwork\CampoGrande_040224_12_force2profs';
if isave==1
	fid=fopen(savedir,'wt');
	for i=1:size(orig(1).p,1)
        fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',999,999,orig(1).p(i,1),orig(1).p(i,2),orig(1).p(i,3),orig(1).p(i,4),999,orig(1).p(i,9),orig(1).p(i,8));
	end
    fclose(fid);
end

 