%one used for e.g. force+0_3th3qv_2t - 9th March,2005.
%forcing for e.g. Campo Grande sounding with top added from Bauru sounding (simpler) 21st Dec 2004

%run readsound with:-
%dire(1).d='c:\documents and settings\Login\my documents\HIBISCUS\baurufield\soundings\ASCII\5\04_02_24_2015_dmi.AED.txt';
%dire(2).d='c:\documents and settings\Login\my documents\HIBISCUS\TrocciBras\radiosondesNetwork\CampoGrande_040224_12_force2profs';

%to get soundings in pr(1) and pr(2)

topheight=9631; %last height (m) for forcing (height after which Bauru sounding was used to fill in RH) 
nqp=54; %no. q fields

outfile='c:/cygwin/home/login/force.dat';

clear force

  

time(1)=(20.25-12)/24; %1st sounding at 12UTC, 2nd at 20:15 (units = day)
time(2)=(24-20.25)/24; 

ktm=size(pr(2).p,1);
force(1).q(1:ktm)=0;
force(1).T(1:ktm)=0;
    
iz=find(pr(2).p(:,1)>topheight+0.01);
iz=iz(1);

    
qv=interp1(pr(1).p(:,1),pr(1).p(:,10),pr(2).p(1:iz,1),'linear','extrap'); %Bauru values interpolated so that are on same grid as CampoGrande grid
th=theta(pr(1).p(:,3),pr(1).p(:,2)); %Bauru theta
th2=theta(pr(2).p(:,3),pr(2).p(:,2)); %CG theta

TH=interp1(pr(1).p(:,1),th,pr(2).p(1:iz,1),'linear','extrap'); %Bauru values interpolated so that are on same grid as CampoGrande grid


force(1).q(1:iz)=(qv-pr(2).p(1:iz,10))/time(1); % ( Bauru - CG ) / time
force(1).T(1:iz)=(TH-th2(1:iz))/time(1);

    
    


force(2)=force(1);


force(3).q(1:ktm(1))=0;
force(3).T(1:ktm(1))=0;



ntm=3;


ftimes=[0 time(1) time(1)+(10/3600/24)].*3600*24; %convert to seconds

fid=fopen(outfile,'w');

fprintf(fid,'%g\n',3600);
fprintf(fid,'%g\n%g\n%s\n',ktm,ntm,'KTMFORntmfor');

fprintf(fid,'%g\n',pr(2).p(1:ktm,1)-pr(2).p(1,1));
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

for i=1:nqp-1
    fprintf(fid,'%g\n%s\n',0,'comment');
end
fclose(fid);
 