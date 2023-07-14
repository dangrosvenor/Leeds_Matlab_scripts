%write a sounding

QnotRH=1; %flag to output q instead of RH - use for soundings where are adding temp and don't want to affect moisture MR
les=0; %flag for LES wind output (v component already taken)
firstline=1; %flag for outputting first line of pr profiles =1 for MilesCity =0 for others

%outfile22='c:\documents and settings\Login\my documents\HIBISCUS\TrocciBras\radiosondesNetwork\DMI_1715_5ppmv_16.6-19km';
%outfile22='c:\documents and settings\Login\my documents\HIBISCUS\TrocciBras\radiosondesNetwork\CampoGrande_040224_12_2';
outfile22='c:\documents and settings\Login\my documents\Leeds_MMOCCA\MileCity.AED';
%outfile22='c:\documents and settings\Login\my documents\HIBISCUS\TrocciBras\radiosondesNetwork\CampoGrande_040224_12_2';
outfile22='c:/documents and settings/login/my documents/hibiscus/troccibras/radiosondesnetwork/lem_ecmwf_hightop'

fid=fopen(outfile22,'wb');

first=2-firstline;

for i=first:size(pr(3).p,1)

    if QnotRH==0
        fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',999,999,pr(3).p(i,1),pr(3).p(i,2),pr(3).p(i,3),pr(3).p(i,4),999,pr(3).p(i,9),pr(3).p(i,8));
    elseif les==1
        fprintf(fid,'%g %g %g %g %g %g %g %g\n',999,999,pr(3).p(i,1),pr(3).p(i,2),pr(3).p(i,3),pr(3).p(i,10),999,pr(3).p(i,8));
        %writes from les output pr(3).p(i,8)=v component of wind 
        %run soundLES first
    else
        fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',999,999,pr(3).p(i,1),pr(3).p(i,2),pr(3).p(i,3),pr(3).p(i,10),999,pr(3).p(i,9),pr(3).p(i,8));
    end  
  
end
fclose(fid);

%run WriteSound.f after with firstline=1;

'done'