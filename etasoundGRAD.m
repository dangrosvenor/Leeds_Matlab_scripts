% takes eta values and makes a sounding and plots it on a tephigram

fid=fopen(strcat('c:/documents and settings/user/my documents/hibiscus/baurufield/soundings/ascii/1/30.01-16eta',direc(ii).d(end-3:end),'.txt'),'w');


ord=[1 2 3 4 5 6]; %order of AED.txt sounding :- blank blank H(m) P(mb) temp(C) rh(%) blank speed(m/s) dir(deg)
%output h p t RH u v - can just use either u or v directly


pr(ii).p(40,ord(1))=18330;
pr(ii).p(40,ord(2))=73.7;
pr(ii).p(40,ord(3))=-76.2;
pr(ii).p(40,ord(4))=pr(ii).p(39,ord(4));
pr(ii).p(40,ord(5))=pr(ii).p(39,ord(5));
pr(ii).p(40,ord(6))=pr(ii).p(39,ord(6));

pr(ii).p(41,ord(1))=19982;
pr(ii).p(41,ord(2))=55.4;
pr(ii).p(41,ord(3))=-69.8;
pr(ii).p(41,ord(4))=pr(ii).p(39,ord(4));
pr(ii).p(41,ord(5))=pr(ii).p(39,ord(5));
pr(ii).p(41,ord(6))=pr(ii).p(39,ord(6));

sp=psfc(ii); %surface pressure in hPa
pind=find(pr(ii).p(:,2)<sp); %find first pressure above surface pressure

%fprintf(fid,'%g %g\n',0.0,0.0); %two blanks
fprintf(fid,'%g %g\n',topo(ii),sp); %print surface height and surface pressure 

for i=pind(1):41
    
for j=1:size(ord,2); 
    fprintf(fid,'%g ',pr(ii).p(i,ord(j)) ); 
end

fprintf(fid,'\n');

end



%fprintf(fid,'%g %g %g %g %g %g\n',18330,73.7,-76.2,pr(ii).p(i,ord(4)),pr(ii).p(i,ord(5)),pr(ii).p(i,ord(6)));
%fprintf(fid,'%g %g %g %g %g %g',19982,55.4,-69.8,pr(ii).p(i,ord(4)),pr(ii).p(i,ord(5)),pr(ii).p(i,ord(6)));
fclose(fid);

% fid=fopen(dire,'rb');
% fscanf(fid,'%g',[1 2]);
istore=ii;
fla=3;
%admulti;
admultiF(pr(ii).p(pind(1):end,3),pr(ii).p(pind(1):end,2),pr(ii).p(pind(1),4),pr(ii).p(pind(1):end,1));
clear ii
ii=istore;
%fclose(fid);