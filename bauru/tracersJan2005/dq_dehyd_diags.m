'dq dehyd'
f=1e6*28.97/18;

totw=sum(TwoDDan(jc).Q(:,:,[1:6]),3); %jc is the number of the eg. GridDan(jc) where the data is stored
vap=sum(TwoDDan(jc).Q(:,:,[1]),3);   %j just loops from 1:n

j
prcs=[0:5:100];
tot_prctiles(j).t(1:size(totw,1),jj,1:length(prcs))=(prctile(totw',prcs))';
vap_prctiles(j).t(1:size(vap,1),jj,1:length(prcs))=(prctile(vap',prcs))';

minpps=[3.67 5 1 2 3 4]/f;
%minpps=3.67/f;

for ipps=1:length(minpps)

    for ikm=1:size(totw,1)
        inon=find(totw(ikm,:)<minpps(ipps));
        dq(j).d(ikm,jj,ipps)=sum(minpps(ipps)-totw(ikm,inon)) * f / size(totw,2);
        nn(j).n(ikm,jj,ipps)=length(inon);

        
        inon2=find(vap(ikm,:)<minpps(ipps));
        dq2(j).d(ikm,jj,ipps)=sum(minpps(ipps)-vap(ikm,inon2)) * f / size(vap,2);
        nn2(j).n(ikm,jj,ipps)=length(inon2);

    end    

end
 
 icediags_5thSept_2005_32;