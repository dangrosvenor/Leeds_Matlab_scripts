function data=getDGAVs(class_str,dgstr,DGAV,jj,fnmin)


area_str=strcat(class_str(1:3),'_A');
dgfind=findhead(area_str,dgstr);
area=DGAV(:,dgfind(1));
areazeroes=find(area==0);
area(areazeroes)=1;

dgfind=findhead(class_str,dgstr);
data=DGAV(:,dgfind(1))./area;
%data=DGAV(:,dgfind(1));