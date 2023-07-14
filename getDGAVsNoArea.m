function data=getDGAVs(class_str,dgstr,DGAV,jj,fnmin)

dgfind=findhead(class_str,dgstr);
data=DGAV(:,dgfind(1));