function vapad=adratef(vapadcum,dumprange) 
%input cumulative sum of data - e.g. vapadcum=changevap + cumsum(microicerate,2)*300;
%and returns the rate of change of it with an extra column added at the beginning to make it the same size as
%other arrays

vapad=(vapadcum(:,2:end)-vapadcum(:,1:end-1))/300; %calculate rate of change of advective source
vapad(:,dumprange(2):dumprange(end))=vapad(:,1:end);
vapad(:,1)=0;