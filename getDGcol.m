function [dgcol,error]=getDGAVs(class_str,dgstr,error)

dgfind=findhead(class_str,dgstr);

if dgfind~=1e99 %this value returned by findhead if not found
	dgcol=dgfind(1);
else
    dgcol=1;    %if not found the item looking for return 1 in case want a dummy variable
    error=error+1; 
end