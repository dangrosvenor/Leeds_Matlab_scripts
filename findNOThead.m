%finds required name from dgstr as made by sortheaderDGS
function j=findhead(stri,dgstr,DGAV)
k=1;
for i=1:size(DGAV,2)
    if length(findstr(stri,dgstr{i}))==0
        j(k)=i;
        k=k+1;
    end
end

        