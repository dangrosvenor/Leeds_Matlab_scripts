%finds required name from dgstr as made by sortheaderDGS
%usgae: j=findhead(stri,dgstr)
function j=findhead(stri,dgstr)
j=1e99;

k=1;
for i=1:length(dgstr)
%    if findstr(lower(stri),lower(dgstr{i}))
    if strmatch(lower(stri),lower(dgstr{i}),'exact')    
        j(k)=i;
        k=k+1;
    end
end


        