function mode=nmode(np,vals)

midvals=(vals(1:end-1)+vals(2:end))/2;
diffs=vals(2:end)-vals(1:end-1);
totn=sum(np);
if totn==0
    mode=0;
else
    prob=np./diffs;
    [maxprob imax]=max(prob);
    mode=midvals(imax); %note, is possible to have more than one mode so make output array at least 2-D
end