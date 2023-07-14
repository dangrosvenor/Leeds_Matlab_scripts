function med=nmedian(nps,vals)

midvals=(vals(1:end-1)+vals(2:end))/2;

for i=1:length(midvals)
    x(i).x=repmat(midvals(i),[1 nps(i)]);
end

xx=[x.x];
med=median(xx);
if size(med)==0;med=NaN;end
