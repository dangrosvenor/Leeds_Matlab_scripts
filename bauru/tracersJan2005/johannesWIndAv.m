n=2; %no. of horizontal point to average over
id=1; %directory number
smat=size(wind(id).W);
for i=1:floor(smat(2)/n)
    ii=(i-1)*n+1;
    av(1:smat(1),i,1:smat(3))=mean(wind(id).W(:,ii:ii+n-1,:),2);
end