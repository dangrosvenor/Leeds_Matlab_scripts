function [rates,midrates,nps,variance,meanval]=meanles(vals,data)

minrate=interp1(reflect,rrates,values(1),'spline');


mini=find(Surf.instant>=minrate);
meles=mean(Surf.instant(mini));

[nples ratevalsles]=count(ratevals,Surf.instant);
rmsles=vari(meles,ratevals,nples);
midvalsles=(ratevals(2:end)+ratevals(1:end-1))/2;
maxins=max(max(Surf.instant));
ratediffs=ratevals(2:end)-ratevals(1:end-1);

function [np,valsnew]=count(vals,dat)
    valsnew=vals;
	for i=1:size(vals,2)-1
          cati=find(dat>vals(i) & dat<vals(i+1));
          np(i)=length(cati);
	end
    maxdat=max(max(dat));
    if maxdat>vals(end)
        cati=find(dat>vals(end));
        np(size(np)+1)=length(cati);
        valsnew(end+1)=maxins;
    end
    
function v=vari(mean,vals,nps)
	tot=sum(nps);
	sqdiff=(mean-vals).^2;
	v=sum(nps.*sqdiff(1:length(nps)))./(tot-1);
