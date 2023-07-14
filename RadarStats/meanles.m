function [rates,midrates,nps,variance,meanval,maxi,tot,diff,means,np,vars,medians,medi,maxs,modes,mode]=meanles(data,type)

rrates=[.15 .3 .6 1.3 2.7 5.6 12 24 49 100 205]; %rain rates corresponding to radar reflectivities
reflect=[10 15 20 25 30 35 40 45 50 55 60]; %reflectivites for above
values=[13.5:1.5:60];



if type==1 %echotops
    vals=[2.0:0.65:20.2 99 999 9999]; %vals for echotops  
elseif type==2 | type==4 %if ppi
    minrate=interp1(reflect,rrates,values(1),'spline');
    vals=interp1(reflect,rrates,values,'spline');
else %cappi
    vals=values;
end

rates=vals;
midrates=(rates(2:end)+rates(1:end-1))/2;


%mini=find(data>=minrate);
%meanval=mean(data(mini));
%median=median(data(mini));
if length(size(data))==3
    [np rates]=count(vals,data);
    maxi=max(max(data));
else
    for i=1:length(data(1).np)-1
        np(i)=0;
        for j=1:length(data)
            np(i)=np(i)+data(j).np(i);
        end
    end
    rates=vals;
    f0=find(np>0);
    maxi=vals(f0(end));
end
meanval=nmean(np,rates);
medi=nmedian(np,rates);
variance=vari(meanval,rates,np);
mode=nmode(np,rates);



tot=sum(np);
diff=rates(2:end)-rates(1:end-1);

[means,nps,vars,medians,maxs,modes]=alltimes(vals,data);





%ratediffs=ratevals(2:end)-ratevals(1:end-1);

function [np,valsnew]=count(vals,dat)
    
     
    valsnew=vals;
	for i=1:size(vals,2)-1
          cati=find(dat>vals(i) & dat<vals(i+1));
          np(i)=length(cati);
	end
    maxdat=max(max(dat));
    if maxdat>vals(end)
        cati=find(dat>vals(end));
        np(length(np)+1)=length(cati);
        valsnew(end+1)=maxdat;
    end
    

    
function  [means,nps,vars,medians,maxs,modes]=alltimes(vals,dat)
	if length(size(dat))==3
        for i=1:size(dat,3)
            for j=1:size(vals,2)-1
                cati=find(dat(:,:,i)>vals(j) & dat(:,:,i)<vals(j+1));
                nps(i,j)=length(cati);
            end
            means(i)=nmean(nps(i,:),vals);
            vars(i)=vari(means(i),vals,nps(i,:));
            medians(i)=nmedian(nps(i,:),vals);
            maxs(i)=max(dat(:,1,i));
%             diffs=vals(2:end)-vals(1:end-1);
%             prob=nps(i,:)./diffs;
%             [maxprob imax]=max(prob);
%             modes(i)=vals(imax);
              modes(i,:)=nmode(nps(i,:),vals);
    end
	else
        for i=1:length(dat)
            for j=1:length(dat(1).np)-1
                nps(i,j)=dat(i).np(j);
            end    
            means(i)=nmean(nps(i,:),vals);
            vars(i)=vari(means(i),vals,nps(i,:));
            medians(i)=nmedian(nps(i,:),vals);
            f0=find(nps(i,:)>0);
            maxs(i)=vals(f0(end));
%             diffs=vals(2:end)-vals(1:end-1);
%             prob=nps(i,:)./diffs;
%             [maxprob imax]=max(prob);
%             modes(i)=vals(imax);
            modes(i,:)=nmode(nps(i,:),vals);
        end
	end
    

