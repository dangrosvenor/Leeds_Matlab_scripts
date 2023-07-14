function rad=getstats(iecho,Surf,np,date,day,month,year,hrs,mins,ia,ib,a,b)

n=length(Surf)+1;
%iecho=2; %1=echotops 2=ppi 3=cappi


for i=1:n-1
    if iecho==1
     [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
            ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanles(Surf(i).echotops(:,1,ia:ib),1);
elseif iecho==2
        [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
            ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanles(Surf(i).instant(:,1,ia:ib),2);
elseif iecho==3 | iecho==5 %for 3.5km and max cappi
    [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
            ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanles(Surf(i).cappi(:,1,ia:ib),3);
elseif iecho==4
    [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
            ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanles(Surf(i).pmax(:,1,ia:ib),3);
    end
    %ydat(i).y=rad(i).means;
    %ydat(i).y=squeeze(mean(Surf(i).instant,1));
%     ydat(i).y=rad(i).medians;
%     xdat(i).x=time(1:length(ydat(i).y));
%     xdat(i).x=datenum(2004,2,24,9+floor(xdat(i).x/3600),60*(xdat(i).x/3600-floor(xdat(i).x/3600)),0); 
end

i=n;



[rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
        ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanles(np(a:b),iecho);

% xdat(n).x=datenum(year,month,day,hrs,mins,0);
% xdat(n).x=xdat(n).x(27:end);
% %ydat(n).y=squeeze(meany2(27:end));
% %ydat(n).y=rms(27:end);
% 
% ydat(n).y=rad(n).vars;
