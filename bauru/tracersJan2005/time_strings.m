function [it1,it2,t1str,t2str]=time_strings(t1,t2,times)

    it1=findheight(times,t1);
    mins=(t1-floor(t1))*60;
    minstr=num2str(mins,'%2.0f');
    
    hrs=mod(floor(t1),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t1str=[hrstr ':' minstr];
    
    it2=findheight(times,t2);
    mins=(t2-floor(t2))*60;
    minstr=num2str(mins,'%2.0f');
    hrs=mod(floor(t2),24);
    hrstr=num2str(hrs,'%2.0f');
    if mins==0; minstr='00';end
    if hrs==0; hrstr='00';end
    t2str=[hrstr ':' minstr];
