function rad=getstatsSlices(npLes,npRad,timeRad,echomax,aa3,Surf2,islice,tstartles,tendles,tstartrad,tendrad)
starttimeles=9; %time at start of les runs = 9am

time=[-300:300:15*3600];
hrs2=starttimeles+time(2:end)/3600;
hrsfloor=floor(hrs2);
mins2=round((hrs2-hrsfloor)*60);

fhrs=find(hrsfloor>=tstartles & hrsfloor<tendles); %indices for required time range
ia=fhrs(1);
ib=fhrs(end);

values=[10:1.5:98.5];

n=length(npLes)+1;

Surf(1)=Surf2(2);
Surf(2)=Surf2(1);

for i=1:n-1
	rad(i).hrs=hrsfloor(fhrs(1):fhrs(end)); %first for LES data
    rad(i).mins=mins2(fhrs(1):fhrs(end));
    
    
    [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
                ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanlesSlices(npLes(i).np(:,fhrs(1):fhrs(end)),values); 
        
    rad(i).maxs=max(echomax(i).prof,[],1); 
    rad(i).echotops=max(Surf(i).echotops(:,fhrs(1)+1:fhrs(end)+1));
    rad(i).maxs2=rad(i).maxs;
    rad(i).hrs2=rad(i).hrs;
    rad(i).mins2=rad(i).mins;

end

    %npLES(1) = normal npLes(2) = ccn720


% np2=np;
% clear np;
% hdiff=hrs(2:end)-hrs(1:end-1);
% fa=find(hdiff>2); %find where step up is large signifying the next slice - usually steps from 00hrs to 14hrs
% fa(end+1)=length(np2);
% fa(2:end+1)=fa(1:end);
% fa(1)=0;

% for i=1:length(fa)-1
%     for j=fa(i)+1:fa(i+1)
%         np(i).np(:,j-fa(i))=np2(j).np; %sorts out np so that is in the same format as with les slice np
% 	end                                %i.e. np(i) i for each different slice, np(i).np(j,k) j=diff categories, k=diff times
% end

values=[10:1.5:56.5];

for i=n:n+length(islice)-1;
    fhrs=find(timeRad(islice(i-n+1)).hrs>=tstartles & timeRad(islice(i-n+1)).hrs<tendles); %indices for required time range
    a=fhrs(1);
	b=fhrs(end);
    rad(i).hrs=timeRad(islice(i-n+1)).hrs(fhrs(1):fhrs(end));
	rad(i).mins=timeRad(islice(i-n+1)).mins(fhrs(1):fhrs(end));    
    
    
    [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
            ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanlesSlices(npRad(islice(i-n+1)).np(:,fhrs(1):fhrs(end)),values);
    
    %work out echotops from slice pixel data in aa3
	
	m=10/(151-450);
	c=-450*m;
	
        rad(i).echotops=m.*aa3(i-n+1).i2max(fhrs(1):fhrs(end))+c;
        
        
        %rad(i).maxs2=;
        rad(i).hrs2=rad(i).hrs;
        rad(i).mins2=rad(i).mins;
        
%     
%     rad(i).hrs=hrs(fa(islice(i-n+1))+1:fa(islice(i-n+1)+1));
%     rad(i).mins=mins(fa(islice(i-n+1))+1:fa(islice(i-n+1)+1));
%     
%     fhrs=find(rad(i).hrs>=tstartles & rad(i).hrs<tendles); %indices for required time range
%     
%     [rad(i).rates rad(i).mid rad(i).nps rad(i).vari rad(i).mean rad(i).max rad(i).tot rad(i).diffs rad(i).means rad(i).np rad(i).vars rad(i).medians rad(i).medi...
%             ,rad(i).maxs,rad(i).modes,rad(i).mode]=meanlesSlices(np(islice(i)).np(:,fhrs(1):fhrs(end)));
%     
%     
end
    