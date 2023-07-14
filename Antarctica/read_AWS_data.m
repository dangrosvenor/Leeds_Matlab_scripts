filename='C:\Documents and Settings\dan\My Documents\WRF\AWS_ant_jan06.txt';
fid=fopen(filename,'rt');
AWS_text=fgetl(fid);


AWSdat=fscanf(fid,'%f'); %make sure that when requesting the data from the BAS website that
                      %select 'type your own' date format and add a space before the year
                      %otherwise it runs the last data point for each pressure into the year
num_cols=8;
AWSdat2=reshape(fliplr(AWSdat),[num_cols length(AWSdat)/num_cols]); 

AWSdat2=fliplr(AWSdat2); %data is written in reverse time order so flip it here

AWShrs = AWSdat2(4,:)+AWSdat2(3,:)*24;

days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];
days_in_month_cum = [0 cumsum(days_in_month(1:11))];

time = days_in_month_cum(AWSdat2(2,:)) + AWSdat2(3,:)+AWSdat2(4,:)/24; %time in days but keeping data from the same hour together
[uni_arr, i_uni]=unique(time);

%AWStime_wis(7).dat = uni_arr;

%data should be every 10 minutes but some of the time slots are missing
clear AWShrs2
inon6=1;
itot=0;
for iav=1:length(uni_arr) %loop through all the unique times
    imean=find(time==uni_arr(iav)); %find times with each unique time
    itot = itot+length(imean);
    Lmean=length(imean)
    if Lmean~=6;  %if there aren't 6 values then store the time and the number we do have
        non6(inon6,1)=uni_arr(iav);
        non6(inon6,2)=Lmean;
        inon6=inon6+1;
    end
    
%    AWSdat3(:,iav)=mean(AWSdat2(:,imean),2); %do mean over the hour
%here are assuming that the missing times all come at the end of the hour - hopefully won't
%matter too much
    AWShrs2(1+itot-Lmean:itot)=uni_arr(iav)*24 + 10/60*[0:Lmean-1];
end


%AWShrs2 = AWSdat2(4,1)+AWSdat2(3,1)*24 + 10/60 * [0:size(AWSdat2,2)-1];

%fields are
%1)Year
%2)Month
%3)Day
%4)Hour
%5)TEMPERATURE 
%6)PRESSURE 
%7)WIND_SPEED 
%8)WIND_DIRECTION 
%HUMIDITY values were all null.....

%Lat : 67.01S  Long :  61.55W  (also says 61.33 on same webiste!) 
%Elev :   17 M  Taken from website - but not
%sure how accurate - 17m abmsl?




'done read of AWS data'                      
                      
