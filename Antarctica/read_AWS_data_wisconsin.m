clear filename AWSdat_wis

days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];
days_in_month_cum = [0 cumsum(days_in_month(1:11))];

%data for Larsen AWS from http://www.antarctica.ac.uk/cgi-bin/metdb-form-2.pl?tabletouse=U_MET.LARSEN_ICE_SHELF_AWS_10&complex=1&idmask=.....&acct=u_met&pass=weather

% (1) Kirkwood Island - Lat : 68.34S  Long :  69.01W    Elev :   30 M  %wind looks ok
% (2) Dismal Island - Lat : 68.09S  Long :  68.82W      Elev :   10 M  %bad wind data
% (3) Bonaparte Point -  Lat : 64.78S  Long :  64.07W   Elev :    8 M  %bad wind data
% (4) Sky Blu   -   Lat : 74.79S  Long :  71.49W        Elev : 1510 M  %wind looks ok - actually says 1700 m at http://ice.ssec.wisc.edu/aws/skyblumain.html
% (5) Limbert   -   Lat : 75.91S  Long :  59.26W        Elev :   40 M  %wind looks ok
% (6) Butler Island  -  Lat : 72.21S  Long :  60.17W    Elev :   91 M  %not sure about wind data - think is ok
% (7) Larsen BAS     -  Lat : 67.01S  Long :  61.55W    Elev :   17 M  (also says 61.33 on same webiste!) 


LAT_AWS(1).dat = -68.34;  LON_AWS(1).dat =  -69.01;    
LAT_AWS(2).dat = -68.09;  LON_AWS(2).dat =  -68.82;    
LAT_AWS(3).dat = -64.78;  LON_AWS(3).dat =  -64.07;   
LAT_AWS(4).dat = -74.79;  LON_AWS(4).dat =  -71.49;   
LAT_AWS(5).dat = -75.91;  LON_AWS(5).dat =  -59.26;   
LAT_AWS(6).dat = -72.21;  LON_AWS(6).dat =  -60.17;
LAT_AWS(7).dat = -67.01;  LON_AWS(7).dat =  -61.55;

ALT_AWS(1).dat = 30;
ALT_AWS(2).dat = 10;
ALT_AWS(3).dat = 8;
ALT_AWS(4).dat = 1597; %get best pressure match if use 1597 m. Says 1700m elevation at http://ice.ssec.wisc.edu/aws/skyblumain.html but altitude of 1510 m - ?
ALT_AWS(5).dat = 40;
ALT_AWS(6).dat = 91;
ALT_AWS(7).dat = 17;

aws_name(1).name='Kirkwood';filename(1).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/Kirkwood_island_jan06.txt';
aws_name(2).name='Dismal Island';filename(2).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/Dismal_island_jan06.txt';
aws_name(3).name='Bonaparte Point';filename(3).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/Bonaparte_jan06.txt';
aws_name(4).name='Sky Blu';filename(4).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/SkyBlu_jan06.txt';
aws_name(5).name='Limbert';filename(5).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/Limbert_jan06.txt';
aws_name(6).name='Butler Island';filename(6).name='C:\Documents and Settings\dan\My Documents\WRF\aws data/ButlerIsland_jan06.txt';
aws_name(7).name='Larsen C';

for i=1:length(filename)

    fid=fopen(filename(i).name,'rt');
    AWS_text=fgetl(fid);
    AWS_text2=fgetl(fid);

    AWSdat=fscanf(fid,'%f'); %make sure that when requesting the data from the BAS website that
    %select 'type your own' date format and add a space before the year
    %otherwise it runs the last data point for each pressure into the year
    num_cols=8;
    AWSdat_wis(i).dat=reshape(fliplr(AWSdat),[num_cols length(AWSdat)/num_cols]);

    AWSdat_wis(i).dat=fliplr(AWSdat_wis(i).dat); %data is written in reverse time order so flip it here


    i444=find(AWSdat_wis(i).dat==444);
    AWSdat_wis(i).dat(i444) = NaN;
    AWStime_wis(i).dat = AWSdat_wis(i).dat(1,:)+AWSdat_wis(i).dat(2,:)*10/60/24;

end

%fields are
%1)Day
%2)Time index (every 10 mins)
%3)Temperature
%4)Pressure (mb)
%5)Wind speed
%6)Wind dir 
%7)RH? 
%8)ALL NULL (=444)


filename='C:\Documents and Settings\dan\My Documents\WRF\AWS_ant_jan06.txt';
filename='C:\Documents and Settings\dan\My Documents\WRF\AWS data\Larsen_whole_jan_feb06.txt';

fid=fopen(filename,'rt');
AWS_text=fgetl(fid);
AWSdat=fscanf(fid,'%f'); %make sure that when requesting the data from the BAS website that
                      %select 'type your own' date format and add a space before the year
                      %otherwise it runs the last data point for each pressure into the year
                      

inan=find(AWSdat==99999);
AWSdat(inan) = NaN;

num_cols=8;
AWSdat_wis_array=reshape(fliplr(AWSdat),[num_cols length(AWSdat)/num_cols]); 


AWSdat_wis_array=fliplr(AWSdat_wis_array); %data is written in reverse time order so flip it here


%calculate hourly averages
time = days_in_month_cum(AWSdat_wis_array(2,:)) + AWSdat_wis_array(3,:)+AWSdat_wis_array(4,:)/24; %time in days but keeping data from the same hour together
[uni_arr, i_uni]=unique(time);

%AWStime_wis(7).dat = uni_arr;

%data should be every 10 minutes but some of the time slots are missing
clear AWSdat3
inon6=1;
for iav=1:length(uni_arr)
    imean=find(time==uni_arr(iav));
%    length(imean)
    if length(imean)~=6;
        non6(inon6,1)=uni_arr(iav);
        non6(inon6,2)=length(imean);
        inon6=inon6+1;
    end
    
    AWSdat3(:,iav)=mean(AWSdat_wis_array(:,imean),2);
    %BUT, will run into problems when doing averages of the wind direction
end


%time2 = days_in_month_cum(AWSdat_wis_array(2,:)) + AWSdat_wis_array(3,:)+AWSdat_wis_array(4,:)/24;

AWStime_wis(7).dat = uni_arr;

AWSdat_wis(7).dat(3,:) = AWSdat3(5,:);
AWSdat_wis(7).dat(4,:) = AWSdat3(6,:);
AWSdat_wis(7).dat(5,:) = AWSdat3(7,:);
AWSdat_wis(7).dat(6,:) = AWSdat3(8,:);

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
                      
